# ============================================================
# BAR-HEADED GOOSE (Anser indicus): relative abundance, India winter
# Hurdle random forest (encounter rate x count given presence),
# eBird data only (no habitat covariates).
# Output = expected count on a STANDARDISED checklist (relative abundance),
# NOT a population size.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(ranger)
  library(mgcv)
  library(pROC)
  library(sf)
  library(rnaturalearth)
  library(ggplot2)
})

# ---- USER INPUTS ----------------------------------------------------------
data_dir <- "D:/EBird project"
smp_file <- file.path(data_dir, "ebd_IN_bahgoo_201511_202511_smp_relOct-2025_sampling.txt")
obs_file <- file.path(data_dir, "ebd_IN_bahgoo_201511_202511_smp_relOct-2025.txt")
out_dir  <- "outputs"

species_common <- "Bar-headed Goose"

# ---- PARAMETERS -----------------------------------------------------------
months_keep      <- c(11, 12, 1, 2, 3)      # Nov-Mar
protocol_keep    <- c("Traveling", "Stationary")
max_duration_min <- 300
max_distance_km  <- 5
max_n_observers  <- 10
cell_km          <- 5          # spatial thinning / prediction grid cell (km)
block_km         <- 100        # spatial blocks for validation and uncertainty
test_frac        <- 0.2        # fraction of BLOCKS held out for validation
support_km       <- 25         # predict only where a checklist exists within this distance
winsor_prob      <- 0.99       # cap extreme flock sizes in the count model
n_trees          <- 500
n_boot           <- 20         # spatial-block subsampling replicates
boot_frac        <- 0.8        # fraction of blocks kept in each replicate
boot_trees       <- 200
n_boot_cells     <- 10000      # cells used to summarise each replicate
predict_dates    <- seq(ymd("2023-11-05"), ymd("2024-03-24"), by = "7 days")
std_start_hour   <- 7
std_hours        <- 1
std_km           <- 1
std_observers    <- 1
nth              <- max(1, parallel::detectCores() - 2)
quick            <- identical(Sys.getenv("QUICK"), "1")   # small test run (subsample, few trees)

if (quick) { n_trees <- 100; n_boot <- 3; boot_trees <- 50 }

dir.create(out_dir, showWarnings = FALSE)
set.seed(2025)
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

# ============================================================
# 1) READ CHECKLISTS, FILTER, DE-DUPLICATE
# ============================================================
cache1 <- file.path(out_dir, "zf_cache.rds")
if (file.exists(cache1)) {
  log_msg("Loading cached zero-filled data ...")
  zf <- readRDS(cache1)
} else {
  log_msg("Reading sampling-event file ...")
  smp_cols <- c("COUNTRY CODE", "LATITUDE", "LONGITUDE", "OBSERVATION DATE",
                "TIME OBSERVATIONS STARTED", "OBSERVER ID", "SAMPLING EVENT IDENTIFIER",
                "PROTOCOL NAME", "DURATION MINUTES", "EFFORT DISTANCE KM",
                "NUMBER OBSERVERS", "ALL SPECIES REPORTED", "GROUP IDENTIFIER")
  smp <- fread(smp_file, sep = "\t", quote = "", na.strings = c("", "NA"),
               select = smp_cols, showProgress = FALSE)
  setnames(smp, c("country_code", "latitude", "longitude", "observation_date",
                  "time_start", "observer_id", "sampling_event_id", "protocol",
                  "duration_min", "distance_km", "n_observers", "all_reported",
                  "group_id"))
  log_msg("Checklists in file: ", format(nrow(smp), big.mark = ","))

  # Observer experience = number of India checklists per observer (before filtering)
  obs_n <- smp[, .(obs_checklists = .N), by = observer_id]
  smp <- merge(smp, obs_n, by = "observer_id", all.x = TRUE)
  smp[, obs_exp := log10(obs_checklists)]

  smp[, date := as.Date(observation_date)]
  smp[, `:=`(month = month(date), year = year(date),
             day_of_year = yday(date), iso_week = isoweek(date))]
  smp[, start_hour := as.integer(substr(time_start, 1, 2)) +
        as.integer(substr(time_start, 4, 5)) / 60]
  # Stationary counts have no travel distance: set to 0 (do NOT drop them)
  smp[protocol == "Stationary", distance_km := 0]

  smp <- smp[country_code == "IN" &
               month %in% months_keep &
               protocol %in% protocol_keep &
               all_reported == 1 &
               !is.na(duration_min) & duration_min > 0 & duration_min <= max_duration_min &
               !is.na(distance_km) & distance_km <= max_distance_km &
               !(protocol == "Traveling" & distance_km <= 0) &
               !is.na(n_observers) & n_observers >= 1 & n_observers <= max_n_observers &
               !is.na(start_hour) & !is.na(latitude) & !is.na(longitude)]
  smp[, effort_hours := duration_min / 60]

  # Shared (group) checklists are copies of one event: keep one per group
  smp[, key := fifelse(is.na(group_id), sampling_event_id, group_id)]
  n_before <- nrow(smp)
  smp <- smp[sample.int(.N)][!duplicated(key)]
  log_msg("Checklists after filters: ", format(n_before, big.mark = ","),
          " | after removing shared-checklist duplicates: ",
          format(nrow(smp), big.mark = ","))

  # ============================================================
  # 2) READ SPECIES OBSERVATIONS (same key), ZERO-FILL
  # ============================================================
  log_msg("Reading observation file ...")
  obs_cols <- c("COMMON NAME", "OBSERVATION COUNT", "SAMPLING EVENT IDENTIFIER",
                "GROUP IDENTIFIER", "COUNTRY CODE")
  obs <- fread(obs_file, sep = "\t", quote = "", na.strings = c("", "NA"),
               select = obs_cols, colClasses = list(character = "OBSERVATION COUNT"),
               showProgress = FALSE)
  setnames(obs, c("common_name", "obs_count", "sampling_event_id", "group_id", "country_code"))
  obs <- obs[common_name == species_common & country_code == "IN"]
  obs[, key := fifelse(is.na(group_id), sampling_event_id, group_id)]
  obs[, count_num := suppressWarnings(as.integer(obs_count))]
  obs <- obs[, .(detected = 1L,
                 count = if (all(is.na(count_num))) NA_integer_ else max(count_num, na.rm = TRUE)),
             by = key]

  zf <- merge(smp, obs, by = "key", all.x = TRUE)
  zf[is.na(detected), detected := 0L]
  zf[detected == 0L, count := NA_integer_]
  zf[, c("observation_date", "time_start", "sampling_event_id", "country_code",
         "all_reported", "duration_min", "group_id", "observer_id",
         "obs_checklists", "month") := NULL]

  # Project to metres (India NSF LCC) for gridding and blocks
  xy <- sf_project("EPSG:4326", "EPSG:7755", cbind(zf$longitude, zf$latitude))
  zf[, `:=`(x = xy[, 1], y = xy[, 2])]
  saveRDS(zf, cache1)
}

n_det  <- sum(zf$detected == 1L)
n_xonly <- sum(zf$detected == 1L & is.na(zf$count))
log_msg("Checklists: ", format(nrow(zf), big.mark = ","),
        " | detections: ", format(n_det, big.mark = ","),
        " (", sprintf("%.2f%%", 100 * n_det / nrow(zf)), ")",
        " | detections with 'X' (no count): ", n_xonly)
if (quick) zf <- zf[sample.int(.N, min(.N, 500000))]

# ============================================================
# 3) SPATIOTEMPORAL THINNING (detection-agnostic) + SPATIAL BLOCKS
# ============================================================
# A winter spans two calendar years, so use winter year and days since 1 Nov
# (calendar year / day-of-year would jump and wrap at 1 January).
zf[, winter_year := fifelse(month(date) >= 11, year(date), year(date) - 1L)]
zf[, season_day := as.integer(date - as.Date(paste0(winter_year, "-11-01")))]
zf[, cell := floor(x / (cell_km * 1000)) * 1e6 + floor(y / (cell_km * 1000))]
zf <- zf[sample.int(.N)]
ss <- unique(zf, by = c("cell", "year", "iso_week"))
log_msg("After thinning (1 checklist / ", cell_km, " km cell / week / year): ",
        format(nrow(ss), big.mark = ","), " checklists, ",
        sum(ss$detected), " detections")

ss[, block := paste(floor(x / (block_km * 1000)), floor(y / (block_km * 1000)))]
blocks <- unique(ss$block)
test_blocks <- sample(blocks, ceiling(test_frac * length(blocks)))
ss[, split := fifelse(block %in% test_blocks, "test", "train")]
log_msg("Blocks: ", length(blocks), " | test blocks: ", length(test_blocks),
        " | train rows: ", sum(ss$split == "train"), " | test rows: ", sum(ss$split == "test"))

enc_vars <- c("winter_year", "season_day", "start_hour", "effort_hours",
              "effort_distance_km", "n_observers", "obs_exp", "latitude", "longitude")
setnames(ss, "distance_km", "effort_distance_km")

# ============================================================
# 4) MODEL PIPELINE: encounter rate -> calibration -> count given presence
# ============================================================
fit_pipeline <- function(d, trees = n_trees) {
  d <- as.data.frame(d)
  y <- d$detected
  er_dat <- d[, enc_vars]
  er_dat$y <- factor(y, levels = c(0, 1))
  # No case weights: up-weighting detections would leave them almost never
  # out-of-bag, so their OOB predictions (needed for calibration) go missing.
  er <- ranger(y ~ ., data = er_dat, probability = TRUE, num.trees = trees,
               num.threads = nth, seed = 1)
  # Out-of-bag predictions (NOT in-sample) -> unbiased calibration
  oob <- er$predictions[, "1"]
  ok <- is.finite(oob)
  cal <- gam(y ~ s(p, k = 6), family = binomial(),
             data = data.frame(y = y[ok], p = oob[ok]))
  er_cal <- rep(NA_real_, length(oob))
  er_cal[ok] <- as.numeric(predict(cal, data.frame(p = oob[ok]), type = "response"))

  # Count model: detections with numeric counts only
  idx <- which(y == 1 & !is.na(d$count) & d$count > 0 & ok)
  cd <- d[idx, enc_vars]
  cd$er_cal <- er_cal[idx]
  upper <- quantile(d$count[idx], winsor_prob)
  cd$count_w <- pmin(d$count[idx], upper)
  cnt <- ranger(count_w ~ ., data = cd, num.trees = trees, num.threads = nth, seed = 1)
  list(er = er, cal = cal, cnt = cnt, upper = upper)
}

predict_pipeline <- function(m, newdata) {
  newdata <- as.data.frame(newdata)
  raw <- predict(m$er, data = newdata[, enc_vars], num.threads = nth)$predictions[, "1"]
  er_c <- as.numeric(predict(m$cal, data.frame(p = raw), type = "response"))
  nd <- newdata[, enc_vars]
  nd$er_cal <- er_c
  cnt <- predict(m$cnt, data = nd, num.threads = nth)$predictions
  data.frame(er = er_c, cnt = cnt, ra = er_c * cnt)
}

# ============================================================
# 5) VALIDATION: (a) spatially held-out blocks (harsh: extrapolation to
#    unsampled regions), (b) random hold-out of thinned checklists
#    (interpolation, optimistic because neighbours are in training)
# ============================================================
validate <- function(train, test, label) {
log_msg("Fitting validation model: ", label)
m_val <- fit_pipeline(train)
pt <- cbind(test[, .(detected, count)], predict_pipeline(m_val, test))

auc   <- as.numeric(auc(roc(pt$detected, pt$er, quiet = TRUE)))
brier <- mean((pt$er - pt$detected)^2)
prev  <- mean(pt$detected)
bss   <- 1 - brier / (prev * (1 - prev))
cal_tab <- pt[, .(mean_pred = mean(er), obs_rate = mean(detected), n = .N),
              by = .(bin = cut(er, c(0, .01, .02, .05, .1, .2, .4, 1), include.lowest = TRUE))][order(bin)]
det_t <- pt[detected == 1 & !is.na(count)]
sp_count <- suppressWarnings(cor(det_t$count, det_t$cnt, method = "spearman"))
all_t <- pt[!(detected == 1 & is.na(count))]           # drop 'X'-only (unknown count)
all_t[, obs_count := fifelse(detected == 1, as.numeric(count), 0)]
sp_ra  <- suppressWarnings(cor(all_t$obs_count, all_t$ra, method = "spearman"))
pe_ra  <- suppressWarnings(cor(log1p(all_t$obs_count), log1p(all_t$ra)))

c(sprintf("=== %s ===", label),
  sprintf("Test checklists: %d | detections: %d (prevalence %.3f%%) | non-finite OOB: %d",
          nrow(pt), sum(pt$detected), 100 * prev, sum(!is.finite(m_val$er$predictions[, "1"]))),
  sprintf("Encounter rate: AUC = %.3f | Brier = %.5f | Brier skill vs prevalence = %.3f", auc, brier, bss),
  sprintf("Count given detection: Spearman = %.3f (n = %d)", sp_count, nrow(det_t)),
  sprintf("Relative abundance vs observed count (all test checklists): Spearman = %.3f | Pearson(log1p) = %.3f",
          sp_ra, pe_ra),
  "Calibration table:",
  capture.output(print(as.data.frame(cal_tab), row.names = FALSE)), "")
}

ss[, split_rand := fifelse(runif(.N) < test_frac, "test", "train")]
val_txt <- c(
  validate(ss[split == "train"], ss[split == "test"], "SPATIAL BLOCK HOLD-OUT"),
  validate(ss[split_rand == "train"], ss[split_rand == "test"], "RANDOM HOLD-OUT"))
writeLines(val_txt, file.path(out_dir, "validation.txt"))
cat(paste(val_txt, collapse = "\n"), "\n")

# ============================================================
# 6) FINAL MODEL ON ALL DATA
# ============================================================
log_msg("Fitting final model on all thinned data ...")
m_all <- fit_pipeline(ss)

# ============================================================
# 7) PREDICTION GRID: standardised checklist, supported cells only
# ============================================================
log_msg("Building India grid ...")
india <- ne_countries(country = "India", scale = 50, returnclass = "sf") |> st_transform(7755)
cent  <- st_make_grid(india, cellsize = cell_km * 1000, what = "centers")
cent  <- cent[lengths(st_intersects(cent, st_union(india))) > 0]
cxy   <- st_coordinates(cent)

# Support mask: nearest checklist location within support_km
loc <- unique(zf[, .(x = round(x / 1000) * 1000, y = round(y / 1000) * 1000)])
loc_sf <- st_as_sf(loc, coords = c("x", "y"), crs = 7755)
nn <- st_nearest_feature(cent, loc_sf)
dist_km <- as.numeric(st_distance(cent, loc_sf[nn, ], by_element = TRUE)) / 1000
sup <- dist_km <= support_km
grid <- data.table(x = cxy[sup, 1], y = cxy[sup, 2])
ll <- sf_project("EPSG:7755", "EPSG:4326", cbind(grid$x, grid$y))
grid[, `:=`(longitude = ll[, 1], latitude = ll[, 2])]
log_msg("Grid cells in India: ", length(cent), " | with data within ", support_km,
        " km: ", nrow(grid))

std_exp <- median(ss$obs_exp)
make_newdata <- function(g, d) {
  wy <- if (month(d) >= 11) year(d) else year(d) - 1L
  data.frame(winter_year = wy, season_day = as.integer(d - as.Date(paste0(wy, "-11-01"))),
             start_hour = std_start_hour,
             effort_hours = std_hours, effort_distance_km = std_km,
             n_observers = std_observers, obs_exp = std_exp,
             latitude = g$latitude, longitude = g$longitude)
}

log_msg("Predicting relative abundance for ", length(predict_dates), " weeks ...")
ra_list <- lapply(predict_dates, function(d) {
  p <- predict_pipeline(m_all, make_newdata(grid, d))
  data.table(date = d, x = grid$x, y = grid$y, longitude = grid$longitude,
             latitude = grid$latitude, p)
})
ra <- rbindlist(ra_list)
curve <- ra[, .(mean_ra = mean(ra), mean_er = mean(er), mean_count = mean(cnt)), by = date]
peak_date <- curve$date[which.max(curve$mean_ra)]
fwrite(curve, file.path(out_dir, "seasonal_curve.csv"))
log_msg("Peak week in the model: ", format(peak_date))

pk <- ra[date == peak_date]
fwrite(pk[order(-ra)][1:50], file.path(out_dir, "top50_cells_peak_week.csv"))
fwrite(pk, file.path(out_dir, "relative_abundance_peak_week.csv"))

# Map
india_b <- st_cast(st_union(india), "MULTILINESTRING")
p_map <- ggplot() +
  geom_tile(data = pk, aes(x, y, fill = ra), width = cell_km * 1000, height = cell_km * 1000) +
  geom_sf(data = india_b, colour = "grey30", linewidth = 0.2) +
  scale_fill_viridis_c(trans = "sqrt", name = "Expected\ncount per\nstandard\nchecklist") +
  coord_sf(crs = 7755) + theme_void() +
  labs(title = sprintf("Bar-headed Goose relative abundance, week of %s", format(peak_date)),
       subtitle = "Standard checklist: 1 h, 1 km, 1 observer, 07:00. White = no checklists within 25 km.")
ggsave(file.path(out_dir, "map_peak_week.png"), p_map, width = 7, height = 7.5, dpi = 150, bg = "white")

# ============================================================
# 8) UNCERTAINTY: spatial-block subsampling (refits both models)
# ============================================================
log_msg("Uncertainty: ", n_boot, " block-subsampling replicates ...")
set.seed(99)
bc <- grid[sample.int(nrow(grid), min(n_boot_cells, nrow(grid)))]
curve_fun <- function(m) sapply(predict_dates, function(d) mean(predict_pipeline(m, make_newdata(bc, d))$ra))
pt_curve <- curve_fun(m_all)
tr_blocks <- unique(ss$block)
boot_mat <- t(sapply(seq_len(n_boot), function(b) {
  keep <- sample(tr_blocks, floor(boot_frac * length(tr_blocks)))
  mb <- fit_pipeline(ss[block %in% keep], trees = boot_trees)
  log_msg("  replicate ", b, "/", n_boot)
  curve_fun(mb)
}))
band <- data.table(date = predict_dates, estimate = pt_curve,
                   lo = apply(boot_mat, 2, quantile, 0.025),
                   hi = apply(boot_mat, 2, quantile, 0.975))
fwrite(band, file.path(out_dir, "seasonal_curve_with_interval.csv"))
p_curve <- ggplot(band, aes(date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) + geom_line(aes(y = estimate)) +
  labs(x = NULL, y = "Mean relative abundance across supported cells",
       title = "Seasonal pattern in relative abundance (India)",
       subtitle = sprintf("Band: 2.5-97.5%% of %d block-subsampling refits", n_boot)) +
  theme_minimal()
ggsave(file.path(out_dir, "seasonal_curve.png"), p_curve, width = 7, height = 4, dpi = 150, bg = "white")

log_msg("Done. Outputs in ", normalizePath(out_dir))

