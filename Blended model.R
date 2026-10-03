#############################################################################################
#
#                          TROPICAL CYCLONES AND HEALTH
#                     Wind Field Reconstruction (Version 4.1)
#
#############################################################################################
#
# Developed by Teng Wang  @ City University of Hong Kong
#              Hanxu Shi  @ Peking University
# Contact: wang.teng19@alumni.imperial.ac.uk
# Version: 4.1
# Last updated: 2025-04-05
#
#############################################################################################


# ================================================================
#  PART 0: PACKAGES
# ================================================================

library(dplyr)
library(lubridate)
library(ggplot2)
library(ncdf4)
library(maps)
library(scales)
library(cowplot)

# Force English locale for date formatting (month abbreviations)
invisible(tryCatch(Sys.setlocale("LC_TIME", "C"),
                   error = function(e)
                     tryCatch(Sys.setlocale("LC_TIME", "English"),
                              error = function(e2) NULL)))

# ================================================================
#  PART 1: PATHS & PARAMETERS
# ================================================================

setwd("C:/Project/Tropical cyclone")

DIR_IBTRACS <- "Src/TC 2016-2024"  # 2016-2024
DIR_ERA5    <- "Src/ERA5/Wind speed"
DIR_RESULT  <- "Result 2016-2024 V5/WindField"
DIR_FIG     <- "Result 2016-2024 V5/Figures/WindField"
DIR_NC      <- "Result 2016-2024 V5/NetCDF"
DIR_VALID   <- "Result 2016-2024 V5/Validation"
DIR_3D      <- "Result 2016-2024 V5/3D"

for (d in c(DIR_RESULT, DIR_FIG, DIR_NC, DIR_VALID, DIR_3D))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

R_EARTH     <- 6371000
P_ENV       <- 101325
KT2MS       <- 0.514444
NM2KM       <- 1.852
GRID_RES    <- 0.25
ALPHA_ASYM  <- 0.55
THETA0_NH   <- pi / 2
N_BLEND     <- 4
INFLOW_ANG  <- 20 * pi / 180
RF_DEFAULT  <- 0.75
TIME_STEP   <- 24
SAVE_ALL_STEPS <- FALSE

# --- Holland Radius Method ---
# Choose one of: "FIXED", "BLEND", "WIND_SPEED"
#
#   FIXED       : fixed MAX_RADIUS_KM for all TCs
#   BLEND       : Holland contributes >= HOLLAND_THRESHOLD fraction
#   WIND_SPEED  : Holland surface wind drops to WIND_SPEED_THRESHOLD m/s
#                 (capped at WIND_SPEED_MAX_RADIUS km)
#
# Literature references for outer radius definitions:
#   Holland (1980)          : parametric profile, outer closed isobar
#   Chavas & Emanuel (2010) : vanishing wind (~0 m/s) from QuikSCAT
#   Chavas et al. (2015)    : physically-based outer wind structure
#   Operational (JTWC/NHC)  : R34 (34kt = 17.5 m/s), R50, R64
#   Recommended thresholds  : 1-5 m/s for total TC extent
#                             17.5 m/s for damaging wind extent
RADIUS_METHOD          <- "WIND_SPEED"
MAX_RADIUS_KM          <- 1000         # for FIXED
HOLLAND_THRESHOLD      <- 0.05         # for BLEND  (Holland >= 5%)
WIND_SPEED_THRESHOLD   <- 5.0          # for WIND_SPEED (m/s)
WIND_SPEED_MAX_RADIUS  <- 1500         # for WIND_SPEED (km cap)

# --- Domain: expanded to cover all of mainland China ---
DOMAIN_LON <- c(70, 150)
DOMAIN_LAT <- c(5, 55)

# --- SSHS classification (full names) ---
SSHS_BREAKS <- c(0, 17.5, 33.0, 43.0, 50.0, 58.0, 70.0, Inf)
SSHS_LABELS <- c("TD", "TS", "Category 1", "Category 2",
                 "Category 3", "Category 4", "Category 5")
SSHS_COLORS <- c("TD"         = "#6BAED6",
                 "TS"         = "#3182BD",
                 "Category 1" = "#FDAE6B",
                 "Category 2" = "#FD8D3C",
                 "Category 3" = "#E6550D",
                 "Category 4" = "#A63603",
                 "Category 5" = "#67000D")

WARNING_LINE <- data.frame(
  lon = c(105, 113, 119, 119, 127, 127),
  lat = c(0,   4.5, 11,  18,  22,  34)
)

# Threshold for Holland "TC-affected" mask in saved NetCDF
HOLLAND_MASK_THRESHOLD <- 0.5   # m/s


# ================================================================
#  PART 2: ESSENTIAL UTILITY FUNCTIONS
# ================================================================

haversine_dist <- function(lon1, lat1, lon2, lat2) {
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  la1  <- lat1 * pi / 180;  la2 <- lat2 * pi / 180
  a <- sin(dlat / 2)^2 + cos(la1) * cos(la2) * sin(dlon / 2)^2
  2 * R_EARTH * asin(pmin(sqrt(a), 1))
}

calc_azimuth <- function(clon, clat, plon, plat) {
  dx <- (plon - clon) * cos(clat * pi / 180)
  dy <- plat - clat
  atan2(dx, dy) %% (2 * pi)
}

gc_bearing <- function(lon1, lat1, lon2, lat2) {
  dl  <- (lon2 - lon1) * pi / 180
  la1 <- lat1 * pi / 180;  la2 <- lat2 * pi / 180
  atan2(sin(dl) * cos(la2),
        cos(la1) * sin(la2) - sin(la1) * cos(la2) * cos(dl)) %% (2 * pi)
}

est_rmax <- function(vmax_ms, lat_deg) {
  46.4 * exp(-0.0155 * vmax_ms + 0.0169 * abs(lat_deg))
}

est_holland_B <- function(vmax_ms, rmax_km, lat_deg) {
  B <- 1.0036 + 0.0173 * vmax_ms -
    0.0313 * log(pmax(rmax_km, 1)) + 0.0087 * abs(lat_deg)
  pmax(1.0, pmin(B, 2.5))
}

holland_Vg <- function(r_m, vmax_grad, rmax_m, B) {
  r_safe <- pmax(r_m, 500)
  ratio  <- (rmax_m / r_safe)^B
  vmax_grad * sqrt(ratio * exp(1 - ratio))
}

blend_xi <- function(r_m, rmax_m, n = N_BLEND) {
  cc  <- r_m / (n * rmax_m)
  cc4 <- cc^4
  cc4 / (1 + cc4)
}

# --- Radius Method 2: BLEND threshold ---
calc_blend_cutoff <- function(rmax_m, n = N_BLEND,
                              threshold = HOLLAND_THRESHOLD) {
  c_val <- ((1 - threshold) / threshold)^(1 / 4)
  min(c_val * n * rmax_m, MAX_RADIUS_KM * 1000)
}

# --- Radius Method 3: WIND_SPEED threshold ---
# Finds r > rmax where Holland surface wind = threshold_ms
calc_wind_cutoff <- function(vmax_grad, rmax_m, B,
                             rf = RF_DEFAULT,
                             threshold_ms = WIND_SPEED_THRESHOLD,
                             max_r_m = WIND_SPEED_MAX_RADIUS * 1000) {
  f <- function(r) {
    rf * holland_Vg(r, vmax_grad, rmax_m, B) - threshold_ms
  }
  # If wind at max_r still exceeds threshold, return max_r
  if (f(max_r_m) >= 0) return(max_r_m)
  # If wind just outside RMW is already below threshold, return small radius
  if (f(rmax_m * 1.5) <= 0) return(rmax_m * 3)
  result <- tryCatch(
    uniroot(f, interval = c(rmax_m * 1.5, max_r_m), tol = 100)$root,
    error = function(e) max_r_m
  )
  return(result)
}

# Unified radius calculation
calc_cutoff_radius <- function(method, vmax_grad, rmax_m, B, has_era5) {
  if (method == "FIXED") {
    r <- MAX_RADIUS_KM * 1000
  } else if (method == "BLEND") {
    r <- calc_blend_cutoff(rmax_m)
    if (!has_era5) r <- max(r, 500000)
  } else if (method == "WIND_SPEED") {
    r <- calc_wind_cutoff(vmax_grad, rmax_m, B)
    if (!has_era5) r <- max(r, 500000)
  } else {
    r <- MAX_RADIUS_KM * 1000
  }
  return(r)
}

safe_spline <- function(x, y, xout) {
  ok <- !is.na(y)
  if (sum(ok) >= 3)
    pmax(spline(x[ok], y[ok], xout = xout, method = "natural")$y, 0)
  else
    rep(mean(y, na.rm = TRUE), length(xout))
}

pal_wind <- colorRampPalette(
  c("white", "#E0F0FF", "#87CEEB", "#4682B4", "#1E5BA8", "#0C2C84", "#081D58")
)

mat2df <- function(mat, glon, glat, varname = "value") {
  df <- expand.grid(lon = glon, lat = glat)
  df[[varname]] <- as.vector(mat)
  df
}

save_both <- function(p, path_base, w = 7, h = 5.5) {
  tryCatch({
    ggsave(paste0(path_base, ".pdf"), p, width = w, height = h, dpi = 300,
           device = cairo_pdf)
  }, error = function(e) {
    tryCatch(ggsave(paste0(path_base, ".pdf"), p, width = w, height = h,
                    dpi = 300),
             error = function(e2) cat("      PDF fail:", e2$message, "\n"))
  })
  tryCatch({
    ggsave(paste0(path_base, ".png"), p, width = w, height = h, dpi = 300,
           bg = "white")
  }, error = function(e) cat("      PNG fail:", e$message, "\n"))
  while (dev.cur() > 1) tryCatch(dev.off(), error = function(e) break)
}

save_data_csv <- function(df, path_base) {
  tryCatch(
    write.csv(df, paste0(path_base, ".csv"), row.names = FALSE),
    error = function(e) cat("      CSV save fail:", e$message, "\n")
  )
}

# --- Publication theme: NO gridlines on any plot ---
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(colour = "black", family = ""),
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 0.5, hjust = 0,
                                      margin = margin(b = 6), colour = "grey30"),
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      legend.title     = element_text(size = base_size - 0.5, face = "bold"),
      legend.text      = element_text(size = base_size - 1.5),
      panel.border     = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.margin      = margin(5, 5, 5, 5, "mm"),
      strip.text       = element_text(size = base_size, face = "bold")
    )
}

build_arrow_df <- function(u_mat, v_mat, grid_lon, grid_lat,
                           nlon, nlat, arrow_skip = 4) {
  idx_i <- seq(1, nlon, by = arrow_skip)
  idx_j <- seq(1, nlat, by = arrow_skip)
  arr_grid <- expand.grid(i = idx_i, j = idx_j)
  arr_df <- data.frame(
    lon = grid_lon[arr_grid$i],
    lat = grid_lat[arr_grid$j],
    u   = mapply(function(ii, jj) u_mat[ii, jj], arr_grid$i, arr_grid$j),
    v   = mapply(function(ii, jj) v_mat[ii, jj], arr_grid$i, arr_grid$j)
  )
  arr_df$spd <- sqrt(arr_df$u^2 + arr_df$v^2)
  arr_df <- arr_df[arr_df$spd > 0.01, ]
  if (nrow(arr_df) > 0) {
    max_spd <- max(arr_df$spd, na.rm = TRUE)
    arr_df$u_norm <- arr_df$u / arr_df$spd
    arr_df$v_norm <- arr_df$v / arr_df$spd
    base_len <- 0.3
    arr_df$total_len <- base_len +
      (arr_df$spd / max(max_spd, 0.1)) * 0.25
    arr_df$end_lon <- arr_df$lon + arr_df$u_norm * arr_df$total_len
    arr_df$end_lat <- arr_df$lat + arr_df$v_norm * arr_df$total_len
    head_frac <- 0.35
    arr_df$wind_angle <- atan2(arr_df$v_norm, arr_df$u_norm)
    arr_df$hl <- arr_df$total_len * head_frac
    arr_df$L_lon <- arr_df$end_lon +
      cos(arr_df$wind_angle + pi - pi / 5) * arr_df$hl
    arr_df$L_lat <- arr_df$end_lat +
      sin(arr_df$wind_angle + pi - pi / 5) * arr_df$hl
    arr_df$R_lon <- arr_df$end_lon +
      cos(arr_df$wind_angle + pi + pi / 5) * arr_df$hl
    arr_df$R_lat <- arr_df$end_lat +
      sin(arr_df$wind_angle + pi + pi / 5) * arr_df$hl
  }
  arr_df
}

add_arrows_to_plot <- function(p, arr_df) {
  if (nrow(arr_df) > 0) {
    p <- p +
      geom_segment(data = arr_df,
                   aes(x = lon, y = lat, xend = end_lon, yend = end_lat),
                   colour = "black", linewidth = 0.2, alpha = 0.7) +
      geom_segment(data = arr_df,
                   aes(x = end_lon, y = end_lat, xend = L_lon, yend = L_lat),
                   colour = "black", linewidth = 0.18, alpha = 0.7) +
      geom_segment(data = arr_df,
                   aes(x = end_lon, y = end_lat, xend = R_lon, yend = R_lat),
                   colour = "black", linewidth = 0.18, alpha = 0.7)
  }
  p
}

make_wf_base <- function(speed_mat, grid_lon, grid_lat,
                         tc_h, tc_lon, tc_lat,
                         world_map, WARNING_LINE,
                         DOMAIN_LON, DOMAIN_LAT,
                         title_str, subtitle_str,
                         xi_mat = NULL, show_blend_contour = FALSE) {
  wdf <- mat2df(speed_mat, grid_lon, grid_lat, "ws")
  wdf <- wdf[wdf$ws > 0.5, ]
  p <- ggplot() +
    geom_raster(data = wdf, aes(lon, lat, fill = ws)) +
    scale_fill_gradientn(
      colours = pal_wind(200), limits = c(0, 70), oob = squish,
      breaks = seq(0, 70, 10),
      name = expression(bold("Wind speed (m s"^{-1}*")"))
    ) +
    geom_polygon(data = world_map, aes(long, lat, group = group),
                 fill = NA, colour = "grey40", linewidth = 0.25) +
    geom_path(data = WARNING_LINE, aes(lon, lat),
              colour = "red", linewidth = 0.5, linetype = "dashed")
  if (show_blend_contour && !is.null(xi_mat)) {
    xi_df <- mat2df(xi_mat, grid_lon, grid_lat, "xi")
    p <- p + geom_contour(data = xi_df, aes(lon, lat, z = xi),
                          breaks = 0.5, colour = "magenta",
                          linewidth = 0.7, linetype = "solid")
  }
  p <- p +
    geom_path(data = tc_h, aes(LON, LAT),
              colour = "darkorange", linewidth = 0.5, alpha = 0.8) +
    geom_point(data = tc_h[tc_h$HOUR == 0, ],
               aes(LON, LAT), colour = "black", fill = "darkorange",
               shape = 21, size = 2.2, stroke = 0.7) +
    geom_text(data = tc_h[!is.na(tc_h$DLABEL), ],
              aes(LON, LAT, label = DLABEL),
              size = 2.2, fontface = "bold", vjust = -1.8, colour = "black") +
    geom_point(aes(x = tc_lon, y = tc_lat),
               shape = 4, colour = "red", size = 3, stroke = 1.2) +
    coord_fixed(xlim = DOMAIN_LON, ylim = DOMAIN_LAT,
                ratio = 1, expand = FALSE) +
    scale_x_continuous(
      breaks = seq(DOMAIN_LON[1], DOMAIN_LON[2], 10),
      labels = paste0(seq(DOMAIN_LON[1], DOMAIN_LON[2], 10), "\u00B0E")) +
    scale_y_continuous(
      breaks = seq(DOMAIN_LAT[1], DOMAIN_LAT[2], 10),
      labels = paste0(seq(DOMAIN_LAT[1], DOMAIN_LAT[2], 10), "\u00B0N")) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Longitude", y = "Latitude") +
    theme_pub()
  p
}

# --- Validation scatter helper ---
make_valid_scatter <- function(df, x_col, y_col, colour_val,
                               title_str, xlab_str = NULL,
                               ylab_str = NULL) {
  ok <- !is.na(df[[x_col]]) & !is.na(df[[y_col]])
  df <- df[ok, ]
  if (nrow(df) < 2) return(NULL)
  mx <- max(c(df[[x_col]], df[[y_col]]), na.rm = TRUE) * 1.08
  r_val    <- cor(df[[x_col]], df[[y_col]], use = "complete")
  rmse_val <- sqrt(mean((df[[x_col]] - df[[y_col]])^2, na.rm = TRUE))
  bias_val <- mean(df[[y_col]] - df[[x_col]], na.rm = TRUE)
  stat_label <- paste0("r = ", sprintf("%.2f", r_val),
                       "\nRMSE = ", sprintf("%.1f", rmse_val), " m/s",
                       "\nbias = ", sprintf("%.1f", bias_val), " m/s",
                       "\nn = ", nrow(df))
  if (is.null(xlab_str))
    xlab_str <- expression("IBTrACS V"[max]*" (m s"^{-1}*")")
  if (is.null(ylab_str))
    ylab_str <- expression("Estimated V"[max]*" (m s"^{-1}*")")
  ggplot(df, aes(.data[[x_col]], .data[[y_col]])) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_point(colour = colour_val, size = 2, alpha = 0.65) +
    geom_smooth(method = "lm", se = FALSE,
                linewidth = 0.6, colour = "black") +
    annotate("text", x = mx * 0.05, y = mx * 0.92, label = stat_label,
             hjust = 0, size = 3, colour = colour_val) +
    coord_fixed(xlim = c(0, mx), ylim = c(0, mx)) +
    labs(x = xlab_str, y = ylab_str, title = title_str) +
    theme_pub()
}

# --- Combined ERA5 vs Blend scatter ---
make_combined_scatter <- function(val_clean, title_str) {
  ok <- !is.na(val_clean$era5_vmax) & !is.na(val_clean$blend_vmax)
  vc <- val_clean[ok, ]
  if (nrow(vc) < 2) return(NULL)
  df_long <- data.frame(
    ibt_vmax = rep(vc$ibt_vmax, 2),
    est_vmax = c(vc$era5_vmax, vc$blend_vmax),
    Source   = rep(c("ERA5", "Holland+ERA5"), each = nrow(vc))
  )
  mx <- max(c(df_long$ibt_vmax, df_long$est_vmax), na.rm = TRUE) * 1.08
  stats_df <- df_long %>%
    group_by(Source) %>%
    summarise(
      r    = cor(ibt_vmax, est_vmax, use = "complete"),
      rmse = sqrt(mean((ibt_vmax - est_vmax)^2, na.rm = TRUE)),
      bias = mean(est_vmax - ibt_vmax, na.rm = TRUE),
      n    = n(), .groups = "drop")
  stats_df$label <- paste0(stats_df$Source,
                           ": r=", sprintf("%.2f", stats_df$r),
                           ", RMSE=", sprintf("%.1f", stats_df$rmse),
                           " m/s")
  ggplot(df_long, aes(ibt_vmax, est_vmax, colour = Source)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    scale_colour_manual(
      values = c("ERA5" = "#E6550D", "Holland+ERA5" = "#0C2C84")) +
    annotate("text", x = mx * 0.05, y = mx * 0.92,
             label = paste(stats_df$label, collapse = "\n"),
             hjust = 0, size = 2.8) +
    coord_fixed(xlim = c(0, mx), ylim = c(0, mx)) +
    labs(x = expression("IBTrACS V"[max]*" (m s"^{-1}*")"),
         y = expression("Estimated V"[max]*" (m s"^{-1}*")"),
         title = title_str) +
    theme_pub()
}

# --- ERA5 file matching by year ---
find_era5_file <- function(era5_files, tc_year) {
  for (f in era5_files) {
    bn <- tolower(basename(f))
    m1 <- regmatches(bn, regexpr("hourly_(\\d{4})\\.nc$", bn))
    if (length(m1) > 0) {
      yr <- as.integer(sub(".*hourly_(\\d{4})\\.nc$", "\\1", m1))
      if (tc_year == yr) return(f)
    }
    m2 <- regmatches(bn, regexpr("daily_(\\d{4})_(\\d{4})\\.nc$", bn))
    if (length(m2) > 0) {
      parts <- regmatches(m2, gregexpr("\\d{4}", m2))[[1]]
      yr1 <- as.integer(parts[1]);  yr2 <- as.integer(parts[2])
      if (tc_year >= yr1 && tc_year <= yr2) return(f)
    }
    yr_matches <- as.integer(regmatches(bn, gregexpr("\\d{4}", bn))[[1]])
    yr_matches <- yr_matches[yr_matches >= 1950 & yr_matches <= 2030]
    if (length(yr_matches) == 1 && yr_matches[1] == tc_year) return(f)
    if (length(yr_matches) == 2 &&
        tc_year >= yr_matches[1] && tc_year <= yr_matches[2])
      return(f)
  }
  return(NULL)
}


# ================================================================
#  PART 3: AUTO-DETECT FILES & BUILD GRID
# ================================================================

csv_files  <- list.files(DIR_IBTRACS, pattern = "\\.csv$", full.names = TRUE)
era5_files <- list.files(DIR_ERA5,    pattern = "\\.nc$",  full.names = TRUE)
n_tc       <- length(csv_files)

cat("\n", paste(rep("=", 65), collapse = ""), "\n")
cat("  TROPICAL CYCLONE WIND FIELD RECONSTRUCTION  (v4.1)\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n")
cat("TC files  :", n_tc, "\n")
cat("ERA5 files:", length(era5_files), "\n")
cat("ERA5 list :", paste(basename(era5_files), collapse = ", "), "\n")
cat("Radius method:", RADIUS_METHOD, "\n")
if (RADIUS_METHOD == "FIXED")
  cat("  Fixed radius:", MAX_RADIUS_KM, "km\n")
if (RADIUS_METHOD == "BLEND")
  cat("  Blend threshold:", HOLLAND_THRESHOLD, "\n")
if (RADIUS_METHOD == "WIND_SPEED")
  cat("  Wind threshold:", WIND_SPEED_THRESHOLD, "m/s | Max:",
      WIND_SPEED_MAX_RADIUS, "km\n")
if (n_tc == 0) stop("No TC CSV files found in ", DIR_IBTRACS)

grid_lon <- seq(DOMAIN_LON[1], DOMAIN_LON[2], by = GRID_RES)
grid_lat <- seq(DOMAIN_LAT[1], DOMAIN_LAT[2], by = GRID_RES)
nlon     <- length(grid_lon)
nlat     <- length(grid_lat)
grd      <- expand.grid(lon = grid_lon, lat = grid_lat)
cat("Grid      :", nlon, "x", nlat, "=", nlon * nlat, "points\n")
cat("Domain    : Lon", DOMAIN_LON[1], "-", DOMAIN_LON[2],
    "| Lat", DOMAIN_LAT[1], "-", DOMAIN_LAT[2], "\n\n")

world_map      <- map_data("world")
all_validation <- data.frame()
tc_summary     <- data.frame()
total_start    <- Sys.time()


# ################################################################
#  PART 4: MAIN LOOP
# ################################################################

for (fi in 1:n_tc) {
  
  tc_csv_path <- csv_files[fi]
  tc_basename <- tools::file_path_sans_ext(basename(tc_csv_path))
  
  # ---- STEP A: Read & clean IBTrACS CSV ----
  tc <- tryCatch(
    read.csv(tc_csv_path, stringsAsFactors = FALSE),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
  if (is.null(tc)) next
  
  names(tc) <- toupper(names(tc))
  tc$ISO_TIME <- as.POSIXct(tc$ISO_TIME,
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  tc <- tc %>%
    filter(!is.na(LON), !is.na(LAT), !is.na(ISO_TIME)) %>%
    arrange(ISO_TIME) %>%
    distinct(ISO_TIME, .keep_all = TRUE)
  tc$LON <- ifelse(tc$LON > 180, tc$LON - 360, tc$LON)
  
  if ("WIND" %in% names(tc)) {
    tc$VMAX_MS <- as.numeric(tc$WIND) * KT2MS
    tc$VMAX_MS[tc$VMAX_MS < 0 | tc$VMAX_MS > 120] <- NA
  } else { cat("  No WIND column, skipping.\n"); next }
  
  tc$PMIN_PA <- NA
  if ("PRES" %in% names(tc)) {
    tc$PMIN_PA <- as.numeric(tc$PRES) * 100
    tc$PMIN_PA[tc$PMIN_PA < 85000 | tc$PMIN_PA > 108000] <- NA
  }
  if (nrow(tc) < 3 || all(is.na(tc$VMAX_MS))) {
    cat("  Insufficient data, skipping.\n"); next
  }
  
  # ---- STEP B: TC identifier ----
  tc_year <- year(min(tc$ISO_TIME))
  tc_sid <- ""
  if ("SID" %in% names(tc) && !is.na(tc$SID[1]) && nchar(tc$SID[1]) > 0) {
    tc_sid <- as.character(tc$SID[1])
  } else {
    sid_match <- regmatches(tc_basename,
                            regexpr("\\d{7}[NS]\\d{5}", tc_basename))
    if (length(sid_match) > 0) tc_sid <- sid_match[1]
  }
  tc_name <- "UNKNOWN"
  if ("NAME" %in% names(tc) && !is.na(tc$NAME[1]) &&
      nchar(tc$NAME[1]) > 0) {
    tc_name <- gsub("[^A-Za-z0-9]", "", tc$NAME[1])
  } else {
    parts <- strsplit(tc_basename, "_")[[1]]
    if (length(parts) >= 3) tc_name <- parts[3]
  }
  date_start <- format(min(tc$ISO_TIME), "%Y%m%d")
  date_end   <- format(max(tc$ISO_TIME), "%Y%m%d")
  tc_id <- paste(tc_year, tc_sid, tc_name, date_start, date_end, sep = "_")
  tc_id <- gsub("__+", "_", tc_id)
  
  elapsed_min <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
  if (fi > 1) {
    avg_per_tc <- elapsed_min / (fi - 1)
    est_remain <- avg_per_tc * (n_tc - fi + 1)
    pct <- round(100 * (fi - 1) / n_tc)
    cat(sprintf("\n[%s%s] %d%% | TC %d/%d: %s | %.1f min | ~%.1f min left\n",
                strrep("#", round(40 * pct / 100)),
                strrep(".", 40 - round(40 * pct / 100)),
                pct, fi, n_tc, tc_id, elapsed_min, est_remain))
  } else {
    cat(sprintf("\n[%s] 0%% | TC %d/%d: %s\n",
                strrep(".", 40), fi, n_tc, tc_id))
  }
  cat(sprintf("  Records: %d | Vmax: %.1f\u2013%.1f m/s\n",
              nrow(tc), min(tc$VMAX_MS, na.rm = TRUE),
              max(tc$VMAX_MS, na.rm = TRUE)))
  
  # ---- STEP C: Hourly interpolation ----
  t0  <- min(tc$ISO_TIME);  t1 <- max(tc$ISO_TIME)
  ht  <- seq(t0, t1, by = "1 hour")
  tn  <- as.numeric(tc$ISO_TIME)
  htn <- as.numeric(ht)
  nh  <- length(ht)
  
  lon_h  <- spline(tn, tc$LON, xout = htn, method = "natural")$y
  lat_h  <- spline(tn, tc$LAT, xout = htn, method = "natural")$y
  vmax_h <- safe_spline(tn, tc$VMAX_MS, htn)
  
  if (sum(!is.na(tc$PMIN_PA)) >= 3) {
    ok_p   <- !is.na(tc$PMIN_PA)
    pmin_h <- spline(tn[ok_p], tc$PMIN_PA[ok_p],
                     xout = htn, method = "natural")$y
    pmin_h <- pmax(85000, pmin(pmin_h, P_ENV))
  } else {
    pmin_h <- P_ENV - (vmax_h / 30)^2 * 2000
    pmin_h <- pmax(85000, pmin(pmin_h, P_ENV))
  }
  
  rmax_h <- est_rmax(vmax_h, lat_h)
  B_h    <- est_holland_B(vmax_h, rmax_h, lat_h)
  
  Vt  <- rep(0, nh);  hdg <- rep(0, nh)
  for (i in 2:nh) {
    dt_s   <- as.numeric(difftime(ht[i], ht[i - 1], units = "secs"))
    Vt[i]  <- haversine_dist(lon_h[i - 1], lat_h[i - 1],
                             lon_h[i], lat_h[i]) / max(dt_s, 1)
    hdg[i] <- gc_bearing(lon_h[i - 1], lat_h[i - 1],
                         lon_h[i], lat_h[i])
  }
  Vt[1] <- Vt[2];  hdg[1] <- hdg[2]
  if (nh > 5) {
    Vt_sm <- stats::filter(Vt, rep(1 / 5, 5), sides = 2)
    Vt_sm[is.na(Vt_sm)] <- Vt[is.na(Vt_sm)]
    Vt <- as.numeric(Vt_sm)
  }
  
  tc_h <- data.frame(
    ISO_TIME = ht, LON = lon_h, LAT = lat_h,
    VMAX_MS = vmax_h, PMIN_PA = pmin_h, RMAX_KM = rmax_h,
    B_PARAM = B_h, VT_MS = Vt, HEADING = hdg
  )
  tc_h$SSHS <- cut(tc_h$VMAX_MS, breaks = SSHS_BREAKS,
                   labels = SSHS_LABELS, right = FALSE)
  tc_h$HOUR <- hour(tc_h$ISO_TIME)
  tc_h$DLABEL <- ifelse(tc_h$HOUR == 0,
                        format(tc_h$ISO_TIME, "%m/%d"), NA)
  
  cat(sprintf("  Hourly records: %d | Duration: %.1f days\n",
              nrow(tc_h), as.numeric(difftime(t1, t0, units = "days"))))
  
  # ---- Output directories ----
  tc_figdir <- file.path(DIR_FIG, tc_id)
  tc_ncdir  <- file.path(DIR_NC, tc_id)
  for (d in c(tc_figdir, tc_ncdir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  
  # ---- STEP D: Pre-load ERA5 ----
  era5_nc   <- NULL;  era5_meta <- NULL
  era5_ei   <- NULL;  era5_ej   <- NULL
  
  era5_path <- find_era5_file(era5_files, tc_year)
  if (!is.null(era5_path)) {
    era5_nc <- tryCatch(nc_open(era5_path), error = function(e) NULL)
    if (!is.null(era5_nc)) {
      era_lon_full <- ncvar_get(era5_nc, "longitude")
      era_lat_full <- ncvar_get(era5_nc, "latitude")
      era_t <- tryCatch(ncvar_get(era5_nc, "valid_time"),
                        error = function(e)
                          ncvar_get(era5_nc, "time"))
      t_units <- tryCatch(
        ncatt_get(era5_nc, "valid_time", "units")$value,
        error = function(e)
          tryCatch(ncatt_get(era5_nc, "time", "units")$value,
                   error = function(e2)
                     "seconds since 1970-01-01"))
      if (grepl("1900", t_units)) {
        era_dt <- as.POSIXct(era_t * 3600,
                             origin = "1900-01-01", tz = "UTC")
      } else {
        era_dt <- as.POSIXct(era_t,
                             origin = "1970-01-01", tz = "UTC")
      }
      if (max(era_lon_full) > 180)
        era_lon_full <- ifelse(era_lon_full > 180,
                               era_lon_full - 360, era_lon_full)
      li <- which(era_lon_full >= DOMAIN_LON[1] &
                    era_lon_full <= DOMAIN_LON[2])
      lj <- which(era_lat_full >= DOMAIN_LAT[1] &
                    era_lat_full <= DOMAIN_LAT[2])
      if (length(li) > 0 && length(lj) > 0) {
        era_lon_sub <- era_lon_full[li]
        era_lat_sub <- era_lat_full[lj]
        era5_ei <- sapply(grid_lon,
                          function(x) which.min(abs(era_lon_sub - x)))
        era5_ej <- sapply(grid_lat,
                          function(x) which.min(abs(era_lat_sub - x)))
        era5_meta <- list(
          lon = era_lon_sub, lat = era_lat_sub, dt = era_dt,
          li = li, lj = lj,
          li_start = min(li), lj_start = min(lj),
          li_count = length(li), lj_count = length(lj)
        )
        t_range <- range(tc_h$ISO_TIME)
        era5_in_range <- sum(
          era5_meta$dt >= t_range[1] - 3600 * 6 &
            era5_meta$dt <= t_range[2] + 3600 * 6)
        cat("  ERA5:", basename(era5_path),
            "| Times:", length(era_dt),
            "| In TC range:", era5_in_range, "\n")
      } else {
        nc_close(era5_nc);  era5_nc <- NULL
        cat("  ERA5: domain mismatch, disabled\n")
      }
    }
  }
  if (is.null(era5_nc))
    cat("  ERA5: not available for year", tc_year, "\n")
  
  # ---- STEP E: Compute wind fields ----
  t_sel   <- seq(1, nrow(tc_h), by = TIME_STEP)
  n_steps <- length(t_sel)
  val_records <- vector("list", n_steps)
  all_blends  <- vector("list", n_steps)
  
  cat(sprintf("  Computing %d wind field steps (method: %s)...\n",
              n_steps, RADIUS_METHOD))
  
  for (k in seq_along(t_sel)) {
    row     <- tc_h[t_sel[k], ]
    clon    <- row$LON;  clat <- row$LAT
    vmax    <- row$VMAX_MS
    rmax_km <- row$RMAX_KM
    Bp      <- row$B_PARAM
    Vt_cur  <- row$VT_MS
    hdg_cur <- row$HEADING
    rmax_m  <- rmax_km * 1000
    vmax_grad <- vmax / RF_DEFAULT
    
    r_m   <- haversine_dist(clon, clat, grd$lon, grd$lat)
    theta <- calc_azimuth(clon, clat, grd$lon, grd$lat)
    
    era5_available <- !is.null(era5_nc)
    cutoff_r <- calc_cutoff_radius(RADIUS_METHOD, vmax_grad,
                                   rmax_m, Bp, era5_available)
    mask <- (r_m <= cutoff_r) & (r_m > 0)
    
    V_sfc <- rep(0, nrow(grd))
    u_sfc <- rep(0, nrow(grd))
    v_sfc <- rep(0, nrow(grd))
    
    Vg <- holland_Vg(r_m[mask], vmax_grad, rmax_m, Bp)
    Vtot <- pmax(Vg + ALPHA_ASYM * Vt_cur *
                   cos(theta[mask] - hdg_cur - THETA0_NH), 0)
    Vs <- RF_DEFAULT * Vtot
    
    V_sfc[mask] <- Vs
    u_sfc[mask] <- Vs * (-cos(theta[mask] - INFLOW_ANG))
    v_sfc[mask] <- Vs * ( sin(theta[mask] - INFLOW_ANG))
    
    hol_speed <- matrix(V_sfc, nrow = nlon, ncol = nlat)
    hol_u     <- matrix(u_sfc, nrow = nlon, ncol = nlat)
    hol_v     <- matrix(v_sfc, nrow = nlon, ncol = nlat)
    hol_dist  <- matrix(r_m,   nrow = nlon, ncol = nlat)
    
    blend_speed    <- hol_speed
    blend_u        <- hol_u
    blend_v        <- hol_v
    xi_mat         <- matrix(0, nlon, nlat)
    era5_speed_mat <- matrix(0, nlon, nlat)
    era5_u_mat     <- matrix(0, nlon, nlat)
    era5_v_mat     <- matrix(0, nlon, nlat)
    era5_vmax_val  <- NA
    era5_loaded    <- FALSE
    
    if (!is.null(era5_nc) && !is.null(era5_meta)) {
      t_idx_era <- which.min(abs(difftime(era5_meta$dt,
                                          row$ISO_TIME, units = "hours")))
      dt_hr <- abs(as.numeric(difftime(era5_meta$dt[t_idx_era],
                                       row$ISO_TIME, units = "hours")))
      if (dt_hr <= 6) {
        era_data <- tryCatch({
          u10 <- ncvar_get(era5_nc, "u10",
                           start = c(era5_meta$li_start,
                                     era5_meta$lj_start, t_idx_era),
                           count = c(era5_meta$li_count,
                                     era5_meta$lj_count, 1))
          v10 <- ncvar_get(era5_nc, "v10",
                           start = c(era5_meta$li_start,
                                     era5_meta$lj_start, t_idx_era),
                           count = c(era5_meta$li_count,
                                     era5_meta$lj_count, 1))
          list(u10 = u10, v10 = v10)
        }, error = function(e) NULL)
        if (!is.null(era_data)) {
          era5_loaded    <- TRUE
          era_u_interp   <- era_data$u10[era5_ei, era5_ej]
          era_v_interp   <- era_data$v10[era5_ei, era5_ej]
          era5_u_mat     <- era_u_interp
          era5_v_mat     <- era_v_interp
          xi_mat         <- blend_xi(hol_dist, rmax_m, N_BLEND)
          blend_u        <- xi_mat * era_u_interp +
            (1 - xi_mat) * hol_u
          blend_v        <- xi_mat * era_v_interp +
            (1 - xi_mat) * hol_v
          blend_speed    <- sqrt(blend_u^2 + blend_v^2)
          era5_speed_mat <- sqrt(era_u_interp^2 + era_v_interp^2)
          grd_e  <- expand.grid(lon = era5_meta$lon,
                                lat = era5_meta$lat)
          d_e    <- haversine_dist(clon, clat,
                                   grd_e$lon, grd_e$lat) / 1000
          era_ws <- sqrt(as.vector(era_data$u10)^2 +
                           as.vector(era_data$v10)^2)
          near   <- era_ws[d_e < 200]
          if (length(near) > 0 && any(!is.na(near)))
            era5_vmax_val <- max(near, na.rm = TRUE)
        }
      }
    }
    
    blend <- list(
      speed = blend_speed, u = blend_u, v = blend_v,
      xi = xi_mat,
      holland_speed = hol_speed, holland_u = hol_u, holland_v = hol_v,
      era5_speed = era5_speed_mat, era5_u = era5_u_mat,
      era5_v = era5_v_mat,
      has_era5 = era5_loaded,
      tc_lon = clon, tc_lat = clat,
      vmax = vmax, rmax_m = rmax_m, B = Bp,
      heading = hdg_cur, Vt = Vt_cur,
      time = row$ISO_TIME,
      cutoff_r = cutoff_r
    )
    all_blends[[k]] <- blend
    
    # --- Validation metrics ---
    r_vec        <- as.vector(hol_dist)
    hol_ws_vec   <- as.vector(hol_speed)
    blend_ws_vec <- as.vector(blend_speed)
    era_ws_vec   <- as.vector(era5_speed_mat)
    
    near_hol   <- hol_ws_vec[r_vec < 300000]
    near_blend <- blend_ws_vec[r_vec < 300000]
    holland_vmax_local <- if (length(near_hol) > 0)
      max(near_hol, na.rm = TRUE) else max(hol_speed, na.rm = TRUE)
    blend_vmax_local <- if (length(near_blend) > 0)
      max(near_blend, na.rm = TRUE) else max(blend_speed, na.rm = TRUE)
    
    # Area-mean at outer radii (200-500 km): shows blending difference
    mask_outer <- r_vec >= 200000 & r_vec <= 500000
    hol_mean_outer   <- mean(hol_ws_vec[mask_outer], na.rm = TRUE)
    blend_mean_outer <- mean(blend_ws_vec[mask_outer], na.rm = TRUE)
    era5_mean_outer  <- if (era5_loaded)
      mean(era_ws_vec[mask_outer], na.rm = TRUE) else NA_real_
    
    val_records[[k]] <- data.frame(
      time             = row$ISO_TIME,
      ibt_vmax         = vmax,
      era5_vmax        = era5_vmax_val,
      holland_vmax     = holland_vmax_local,
      blend_vmax       = blend_vmax_local,
      holland_mean_200_500 = hol_mean_outer,
      blend_mean_200_500   = blend_mean_outer,
      era5_mean_200_500    = era5_mean_outer,
      holland_radius_km    = round(cutoff_r / 1000, 1)
    )
    
    # --- Save BLENDED NetCDF ---
    ts_str   <- format(row$ISO_TIME, "%Y%m%d_%H%M")
    nc_fname <- file.path(tc_ncdir,
                          paste0("WF_", tc_id, "_", ts_str, ".nc"))
    dim_lon <- ncdim_def("longitude", "degrees_east",  grid_lon)
    dim_lat <- ncdim_def("latitude",  "degrees_north", grid_lat)
    var_spd <- ncvar_def("wind_speed",      "m s-1",
                         list(dim_lon, dim_lat), missval = -999)
    var_u   <- ncvar_def("u10",             "m s-1",
                         list(dim_lon, dim_lat), missval = -999)
    var_v   <- ncvar_def("v10",             "m s-1",
                         list(dim_lon, dim_lat), missval = -999)
    var_xi  <- ncvar_def("blending_weight", "1",
                         list(dim_lon, dim_lat), missval = -999)
    ncf <- nc_create(nc_fname, list(var_spd, var_u, var_v, var_xi))
    ncvar_put(ncf, var_spd, blend_speed)
    ncvar_put(ncf, var_u,   blend_u)
    ncvar_put(ncf, var_v,   blend_v)
    ncvar_put(ncf, var_xi,  xi_mat)
    ncatt_put(ncf, 0, "title",    "TC wind field (Holland+ERA5)")
    ncatt_put(ncf, 0, "time",     as.character(row$ISO_TIME))
    ncatt_put(ncf, 0, "tc_lon",   clon)
    ncatt_put(ncf, 0, "tc_lat",   clat)
    ncatt_put(ncf, 0, "vmax_ms",  vmax)
    ncatt_put(ncf, 0, "tc_id",    tc_id)
    ncatt_put(ncf, 0, "holland_radius_km", round(cutoff_r / 1000, 1))
    ncatt_put(ncf, 0, "radius_method", RADIUS_METHOD)
    nc_close(ncf)
    
    # --- Save HOLLAND-ONLY NetCDF (for Step 2 county analysis) ---
    nc_hol_fname <- file.path(tc_ncdir,
                              paste0("WF_Holland_", tc_id, "_",
                                     ts_str, ".nc"))
    var_hspd <- ncvar_def("holland_wind_speed", "m s-1",
                          list(dim_lon, dim_lat), missval = -999)
    var_hu   <- ncvar_def("holland_u10",        "m s-1",
                          list(dim_lon, dim_lat), missval = -999)
    var_hv   <- ncvar_def("holland_v10",        "m s-1",
                          list(dim_lon, dim_lat), missval = -999)
    var_mask <- ncvar_def("tc_affected",        "1",
                          list(dim_lon, dim_lat), missval = -999)
    ncf_h <- nc_create(nc_hol_fname,
                       list(var_hspd, var_hu, var_hv, var_mask))
    ncvar_put(ncf_h, var_hspd, hol_speed)
    ncvar_put(ncf_h, var_hu,   hol_u)
    ncvar_put(ncf_h, var_hv,   hol_v)
    tc_affected_mask <- ifelse(hol_speed >= HOLLAND_MASK_THRESHOLD,
                               1, 0)
    ncvar_put(ncf_h, var_mask, tc_affected_mask)
    ncatt_put(ncf_h, 0, "title",
              "TC wind field (Holland model only)")
    ncatt_put(ncf_h, 0, "description",
              paste0("Holland parametric model wind field. ",
                     "tc_affected = 1 where holland_wind_speed >= ",
                     HOLLAND_MASK_THRESHOLD, " m/s. ",
                     "Use this to identify TC-influenced grid cells."))
    ncatt_put(ncf_h, 0, "time",     as.character(row$ISO_TIME))
    ncatt_put(ncf_h, 0, "tc_lon",   clon)
    ncatt_put(ncf_h, 0, "tc_lat",   clat)
    ncatt_put(ncf_h, 0, "vmax_ms",  vmax)
    ncatt_put(ncf_h, 0, "tc_id",    tc_id)
    ncatt_put(ncf_h, 0, "holland_radius_km",
              round(cutoff_r / 1000, 1))
    ncatt_put(ncf_h, 0, "radius_method", RADIUS_METHOD)
    ncatt_put(ncf_h, 0, "mask_threshold_ms", HOLLAND_MASK_THRESHOLD)
    nc_close(ncf_h)
    
    if (k %% 10 == 0 || k == n_steps)
      cat(sprintf("    Step %d/%d done (R=%.0f km)\n",
                  k, n_steps, cutoff_r / 1000))
  }
  
  if (!is.null(era5_nc)) { nc_close(era5_nc); era5_nc <- NULL }
  
  
  # ===========================================================
  #  STEP F: GENERATE ALL PLOTS
  # ===========================================================
  
  cat("  Generating plots...\n")
  tryCatch(graphics.off(), error = function(e) NULL)
  
  peak_k <- which.max(sapply(all_blends, function(b) b$vmax))
  pk     <- all_blends[[peak_k]]
  ts_pk  <- format(pk$time, "%Y%m%d_%H%M")
  
  max_vmax_tc <- max(tc_h$VMAX_MS, na.rm = TRUE)
  actual_cats <- levels(tc_h$SSHS)[levels(tc_h$SSHS) %in%
                                     unique(as.character(tc_h$SSHS))]
  sshs_present <- SSHS_LABELS[SSHS_LABELS %in% actual_cats]
  sshs_col_sub <- SSHS_COLORS[sshs_present]
  
  y_display_max <- max_vmax_tc * 1.15
  band_visible  <- SSHS_BREAKS[-length(SSHS_BREAKS)] < y_display_max
  band_labels   <- SSHS_LABELS[band_visible]
  band_ymin     <- SSHS_BREAKS[1:length(SSHS_LABELS)][band_visible]
  band_ymax     <- SSHS_BREAKS[2:(length(SSHS_LABELS) + 1)][band_visible]
  band_ymax[band_ymax > y_display_max] <- y_display_max
  
  
  # F.1  TC Track Map
  tryCatch({
    tc_h$SSHS <- factor(tc_h$SSHS, levels = sshs_present)
    p_track <- ggplot() +
      geom_polygon(data = world_map, aes(long, lat, group = group),
                   fill = "grey92", colour = "grey60", linewidth = 0.2) +
      geom_path(data = WARNING_LINE, aes(lon, lat),
                colour = "red", linewidth = 0.6, linetype = "dashed") +
      geom_path(data = tc_h, aes(LON, LAT),
                colour = "grey40", linewidth = 0.4) +
      geom_point(data = tc_h, aes(LON, LAT, colour = SSHS),
                 size = 1.2, shape = 16) +
      scale_colour_manual(values = sshs_col_sub,
                          name = "SSHS Category", drop = TRUE) +
      geom_point(data = tc_h[tc_h$HOUR == 0, ],
                 aes(LON, LAT), colour = "black", size = 3,
                 shape = 21, fill = "white", stroke = 0.9) +
      geom_text(data = tc_h[!is.na(tc_h$DLABEL), ],
                aes(LON, LAT, label = DLABEL),
                size = 2.8, fontface = "bold", vjust = -2.0,
                colour = "black") +
      coord_fixed(xlim = DOMAIN_LON, ylim = DOMAIN_LAT, ratio = 1) +
      scale_x_continuous(
        breaks = seq(DOMAIN_LON[1], DOMAIN_LON[2], 10),
        labels = paste0(seq(DOMAIN_LON[1], DOMAIN_LON[2], 10),
                        "\u00B0E")) +
      scale_y_continuous(
        breaks = seq(DOMAIN_LAT[1], DOMAIN_LAT[2], 10),
        labels = paste0(seq(DOMAIN_LAT[1], DOMAIN_LAT[2], 10),
                        "\u00B0N")) +
      labs(title = paste0("TC Track: ", tc_name, " (", tc_year, ")"),
           subtitle = paste(date_start, "\u2013", date_end,
                            "| Peak Vmax =",
                            round(max_vmax_tc, 1), "m/s"),
           x = "Longitude", y = "Latitude") +
      theme_pub()
    save_both(p_track,
              file.path(tc_figdir, paste0("Track_", tc_id)),
              w = 8, h = 5.5)
    cat("    Track map saved\n")
  }, error = function(e) cat("    Track error:", e$message, "\n"))
  
  
  # F.2  Wind Field Map (clean, blended)
  tryCatch({
    blend_note <- if (pk$has_era5)
      "| Magenta contour: 50% blend boundary (\u03be=0.5)"
    else "| ERA5 not available"
    p_wf_clean <- make_wf_base(
      pk$speed, grid_lon, grid_lat, tc_h, pk$tc_lon, pk$tc_lat,
      world_map, WARNING_LINE, DOMAIN_LON, DOMAIN_LAT,
      title_str = paste0("TC Wind Field (Blended): ",
                         tc_name, " (", tc_year, ")"),
      subtitle_str = paste(format(pk$time, "%Y-%m-%d %H:%M UTC"),
                           "| Vmax =", round(pk$vmax, 1), "m/s",
                           "| R =", round(pk$cutoff_r / 1000), "km",
                           blend_note),
      xi_mat = pk$xi, show_blend_contour = pk$has_era5
    )
    save_both(p_wf_clean,
              file.path(tc_figdir,
                        paste0("WF_Clean_", tc_id, "_", ts_pk)),
              w = 8, h = 5.5)
    cat("    Wind field (clean) saved\n")
  }, error = function(e) cat("    WF clean error:", e$message, "\n"))
  
  
  # F.3  Arrow maps
  tryCatch({
    p_hol_base <- make_wf_base(
      pk$holland_speed, grid_lon, grid_lat, tc_h,
      pk$tc_lon, pk$tc_lat,
      world_map, WARNING_LINE, DOMAIN_LON, DOMAIN_LAT,
      title_str = paste0("Holland Model Wind Field: ",
                         tc_name, " (", tc_year, ")"),
      subtitle_str = paste(format(pk$time, "%Y-%m-%d %H:%M UTC"),
                           "| Vmax =", round(pk$vmax, 1), "m/s")
    )
    arr_hol <- build_arrow_df(pk$holland_u, pk$holland_v,
                              grid_lon, grid_lat, nlon, nlat)
    p_hol_arrows <- add_arrows_to_plot(p_hol_base, arr_hol)
    save_both(p_hol_arrows,
              file.path(tc_figdir,
                        paste0("WF_Arrows_Holland_", tc_id, "_", ts_pk)),
              w = 8, h = 5.5)
    cat("    Holland arrow map saved\n")
  }, error = function(e) cat("    Holland arrows error:", e$message, "\n"))
  
  if (pk$has_era5) {
    tryCatch({
      p_era_base <- make_wf_base(
        pk$era5_speed, grid_lon, grid_lat, tc_h,
        pk$tc_lon, pk$tc_lat,
        world_map, WARNING_LINE, DOMAIN_LON, DOMAIN_LAT,
        title_str = paste0("ERA5 Reanalysis Wind Field: ",
                           tc_name, " (", tc_year, ")"),
        subtitle_str = paste(format(pk$time, "%Y-%m-%d %H:%M UTC"),
                             "| ERA5 10-m wind")
      )
      arr_era <- build_arrow_df(pk$era5_u, pk$era5_v,
                                grid_lon, grid_lat, nlon, nlat)
      p_era_arrows <- add_arrows_to_plot(p_era_base, arr_era)
      save_both(p_era_arrows,
                file.path(tc_figdir,
                          paste0("WF_Arrows_ERA5_", tc_id, "_", ts_pk)),
                w = 8, h = 5.5)
      cat("    ERA5 arrow map saved\n")
    }, error = function(e) cat("    ERA5 arrows error:", e$message, "\n"))
  }
  
  tryCatch({
    blend_note2 <- if (pk$has_era5)
      "| Magenta contour: 50% blend boundary"
    else "| Holland only (no ERA5)"
    p_bld_base <- make_wf_base(
      pk$speed, grid_lon, grid_lat, tc_h,
      pk$tc_lon, pk$tc_lat,
      world_map, WARNING_LINE, DOMAIN_LON, DOMAIN_LAT,
      title_str = paste0("Blended Wind Field (Holland+ERA5): ",
                         tc_name, " (", tc_year, ")"),
      subtitle_str = paste(format(pk$time, "%Y-%m-%d %H:%M UTC"),
                           "| Vmax =", round(pk$vmax, 1), "m/s",
                           blend_note2),
      xi_mat = pk$xi, show_blend_contour = pk$has_era5
    )
    arr_bld <- build_arrow_df(pk$u, pk$v,
                              grid_lon, grid_lat, nlon, nlat)
    p_bld_arrows <- add_arrows_to_plot(p_bld_base, arr_bld)
    save_both(p_bld_arrows,
              file.path(tc_figdir,
                        paste0("WF_Arrows_Blend_", tc_id, "_", ts_pk)),
              w = 8, h = 5.5)
    cat("    Blended arrow map saved\n")
  }, error = function(e) cat("    Blend arrows error:", e$message, "\n"))
  
  
  # F.4  3D Wind Field Surface
  tryCatch({
    buf3d <- 5
    sub_i <- which(grid_lon >= pk$tc_lon - buf3d &
                     grid_lon <= pk$tc_lon + buf3d)
    sub_j <- which(grid_lat >= pk$tc_lat - buf3d &
                     grid_lat <= pk$tc_lat + buf3d)
    if (length(sub_i) >= 5 && length(sub_j) >= 5) {
      x_km <- (grid_lon[sub_i] - pk$tc_lon) *
        111.32 * cos(pk$tc_lat * pi / 180)
      y_km <- (grid_lat[sub_j] - pk$tc_lat) * 111.32
      z_3d <- pk$holland_speed[sub_i, sub_j]
      z_range  <- range(z_3d, na.rm = TRUE)
      n_colors <- 100
      col_pal  <- pal_wind(n_colors)
      z_facet <- (z_3d[-1, -1] + z_3d[-1, -ncol(z_3d)] +
                    z_3d[-nrow(z_3d), -1] +
                    z_3d[-nrow(z_3d), -ncol(z_3d)]) / 4
      z_norm  <- (z_facet - z_range[1]) /
        max(z_range[2] - z_range[1], 0.1)
      z_norm  <- pmin(pmax(z_norm, 0), 1)
      col_idx <- floor(z_norm * (n_colors - 1)) + 1
      facet_colors <- col_pal[col_idx]
      fn_3d <- file.path(tc_figdir,
                         paste0("WF3D_", tc_id, "_", ts_pk, ".png"))
      png(fn_3d, width = 2400, height = 1800, res = 300)
      par(mar = c(1.5, 1.5, 3, 1.5), bg = "white")
      pmat <- persp(
        x = x_km, y = y_km, z = z_3d,
        xlab = "East\u2013West (km)",
        ylab = "North\u2013South (km)",
        zlab = "Wind speed (m/s)",
        main = paste0("TC ", tc_name, " (", tc_year, ")\n",
                      format(pk$time, "%Y-%m-%d %H:%M UTC"),
                      "  |  Vmax = ", round(pk$vmax, 1), " m/s"),
        col = facet_colors, border = NA,
        theta = 315, phi = 25, expand = 0.55,
        shade = 0.35, ltheta = 120,
        ticktype = "detailed",
        cex.main = 0.85, cex.lab = 0.7, cex.axis = 0.6,
        box = TRUE, nticks = 5,
        zlim = c(0, max(z_3d, na.rm = TRUE) * 1.05)
      )
      par(fig = c(0.88, 0.92, 0.15, 0.85), new = TRUE,
          mar = c(0, 0, 0, 0))
      legend_z <- seq(z_range[1], z_range[2], length.out = n_colors)
      image(1, legend_z, t(as.matrix(legend_z)),
            col = col_pal, axes = FALSE, xlab = "", ylab = "")
      axis(4, at = pretty(legend_z, n = 5),
           labels = round(pretty(legend_z, n = 5), 0),
           cex.axis = 0.6, las = 1, tck = -0.3)
      mtext(expression("m s"^{-1}), side = 4, line = 2, cex = 0.55)
      dev.off()
      cat("    3D surface saved\n")
    }
  }, error = function(e) {
    tryCatch(dev.off(), error = function(e2) NULL)
    cat("    3D plot error:", e$message, "\n")
  })
  
  
  # F.5  Radial Wind Profile with asymmetry decomposition
  tryCatch({
    r_km_vec   <- haversine_dist(pk$tc_lon, pk$tc_lat,
                                 grd$lon, grd$lat) / 1000
    theta_vec  <- calc_azimuth(pk$tc_lon, pk$tc_lat,
                               grd$lon, grd$lat)
    hol_ws_vec <- as.vector(pk$holland_speed)
    era_ws_vec <- as.vector(pk$era5_speed)
    bld_ws_vec <- as.vector(pk$speed)
    
    r_bins <- seq(0, 500, by = 5)
    r_mid  <- (r_bins[-length(r_bins)] + r_bins[-1]) / 2
    n_bins <- length(r_mid)
    r_cut  <- cut(r_km_vec, breaks = r_bins, labels = FALSE)
    
    hol_mean_raw <- tapply(hol_ws_vec, r_cut, mean, na.rm = TRUE)
    bld_mean_raw <- tapply(bld_ws_vec, r_cut, mean, na.rm = TRUE)
    hol_mean_full <- rep(NA_real_, n_bins)
    bld_mean_full <- rep(NA_real_, n_bins)
    hol_mean_full[as.integer(names(hol_mean_raw))] <- hol_mean_raw
    bld_mean_full[as.integer(names(bld_mean_raw))] <- bld_mean_raw
    
    peak_azi <- (pk$heading + THETA0_NH) %% (2 * pi)
    azi_diff <- (theta_vec - peak_azi + pi) %% (2 * pi) - pi
    right_mask <- abs(azi_diff) < pi / 3
    left_mask  <- abs(azi_diff) > 2 * pi / 3
    
    hol_right_raw <- tapply(hol_ws_vec[right_mask],
                            r_cut[right_mask], mean, na.rm = TRUE)
    hol_left_raw  <- tapply(hol_ws_vec[left_mask],
                            r_cut[left_mask], mean, na.rm = TRUE)
    hol_right_full <- rep(NA_real_, n_bins)
    hol_left_full  <- rep(NA_real_, n_bins)
    if (length(hol_right_raw) > 0)
      hol_right_full[as.integer(names(hol_right_raw))] <- hol_right_raw
    if (length(hol_left_raw) > 0)
      hol_left_full[as.integer(names(hol_left_raw))] <- hol_left_raw
    
    df_radial <- data.frame(
      r  = rep(r_mid, 4),
      ws = c(hol_mean_full, bld_mean_full,
             hol_right_full, hol_left_full),
      Source = rep(c("Holland (azimuthal mean)",
                     "Blended (Holland+ERA5)",
                     "Holland (right-of-track)",
                     "Holland (left-of-track)"),
                   each = n_bins)
    )
    if (pk$has_era5) {
      era_mean_raw  <- tapply(era_ws_vec, r_cut, mean, na.rm = TRUE)
      era_mean_full <- rep(NA_real_, n_bins)
      era_mean_full[as.integer(names(era_mean_raw))] <- era_mean_raw
      df_radial <- rbind(df_radial,
                         data.frame(r = r_mid, ws = era_mean_full,
                                    Source = "ERA5"))
    }
    df_radial <- df_radial[!is.na(df_radial$ws), ]
    
    r_th_m <- seq(1, 500, by = 1) * 1000
    vg_th  <- holland_Vg(r_th_m, pk$vmax / RF_DEFAULT, pk$rmax_m,
                         est_holland_B(pk$vmax, pk$rmax_m / 1000,
                                       pk$tc_lat))
    vs_th  <- RF_DEFAULT * vg_th
    df_th  <- data.frame(r = r_th_m / 1000, ws = vs_th)
    
    source_cols <- c("Holland (azimuthal mean)" = "#3182BD",
                     "ERA5"                      = "#E6550D",
                     "Blended (Holland+ERA5)"    = "#0C2C84",
                     "Holland (right-of-track)"  = "#D62728",
                     "Holland (left-of-track)"   = "#2CA02C")
    source_lty <- c("Holland (azimuthal mean)" = "solid",
                    "ERA5"                      = "solid",
                    "Blended (Holland+ERA5)"    = "solid",
                    "Holland (right-of-track)"  = "dashed",
                    "Holland (left-of-track)"   = "dashed")
    
    p_rad <- ggplot() +
      geom_line(data = df_th, aes(r, ws),
                colour = "grey50", linewidth = 0.6,
                linetype = "dotted") +
      geom_line(data = df_radial,
                aes(r, ws, colour = Source, linetype = Source),
                linewidth = 0.7) +
      scale_colour_manual(values = source_cols, name = "Source") +
      scale_linetype_manual(values = source_lty, name = "Source") +
      geom_vline(xintercept = pk$rmax_m / 1000,
                 linetype = "dotted", colour = "red",
                 linewidth = 0.5) +
      annotate("text", x = pk$rmax_m / 1000 + 12,
               y = pk$vmax * 0.85,
               label = paste0("RMW = ",
                              round(pk$rmax_m / 1000), " km"),
               colour = "red", size = 3.2, hjust = 0) +
      labs(x = "Distance from TC centre (km)",
           y = expression("Surface wind speed (m s"^{-1}*")"),
           title = paste0("Radial Wind Profile: ", tc_name,
                          " (", tc_year, ")"),
           subtitle = format(pk$time, "%Y-%m-%d %H:%M UTC")) +
      theme_pub() +
      theme(legend.position = c(0.72, 0.72),
            legend.background = element_rect(
              fill = alpha("white", 0.85), colour = NA))
    save_both(p_rad,
              file.path(tc_figdir, paste0("Radial_", tc_id)),
              w = 7.5, h = 5)
    cat("    Radial profile saved\n")
  }, error = function(e) cat("    Radial error:", e$message, "\n"))
  
  
  # F.6  Three-panel comparison
  tryCatch({
    buf_cmp <- 8
    sub_lon <- c(pk$tc_lon - buf_cmp, pk$tc_lon + buf_cmp)
    sub_lat <- c(pk$tc_lat - buf_cmp, pk$tc_lat + buf_cmp)
    make_panel <- function(mat, lab, show_legend = FALSE) {
      df <- mat2df(mat, grid_lon, grid_lat, "ws")
      p <- ggplot() +
        geom_raster(data = df[df$ws > 0.5, ],
                    aes(lon, lat, fill = ws)) +
        scale_fill_gradientn(
          colours = pal_wind(200), limits = c(0, 70),
          oob = squish,
          name = expression("Wind (m s"^{-1}*")"),
          breaks = seq(0, 70, 10)
        ) +
        geom_polygon(data = world_map,
                     aes(long, lat, group = group),
                     fill = NA, colour = "grey50",
                     linewidth = 0.2) +
        geom_point(aes(x = pk$tc_lon, y = pk$tc_lat),
                   shape = 4, colour = "red", size = 2.5,
                   stroke = 1) +
        coord_fixed(xlim = sub_lon, ylim = sub_lat,
                    expand = FALSE) +
        labs(title = lab, x = "Longitude", y = "Latitude") +
        theme_pub()
      if (!show_legend) p <- p + theme(legend.position = "none")
      p
    }
    p1 <- make_panel(pk$holland_speed, "a  Holland model")
    p2 <- make_panel(pk$era5_speed,    "b  ERA5 reanalysis")
    p3 <- make_panel(pk$speed,         "c  Blended (Holland+ERA5)")
    leg_df <- data.frame(ws = seq(0, 70, 1))
    pleg <- ggplot(leg_df, aes(ws, 1, fill = ws)) + geom_tile() +
      scale_fill_gradientn(
        colours = pal_wind(200), limits = c(0, 70),
        name = expression("Wind (m s"^{-1}*")"),
        breaks = seq(0, 70, 10)
      ) + theme_void() +
      theme(legend.position = "right",
            legend.title = element_text(size = 9, face = "bold"),
            legend.text  = element_text(size = 8))
    shared_leg <- cowplot::get_legend(pleg)
    p_cmp <- cowplot::plot_grid(
      cowplot::plot_grid(p1, p2, p3, nrow = 1),
      shared_leg, nrow = 1, rel_widths = c(3, 0.35)
    )
    save_both(p_cmp,
              file.path(tc_figdir,
                        paste0("Comparison_", tc_id, "_", ts_pk)),
              w = 14, h = 5)
    p1_ind <- make_panel(pk$holland_speed,
                         paste0("Holland Model: ", tc_name,
                                " (", tc_year, ")"),
                         show_legend = TRUE)
    p2_ind <- make_panel(pk$era5_speed,
                         paste0("ERA5 Reanalysis: ", tc_name,
                                " (", tc_year, ")"),
                         show_legend = TRUE)
    p3_ind <- make_panel(pk$speed,
                         paste0("Blended (Holland+ERA5): ", tc_name,
                                " (", tc_year, ")"),
                         show_legend = TRUE)
    save_both(p1_ind,
              file.path(tc_figdir,
                        paste0("Comparison_Holland_", tc_id, "_",
                               ts_pk)), w = 6, h = 5)
    save_both(p2_ind,
              file.path(tc_figdir,
                        paste0("Comparison_ERA5_", tc_id, "_",
                               ts_pk)), w = 6, h = 5)
    save_both(p3_ind,
              file.path(tc_figdir,
                        paste0("Comparison_Blend_", tc_id, "_",
                               ts_pk)), w = 6, h = 5)
    cat("    Comparison panels saved\n")
  }, error = function(e) cat("    Comparison error:", e$message, "\n"))
  
  
  # F.7  Vmax Time Series (English dates, SSHS labels on left)
  tryCatch({
    band_df <- data.frame(
      ymin = band_ymin, ymax = band_ymax, label = band_labels,
      stringsAsFactors = FALSE
    )
    band_df$ymid <- (band_df$ymin + band_df$ymax) / 2
    
    p_ts <- ggplot() +
      geom_rect(data = band_df,
                aes(xmin = min(tc_h$ISO_TIME),
                    xmax = max(tc_h$ISO_TIME),
                    ymin = ymin, ymax = ymax, fill = label),
                alpha = 0.15) +
      scale_fill_manual(values = SSHS_COLORS[band_labels],
                        name = "SSHS", guide = "none") +
      geom_line(data = tc_h, aes(ISO_TIME, VMAX_MS),
                colour = "#0C2C84", linewidth = 0.8) +
      geom_point(data = tc_h[tc_h$HOUR == 0, ],
                 aes(ISO_TIME, VMAX_MS),
                 colour = "black", fill = "white",
                 shape = 21, size = 2, stroke = 0.7) +
      # SSHS labels on the LEFT side of the plot
      annotate("text",
               x = min(tc_h$ISO_TIME) + 3600 * 3,
               y = band_df$ymid, label = band_df$label,
               size = 2.4, hjust = 0, fontface = "bold",
               colour = SSHS_COLORS[band_df$label]) +
      scale_x_datetime(date_labels = "%b %d",
                       date_breaks = "2 days",
                       expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(limits = c(0, y_display_max),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Date (UTC)",
           y = expression("Maximum sustained wind (m s"^{-1}*")"),
           title = paste0("Intensity Evolution: ", tc_name,
                          " (", tc_year, ")")) +
      theme_pub()
    save_both(p_ts,
              file.path(tc_figdir,
                        paste0("VmaxTimeSeries_", tc_id)),
              w = 8, h = 4.5)
    cat("    Vmax time series saved\n")
  }, error = function(e) cat("    Vmax TS error:", e$message, "\n"))
  
  
  # F.8  Structure Close-up
  tryCatch({
    buf_s <- 5
    si <- which(grid_lon >= pk$tc_lon - buf_s &
                  grid_lon <= pk$tc_lon + buf_s)
    sj <- which(grid_lat >= pk$tc_lat - buf_s &
                  grid_lat <= pk$tc_lat + buf_s)
    if (length(si) >= 5 && length(sj) >= 5) {
      x_km <- (grid_lon[si] - pk$tc_lon) *
        111.32 * cos(pk$tc_lat * pi / 180)
      y_km <- (grid_lat[sj] - pk$tc_lat) * 111.32
      z <- pk$speed[si, sj]
      df_str <- expand.grid(x = x_km, y = y_km)
      df_str$z <- as.vector(z)
      p_str <- ggplot(df_str, aes(x, y)) +
        geom_raster(aes(fill = z)) +
        geom_contour(aes(z = z), colour = "grey20",
                     linewidth = 0.3, bins = 12) +
        scale_fill_gradientn(
          colours = pal_wind(200),
          limits = c(0, max(z, na.rm = TRUE) * 1.05),
          name = expression("Wind (m s"^{-1}*")")
        ) +
        geom_point(aes(x = 0, y = 0), shape = 3, colour = "red",
                   size = 3, stroke = 1.2) +
        annotate("path",
                 x = (pk$rmax_m / 1000) *
                   cos(seq(0, 2 * pi, len = 100)),
                 y = (pk$rmax_m / 1000) *
                   sin(seq(0, 2 * pi, len = 100)),
                 colour = "red", linetype = "dashed",
                 linewidth = 0.5) +
        coord_fixed() +
        labs(x = "East\u2013West distance (km)",
             y = "North\u2013South distance (km)",
             title = paste0("Wind Field Structure (Blended): ",
                            tc_name, " (", tc_year, ")"),
             subtitle = paste(
               format(pk$time, "%Y-%m-%d %H:%M UTC"),
               "| Dashed red circle: RMW =",
               round(pk$rmax_m / 1000), "km")) +
        theme_pub()
      save_both(p_str,
                file.path(tc_figdir,
                          paste0("Structure_", tc_id, "_", ts_pk)),
                w = 6.5, h = 5.5)
      cat("    Structure close-up saved\n")
    }
  }, error = function(e) cat("    Structure error:", e$message, "\n"))
  
  
  # F.9  Pressure-Wind
  tryCatch({
    df_pw <- tc_h[!is.na(tc_h$PMIN_PA) & !is.na(tc_h$VMAX_MS), ]
    if (nrow(df_pw) > 3) {
      df_pw$PMIN_HPA <- df_pw$PMIN_PA / 100
      df_pw$SSHS <- factor(df_pw$SSHS, levels = sshs_present)
      p_pw <- ggplot(df_pw, aes(PMIN_HPA, VMAX_MS)) +
        geom_path(colour = "grey60", linewidth = 0.3) +
        geom_point(aes(colour = SSHS), size = 1.5) +
        scale_colour_manual(values = sshs_col_sub,
                            name = "SSHS Category", drop = TRUE) +
        labs(x = "Minimum central pressure (hPa)",
             y = expression("V"[max]*" (m s"^{-1}*")"),
             title = paste0("Pressure\u2013Wind Relationship: ",
                            tc_name, " (", tc_year, ")")) +
        theme_pub()
      save_both(p_pw,
                file.path(tc_figdir,
                          paste0("PressureWind_", tc_id)),
                w = 6, h = 4.5)
      cat("    Pressure-wind saved\n")
    }
  }, error = function(e) cat("    P-W error:", e$message, "\n"))
  
  
  # F.10  Holland Asymmetry Close-up
  tryCatch({
    buf_a <- 5
    ai <- which(grid_lon >= pk$tc_lon - buf_a &
                  grid_lon <= pk$tc_lon + buf_a)
    aj <- which(grid_lat >= pk$tc_lat - buf_a &
                  grid_lat <= pk$tc_lat + buf_a)
    if (length(ai) >= 5 && length(aj) >= 5) {
      x_km <- (grid_lon[ai] - pk$tc_lon) *
        111.32 * cos(pk$tc_lat * pi / 180)
      y_km <- (grid_lat[aj] - pk$tc_lat) * 111.32
      z_hol <- pk$holland_speed[ai, aj]
      df_hol <- expand.grid(x = x_km, y = y_km)
      df_hol$z <- as.vector(z_hol)
      hdg_rad   <- pk$heading
      arrow_len <- 100
      arrow_dx  <- arrow_len * sin(hdg_rad)
      arrow_dy  <- arrow_len * cos(hdg_rad)
      
      p_asym_base <- ggplot(df_hol, aes(x, y)) +
        geom_raster(aes(fill = z)) +
        geom_contour(aes(z = z), colour = "grey20",
                     linewidth = 0.3, bins = 12) +
        scale_fill_gradientn(
          colours = pal_wind(200),
          limits = c(0, max(z_hol, na.rm = TRUE) * 1.05),
          name = expression("Wind (m s"^{-1}*")")
        ) +
        geom_point(aes(x = 0, y = 0), shape = 3, colour = "red",
                   size = 3, stroke = 1.2) +
        annotate("path",
                 x = (pk$rmax_m / 1000) *
                   cos(seq(0, 2 * pi, len = 100)),
                 y = (pk$rmax_m / 1000) *
                   sin(seq(0, 2 * pi, len = 100)),
                 colour = "red", linetype = "dashed",
                 linewidth = 0.5) +
        coord_fixed() +
        labs(x = "East\u2013West distance (km)",
             y = "North\u2013South distance (km)",
             title = paste0("Holland Model Asymmetry: ", tc_name,
                            " (", tc_year, ")"),
             subtitle = paste(
               format(pk$time, "%Y-%m-%d %H:%M UTC"),
               "| Vt =", round(pk$Vt, 1), "m/s",
               "| Dashed red circle: RMW =",
               round(pk$rmax_m / 1000), "km")) +
        theme_pub()
      
      p_asym_arrow <- p_asym_base +
        annotate("segment",
                 x = 0, y = 0, xend = arrow_dx, yend = arrow_dy,
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "black", linewidth = 1) +
        annotate("text", x = arrow_dx * 1.15, y = arrow_dy * 1.15,
                 label = "Motion", size = 3, fontface = "bold")
      save_both(p_asym_arrow,
                file.path(tc_figdir,
                          paste0("Holland_Asymmetry_Arrow_", tc_id,
                                 "_", ts_pk)),
                w = 6.5, h = 5.5)
      save_both(p_asym_base,
                file.path(tc_figdir,
                          paste0("Holland_Asymmetry_", tc_id,
                                 "_", ts_pk)),
                w = 6.5, h = 5.5)
      cat("    Holland asymmetry saved\n")
    }
  }, error = function(e) cat("    Asymmetry error:", e$message, "\n"))
  
  
  # F.11  Blended Close-up
  tryCatch({
    buf_b <- 5
    bi <- which(grid_lon >= pk$tc_lon - buf_b &
                  grid_lon <= pk$tc_lon + buf_b)
    bj <- which(grid_lat >= pk$tc_lat - buf_b &
                  grid_lat <= pk$tc_lat + buf_b)
    if (length(bi) >= 5 && length(bj) >= 5) {
      x_km <- (grid_lon[bi] - pk$tc_lon) *
        111.32 * cos(pk$tc_lat * pi / 180)
      y_km <- (grid_lat[bj] - pk$tc_lat) * 111.32
      z_bld <- pk$speed[bi, bj]
      df_bld <- expand.grid(x = x_km, y = y_km)
      df_bld$z <- as.vector(z_bld)
      p_bld_close <- ggplot(df_bld, aes(x, y)) +
        geom_raster(aes(fill = z)) +
        geom_contour(aes(z = z), colour = "grey20",
                     linewidth = 0.3, bins = 12) +
        scale_fill_gradientn(
          colours = pal_wind(200),
          limits = c(0, max(z_bld, na.rm = TRUE) * 1.05),
          name = expression("Wind (m s"^{-1}*")")
        ) +
        geom_point(aes(x = 0, y = 0), shape = 3, colour = "red",
                   size = 3, stroke = 1.2) +
        annotate("path",
                 x = (pk$rmax_m / 1000) *
                   cos(seq(0, 2 * pi, len = 100)),
                 y = (pk$rmax_m / 1000) *
                   sin(seq(0, 2 * pi, len = 100)),
                 colour = "red", linetype = "dashed",
                 linewidth = 0.5) +
        coord_fixed() +
        labs(x = "East\u2013West distance (km)",
             y = "North\u2013South distance (km)",
             title = paste0("Blended Wind Field (Holland+ERA5): ",
                            tc_name, " (", tc_year, ")"),
             subtitle = paste(
               format(pk$time, "%Y-%m-%d %H:%M UTC"),
               "| Vmax =", round(pk$vmax, 1), "m/s",
               "| Dashed red circle: RMW")) +
        theme_pub()
      save_both(p_bld_close,
                file.path(tc_figdir,
                          paste0("Blend_Closeup_", tc_id,
                                 "_", ts_pk)),
                w = 6.5, h = 5.5)
      cat("    Blended close-up saved\n")
    }
  }, error = function(e) cat("    Blend close-up error:", e$message, "\n"))
  
  
  # ===========================================================
  #  STEP G: Cross-validation (SEPARATE plots + CSV + combined)
  # ===========================================================
  val_df    <- do.call(rbind, val_records)
  val_clean <- val_df[!is.na(val_df$era5_vmax), ]
  
  # Save per-TC validation data as CSV
  save_data_csv(val_df,
                file.path(DIR_VALID,
                          paste0("Valid_data_", tc_id)))
  cat("    Validation data CSV saved\n")
  
  if (nrow(val_clean) >= 3) {
    tryCatch({
      # --- Individual scatter: ERA5 vs IBTrACS ---
      p_v1 <- make_valid_scatter(
        val_clean, "ibt_vmax", "era5_vmax", "#E6550D",
        paste0("ERA5 vs IBTrACS: ", tc_name, " (", tc_year, ")"))
      if (!is.null(p_v1))
        save_both(p_v1,
                  file.path(DIR_VALID,
                            paste0("Valid_ERA5_", tc_id)),
                  w = 5.5, h = 5)
      
      # --- Individual scatter: Holland vs IBTrACS ---
      p_v2 <- make_valid_scatter(
        val_clean, "ibt_vmax", "holland_vmax", "#3182BD",
        paste0("Holland vs IBTrACS: ", tc_name, " (", tc_year, ")"))
      if (!is.null(p_v2))
        save_both(p_v2,
                  file.path(DIR_VALID,
                            paste0("Valid_Holland_", tc_id)),
                  w = 5.5, h = 5)
      
      # --- Individual scatter: Blend vs IBTrACS ---
      # NOTE: Holland Vmax ~ Blend Vmax is EXPECTED because peak wind
      # occurs at RMW where blending weight xi ~ 0 (Holland dominates).
      # The blending improves outer radii, not the peak.
      p_v3 <- make_valid_scatter(
        val_clean, "ibt_vmax", "blend_vmax", "#0C2C84",
        paste0("Holland+ERA5 vs IBTrACS: ", tc_name,
               " (", tc_year, ")"))
      if (!is.null(p_v3))
        save_both(p_v3,
                  file.path(DIR_VALID,
                            paste0("Valid_Blend_", tc_id)),
                  w = 5.5, h = 5)
      
      # --- Combined: ERA5 & Holland+ERA5 on same plot ---
      p_v4 <- make_combined_scatter(
        val_clean,
        paste0("ERA5 vs Holland+ERA5: ", tc_name,
               " (", tc_year, ")"))
      if (!is.null(p_v4))
        save_both(p_v4,
                  file.path(DIR_VALID,
                            paste0("Valid_ERA5_vs_Blend_", tc_id)),
                  w = 6, h = 5)
      
      # --- Area-mean at outer radii (200-500 km) ---
      val_outer <- val_clean[!is.na(val_clean$era5_mean_200_500), ]
      if (nrow(val_outer) >= 3) {
        p_v5_era <- make_valid_scatter(
          val_outer, "era5_mean_200_500", "holland_mean_200_500",
          "#3182BD",
          paste0("Outer wind (200-500 km): Holland vs ERA5\n",
                 tc_name, " (", tc_year, ")"),
          xlab_str = expression("ERA5 mean wind 200-500 km (m s"^{-1}*")"),
          ylab_str = expression("Holland mean wind 200-500 km (m s"^{-1}*")"))
        if (!is.null(p_v5_era))
          save_both(p_v5_era,
                    file.path(DIR_VALID,
                              paste0("Valid_OuterWind_HolVsERA5_",
                                     tc_id)),
                    w = 5.5, h = 5)
        
        p_v5_bld <- make_valid_scatter(
          val_outer, "era5_mean_200_500", "blend_mean_200_500",
          "#0C2C84",
          paste0("Outer wind (200-500 km): Blend vs ERA5\n",
                 tc_name, " (", tc_year, ")"),
          xlab_str = expression("ERA5 mean wind 200-500 km (m s"^{-1}*")"),
          ylab_str = expression("Blend mean wind 200-500 km (m s"^{-1}*")"))
        if (!is.null(p_v5_bld))
          save_both(p_v5_bld,
                    file.path(DIR_VALID,
                              paste0("Valid_OuterWind_BldVsERA5_",
                                     tc_id)),
                    w = 5.5, h = 5)
      }
      
      cat("    Validation plots saved\n")
    }, error = function(e) cat("    Validation error:", e$message, "\n"))
  }
  
  # ---- STEP H: Accumulate ----
  all_validation <- rbind(all_validation, val_df)
  peak_vmax <- max(tc_h$VMAX_MS, na.rm = TRUE)
  sshs_cat  <- as.character(cut(peak_vmax, breaks = SSHS_BREAKS,
                                labels = SSHS_LABELS, right = FALSE))
  tc_summary <- rbind(tc_summary, data.frame(
    TC_ID = tc_id, Name = tc_name, Year = tc_year, SID = tc_sid,
    N_Records = nrow(tc), Vmax_ms = round(peak_vmax, 1),
    SSHS = sshs_cat,
    Start = format(t0, "%Y-%m-%d"), End = format(t1, "%Y-%m-%d"),
    N_Steps = n_steps, stringsAsFactors = FALSE
  ))
  cat(sprintf("  >> %s completed (%s, Vmax = %.1f m/s)\n",
              tc_id, sshs_cat, peak_vmax))
  
}  # === end TC loop ===


# ################################################################
#  PART 5: GLOBAL SUMMARY
# ################################################################

cat("\n", paste(rep("=", 65), collapse = ""), "\n")
cat("  ALL PROCESSING COMPLETE\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n")

total_min <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
cat(sprintf("Total: %.1f min for %d TCs (%.1f min/TC)\n\n",
            total_min, nrow(tc_summary),
            total_min / max(nrow(tc_summary), 1)))

print(tc_summary)
write.csv(tc_summary,
          file.path(DIR_RESULT, "TC_Processing_Summary.csv"),
          row.names = FALSE)

# --- Combined cross-validation ---
val_all <- all_validation[!is.na(all_validation$era5_vmax), ]

# Save combined validation data
save_data_csv(all_validation,
              file.path(DIR_VALID, "Validation_AllTCs_data"))

if (nrow(val_all) > 5) {
  tryCatch({
    mx <- max(c(val_all$ibt_vmax, val_all$era5_vmax,
                val_all$holland_vmax, val_all$blend_vmax),
              na.rm = TRUE) * 1.08
    
    # --- ERA5 vs IBTrACS (all TCs) ---
    p_all_era <- make_valid_scatter(
      val_all, "ibt_vmax", "era5_vmax", "#E6550D",
      paste0("ERA5 vs IBTrACS (All TCs, N=", nrow(val_all), ")"))
    if (!is.null(p_all_era))
      save_both(p_all_era,
                file.path(DIR_VALID, "Validation_AllTCs_ERA5"),
                w = 5.5, h = 5)
    
    # --- Holland vs IBTrACS (all TCs) ---
    p_all_hol <- make_valid_scatter(
      val_all, "ibt_vmax", "holland_vmax", "#3182BD",
      paste0("Holland vs IBTrACS (All TCs, N=", nrow(val_all), ")"))
    if (!is.null(p_all_hol))
      save_both(p_all_hol,
                file.path(DIR_VALID, "Validation_AllTCs_Holland"),
                w = 5.5, h = 5)
    
    # --- Blend vs IBTrACS (all TCs) ---
    p_all_bld <- make_valid_scatter(
      val_all, "ibt_vmax", "blend_vmax", "#0C2C84",
      paste0("Holland+ERA5 vs IBTrACS (All TCs, N=",
             nrow(val_all), ")"))
    if (!is.null(p_all_bld))
      save_both(p_all_bld,
                file.path(DIR_VALID, "Validation_AllTCs_Blend"),
                w = 5.5, h = 5)
    
    # --- Combined: ERA5 & Blend on same plot (all TCs) ---
    p_all_comb <- make_combined_scatter(
      val_all,
      paste0("ERA5 vs Holland+ERA5 (All TCs, N=",
             nrow(val_all), ", ", nrow(tc_summary), " TCs)"))
    if (!is.null(p_all_comb))
      save_both(p_all_comb,
                file.path(DIR_VALID, "Validation_AllTCs_ERA5_vs_Blend"),
                w = 6, h = 5)
    
    # --- Outer-radii area-mean validation (all TCs) ---
    val_outer_all <- val_all[!is.na(val_all$era5_mean_200_500), ]
    if (nrow(val_outer_all) >= 5) {
      p_outer_hol <- make_valid_scatter(
        val_outer_all, "era5_mean_200_500", "holland_mean_200_500",
        "#3182BD",
        paste0("Outer wind (200-500 km): Holland vs ERA5 (All TCs)"),
        xlab_str = expression("ERA5 mean wind 200-500 km (m s"^{-1}*")"),
        ylab_str = expression("Holland mean wind 200-500 km (m s"^{-1}*")"))
      if (!is.null(p_outer_hol))
        save_both(p_outer_hol,
                  file.path(DIR_VALID,
                            "Validation_AllTCs_OuterWind_HolVsERA5"),
                  w = 5.5, h = 5)
      
      p_outer_bld <- make_valid_scatter(
        val_outer_all, "era5_mean_200_500", "blend_mean_200_500",
        "#0C2C84",
        paste0("Outer wind (200-500 km): Blend vs ERA5 (All TCs)"),
        xlab_str = expression("ERA5 mean wind 200-500 km (m s"^{-1}*")"),
        ylab_str = expression("Blend mean wind 200-500 km (m s"^{-1}*")"))
      if (!is.null(p_outer_bld))
        save_both(p_outer_bld,
                  file.path(DIR_VALID,
                            "Validation_AllTCs_OuterWind_BldVsERA5"),
                  w = 5.5, h = 5)
    }
    
    # --- Summary statistics ---
    df_stats <- data.frame(
      Source = c("ERA5", "Holland", "Holland+ERA5"),
      r = c(cor(val_all$ibt_vmax, val_all$era5_vmax, use = "complete"),
            cor(val_all$ibt_vmax, val_all$holland_vmax, use = "complete"),
            cor(val_all$ibt_vmax, val_all$blend_vmax, use = "complete")),
      RMSE = c(sqrt(mean((val_all$ibt_vmax - val_all$era5_vmax)^2,
                         na.rm = TRUE)),
               sqrt(mean((val_all$ibt_vmax - val_all$holland_vmax)^2,
                         na.rm = TRUE)),
               sqrt(mean((val_all$ibt_vmax - val_all$blend_vmax)^2,
                         na.rm = TRUE))),
      Bias = c(mean(val_all$era5_vmax - val_all$ibt_vmax, na.rm = TRUE),
               mean(val_all$holland_vmax - val_all$ibt_vmax, na.rm = TRUE),
               mean(val_all$blend_vmax - val_all$ibt_vmax, na.rm = TRUE)),
      N = nrow(val_all)
    )
    cat("\n--- Vmax Validation Statistics ---\n")
    print(df_stats)
    save_data_csv(df_stats,
                  file.path(DIR_VALID, "Validation_AllTCs_stats"))
    
    cat("\nNOTE: Holland Vmax ~ Blend Vmax is EXPECTED.\n")
    cat("  Peak wind occurs at RMW (~30-60 km) where blending\n")
    cat("  weight xi ~ 0 (Holland contributes >99%).\n")
    cat("  Blending improves outer-radii winds, not Vmax.\n")
    cat("  See outer-radii (200-500 km) plots for differences.\n\n")
    
    cat("Combined validation saved\n")
  }, error = function(e) cat("Combined validation error:",
                             e$message, "\n"))
}

# --- Holland B Sensitivity ---
tryCatch({
  vmax_test <- 40;  lat_test <- 22
  r_seq_m   <- seq(500, 500000, by = 500)
  df_sens   <- data.frame()
  for (rmax_km in c(20, 35, 50, 70)) {
    B_val  <- est_holland_B(vmax_test, rmax_km, lat_test)
    rmax_m <- rmax_km * 1000
    vg     <- holland_Vg(r_seq_m, vmax_test / RF_DEFAULT,
                         rmax_m, B_val)
    vs     <- RF_DEFAULT * vg
    df_sens <- rbind(df_sens, data.frame(
      r_km = r_seq_m / 1000, ws = vs,
      label = paste0("Rmax=", rmax_km,
                     " km, B=", round(B_val, 2))
    ))
  }
  p_sens <- ggplot(df_sens, aes(r_km, ws, colour = label)) +
    geom_line(linewidth = 0.7) +
    labs(x = "Radius (km)",
         y = expression("Surface wind speed (m s"^{-1}*")"),
         title = paste0("Holland Profile Sensitivity (Vmax = ",
                        vmax_test, " m/s, lat = ",
                        lat_test, "\u00B0N)"),
         colour = NULL) +
    theme_pub() +
    theme(legend.position = c(0.7, 0.8),
          legend.background = element_rect(
            fill = alpha("white", 0.85), colour = NA))
  save_both(p_sens, file.path(DIR_FIG, "Holland_B_Sensitivity"),
            w = 6.5, h = 4.5)
  cat("B sensitivity plot saved\n")
}, error = function(e) cat("B sensitivity error:", e$message, "\n"))

# --- Blending Weight ---
tryCatch({
  rmax_demo <- 40
  r_demo    <- seq(0, 800, by = 1)
  xi_demo   <- blend_xi(r_demo * 1000, rmax_demo * 1000)
  tr_r      <- N_BLEND * rmax_demo
  
  # Compute dynamic cutoff for all three methods
  vmax_demo <- 40
  vmax_grad_demo <- vmax_demo / RF_DEFAULT
  B_demo <- est_holland_B(vmax_demo, rmax_demo, 22)
  
  blend_cut <- calc_blend_cutoff(rmax_demo * 1000) / 1000
  wind_cut  <- calc_wind_cutoff(vmax_grad_demo, rmax_demo * 1000,
                                B_demo) / 1000
  
  p_blend <- ggplot(data.frame(r = r_demo, xi = xi_demo),
                    aes(r, xi)) +
    geom_line(colour = "#0C2C84", linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey50") +
    geom_vline(xintercept = tr_r, linetype = "dashed",
               colour = "red") +
    geom_vline(xintercept = blend_cut, linetype = "dotted",
               colour = "darkgreen", linewidth = 0.6) +
    geom_vline(xintercept = wind_cut, linetype = "dotdash",
               colour = "darkorange", linewidth = 0.6) +
    annotate("text", x = tr_r + 8, y = 0.55,
             label = paste0("r = ", N_BLEND,
                            " \u00D7 Rmax\n= ", tr_r, " km"),
             colour = "red", size = 3, hjust = 0) +
    annotate("text", x = blend_cut + 8, y = 0.98,
             label = paste0("BLEND cutoff\n= ",
                            round(blend_cut), " km"),
             colour = "darkgreen", size = 2.6, hjust = 0) +
    annotate("text", x = wind_cut + 8, y = 0.85,
             label = paste0("WIND_SPEED cutoff\n= ",
                            round(wind_cut), " km\n(V=",
                            WIND_SPEED_THRESHOLD, " m/s)"),
             colour = "darkorange", size = 2.6, hjust = 0) +
    annotate("text", x = 25, y = 0.12,
             label = "Holland\ndominates",
             colour = "#0C2C84", size = 3,
             fontface = "italic") +
    annotate("text", x = 600, y = 0.88,
             label = "ERA5\ndominates",
             colour = "#E6550D", size = 3,
             fontface = "italic") +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    labs(x = "Distance from TC centre (km)",
         y = expression("Blending weight (" * xi * ")"),
         title = paste0("Blending Weight (Rmax = ",
                        rmax_demo, " km, n = ", N_BLEND,
                        ", Vmax = ", vmax_demo, " m/s)")) +
    theme_pub()
  save_both(p_blend, file.path(DIR_FIG, "Blending_Weight"),
            w = 7, h = 4.5)
  cat("Blending weight plot saved\n")
}, error = function(e) cat("Blend weight error:", e$message, "\n"))

cat("\n\u2714 Step 1 complete.\n")
















