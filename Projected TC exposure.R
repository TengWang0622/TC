#############################################################################################
#
#                          TROPICAL CYCLONES AND HEALTH
#         Population Exposure, Panel Construction  (Step 2.1 v2.3 — Accelerated, FIXED)
#
#############################################################################################
#
# Developed by Teng Wang  @ The University of Hong Kong
#              Hanxu Shi  @ Peking University
# Version: 2.3 (Accelerated future projection — coastal-exposure bug fixed)
# Last updated: 2025-04-28
#
#############################################################################################


# ================================================================
#  PART 0: PACKAGES
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(ggplot2); library(sf); library(raster)
  library(ncdf4); library(exactextractr); library(openxlsx)
  library(scales); library(cowplot); library(maps)
  library(readxl); library(stringr); library(ggrepel)
})
options(scipen = 6)


# ================================================================
#  PART 1: CONFIGURATION
# ================================================================

setwd("C:/Project/Tropical cyclone")

YEARS      <- 2025:2050
N_YEARS    <- length(YEARS)
SSP_TAG    <- "SSP585"

PERIOD_LABEL   <- "CMCC-CM2-VHR4-highres-TRACK"      # ---------------------- Input

DIR_NC      <- paste0("Result HighResMIP/", PERIOD_LABEL, "/NetCDF")
#DIR_TRACK      <- "Src/TC future/CMCC-CM2-VHR4-highresSST-TRACK"
DIR_TRACK      <- paste0("Src/TC future/", PERIOD_LABEL)

DIR_SHP       <- "C:/Project/Database/MapExtract/Src/ChinaMap"
FILE_CMIP6_POP_NC <- "Src/Pop_CMIP6/grid_pop_age_gender_ssp585.nc"

FILE_POP_CNTY_AGE_FUTURE <- file.path(
  "Src/Pop_CMIP6/Processed/Extract",
  paste0("Pop_County_", SSP_TAG, "_long.xlsx"))

FILE_POP_CITY_AGE_FUTURE <- file.path(
  "Src/Pop_CMIP6/Processed/Extract",
  paste0("Pop_City_", SSP_TAG, "_long.xlsx"))

DIR_OUT        <- paste0("Result HighResMIP/", PERIOD_LABEL, "/Exposure")
DIR_FIG_EXP    <- paste0("Result HighResMIP/", PERIOD_LABEL, "/Exposure")
DIR_MAP        <- paste0("Result HighResMIP/", PERIOD_LABEL, "/Maps")
DIR_MAP_TC     <- file.path(DIR_MAP, "TC_Individual")
DIR_TABLE      <- paste0("Result HighResMIP/", PERIOD_LABEL, "/Tables")
DIR_DATA       <- paste0("Result HighResMIP/", PERIOD_LABEL, "/PlotData")
DIR_PANEL      <- paste0("Result HighResMIP/", PERIOD_LABEL, "/Panel")

for (d in c(DIR_OUT, DIR_FIG_EXP, DIR_MAP, DIR_MAP_TC, DIR_TABLE,
            DIR_DATA, DIR_PANEL))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

GRID_RES   <- 0.25
DOMAIN_LON <- c(70, 150)
DOMAIN_LAT <- c(5, 55)
WS_THRESHOLD <- 17.5
KT2MS      <- 0.514444
HALF_RES   <- GRID_RES / 2

grid_lon <- seq(DOMAIN_LON[1], DOMAIN_LON[2], by = GRID_RES)
grid_lat <- seq(DOMAIN_LAT[1], DOMAIN_LAT[2], by = GRID_RES)
nlon     <- length(grid_lon)
nlat     <- length(grid_lat)

SSHS_BREAKS <- c(0, 17.5, 33.0, 43.0, 50.0, 58.0, 70.0, Inf)
SSHS_LABELS <- c("TD", "TS", "Category 1", "Category 2",
                 "Category 3", "Category 4", "Category 5")
SSHS_COLORS <- c("TD"         = "#6BAED6", "TS"         = "#3182BD",
                 "Category 1" = "#FDAE6B", "Category 2" = "#FD8D3C",
                 "Category 3" = "#E6550D", "Category 4" = "#A63603",
                 "Category 5" = "#67000D")

COL_CNTY  <- "grey88";  LW_CNTY <- 0.06
COL_CITY  <- "grey55";  LW_CITY <- 0.12
COL_NATL  <- "grey55";  LW_NATL <- 0.12

CHINA_LON <- c(72, 140)
CHINA_LAT <- c(15, 55)
FULL_LON  <- DOMAIN_LON
FULL_LAT  <- DOMAIN_LAT

EXCLUDE_PROV <- c("香港特别行政区", "澳门特别行政区", "台湾省",
                  "香港", "澳门", "台湾")

TC_SEASON_START <- 1
TC_SEASON_END   <- 12
MAX_LAG  <- 21
MAX_LEAD <- 3

CMIP6_LON   <- seq(73.75, 135.25, by = 0.5)
CMIP6_LAT   <- seq(16.25, 53.75,  by = 0.5)
CMIP6_YEARS <- 2010:2100
CMIP6_RES   <- 0.5


# ================================================================
#  PART 2: UTILITY FUNCTIONS
# ================================================================

pal_wind <- colorRampPalette(
  c("white","#E0F0FF","#87CEEB","#4682B4","#1E5BA8","#0C2C84","#081D58"))

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text          = element_text(colour = "black"),
      plot.title    = element_text(size = base_size + 2, face = "bold",
                                   hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 0.5, hjust = 0,
                                   margin = margin(b = 6), colour = "grey30"),
      axis.title    = element_text(size = base_size, face = "bold"),
      axis.text     = element_text(size = base_size - 1, colour = "black"),
      legend.title  = element_text(size = base_size - 0.5, face = "bold"),
      legend.text   = element_text(size = base_size - 1.5),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.margin   = margin(5, 5, 5, 5, "mm"),
      strip.text    = element_text(size = base_size, face = "bold"))
}

save_both <- function(p, path_base, w = 7, h = 5.5) {
  tryCatch(ggsave(paste0(path_base, ".pdf"), p, width = w, height = h,
                  dpi = 300, device = cairo_pdf),
           error = function(e)
             tryCatch(ggsave(paste0(path_base, ".pdf"), p, width = w,
                             height = h, dpi = 300),
                      error = function(e2) cat("  PDF fail\n")))
  tryCatch(ggsave(paste0(path_base, ".png"), p, width = w, height = h,
                  dpi = 300, bg = "white"),
           error = function(e) cat("  PNG fail\n"))
  while (dev.cur() > 1) tryCatch(dev.off(), error = function(e) break)
}

save_data <- function(df, name, dir = DIR_DATA) {
  tryCatch(write.xlsx(as.data.frame(df),
                      file.path(dir, paste0(name, ".xlsx"))),
           error = function(e)
             tryCatch(write.csv(df, file.path(dir, paste0(name, ".csv")),
                                row.names = FALSE),
                      error = function(e2)
                        cat("  Data save fail:", name, "\n")))
}

haversine_dist <- function(lon1, lat1, lon2, lat2) {
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  la1  <- lat1 * pi / 180; la2 <- lat2 * pi / 180
  a <- sin(dlat / 2)^2 + cos(la1) * cos(la2) * sin(dlon / 2)^2
  2 * 6371000 * asin(pmin(sqrt(a), 1))
}

mat_to_raster <- function(mat, glon, glat) {
  if (!is.matrix(mat))
    mat <- matrix(mat, nrow = length(glon), ncol = length(glat))
  mat_t <- t(mat)[length(glat):1, ]
  raster(mat_t,
         xmn = min(glon) - HALF_RES, xmx = max(glon) + HALF_RES,
         ymn = min(glat) - HALF_RES, ymx = max(glat) + HALF_RES,
         crs = "+proj=longlat +datum=WGS84")
}

mat2df <- function(mat, glon, glat, varname = "value") {
  df <- expand.grid(lon = glon, lat = glat)
  df[[varname]] <- as.vector(mat)
  df
}

sshs_classify <- function(ws) {
  as.character(cut(ws, breaks = SSHS_BREAKS, labels = SSHS_LABELS,
                   right = FALSE))
}

nc_get_var <- function(ncf, var_names) {
  for (vn in var_names) {
    val <- tryCatch(ncvar_get(ncf, vn), error = function(e) NULL)
    if (!is.null(val)) return(val)
  }
  NULL
}

ensure_matrix <- function(x, nlon, nlat) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x) && nrow(x) == nlon && ncol(x) == nlat) return(x)
  if (length(x) == nlon * nlat) return(matrix(x, nrow = nlon, ncol = nlat))
  return(NULL)
}

check_landfalling <- function(lon_vec, lat_vec, boundary) {
  if (is.null(boundary) || length(lon_vec) == 0)
    return(list(lf = FALSE, idx = NA_integer_))
  pts <- st_as_sf(data.frame(lon = lon_vec, lat = lat_vec),
                  coords = c("lon", "lat"), crs = 4326)
  inside <- tryCatch(
    st_intersects(pts, boundary, sparse = FALSE)[, 1],
    error = function(e) rep(FALSE, length(lon_vec)))
  if (any(inside)) list(lf = TRUE, idx = which(inside)[1])
  else             list(lf = FALSE, idx = NA_integer_)
}

make_track_segs <- function(tk) {
  segs <- list()
  for (sg in unique(tk$segment)) {
    s <- tk[tk$segment == sg, ]
    if (nrow(s) < 2) next
    for (j in 2:nrow(s)) {
      segs[[length(segs) + 1]] <- data.frame(
        x = s$lon[j-1], y = s$lat[j-1],
        xend = s$lon[j], yend = s$lat[j],
        sshs = s$sshs[j], tc_id = s$tc_id[j],
        tc_name = s$tc_name[j], tc_year = s$tc_year[j],
        stringsAsFactors = FALSE)
    }
  }
  if (length(segs) == 0) return(NULL)
  do.call(rbind, segs)
}

fmt_lon <- function(x) paste0(x, "\u00B0E")
fmt_lat <- function(y) ifelse(y >= 0, paste0(y, "\u00B0N"), paste0(-y, "\u00B0S"))

make_tc_labels <- function(tc_name, tc_year) {
  lab <- paste0(tc_name, " (", tc_year, ")")
  make.unique(lab, sep = " #")
}


# ================================================================
#  PART 3: LOAD ADMINISTRATIVE BOUNDARIES & POPULATION AGE DATA
# ================================================================

cat("Loading shapefiles...\n")

county_sf <- st_read(file.path(DIR_SHP, "县.shp"), quiet = TRUE)
city_sf   <- st_read(file.path(DIR_SHP, "市.shp"), quiet = TRUE)

county_sf <- st_transform(county_sf, 4326)
city_sf   <- st_transform(city_sf,   4326)

county_sf$省代码 <- as.character(county_sf$省代码)
county_sf$市代码 <- as.character(county_sf$市代码)
county_sf$县代码 <- as.character(county_sf$县代码)
city_sf$省代码   <- as.character(city_sf$省代码)
city_sf$市代码   <- as.character(city_sf$市代码)

cat(sprintf("  Counties: %d | Cities: %d\n", nrow(county_sf), nrow(city_sf)))

china_land <- tryCatch({
  cv <- st_make_valid(county_sf)
  cl <- st_union(cv)
  if (!is.null(cl) && !st_is_empty(cl)) cl else NULL
}, error = function(e) {
  cat("  Warning: st_union failed, trying simplified...\n")
  tryCatch({
    cs <- st_simplify(st_make_valid(county_sf), dTolerance = 0.01)
    st_union(cs)
  }, error = function(e2) NULL)
})

china_land_buf <- NULL
if (!is.null(china_land)) {
  china_land_buf <- tryCatch(st_buffer(china_land, 0.15),
                             error = function(e) china_land)
  cat("  China land boundary OK (with 0.15-deg buffer)\n")
} else {
  cat("  WARNING: China land boundary creation FAILED\n")
}

age_name_map <- c("0-19"  = "<20 years",
                  "20-44" = "20 to 44 years",
                  "45-64" = "45 to 64 years",
                  "65+"   = ">64 years")

cat("\nLoading future population data (SSP5-8.5)...\n")

pop_county_age <- NULL
pop_city_age   <- NULL

if (file.exists(FILE_POP_CNTY_AGE_FUTURE)) {
  
  pop_county_long <- read.xlsx(FILE_POP_CNTY_AGE_FUTURE)
  pop_county_long$CountyCode <- as.character(pop_county_long$CountyCode)
  pop_county_long$CityCode   <- as.character(pop_county_long$CityCode)
  pop_county_long$Year       <- as.integer(pop_county_long$Year)
  pop_county_long$Agegroup   <- as.character(pop_county_long$Agegroup)
  pop_county_long$Population <- as.numeric(pop_county_long$Population)
  
  cat(sprintf("  County pop long: %d rows  |  Years: %d-%d\n",
              nrow(pop_county_long),
              min(pop_county_long$Year), max(pop_county_long$Year)))
  
  pop_county_long$AgeLabel <- age_name_map[pop_county_long$Agegroup]
  
  pop_county_age <- pop_county_long %>%
    mutate(ColName = paste0(AgeLabel, "_", Year)) %>%
    dplyr::select(CountyCode, ColName, Population) %>%
    pivot_wider(id_cols = CountyCode,
                names_from = ColName, values_from = Population,
                values_fn = sum)
  pop_county_age$CountyCode <- as.character(pop_county_age$CountyCode)
  
  cat(sprintf("  County pop wide: %d units x %d columns\n",
              nrow(pop_county_age), ncol(pop_county_age)))
  
  rm(pop_county_long); gc()
  
} else {
  cat("  WARNING: Future county pop file NOT found:\n")
  cat("    ", FILE_POP_CNTY_AGE_FUTURE, "\n")
  cat("  Age-specific population will NOT be available.\n")
}

if (file.exists(FILE_POP_CITY_AGE_FUTURE)) {
  
  cat("  Reading city pop directly from:", basename(FILE_POP_CITY_AGE_FUTURE), "\n")
  
  pop_city_long <- read.xlsx(FILE_POP_CITY_AGE_FUTURE)
  pop_city_long$CityCode   <- as.character(pop_city_long$CityCode)
  pop_city_long$Year       <- as.integer(pop_city_long$Year)
  pop_city_long$Agegroup   <- as.character(pop_city_long$Agegroup)
  pop_city_long$Population <- as.numeric(pop_city_long$Population)
  
  cat(sprintf("  City pop long: %d rows  |  Years: %d-%d\n",
              nrow(pop_city_long),
              min(pop_city_long$Year), max(pop_city_long$Year)))
  
  pop_city_long$AgeLabel <- age_name_map[pop_city_long$Agegroup]
  
  pop_city_age <- pop_city_long %>%
    mutate(ColName = paste0(AgeLabel, "_", Year)) %>%
    dplyr::select(CityCode, ColName, Population) %>%
    pivot_wider(id_cols = CityCode,
                names_from = ColName, values_from = Population,
                values_fn = sum)
  pop_city_age$CityCode <- as.character(pop_city_age$CityCode)
  
  cat(sprintf("  City pop wide (direct file): %d units x %d columns\n",
              nrow(pop_city_age), ncol(pop_city_age)))
  
  rm(pop_city_long); gc()
  
} else if (!is.null(pop_county_age) && file.exists(FILE_POP_CNTY_AGE_FUTURE)) {
  
  cat("  City pop file NOT found, aggregating from county data...\n")
  
  pop_county_long_tmp <- read.xlsx(FILE_POP_CNTY_AGE_FUTURE)
  pop_county_long_tmp$CityCode   <- as.character(pop_county_long_tmp$CityCode)
  pop_county_long_tmp$Year       <- as.integer(pop_county_long_tmp$Year)
  pop_county_long_tmp$Agegroup   <- as.character(pop_county_long_tmp$Agegroup)
  pop_county_long_tmp$Population <- as.numeric(pop_county_long_tmp$Population)
  pop_county_long_tmp$AgeLabel   <- age_name_map[pop_county_long_tmp$Agegroup]
  
  pop_city_long_agg <- pop_county_long_tmp %>%
    group_by(CityCode, Year, AgeLabel) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop")
  
  pop_city_age <- pop_city_long_agg %>%
    mutate(ColName = paste0(AgeLabel, "_", Year)) %>%
    dplyr::select(CityCode, ColName, Population) %>%
    pivot_wider(id_cols = CityCode,
                names_from = ColName, values_from = Population,
                values_fn = sum)
  pop_city_age$CityCode <- as.character(pop_city_age$CityCode)
  
  cat(sprintf("  City pop wide (aggregated): %d units x %d columns\n",
              nrow(pop_city_age), ncol(pop_city_age)))
  
  rm(pop_county_long_tmp, pop_city_long_agg); gc()
  
} else {
  cat("  WARNING: Neither city nor county pop file found.\n")
}


# ================================================================
#  PART 4: POPULATION WEIGHT RASTERS
#          ★ FIX v2.3: NA cells filled with 0 after resample;
#                       cache filename bumped to force rebuild
# ================================================================

pop_cache_dir <- file.path(DIR_OUT, "pop_cache")
if (!dir.exists(pop_cache_dir)) dir.create(pop_cache_dir, recursive = TRUE)

ref_raster <- raster(nrows = nlat, ncols = nlon,
                     xmn = min(grid_lon) - HALF_RES,
                     xmx = max(grid_lon) + HALF_RES,
                     ymn = min(grid_lat) - HALF_RES,
                     ymx = max(grid_lat) + HALF_RES,
                     crs = "+proj=longlat +datum=WGS84")

cat("\nReading CMIP6 SSP5-8.5 gridded population...\n")

pop_total_cmip6 <- NULL

if (file.exists(FILE_CMIP6_POP_NC)) {
  
  nc_pop <- nc_open(FILE_CMIP6_POP_NC)
  pop_raw <- ncvar_get(nc_pop, "grid_pop_age_gender")
  nc_close(nc_pop)
  
  cat(sprintf("  NC dims: %s\n", paste(dim(pop_raw), collapse = " x ")))
  
  pop_both <- pop_raw[,,,, 1] + pop_raw[,,,, 2]
  rm(pop_raw); gc()
  
  pop_total_cmip6 <- apply(pop_both, c(1, 2, 4), sum, na.rm = TRUE)
  rm(pop_both); gc()
  
  cat(sprintf("  Total pop grid: %d x %d x %d years\n",
              dim(pop_total_cmip6)[1], dim(pop_total_cmip6)[2],
              dim(pop_total_cmip6)[3]))
  
} else {
  cat("  WARNING: CMIP6 pop NC not found\n")
}

build_pop_weight_cmip6 <- function(yr) {
  # ★ FIX v2.3: new cache filename — old cached tifs (which may contain NA
  #   cells outside the CMIP6 73.5-135.5E / 16-54N box) are bypassed.
  cache_tif <- file.path(pop_cache_dir,
                         paste0("pop_weight_025deg_cmip6_v23_", yr, ".tif"))
  if (file.exists(cache_tif)) return(raster(cache_tif))
  if (is.null(pop_total_cmip6)) return(NULL)
  
  yr_idx <- which(CMIP6_YEARS == yr)
  if (length(yr_idx) == 0) {
    cat(sprintf("  Pop: year %d not in CMIP6 range, using nearest\n", yr))
    yr_idx <- which.min(abs(CMIP6_YEARS - yr))
  }
  
  mat <- pop_total_cmip6[,, yr_idx]
  mat[is.na(mat)] <- 0
  
  mat_t <- t(mat)[length(CMIP6_LAT):1, ]
  r05 <- raster(mat_t,
                xmn = min(CMIP6_LON) - CMIP6_RES / 2,
                xmx = max(CMIP6_LON) + CMIP6_RES / 2,
                ymn = min(CMIP6_LAT) - CMIP6_RES / 2,
                ymx = max(CMIP6_LAT) + CMIP6_RES / 2,
                crs = "+proj=longlat +datum=WGS84")
  
  pop_aligned <- resample(r05, ref_raster, method = "bilinear")
  
  # ★ FIX v2.3: cells outside the CMIP6 source extent become NA after
  #   resample. NA weights poison exact_extract 'weighted_mean' (returns NA)
  #   and coastal counties are then silently dropped by the keep-filter.
  #   Fill all NA with 0 so the weight raster is defined over the FULL domain.
  pop_aligned[is.na(pop_aligned)] <- 0
  
  writeRaster(pop_aligned, cache_tif, format = "GTiff", overwrite = TRUE)
  pop_aligned
}

cat("\nBuilding population weight rasters (CMIP6 SSP5-8.5, v2.3 fixed)...\n")
pop_weights <- list()
for (yr in YEARS) {
  pop_weights[[as.character(yr)]] <- build_pop_weight_cmip6(yr)
  cat(sprintf("  %d: %s\n", yr,
              ifelse(is.null(pop_weights[[as.character(yr)]]), "MISSING", "OK")))
}

avail_pop_yrs <- names(pop_weights)[!sapply(pop_weights, is.null)]
if (length(avail_pop_yrs) > 0) {
  for (yr in YEARS) {
    yr_c <- as.character(yr)
    if (is.null(pop_weights[[yr_c]])) {
      nearest <- avail_pop_yrs[which.min(abs(as.integer(avail_pop_yrs) - yr))]
      pop_weights[[yr_c]] <- pop_weights[[nearest]]
      cat(sprintf("  Pop: no %d data, using %s\n", yr, nearest))
    }
  }
}

rm(pop_total_cmip6); gc()
cat("Population weight rasters ready.\n\n")


# ================================================================
#  PART 5: INVENTORY NC FILES
# ================================================================

cat("=== Inventorying NC files ===\n")

tc_dirs <- list.dirs(DIR_NC, recursive = FALSE)
tc_dirs <- tc_dirs[!grepl("pop_cache", basename(tc_dirs))]

nc_inventory <- do.call(rbind, lapply(tc_dirs, function(td) {
  tc_id <- basename(td)
  all_nc <- list.files(td, pattern = "\\.nc$", full.names = TRUE)
  if (length(all_nc) == 0) return(NULL)
  yr <- as.integer(substr(tc_id, 1, 4))
  parts <- strsplit(tc_id, "_")[[1]]
  tc_name <- if (length(parts) >= 3) parts[3] else "UNKNOWN"
  is_holland <- grepl("^WF_Holland_", basename(all_nc))
  data.frame(tc_id = tc_id, tc_year = yr, tc_name = tc_name,
             nc_path = all_nc, is_holland = is_holland,
             stringsAsFactors = FALSE)
}))

nc_inventory$ts_str <- sub(".*_(\\d{8}_\\d{4})\\.nc$", "\\1",
                           basename(nc_inventory$nc_path))
nc_inventory$time <- as.POSIXct(nc_inventory$ts_str,
                                format = "%Y%m%d_%H%M", tz = "UTC")

holland_inv <- nc_inventory %>% filter(is_holland)
blend_inv   <- nc_inventory %>% filter(!is_holland)

cat(sprintf("  Holland NC: %d | Blended NC: %d | TCs: %d\n",
            nrow(holland_inv), nrow(blend_inv),
            n_distinct(nc_inventory$tc_id)))

if (nrow(blend_inv) == 0) {
  cat("  >>> Future mode: NO blended (ERA5) files detected.\n")
  cat("      Holland wind will be used directly for exposure.\n")
}


# ================================================================
#  PART 6: LOAD TRACK CSV
# ================================================================

cat("\n=== Loading HighResMIP Track CSV ===\n")

csv_files <- list.files(DIR_TRACK, pattern = "\\.csv$", full.names = TRUE)
world_map <- map_data("world")

world_sf <- st_as_sf(maps::map("world", fill = TRUE, plot = FALSE))
st_crs(world_sf) <- 4326

ibtracs_tracks <- do.call(rbind, lapply(csv_files, function(f) {
  tc <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(tc)) return(NULL)
  names(tc) <- toupper(names(tc))
  tc$ISO_TIME <- as.POSIXct(tc$ISO_TIME, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  tc <- tc %>% filter(!is.na(LON), !is.na(LAT), !is.na(ISO_TIME)) %>%
    arrange(ISO_TIME) %>% distinct(ISO_TIME, .keep_all = TRUE)
  tc$LON <- ifelse(tc$LON > 180, tc$LON - 360, tc$LON)
  if (!"WIND" %in% names(tc)) return(NULL)
  tc$VMAX_MS <- as.numeric(tc$WIND) * KT2MS
  tc_year <- year(min(tc$ISO_TIME, na.rm = TRUE))
  tc_name <- "UNKNOWN"
  if ("NAME" %in% names(tc) && !is.na(tc$NAME[1]) && nchar(tc$NAME[1]) > 0)
    tc_name <- gsub("[^A-Za-z0-9]", "", tc$NAME[1])
  tc_sid <- ifelse("SID" %in% names(tc) && !is.na(tc$SID[1]),
                   as.character(tc$SID[1]), "")
  date_s <- format(min(tc$ISO_TIME), "%Y%m%d")
  date_e <- format(max(tc$ISO_TIME), "%Y%m%d")
  
  if (nchar(tc_name) > 15) {
    t_num <- regmatches(tc_sid, regexpr("T\\d+", tc_sid))
    if (length(t_num) > 0) {
      tc_name <- gsub("^T", "TC", t_num[1])
    } else {
      tc_name <- "TC0000"
    }
  }
  if (nchar(tc_sid) > 25) {
    t_num <- regmatches(tc_sid, regexpr("T\\d+", tc_sid))
    if (length(t_num) > 0) {
      tc_sid <- t_num[1]
    } else {
      tc_sid <- substr(tc_sid, 1, 25)
    }
  }
  
  tc_id  <- gsub("__+", "_", paste(tc_year, tc_sid, tc_name,
                                   date_s, date_e, sep = "_"))
  
  data.frame(tc_id = tc_id, tc_name = tc_name, tc_year = tc_year,
             lon = tc$LON, lat = tc$LAT, vmax_ms = tc$VMAX_MS,
             time = tc$ISO_TIME, stringsAsFactors = FALSE)
}))

ibtracs_tracks$sshs <- sshs_classify(ibtracs_tracks$vmax_ms)
ibtracs_tracks <- ibtracs_tracks %>% arrange(tc_id, time)
ibtracs_tracks$segment <- ibtracs_tracks$tc_id
for (tid_seg in unique(ibtracs_tracks$tc_id)) {
  idx_seg <- which(ibtracs_tracks$tc_id == tid_seg)
  if (length(idx_seg) < 2) next
  dl <- c(0, abs(diff(ibtracs_tracks$lon[idx_seg])))
  da <- c(0, abs(diff(ibtracs_tracks$lat[idx_seg])))
  ibtracs_tracks$segment[idx_seg] <- paste0(tid_seg, "_s",
                                            cumsum(dl > 10 | da > 10))
}

cat(sprintf("  Loaded %d track points for %d TCs\n",
            nrow(ibtracs_tracks), n_distinct(ibtracs_tracks$tc_id)))


# ================================================================
#  PART 6.5: PRE-COMPUTE SPATIAL INDICES  (SPEED)
#           ★ FIX v2.3: safe centroids — guaranteed 1 row per feature
# ================================================================

cat("\n=== Pre-computing spatial indices for acceleration ===\n")

# ★ FIX v2.3: st_coordinates(st_centroid()) silently DROPS rows for empty /
#   invalid geometries, which would misalign the spatial pre-filter indices
#   with the sf rows and select the WRONG counties/cities. This function
#   guarantees exactly one coordinate row per sf feature (bbox-centre fallback).
get_centroids_safe <- function(sf_obj) {
  n  <- nrow(sf_obj)
  xy <- matrix(NA_real_, nrow = n, ncol = 2)
  cent <- suppressWarnings(tryCatch(
    st_centroid(st_geometry(sf_obj), of_largest_polygon = TRUE),
    error = function(e) NULL))
  if (!is.null(cent) && length(cent) == n) {
    ok <- !st_is_empty(cent)
    if (any(ok)) xy[ok, ] <- st_coordinates(cent[ok])
  }
  bad <- is.na(xy[, 1]) | is.na(xy[, 2])
  if (any(bad)) {
    cat(sprintf("  NOTE: %d feature(s) had no valid centroid; using bbox centre\n",
                sum(bad)))
    for (i in which(bad)) {
      bb <- suppressWarnings(st_bbox(st_geometry(sf_obj)[i]))
      xy[i, ] <- c(mean(bb[c("xmin", "xmax")]), mean(bb[c("ymin", "ymax")]))
    }
  }
  xy
}

cnty_cent <- get_centroids_safe(county_sf)
city_cent <- get_centroids_safe(city_sf)

# ★ FIX v2.3: hard guarantee of row alignment
stopifnot(nrow(cnty_cent) == nrow(county_sf),
          nrow(city_cent) == nrow(city_sf))

SPATIAL_FILTER_KM <- 2000

cat(sprintf("  County centroids: %d | City centroids: %d\n",
            nrow(cnty_cent), nrow(city_cent)))
cat(sprintf("  Spatial filter radius: %d km\n", SPATIAL_FILTER_KM))


# ================================================================
#  PART 7: MAIN PROCESSING LOOP  (OPTIMIZED + v2.3 FIXED)
# ================================================================

cat("\n=== Processing all TCs (v2.3 accelerated, coastal-exposure fixed) ===\n")
step2_start <- Sys.time()

tc_ids_all <- unique(holland_inv$tc_id)
n_tcs      <- length(tc_ids_all)

grid_freq_all   <- matrix(0L, nlon, nlat)
grid_freq_tdts  <- matrix(0L, nlon, nlat)
grid_freq_cat12 <- matrix(0L, nlon, nlat)
grid_freq_cat35 <- matrix(0L, nlon, nlat)
grid_max_ws_bld <- matrix(0,  nlon, nlat)

grid_freq_td   <- matrix(0L, nlon, nlat)
grid_freq_ts   <- matrix(0L, nlon, nlat)
grid_freq_cat1 <- matrix(0L, nlon, nlat)
grid_freq_cat2 <- matrix(0L, nlon, nlat)
grid_freq_cat3 <- matrix(0L, nlon, nlat)
grid_freq_cat4 <- matrix(0L, nlon, nlat)
grid_freq_cat5 <- matrix(0L, nlon, nlat)

tc_meta_list    <- vector("list", n_tcs)
tc_centers_all  <- vector("list", n_tcs)

county_raw_list <- vector("list", 100000L)
city_raw_list   <- vector("list", 100000L)
cnty_counter    <- 0L
city_counter    <- 0L

# split FIRST (correct grouping), THEN sort within each group (v2.2 fix kept)
hol_by_tc <- split(holland_inv, holland_inv$tc_id)
hol_by_tc <- lapply(hol_by_tc, function(x) x[order(x$time), ])

if (nrow(blend_inv) > 0) {
  bld_by_tc <- split(blend_inv, blend_inv$tc_id)
  bld_by_tc <- lapply(bld_by_tc, function(x) x[order(x$time), ])
} else {
  bld_by_tc <- list()
}

# ★ FIX v2.3: 18 deg lon at 50N is only ~1420 km, less than the 1500 km
#   WIND_SPEED_MAX_RADIUS used in Step 1 — parts of the wind field could be
#   cropped away. 21 deg safely covers 1500 km everywhere in the domain.
CROP_BUF_DEG <- 21

cat(sprintf("  TCs to process: %d | Optimizations: spatial filter + crop + safe extract\n",
            n_tcs))

for (ti in seq_along(tc_ids_all)) {
  
  tid <- tc_ids_all[ti]
  
  hol_f <- hol_by_tc[[tid]]
  if (is.null(hol_f) || nrow(hol_f) == 0) next
  
  bld_f <- bld_by_tc[[tid]]
  if (is.null(bld_f)) bld_f <- blend_inv[0, ]
  
  tc_yr <- hol_f$tc_year[1]
  tc_nm <- hol_f$tc_name[1]
  pw    <- pop_weights[[as.character(tc_yr)]]
  
  elapsed <- as.numeric(difftime(Sys.time(), step2_start, units = "mins"))
  if (ti > 1) {
    avg_per_tc <- elapsed / (ti - 1)
    remaining  <- avg_per_tc * (n_tcs - ti + 1)
  } else {
    remaining <- 0
  }
  cat(sprintf("  [%d/%d] %s (%s) | elapsed %.1f min | ETA %.1f min\n",
              ti, n_tcs, tid, tc_nm, elapsed, remaining))
  
  tc_max_hol  <- matrix(0, nlon, nlat)
  tc_max_bld  <- matrix(0, nlon, nlat)
  tc_affected <- matrix(FALSE, nlon, nlat)
  
  n_steps_tc  <- nrow(hol_f)
  tc_pts_list <- vector("list", n_steps_tc)
  
  affects_china_flag <- FALSE
  
  for (fi in seq_len(n_steps_tc)) {
    ncf_h <- tryCatch(nc_open(hol_f$nc_path[fi]), error = function(e) NULL)
    if (is.null(ncf_h)) next
    hol_ws_raw <- nc_get_var(ncf_h, c("holland_wind_speed", "wind_speed"))
    tc_aff_raw <- nc_get_var(ncf_h, c("tc_affected", "tc_mask"))
    tc_lon_i <- ncatt_get(ncf_h, 0, "tc_lon")$value
    tc_lat_i <- ncatt_get(ncf_h, 0, "tc_lat")$value
    vmax_i   <- ncatt_get(ncf_h, 0, "vmax_ms")$value
    time_i   <- ncatt_get(ncf_h, 0, "time")$value
    nc_close(ncf_h)
    
    hol_ws <- ensure_matrix(hol_ws_raw, nlon, nlat)
    if (is.null(hol_ws)) next
    
    tc_aff <- ensure_matrix(tc_aff_raw, nlon, nlat)
    if (is.null(tc_aff))
      tc_aff <- ifelse(!is.na(hol_ws) & hol_ws >= 0.5, 1, 0)
    
    aff_mask <- (tc_aff >= 0.5)
    tc_max_hol  <- pmax(tc_max_hol, hol_ws * aff_mask)
    tc_affected <- tc_affected | aff_mask
    
    tc_pts_list[[fi]] <- data.frame(
      tc_id = tid, tc_name = tc_nm, tc_year = tc_yr,
      tc_lon = tc_lon_i, tc_lat = tc_lat_i, vmax_ms = vmax_i,
      time = as.POSIXct(time_i, tz = "UTC"),
      stringsAsFactors = FALSE)
    
    # ---- Determine exposure wind speed ----
    exp_ws <- NULL
    if (nrow(bld_f) > 0) {
      bld_idx <- which(bld_f$ts_str == hol_f$ts_str[fi])
      if (length(bld_idx) > 0) {
        ncf_b <- tryCatch(nc_open(bld_f$nc_path[bld_idx[1]]),
                          error = function(e) NULL)
        if (!is.null(ncf_b)) {
          bld_ws_raw <- nc_get_var(ncf_b, c("wind_speed"))
          nc_close(ncf_b)
          exp_ws <- ensure_matrix(bld_ws_raw, nlon, nlat)
        }
      }
    }
    if (is.null(exp_ws)) exp_ws <- hol_ws
    
    masked_exp <- exp_ws * aff_mask
    tc_max_bld <- pmax(tc_max_bld, masked_exp)
    
    # SPEED: Skip extraction if wind field is essentially zero
    if (max(masked_exp, na.rm = TRUE) < 0.1) next
    
    # SPEED: Spatial pre-filtering (safe centroids guarantee alignment)
    cnty_dist_km <- haversine_dist(tc_lon_i, tc_lat_i,
                                   cnty_cent[, 1], cnty_cent[, 2]) / 1000
    city_dist_km <- haversine_dist(tc_lon_i, tc_lat_i,
                                   city_cent[, 1], city_cent[, 2]) / 1000
    cnty_near <- which(cnty_dist_km < SPATIAL_FILTER_KM)
    city_near <- which(city_dist_km < SPATIAL_FILTER_KM)
    
    if (length(cnty_near) == 0 && length(city_near) == 0) next
    
    # SPEED: Create raster and crop to TC vicinity
    exp_r <- mat_to_raster(masked_exp, grid_lon, grid_lat)
    
    crop_ext <- extent(
      max(tc_lon_i - CROP_BUF_DEG, DOMAIN_LON[1] - HALF_RES),
      min(tc_lon_i + CROP_BUF_DEG, DOMAIN_LON[2] + HALF_RES),
      max(tc_lat_i - CROP_BUF_DEG, DOMAIN_LAT[1] - HALF_RES),
      min(tc_lat_i + CROP_BUF_DEG, DOMAIN_LAT[2] + HALF_RES))
    exp_r_crop <- tryCatch(crop(exp_r, crop_ext), error = function(e) exp_r)
    
    # ★ FIX v2.3: crop pw to EXACTLY the extent of exp_r_crop and verify
    #   alignment; on any mismatch fall back to NULL (-> unweighted mean),
    #   never to a silently misaligned weight raster.
    pw_crop <- NULL
    if (!is.null(pw)) {
      pw_crop <- tryCatch(crop(pw, extent(exp_r_crop)),
                          error = function(e) NULL)
      if (!is.null(pw_crop)) {
        same_grid <- tryCatch(
          compareRaster(exp_r_crop, pw_crop, extent = TRUE, rowcol = TRUE,
                        crs = FALSE, stopiffalse = FALSE),
          error = function(e) FALSE)
        if (!isTRUE(same_grid)) {
          pw_crop <- tryCatch(resample(pw, exp_r_crop, method = "ngb"),
                              error = function(e) NULL)
        }
      }
    }
    
    time_b <- hol_f$time[fi]
    date_b <- as.Date(time_b)
    
    # ---- County extraction (nearby only) ----
    # ★ FIX v2.3: separate exact_extract calls (identical semantics to the
    #   validated v2.1/History pipeline) + row-wise NA/NaN fallback so that
    #   coastal counties with real wind can NEVER be silently dropped just
    #   because the weight raster is NA/0 inside the cropped window.
    if (length(cnty_near) > 0) {
      county_sub <- county_sf[cnty_near, ]
      
      cnty_ws_max  <- exact_extract(exp_r_crop, county_sub, fun = "max")
      cnty_ws_mean <- exact_extract(exp_r_crop, county_sub, fun = "mean")
      
      if (!is.null(pw_crop)) {
        cnty_ws_pw <- tryCatch(
          exact_extract(exp_r_crop, county_sub, fun = "weighted_mean",
                        weights = pw_crop),
          error = function(e) rep(NA_real_, nrow(county_sub)))
      } else {
        cnty_ws_pw <- cnty_ws_mean
      }
      
      # ★ FIX v2.3 (the critical line): weighted_mean = NA/NaN but the
      #   county DOES have valid wind -> use the unweighted area mean.
      bad_pw <- !is.finite(cnty_ws_pw) & is.finite(cnty_ws_mean)
      if (any(bad_pw)) cnty_ws_pw[bad_pw] <- cnty_ws_mean[bad_pw]
      
      keep <- !is.na(cnty_ws_pw) & cnty_ws_pw > 0.5
      if (sum(keep) > 0) {
        affects_china_flag <- TRUE
        ki <- cnty_near[keep]
        cnty_counter <- cnty_counter + 1L
        if (cnty_counter > length(county_raw_list))
          county_raw_list <- c(county_raw_list, vector("list", 50000L))
        county_raw_list[[cnty_counter]] <- data.frame(
          tc_id = tid, tc_name = tc_nm, tc_year = tc_yr,
          time = time_b, date = date_b,
          tc_lon = tc_lon_i, tc_lat = tc_lat_i, tc_vmax = vmax_i,
          province = county_sf$省[ki], city = county_sf$市[ki],
          county = county_sf$县[ki],
          province_code = county_sf$省代码[ki],
          city_code     = county_sf$市代码[ki],
          county_code   = county_sf$县代码[ki],
          ws_popweighted = round(cnty_ws_pw[keep], 2),
          ws_max         = round(cnty_ws_max[keep], 2),
          ws_area_mean   = round(cnty_ws_mean[keep], 2),
          dist_km        = round(cnty_dist_km[ki], 1),
          stringsAsFactors = FALSE)
      }
    }
    
    # ---- City extraction (nearby only) ----
    if (length(city_near) > 0) {
      city_sub <- city_sf[city_near, ]
      
      city_ws_max  <- exact_extract(exp_r_crop, city_sub, fun = "max")
      city_ws_mean <- exact_extract(exp_r_crop, city_sub, fun = "mean")
      
      if (!is.null(pw_crop)) {
        city_ws_pw <- tryCatch(
          exact_extract(exp_r_crop, city_sub, fun = "weighted_mean",
                        weights = pw_crop),
          error = function(e) rep(NA_real_, nrow(city_sub)))
      } else {
        city_ws_pw <- city_ws_mean
      }
      
      # ★ FIX v2.3: same NA/NaN safety net for cities
      bad_pw_c <- !is.finite(city_ws_pw) & is.finite(city_ws_mean)
      if (any(bad_pw_c)) city_ws_pw[bad_pw_c] <- city_ws_mean[bad_pw_c]
      
      keep_c <- !is.na(city_ws_pw) & city_ws_pw > 0.5
      if (sum(keep_c) > 0) {
        affects_china_flag <- TRUE
        kic <- city_near[keep_c]
        city_counter <- city_counter + 1L
        if (city_counter > length(city_raw_list))
          city_raw_list <- c(city_raw_list, vector("list", 50000L))
        city_raw_list[[city_counter]] <- data.frame(
          tc_id = tid, tc_name = tc_nm, tc_year = tc_yr,
          time = time_b, date = date_b,
          tc_lon = tc_lon_i, tc_lat = tc_lat_i, tc_vmax = vmax_i,
          province = city_sf$省[kic], city = city_sf$市[kic],
          province_code = city_sf$省代码[kic],
          city_code     = city_sf$市代码[kic],
          ws_popweighted = round(city_ws_pw[keep_c], 2),
          ws_max         = round(city_ws_max[keep_c], 2),
          ws_area_mean   = round(city_ws_mean[keep_c], 2),
          dist_km        = round(city_dist_km[kic], 1),
          stringsAsFactors = FALSE)
      }
    }
  }  # end inner loop (timesteps)
  
  tc_pts <- do.call(rbind, tc_pts_list[!sapply(tc_pts_list, is.null)])
  if (is.null(tc_pts)) tc_pts <- data.frame()
  
  # Grid accumulators
  if (any(tc_affected)) {
    tc_aff_int <- (tc_affected * 1L)
    grid_freq_all   <- grid_freq_all + tc_aff_int
    grid_max_ws_bld <- pmax(grid_max_ws_bld, tc_max_bld)
    hv <- tc_max_hol
    grid_freq_tdts  <- grid_freq_tdts  + ((tc_affected & hv < 33.0) * 1L)
    grid_freq_cat12 <- grid_freq_cat12 + ((tc_affected & hv >= 33.0 & hv < 50.0) * 1L)
    grid_freq_cat35 <- grid_freq_cat35 + ((tc_affected & hv >= 50.0) * 1L)
    grid_freq_td   <- grid_freq_td   + ((tc_affected & hv < 17.5) * 1L)
    grid_freq_ts   <- grid_freq_ts   + ((tc_affected & hv >= 17.5 & hv < 33.0) * 1L)
    grid_freq_cat1 <- grid_freq_cat1 + ((tc_affected & hv >= 33.0 & hv < 43.0) * 1L)
    grid_freq_cat2 <- grid_freq_cat2 + ((tc_affected & hv >= 43.0 & hv < 50.0) * 1L)
    grid_freq_cat3 <- grid_freq_cat3 + ((tc_affected & hv >= 50.0 & hv < 58.0) * 1L)
    grid_freq_cat4 <- grid_freq_cat4 + ((tc_affected & hv >= 58.0 & hv < 70.0) * 1L)
    grid_freq_cat5 <- grid_freq_cat5 + ((tc_affected & hv >= 70.0) * 1L)
  }
  
  # Landfalling detection
  landfalling <- FALSE; landfall_time <- NA_character_
  
  if (!is.null(china_land_buf) && nrow(tc_pts) > 0) {
    lf_nc <- check_landfalling(tc_pts$tc_lon, tc_pts$tc_lat, china_land_buf)
    if (lf_nc$lf) {
      landfalling <- TRUE
      landfall_time <- as.character(tc_pts$time[lf_nc$idx])
    }
  }
  
  if (!landfalling && !is.null(china_land_buf)) {
    ibt_sub <- ibtracs_tracks[ibtracs_tracks$tc_id == tid, ]
    if (nrow(ibt_sub) > 0) {
      lf_ibt <- check_landfalling(ibt_sub$lon, ibt_sub$lat, china_land_buf)
      if (lf_ibt$lf) {
        landfalling <- TRUE
        landfall_time <- as.character(ibt_sub$time[lf_ibt$idx])
      }
    }
  }
  
  # ★ FIX v2.3: extraction flag (fast path) + ORG-style grid check as backup.
  #   The flag alone requires pop-weighted wind > 0.5 somewhere; the ORG
  #   definition only requires the Holland mask to touch a county. Restore
  #   full ORG semantics when the flag is FALSE.
  affects_china <- affects_china_flag
  if (!affects_china && any(tc_affected)) {
    hol_r <- mat_to_raster(tc_affected * 1, grid_lon, grid_lat)
    ca <- tryCatch(exact_extract(hol_r, county_sf, fun = "max"),
                   error = function(e) rep(0, nrow(county_sf)))
    affects_china <- any(ca >= 0.5, na.rm = TRUE)
  }
  
  peak_vmax <- if (nrow(tc_pts) > 0) max(tc_pts$vmax_ms, na.rm = TRUE) else 0
  peak_sshs <- sshs_classify(peak_vmax)
  
  tc_meta_list[[ti]] <- data.frame(
    tc_id = tid, tc_name = tc_nm, tc_year = tc_yr,
    peak_vmax = round(peak_vmax, 1), peak_sshs = peak_sshs,
    landfalling = landfalling, landfall_time = landfall_time,
    affects_china = affects_china, n_timesteps = nrow(hol_f),
    date_start = if (nrow(tc_pts) > 0) format(min(tc_pts$time), "%Y-%m-%d") else NA,
    date_end   = if (nrow(tc_pts) > 0) format(max(tc_pts$time), "%Y-%m-%d") else NA,
    stringsAsFactors = FALSE)
  
  tc_centers_all[[ti]] <- tc_pts
  
  cat(sprintf("    %s: %d pts | aff=%s | LF=%s\n",
              tid, nrow(tc_pts), affects_china, landfalling))
  
  if (ti %% 25 == 0) {
    gc(verbose = FALSE)
    cat(sprintf("    [GC] Memory cleaned after %d TCs\n", ti))
  }
  
}  # === end TC loop ===

tc_meta    <- do.call(rbind, tc_meta_list[!sapply(tc_meta_list, is.null)])
tc_centers <- do.call(rbind, tc_centers_all[!sapply(tc_centers_all, is.null)])
county_raw <- if (cnty_counter > 0)
  do.call(rbind, county_raw_list[seq_len(cnty_counter)]) else data.frame()
city_raw   <- if (city_counter > 0)
  do.call(rbind, city_raw_list[seq_len(city_counter)]) else data.frame()

rm(county_raw_list, city_raw_list, hol_by_tc, bld_by_tc); gc()

n_lf  <- sum(tc_meta$landfalling, na.rm = TRUE)
n_aff <- sum(!tc_meta$landfalling & tc_meta$affects_china, na.rm = TRUE)
n_tot <- nrow(tc_meta)

cat(sprintf("\n  Total TCs: %d | Landfalling: %d | Non-LF affecting: %d\n",
            n_tot, n_lf, n_aff))

# ★ FIX v2.3: quick sanity check — coastal provinces MUST appear
if (nrow(county_raw) > 0) {
  coastal_check <- c("福建省", "浙江省", "江苏省", "上海市", "广东省", "海南省")
  n_coastal <- sum(county_raw$province %in% coastal_check)
  cat(sprintf("  [CHECK] Coastal-province county-timestep records: %d %s\n",
              n_coastal,
              ifelse(n_coastal > 0, "(OK)", "(!! STILL ZERO — INVESTIGATE !!)")))
}

ibtracs_tracks <- ibtracs_tracks %>%
  left_join(tc_meta %>% dplyr::select(tc_id, landfalling, affects_china,
                                      peak_sshs, peak_vmax),
            by = "tc_id")


# ================================================================
#  PART 8: AGGREGATE TO DAILY EXPOSURE
# ================================================================

cat("\n=== Aggregating to daily ===\n")

agg_daily <- function(raw_df, geo_cols) {
  if (nrow(raw_df) == 0) return(raw_df)
  raw_df %>%
    group_by(across(all_of(c("tc_id","tc_name","tc_year","date", geo_cols)))) %>%
    summarise(
      ws_daily_max_pw  = max(ws_popweighted, na.rm = TRUE),
      ws_daily_mean_pw = mean(ws_popweighted, na.rm = TRUE),
      ws_daily_max     = max(ws_max, na.rm = TRUE),
      ws_daily_mean    = mean(ws_area_mean, na.rm = TRUE),
      dist_min_km      = min(dist_km, na.rm = TRUE),
      tc_vmax_max      = max(tc_vmax, na.rm = TRUE),
      n_timesteps      = n(), .groups = "drop") %>%
    mutate(exposed_ts = ws_daily_max_pw >= WS_THRESHOLD,
           sshs_cat   = sshs_classify(ws_daily_max_pw))
}

county_geo <- c("province","city","county","province_code","city_code","county_code")
city_geo   <- c("province","city","province_code","city_code")

county_daily <- agg_daily(county_raw, county_geo)
city_daily   <- agg_daily(city_raw,   city_geo)

if (nrow(county_daily) > 0)
  county_daily <- county_daily %>%
  left_join(tc_meta %>% dplyr::select(tc_id, landfalling, peak_vmax,
                                      peak_sshs, landfall_time),
            by = "tc_id")
if (nrow(city_daily) > 0)
  city_daily <- city_daily %>%
  left_join(tc_meta %>% dplyr::select(tc_id, landfalling, peak_vmax,
                                      peak_sshs, landfall_time),
            by = "tc_id")

cat(sprintf("  County daily: %d (exposed TS+: %d)\n",
            nrow(county_daily),
            if (nrow(county_daily) > 0) sum(county_daily$exposed_ts) else 0))
cat(sprintf("  City daily:   %d (exposed TS+: %d)\n",
            nrow(city_daily),
            if (nrow(city_daily) > 0) sum(city_daily$exposed_ts) else 0))


# ================================================================
#  PART 9: MERGE AGE-SPECIFIC POPULATION
# ================================================================

cat("\n=== Merging age-specific population ===\n")

merge_pop <- function(df, pop_df, code_col_df, code_col_pop) {
  if (nrow(df) == 0) return(df)
  df$pop_lt20 <- df$pop_20to44 <- df$pop_45to64 <-
    df$pop_ge65 <- df$pop_total <- NA_real_
  if (is.null(pop_df)) return(df)
  df[[code_col_df]]    <- as.character(df[[code_col_df]])
  pop_df[[code_col_pop]] <- as.character(pop_df[[code_col_pop]])
  for (yr in unique(df$tc_year)) {
    idx <- which(df$tc_year == yr)
    if (length(idx) == 0) next
    cols <- c(paste0("<20 years_", yr), paste0("20 to 44 years_", yr),
              paste0("45 to 64 years_", yr), paste0(">64 years_", yr))
    if (!all(cols %in% names(pop_df))) {
      avail_yrs <- as.integer(
        gsub(".*_(\\d{4})$", "\\1",
             grep("^<20 years_\\d{4}$", names(pop_df), value = TRUE)))
      if (length(avail_yrs) > 0) {
        nearest_yr <- avail_yrs[which.min(abs(avail_yrs - yr))]
        cols <- c(paste0("<20 years_", nearest_yr),
                  paste0("20 to 44 years_", nearest_yr),
                  paste0("45 to 64 years_", nearest_yr),
                  paste0(">64 years_", nearest_yr))
        if (!all(cols %in% names(pop_df))) next
      } else next
    }
    m <- match(df[[code_col_df]][idx], pop_df[[code_col_pop]])
    df$pop_lt20[idx]   <- pop_df[[cols[1]]][m]
    df$pop_20to44[idx] <- pop_df[[cols[2]]][m]
    df$pop_45to64[idx] <- pop_df[[cols[3]]][m]
    df$pop_ge65[idx]   <- pop_df[[cols[4]]][m]
  }
  df$pop_total <- rowSums(df[, c("pop_lt20","pop_20to44",
                                 "pop_45to64","pop_ge65")], na.rm = TRUE)
  df$pop_total[df$pop_total == 0] <- NA
  df
}

if (nrow(county_daily) > 0)
  county_daily <- merge_pop(county_daily, pop_county_age,
                            "county_code", "CountyCode")

if (!is.null(pop_city_age) && nrow(city_daily) > 0)
  city_daily <- merge_pop(city_daily, pop_city_age,
                          "city_code", "CityCode")

cat("  Population merge done\n")

if (nrow(city_daily) > 0) {
  city_pop_cov <- mean(!is.na(city_daily$pop_total))
  cat(sprintf("  City pop coverage: %.1f%%\n", city_pop_cov * 100))
}


# ================================================================
#  PART 10: COUNTY/CITY FREQUENCY & MAX WIND FROM GRID
# ================================================================

cat("\n=== County/city freq & max wind from grid ===\n")

freq_r       <- mat_to_raster(grid_freq_all,   grid_lon, grid_lat)
freq_tdts_r  <- mat_to_raster(grid_freq_tdts,  grid_lon, grid_lat)
freq_cat12_r <- mat_to_raster(grid_freq_cat12, grid_lon, grid_lat)
freq_cat35_r <- mat_to_raster(grid_freq_cat35, grid_lon, grid_lat)
maxws_r      <- mat_to_raster(grid_max_ws_bld, grid_lon, grid_lat)

freq_td_r   <- mat_to_raster(grid_freq_td,   grid_lon, grid_lat)
freq_ts_r   <- mat_to_raster(grid_freq_ts,   grid_lon, grid_lat)
freq_cat1_r <- mat_to_raster(grid_freq_cat1, grid_lon, grid_lat)
freq_cat2_r <- mat_to_raster(grid_freq_cat2, grid_lon, grid_lat)
freq_cat3_r <- mat_to_raster(grid_freq_cat3, grid_lon, grid_lat)
freq_cat4_r <- mat_to_raster(grid_freq_cat4, grid_lon, grid_lat)
freq_cat5_r <- mat_to_raster(grid_freq_cat5, grid_lon, grid_lat)

ee <- function(ras, sf_obj, fun = "max")
  exact_extract(ras, sf_obj, fun = fun)

county_grid_stats <- data.frame(
  county_code = county_sf$县代码,
  province = county_sf$省, city = county_sf$市, county = county_sf$县,
  freq_all = ee(freq_r, county_sf), freq_tdts = ee(freq_tdts_r, county_sf),
  freq_cat12 = ee(freq_cat12_r, county_sf),
  freq_cat35 = ee(freq_cat35_r, county_sf),
  freq_td   = ee(freq_td_r, county_sf),
  freq_ts   = ee(freq_ts_r, county_sf),
  freq_cat1 = ee(freq_cat1_r, county_sf),
  freq_cat2 = ee(freq_cat2_r, county_sf),
  freq_cat3 = ee(freq_cat3_r, county_sf),
  freq_cat4 = ee(freq_cat4_r, county_sf),
  freq_cat5 = ee(freq_cat5_r, county_sf),
  max_ws = round(ee(maxws_r, county_sf), 2),
  mean_ws = round(ee(maxws_r, county_sf, "mean"), 2),
  stringsAsFactors = FALSE)
county_grid_stats$ann_freq_all   <- county_grid_stats$freq_all   / N_YEARS
county_grid_stats$ann_freq_tdts  <- county_grid_stats$freq_tdts  / N_YEARS
county_grid_stats$ann_freq_cat12 <- county_grid_stats$freq_cat12 / N_YEARS
county_grid_stats$ann_freq_cat35 <- county_grid_stats$freq_cat35 / N_YEARS
county_grid_stats$ann_freq_td    <- county_grid_stats$freq_td    / N_YEARS
county_grid_stats$ann_freq_ts    <- county_grid_stats$freq_ts    / N_YEARS
county_grid_stats$ann_freq_cat1  <- county_grid_stats$freq_cat1  / N_YEARS
county_grid_stats$ann_freq_cat2  <- county_grid_stats$freq_cat2  / N_YEARS
county_grid_stats$ann_freq_cat3  <- county_grid_stats$freq_cat3  / N_YEARS
county_grid_stats$ann_freq_cat4  <- county_grid_stats$freq_cat4  / N_YEARS
county_grid_stats$ann_freq_cat5  <- county_grid_stats$freq_cat5  / N_YEARS

city_grid_stats <- data.frame(
  city_code = city_sf$市代码,
  province = city_sf$省, city = city_sf$市,
  freq_all = ee(freq_r, city_sf), freq_tdts = ee(freq_tdts_r, city_sf),
  freq_cat12 = ee(freq_cat12_r, city_sf),
  freq_cat35 = ee(freq_cat35_r, city_sf),
  freq_td   = ee(freq_td_r, city_sf),
  freq_ts   = ee(freq_ts_r, city_sf),
  freq_cat1 = ee(freq_cat1_r, city_sf),
  freq_cat2 = ee(freq_cat2_r, city_sf),
  freq_cat3 = ee(freq_cat3_r, city_sf),
  freq_cat4 = ee(freq_cat4_r, city_sf),
  freq_cat5 = ee(freq_cat5_r, city_sf),
  max_ws = round(ee(maxws_r, city_sf), 2),
  mean_ws = round(ee(maxws_r, city_sf, "mean"), 2),
  stringsAsFactors = FALSE)
city_grid_stats$ann_freq_all   <- city_grid_stats$freq_all   / N_YEARS
city_grid_stats$ann_freq_tdts  <- city_grid_stats$freq_tdts  / N_YEARS
city_grid_stats$ann_freq_cat12 <- city_grid_stats$freq_cat12 / N_YEARS
city_grid_stats$ann_freq_cat35 <- city_grid_stats$freq_cat35 / N_YEARS
city_grid_stats$ann_freq_td    <- city_grid_stats$freq_td    / N_YEARS
city_grid_stats$ann_freq_ts    <- city_grid_stats$freq_ts    / N_YEARS
city_grid_stats$ann_freq_cat1  <- city_grid_stats$freq_cat1  / N_YEARS
city_grid_stats$ann_freq_cat2  <- city_grid_stats$freq_cat2  / N_YEARS
city_grid_stats$ann_freq_cat3  <- city_grid_stats$freq_cat3  / N_YEARS
city_grid_stats$ann_freq_cat4  <- city_grid_stats$freq_cat4  / N_YEARS
city_grid_stats$ann_freq_cat5  <- city_grid_stats$freq_cat5  / N_YEARS

save_data(county_grid_stats, "County_Frequency_MaxWS", DIR_TABLE)
save_data(city_grid_stats,   "City_Frequency_MaxWS",   DIR_TABLE)


# ================================================================
#  PART 11: EVENT-LEVEL TABLES
# ================================================================

cat("\n=== Building event-level tables ===\n")

build_event <- function(daily_df, geo_cols) {
  if (nrow(daily_df) == 0) return(data.frame())
  daily_df %>%
    group_by(across(all_of(c(geo_cols, "tc_id","tc_name","tc_year")))) %>%
    summarise(
      event_start = min(date), event_end = max(date),
      max_ws_pw = max(ws_daily_max_pw, na.rm = TRUE),
      mean_ws_pw = mean(ws_daily_mean_pw, na.rm = TRUE),
      max_ws = max(ws_daily_max, na.rm = TRUE),
      mean_ws = mean(ws_daily_mean, na.rm = TRUE),
      dist_min_km = min(dist_min_km, na.rm = TRUE),
      n_days = n(),
      pop_total = first(pop_total), pop_lt20 = first(pop_lt20),
      pop_20to44 = first(pop_20to44), pop_45to64 = first(pop_45to64),
      pop_ge65 = first(pop_ge65), .groups = "drop") %>%
    mutate(event_sshs = sshs_classify(max_ws_pw)) %>%
    left_join(tc_meta %>% dplyr::select(tc_id, peak_sshs, landfalling),
              by = "tc_id")
}

county_event <- build_event(county_daily, county_geo)
city_event   <- build_event(city_daily, city_geo)

if (nrow(county_event) > 0) {
  county_event <- county_event %>%
    group_by(county_code) %>% mutate(event_seq = row_number()) %>%
    ungroup() %>% arrange(county_code, event_start)
  save_data(county_event, "TC_County_Event_Summary", DIR_TABLE)
}
if (nrow(city_event) > 0) {
  city_event <- city_event %>%
    group_by(city_code) %>% mutate(event_seq = row_number()) %>%
    ungroup() %>% arrange(city_code, event_start)
  save_data(city_event, "TC_City_Event_Summary", DIR_TABLE)
}

tc_summary_exp <- data.frame()
if (nrow(county_daily) > 0 && sum(county_daily$exposed_ts) > 0) {
  tc_summary_exp <- county_daily %>%
    filter(exposed_ts) %>%
    group_by(tc_id, tc_name, tc_year, landfalling, peak_sshs) %>%
    summarise(
      n_days = n_distinct(date), n_counties = n_distinct(county_code),
      n_cities = n_distinct(city_code),
      max_ws_pw = max(ws_daily_max_pw, na.rm = TRUE),
      pop_total  = sum(pop_total[!duplicated(paste0(county_code,date))],  na.rm = TRUE),
      pop_lt20   = sum(pop_lt20[!duplicated(paste0(county_code,date))],   na.rm = TRUE),
      pop_20to44 = sum(pop_20to44[!duplicated(paste0(county_code,date))], na.rm = TRUE),
      pop_45to64 = sum(pop_45to64[!duplicated(paste0(county_code,date))], na.rm = TRUE),
      pop_ge65   = sum(pop_ge65[!duplicated(paste0(county_code,date))],   na.rm = TRUE),
      .groups = "drop") %>%
    left_join(tc_meta %>% dplyr::select(tc_id, date_start, date_end,
                                        landfall_time, affects_china),
              by = "tc_id") %>%
    arrange(tc_year, tc_id)
  save_data(tc_summary_exp, "TC_Summary_Exposure_TSplus", DIR_TABLE)
}

tc_summary_all <- data.frame()
if (nrow(county_daily) > 0) {
  tc_summary_all <- county_daily %>%
    group_by(tc_id, tc_name, tc_year, landfalling, peak_sshs) %>%
    summarise(
      n_days = n_distinct(date), n_counties = n_distinct(county_code),
      n_cities = n_distinct(city_code),
      max_ws_pw = max(ws_daily_max_pw, na.rm = TRUE),
      pop_total  = sum(pop_total[!duplicated(paste0(county_code,date))],  na.rm = TRUE),
      pop_lt20   = sum(pop_lt20[!duplicated(paste0(county_code,date))],   na.rm = TRUE),
      pop_20to44 = sum(pop_20to44[!duplicated(paste0(county_code,date))], na.rm = TRUE),
      pop_45to64 = sum(pop_45to64[!duplicated(paste0(county_code,date))], na.rm = TRUE),
      pop_ge65   = sum(pop_ge65[!duplicated(paste0(county_code,date))],   na.rm = TRUE),
      .groups = "drop") %>%
    left_join(tc_meta %>% dplyr::select(tc_id, date_start, date_end,
                                        landfall_time, affects_china),
              by = "tc_id") %>%
    arrange(tc_year, tc_id)
  save_data(tc_summary_all, "TC_Summary_Exposure_All", DIR_TABLE)
}


# ================================================================
#  PART 12: SAVE CORE TABLES
# ================================================================

cat("\n=== Saving tables ===\n")

save_data(county_daily, "TC_County_Daily_Exposure", DIR_TABLE)
save_data(city_daily,   "TC_City_Daily_Exposure",   DIR_TABLE)
save_data(tc_meta,      "TC_Metadata",              DIR_TABLE)

annual_all <- data.frame()
if (nrow(county_daily) > 0 && sum(county_daily$exposed_ts) > 0) {
  annual_all <- county_daily %>%
    filter(exposed_ts) %>%
    group_by(tc_year) %>%
    summarise(
      n_tcs = n_distinct(tc_id),
      n_tcs_lf = n_distinct(tc_id[landfalling == TRUE]),
      n_counties = n_distinct(county_code),
      n_cities   = n_distinct(city_code),
      pop_total = sum(pop_total[!duplicated(paste0(county_code,tc_id,date))],
                      na.rm = TRUE),
      pop_lt20  = sum(pop_lt20[!duplicated(paste0(county_code,tc_id,date))],
                      na.rm = TRUE),
      pop_ge65  = sum(pop_ge65[!duplicated(paste0(county_code,tc_id,date))],
                      na.rm = TRUE),
      .groups = "drop")
  save_data(annual_all, "TC_Annual_Summary_TSplus", DIR_TABLE)
}

prov_summary <- data.frame()
if (nrow(county_daily) > 0 && sum(county_daily$exposed_ts) > 0) {
  prov_summary <- county_daily %>%
    filter(exposed_ts, !province %in% EXCLUDE_PROV) %>%
    group_by(province, province_code) %>%
    summarise(
      n_tc_events = n_distinct(tc_id), n_counties = n_distinct(county_code),
      total_days = n(),
      pop_total  = sum(pop_total[!duplicated(paste0(county_code,tc_id,date))],
                       na.rm = TRUE),
      pop_lt20   = sum(pop_lt20[!duplicated(paste0(county_code,tc_id,date))],
                       na.rm = TRUE),
      pop_20to44 = sum(pop_20to44[!duplicated(paste0(county_code,tc_id,date))],
                       na.rm = TRUE),
      pop_45to64 = sum(pop_45to64[!duplicated(paste0(county_code,tc_id,date))],
                       na.rm = TRUE),
      pop_ge65   = sum(pop_ge65[!duplicated(paste0(county_code,tc_id,date))],
                       na.rm = TRUE),
      mean_ws = mean(ws_daily_max_pw, na.rm = TRUE),
      max_ws  = max(ws_daily_max_pw, na.rm = TRUE),
      .groups = "drop") %>%
    arrange(desc(pop_total))
  save_data(prov_summary, "TC_Provincial_Summary", DIR_TABLE)
}
cat("  All tables saved\n")


# ================================================================
#  PART 13: MAP01 - TC TRACK OVERVIEW
# ================================================================

cat("\n=== Map01: TC Track Overview ===\n")

track_segs <- make_track_segs(ibtracs_tracks)
if (!is.null(track_segs)) {
  track_segs$sshs <- factor(track_segs$sshs, levels = SSHS_LABELS)
  track_segs <- track_segs %>%
    left_join(tc_meta %>% dplyr::select(tc_id, landfalling, affects_china,
                                        peak_sshs, peak_vmax), by = "tc_id")
}

if (!is.null(track_segs) && nrow(track_segs) > 0) {
  
  track_xlim <- range(c(track_segs$x, track_segs$xend), na.rm = TRUE) + c(-3, 3)
  track_ylim <- range(c(track_segs$y, track_segs$yend), na.rm = TRUE) + c(-2, 2)
  track_xlim <- c(max(track_xlim[1], 95),  min(track_xlim[2], 180))
  track_ylim <- c(max(track_ylim[1], -5),  min(track_ylim[2], 55))
  
  ax_brk <- seq(10 * floor(track_xlim[1]/10), 10 * ceiling(track_xlim[2]/10), by = 10)
  ay_brk <- seq(10 * floor(track_ylim[1]/10), 10 * ceiling(track_ylim[2]/10), by = 10)
  
  label_pts <- track_segs %>%
    filter(peak_sshs %in% c("Category 4","Category 5")) %>%
    group_by(tc_id) %>% slice(which.min(y)) %>% ungroup() %>%
    mutate(label = paste0(tc_name, " (", tc_year, ")"))
  
  tryCatch({
    present_sshs <- levels(track_segs$sshs)[levels(track_segs$sshs) %in%
                                              unique(as.character(track_segs$sshs))]
    
    p <- ggplot() +
      geom_polygon(data = world_map, aes(long, lat, group = group),
                   fill = "grey92", colour = "grey70", linewidth = 0.15) +
      geom_segment(data = track_segs,
                   aes(x = x, y = y, xend = xend, yend = yend,
                       colour = sshs), linewidth = 0.6, alpha = 0.7) +
      scale_colour_manual(values = SSHS_COLORS[present_sshs],
                          name = "SSHS", drop = TRUE)
    if (nrow(label_pts) > 0)
      p <- p + geom_text_repel(data = label_pts,
                               aes(x = x, y = y, label = label),
                               size = 2.2, fontface = "bold",
                               max.overlaps = 20, segment.size = 0.3,
                               box.padding = 0.4, seed = 42)
    p <- p + coord_sf(xlim = track_xlim, ylim = track_ylim,
                      datum = st_crs(4326)) +
      scale_x_continuous(breaks = ax_brk, labels = fmt_lon) +
      scale_y_continuous(breaks = ay_brk, labels = fmt_lat) +
      labs(title = paste0("TC Track Overview (", min(YEARS), "-", max(YEARS),
                          ") \u2014 SSP5-8.5 Projection"),
           subtitle = sprintf("Total: %d | Landfalling: %d | Non-LF affecting: %d",
                              n_tot, n_lf, n_aff),
           x = "Longitude", y = "Latitude") +
      theme_pub()
    save_both(p, file.path(DIR_MAP, "Map01a_AllTracks_SSHS"), w = 12, h = 7)
    save_data(track_segs %>% dplyr::select(tc_id, tc_name, tc_year,
                                           x, y, xend, yend, sshs),
              "Map01a_AllTracks_SSHS")
    cat("  Map01a saved\n")
    
    lf_segs <- track_segs %>% filter(landfalling == TRUE)
    if (nrow(lf_segs) > 0) {
      lf_label <- label_pts %>% filter(tc_id %in% lf_segs$tc_id)
      ps_lf <- levels(lf_segs$sshs)[levels(lf_segs$sshs) %in%
                                      unique(as.character(lf_segs$sshs))]
      p2 <- ggplot() +
        geom_polygon(data = world_map, aes(long, lat, group = group),
                     fill = "grey92", colour = "grey70", linewidth = 0.15) +
        geom_segment(data = lf_segs,
                     aes(x = x, y = y, xend = xend, yend = yend,
                         colour = sshs), linewidth = 0.7, alpha = 0.7) +
        scale_colour_manual(values = SSHS_COLORS[ps_lf], name = "SSHS",
                            drop = TRUE)
      if (nrow(lf_label) > 0)
        p2 <- p2 + geom_text_repel(data = lf_label,
                                   aes(x = x, y = y, label = label),
                                   size = 2.2, fontface = "bold",
                                   max.overlaps = 20, segment.size = 0.3,
                                   seed = 42)
      p2 <- p2 + coord_sf(xlim = track_xlim, ylim = track_ylim,
                          datum = st_crs(4326)) +
        scale_x_continuous(breaks = ax_brk, labels = fmt_lon) +
        scale_y_continuous(breaks = ay_brk, labels = fmt_lat) +
        labs(title = "Landfalling TC Tracks (SSP5-8.5)",
             subtitle = sprintf("N = %d", n_lf),
             x = "Longitude", y = "Latitude") + theme_pub()
      save_both(p2, file.path(DIR_MAP, "Map01b_Landfalling"), w = 12, h = 7)
      cat("  Map01b saved\n")
    }
    
    nf_segs <- track_segs %>% filter(landfalling == FALSE, affects_china == TRUE)
    if (nrow(nf_segs) > 0) {
      ps_nf <- levels(nf_segs$sshs)[levels(nf_segs$sshs) %in%
                                      unique(as.character(nf_segs$sshs))]
      p3 <- ggplot() +
        geom_polygon(data = world_map, aes(long, lat, group = group),
                     fill = "grey92", colour = "grey70", linewidth = 0.15) +
        geom_segment(data = nf_segs,
                     aes(x = x, y = y, xend = xend, yend = yend,
                         colour = sshs), linewidth = 0.7, alpha = 0.7) +
        scale_colour_manual(values = SSHS_COLORS[ps_nf], name = "SSHS",
                            drop = TRUE) +
        coord_sf(xlim = track_xlim, ylim = track_ylim,
                 datum = st_crs(4326)) +
        scale_x_continuous(breaks = ax_brk, labels = fmt_lon) +
        scale_y_continuous(breaks = ay_brk, labels = fmt_lat) +
        labs(title = "Non-landfalling TCs Affecting China (SSP5-8.5)",
             subtitle = sprintf("N = %d", n_aff),
             x = "Longitude", y = "Latitude") + theme_pub()
      save_both(p3, file.path(DIR_MAP, "Map01c_NonLF_Affecting"), w = 12, h = 7)
      cat("  Map01c saved\n")
    }
    
    ts_sub <- track_segs %>% filter(landfalling | affects_china)
    if (nrow(ts_sub) > 0) {
      ts_sub$lf_label <- ifelse(ts_sub$landfalling,
                                "Landfalling", "Non-LF affecting")
      p4 <- ggplot() +
        geom_polygon(data = world_map, aes(long, lat, group = group),
                     fill = "grey92", colour = "grey70", linewidth = 0.15) +
        geom_segment(data = ts_sub,
                     aes(x = x, y = y, xend = xend, yend = yend,
                         colour = lf_label), linewidth = 0.5, alpha = 0.6) +
        scale_colour_manual(values = c("Landfalling"="#D7191C",
                                       "Non-LF affecting"="#2C7BB6"),
                            name = NULL) +
        coord_sf(xlim = track_xlim, ylim = track_ylim,
                 datum = st_crs(4326)) +
        scale_x_continuous(breaks = ax_brk, labels = fmt_lon) +
        scale_y_continuous(breaks = ay_brk, labels = fmt_lat) +
        labs(title = "Landfalling vs. Non-landfalling Affecting (SSP5-8.5)",
             x = "Longitude", y = "Latitude") + theme_pub()
      save_both(p4, file.path(DIR_MAP, "Map01d_LF_vs_NonLF"), w = 12, h = 7)
      cat("  Map01d saved\n")
    }
  }, error = function(e) cat("  Map01 error:", e$message, "\n"))
}


# ================================================================
#  PART 14: MAP02 - FREQUENCY MAPS
# ================================================================

cat("\n=== Map02: Frequency maps ===\n")

grid_ann_all   <- grid_freq_all   / N_YEARS
grid_ann_tdts  <- grid_freq_tdts  / N_YEARS
grid_ann_cat12 <- grid_freq_cat12 / N_YEARS
grid_ann_cat35 <- grid_freq_cat35 / N_YEARS

plot_freq_grid <- function(mat, ttl, sub, fname, xlm, ylm, show_cn = FALSE) {
  df <- mat2df(mat, grid_lon, grid_lat, "freq")
  df <- df[df$freq > 0, ]
  if (nrow(df) == 0) { cat("    Skip:", fname, "\n"); return() }
  p <- ggplot()
  if (show_cn) {
    p <- p + geom_sf(data = county_sf, fill = "grey97",
                     colour = COL_CNTY, linewidth = LW_CNTY) +
      geom_sf(data = city_sf, fill = NA, colour = COL_CITY, linewidth = LW_CITY)
  } else {
    p <- p + geom_polygon(data = world_map, aes(long, lat, group = group),
                          fill = "grey95", colour = "grey80", linewidth = 0.1)
  }
  p <- p + geom_raster(data = df, aes(lon, lat, fill = freq)) +
    scale_fill_viridis_c(option = "inferno", direction = -1,
                         name = "Annual\nfrequency", trans = "sqrt") +
    coord_sf(xlim = xlm, ylim = ylm, datum = st_crs(4326)) +
    labs(title = ttl, subtitle = sub, x = "Longitude", y = "Latitude") +
    theme_pub()
  save_both(p, file.path(DIR_MAP, fname), w = 10, h = 8)
}

plot_freq_admin <- function(sf_obj, code_col_name, stats_df, freq_col_name,
                            ttl, sub, fname, xlm, ylm) {
  sf_j <- sf_obj
  sf_j$freq_val <- stats_df[[freq_col_name]][
    match(sf_obj[[code_col_name]], stats_df[[1]])]
  p <- ggplot() +
    geom_sf(data = sf_j, aes(fill = freq_val),
            colour = COL_CNTY, linewidth = LW_CNTY) +
    scale_fill_viridis_c(option = "inferno", direction = -1,
                         na.value = "grey97", name = "Annual\nfrequency",
                         trans = "sqrt") +
    geom_sf(data = city_sf, fill = NA, colour = COL_CITY, linewidth = LW_CITY) +
    coord_sf(xlim = xlm, ylim = ylm, datum = st_crs(4326)) +
    labs(title = ttl, subtitle = sub, x = "Longitude", y = "Latitude") +
    theme_pub()
  save_both(p, file.path(DIR_MAP, fname), w = 10, h = 8)
}

int_sets <- list(
  list(mat = grid_ann_all,   lab = "All",   sub = "All intensities",
       cnty_col = "ann_freq_all",   city_col = "ann_freq_all"),
  list(mat = grid_ann_tdts,  lab = "TD_TS", sub = "TD + TS",
       cnty_col = "ann_freq_tdts",  city_col = "ann_freq_tdts"),
  list(mat = grid_ann_cat12, lab = "Cat12", sub = "Category 1-2",
       cnty_col = "ann_freq_cat12", city_col = "ann_freq_cat12"),
  list(mat = grid_ann_cat35, lab = "Cat35", sub = "Category 3-5",
       cnty_col = "ann_freq_cat35", city_col = "ann_freq_cat35"))

period_label <- paste0(min(YEARS), "-", max(YEARS), " SSP5-8.5")

for (iset in int_sets) {
  tag <- iset$lab
  plot_freq_grid(iset$mat, paste0("TC Frequency (grid): ", iset$sub),
                 paste0(period_label, " annual avg"),
                 paste0("Map02_Grid_Full_", tag), FULL_LON, FULL_LAT, FALSE)
  plot_freq_grid(iset$mat, paste0("TC Frequency (grid, China): ", iset$sub),
                 period_label,
                 paste0("Map02_Grid_China_", tag), CHINA_LON, CHINA_LAT, TRUE)
  plot_freq_admin(county_sf, "县代码", county_grid_stats, iset$cnty_col,
                  paste0("TC Frequency (county): ", iset$sub),
                  period_label,
                  paste0("Map02_County_China_", tag), CHINA_LON, CHINA_LAT)
  plot_freq_admin(city_sf, "市代码", city_grid_stats, iset$city_col,
                  paste0("TC Frequency (city): ", iset$sub),
                  period_label,
                  paste0("Map02_City_China_", tag), CHINA_LON, CHINA_LAT)
}
cat("  Map02 done\n")


# ================================================================
#  PART 15: MAP03 - MAX WIND SPEED MAPS
# ================================================================

cat("\n=== Map03: Max Wind Speed maps ===\n")

plot_ws_grid <- function(mat, ttl, sub, fname, xlm, ylm, show_cn = FALSE) {
  df <- mat2df(mat, grid_lon, grid_lat, "ws")
  df <- df[df$ws > 0.5, ]
  if (nrow(df) == 0) return()
  p <- ggplot()
  if (show_cn) {
    p <- p + geom_sf(data = county_sf, fill = "grey97",
                     colour = COL_CNTY, linewidth = LW_CNTY) +
      geom_sf(data = city_sf, fill = NA, colour = COL_CITY, linewidth = LW_CITY)
  } else {
    p <- p + geom_polygon(data = world_map, aes(long, lat, group = group),
                          fill = "grey95", colour = "grey80", linewidth = 0.1)
  }
  p <- p + geom_raster(data = df, aes(lon, lat, fill = ws)) +
    scale_fill_gradientn(colours = pal_wind(200), limits = c(0, 70),
                         oob = squish,
                         name = expression("Max wind (m s"^{-1}*")")) +
    coord_sf(xlim = xlm, ylim = ylm, datum = st_crs(4326)) +
    labs(title = ttl, subtitle = sub, x = "Longitude", y = "Latitude") +
    theme_pub()
  save_both(p, file.path(DIR_MAP, fname), w = 10, h = 8)
}

plot_ws_admin <- function(sf_obj, code_col_name, stats_df,
                          ttl, fname, xlm, ylm) {
  sf_j <- sf_obj
  sf_j$mws <- stats_df$max_ws[match(sf_obj[[code_col_name]], stats_df[[1]])]
  p <- ggplot() +
    geom_sf(data = sf_j, aes(fill = mws),
            colour = COL_CNTY, linewidth = LW_CNTY) +
    scale_fill_gradientn(colours = pal_wind(200), limits = c(0, 70),
                         oob = squish, na.value = "grey97",
                         name = expression("Max wind (m s"^{-1}*")")) +
    geom_sf(data = city_sf, fill = NA, colour = COL_CITY, linewidth = LW_CITY) +
    coord_sf(xlim = xlm, ylim = ylm, datum = st_crs(4326)) +
    labs(title = ttl, x = "Longitude", y = "Latitude") + theme_pub()
  save_both(p, file.path(DIR_MAP, fname), w = 10, h = 8)
}

plot_ws_grid(grid_max_ws_bld, "Max TC Wind (grid)", period_label,
             "Map03_Grid_Full_MaxWS", FULL_LON, FULL_LAT, FALSE)
plot_ws_grid(grid_max_ws_bld, "Max TC Wind (grid, China)", "",
             "Map03_Grid_China_MaxWS", CHINA_LON, CHINA_LAT, TRUE)
plot_ws_admin(county_sf, "县代码", county_grid_stats,
              "Max TC Wind (county)", "Map03_County_China_MaxWS",
              CHINA_LON, CHINA_LAT)
plot_ws_admin(city_sf, "市代码", city_grid_stats,
              "Max TC Wind (city)", "Map03_City_China_MaxWS",
              CHINA_LON, CHINA_LAT)
cat("  Map03 done\n")


# ================================================================
#  PART 16: MAP04 - POPULATION + TC OVERLAY
# ================================================================

cat("\n=== Map04: Pop overlay maps ===\n")

pop_ref_yr <- as.character(min(YEARS))
pop_r_ref  <- pop_weights[[pop_ref_yr]]
if (is.null(pop_r_ref)) pop_r_ref <- pop_weights[[as.character(max(YEARS))]]

make_pop_tc_map <- function(pop_r, ov_mat, ov_name, ttl, fname,
                            xlm, ylm, show_cn = TRUE, use_log = TRUE) {
  if (is.null(pop_r)) { cat("    No pop raster\n"); return() }
  pop_df <- as.data.frame(rasterToPoints(pop_r))
  names(pop_df) <- c("lon", "lat", "pop")
  pop_df <- pop_df[pop_df$pop > 10, ]
  if (use_log) {
    pop_df$pop_val <- log10(pmax(pop_df$pop, 1))
    pop_lab <- expression("log"[10]*"(Pop)")
  } else {
    pop_df$pop_val <- pop_df$pop / 1000
    pop_lab <- "Pop (thousands)"
  }
  
  ov_df <- mat2df(ov_mat, grid_lon, grid_lat, "freq")
  ov_df <- ov_df[ov_df$freq > 0, ]
  
  p <- ggplot()
  if (show_cn) {
    p <- p + geom_sf(data = county_sf, fill = "grey97",
                     colour = COL_CNTY, linewidth = LW_CNTY)
  } else {
    p <- p + geom_polygon(data = world_map, aes(long, lat, group = group),
                          fill = "grey97", colour = "grey85", linewidth = 0.1)
  }
  
  p <- p + geom_raster(data = pop_df, aes(lon, lat, fill = pop_val)) +
    scale_fill_gradient(low = "#FFFFFF", high = "#FFB7C5", name = pop_lab)
  
  if (nrow(ov_df) > 0) {
    p <- p + geom_raster(data = ov_df, aes(lon, lat, alpha = freq),
                         fill = "#08519C") +
      scale_alpha_continuous(range = c(0.05, 0.55), name = ov_name)
  }
  
  if (!is.null(china_land) && show_cn)
    p <- p + geom_sf(data = china_land, fill = NA,
                     colour = COL_NATL, linewidth = LW_NATL)
  if (show_cn)
    p <- p + geom_sf(data = city_sf, fill = NA, colour = COL_CITY,
                     linewidth = LW_CITY)
  
  p <- p + coord_sf(xlim = xlm, ylim = ylm, datum = st_crs(4326)) +
    labs(title = ttl,
         subtitle = paste0("Pink: population | Blue overlay: ", ov_name),
         x = "Longitude", y = "Latitude") + theme_pub()
  save_both(p, file.path(DIR_MAP, fname),
            w = ifelse(diff(xlm) > 60, 13, 11),
            h = ifelse(diff(ylm) > 35, 9, 8))
}

make_pop_tc_map(pop_r_ref, grid_ann_all, "TC frequency",
                paste0("Population x TC Frequency (log) \u2014 ", pop_ref_yr),
                "Map04a_Pop_Freq_China",
                CHINA_LON, CHINA_LAT, TRUE, TRUE)
make_pop_tc_map(pop_r_ref, grid_ann_all, "TC frequency",
                paste0("Population x TC Frequency (log, full) \u2014 ", pop_ref_yr),
                "Map04b_Pop_Freq_Full", FULL_LON, FULL_LAT, FALSE, TRUE)
make_pop_tc_map(pop_r_ref, grid_max_ws_bld, "Max TC wind",
                paste0("Population x Max TC Wind (log) \u2014 ", pop_ref_yr),
                "Map04c_Pop_MaxWS_China",
                CHINA_LON, CHINA_LAT, TRUE, TRUE)
make_pop_tc_map(pop_r_ref, grid_max_ws_bld, "Max TC wind",
                paste0("Population x Max TC Wind (log, full) \u2014 ", pop_ref_yr),
                "Map04d_Pop_MaxWS_Full", FULL_LON, FULL_LAT, FALSE, TRUE)
make_pop_tc_map(pop_r_ref, grid_ann_all, "TC frequency",
                paste0("Population x TC Frequency (linear) \u2014 ", pop_ref_yr),
                "Map04e_Pop_Freq_China_Linear",
                CHINA_LON, CHINA_LAT, TRUE, FALSE)
make_pop_tc_map(pop_r_ref, grid_ann_all, "TC frequency",
                paste0("Population x TC Frequency (linear, full) \u2014 ", pop_ref_yr),
                "Map04f_Pop_Freq_Full_Linear", FULL_LON, FULL_LAT, FALSE, FALSE)
cat("  Map04 done\n")


# ================================================================
#  PART 17: TC INDIVIDUAL TRACK + FOOTPRINT MAPS
# ================================================================

cat("\n=== TC Individual maps ===\n")
tryCatch({
  tc_ids_map <- unique(tc_meta$tc_id[tc_meta$affects_china | tc_meta$landfalling])
  tc_map_xlim <- c(100, 145)
  tc_map_ylim <- c(5, 50)
  
  tc_cnty_foot <- if (nrow(county_daily) > 0) {
    county_daily %>%
      group_by(tc_id, tc_name, tc_year, county_code) %>%
      summarise(max_ws = max(ws_daily_max_pw, na.rm = TRUE), .groups = "drop")
  } else data.frame()
  
  for (tid in tc_ids_map) {
    tryCatch({
      tm <- tc_meta[tc_meta$tc_id == tid, ]
      tk <- ibtracs_tracks[ibtracs_tracks$tc_id == tid, ]
      foot <- if (nrow(tc_cnty_foot) > 0)
        tc_cnty_foot %>% filter(tc_id == tid) else data.frame()
      if (nrow(tk) < 2 || nrow(foot) == 0) next
      
      tk$sshs_f <- factor(sshs_classify(tk$vmax_ms), levels = SSHS_LABELS)
      present_cats <- levels(tk$sshs_f)[levels(tk$sshs_f) %in%
                                          unique(as.character(tk$sshs_f))]
      
      county_foot <- county_sf
      county_foot$foot_ws <- foot$max_ws[
        match(county_sf$县代码, foot$county_code)]
      
      tk_segs <- make_track_segs(tk)
      if (is.null(tk_segs) || nrow(tk_segs) == 0) next
      tk_segs$sshs <- factor(tk_segs$sshs, levels = present_cats)
      
      p <- ggplot() +
        geom_sf(data = world_sf, fill = "grey93",
                colour = "grey75", linewidth = 0.1) +
        geom_sf(data = county_foot, aes(fill = foot_ws),
                colour = COL_CNTY, linewidth = LW_CNTY) +
        scale_fill_gradientn(
          colours = pal_wind(200), limits = c(0, 70), oob = squish,
          na.value = NA,
          name = expression("Max wind (m s"^{-1}*")")) +
        geom_sf(data = city_sf, fill = NA, colour = COL_CITY,
                linewidth = LW_CITY) +
        geom_segment(data = tk_segs,
                     aes(x = x, y = y, xend = xend, yend = yend,
                         colour = sshs), linewidth = 1.2, alpha = 0.85) +
        scale_colour_manual(values = SSHS_COLORS[present_cats],
                            name = "Intensity", drop = TRUE) +
        geom_point(data = tk[1, ], aes(lon, lat),
                   colour = "green3", size = 3.5, shape = 17) +
        geom_point(data = tk[which.max(tk$vmax_ms), ], aes(lon, lat),
                   colour = "red", size = 3.5, shape = 18) +
        coord_sf(xlim = tc_map_xlim, ylim = tc_map_ylim,
                 datum = st_crs(4326)) +
        labs(title = sprintf("%s (%d) \u2014 SSP5-8.5", tm$tc_name, tm$tc_year),
             subtitle = sprintf("Peak: %s | %s | %d counties",
                                tm$peak_sshs,
                                ifelse(tm$landfalling, "Landfalling",
                                       "Non-landfalling"),
                                nrow(foot)),
             x = "Longitude", y = "Latitude") +
        theme_pub(base_size = 9)
      save_both(p, file.path(DIR_MAP_TC,
                             sprintf("TC_%04d_%s", tm$tc_year, tm$tc_name)),
                w = 10, h = 8)
    }, error = function(e) NULL)
  }
  cat(sprintf("  Done (%d TCs)\n", length(tc_ids_map)))
}, error = function(e) cat("  TC individual error:", e$message, "\n"))


# ================================================================
#  PART 18: FIGURES
# ================================================================

cat("\n=== Generating figures ===\n")

tc_year_sshs <- tc_meta %>%
  filter(affects_china | landfalling) %>%
  count(tc_year, peak_sshs, name = "n")

tc_year_sshs_all <- tc_meta %>%
  count(tc_year, peak_sshs, name = "n")

tryCatch({
  d1 <- tc_year_sshs
  if (nrow(d1) > 0) {
    d1$peak_sshs <- factor(d1$peak_sshs, levels = SSHS_LABELS)
    yr_tot <- d1 %>% group_by(tc_year) %>% summarise(total = sum(n))
    p <- ggplot(d1, aes(tc_year, n, fill = peak_sshs)) +
      geom_col(width = 0.7) +
      geom_text(data = yr_tot, aes(tc_year, total, label = total),
                inherit.aes = FALSE, vjust = -0.5, size = 3) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = seq(min(YEARS), max(YEARS), by = 5)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Year", y = "Number of TCs",
           title = "Annual TC Count (Affecting China) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig01a_AnnualTCs_SSHS"), w = 12, h = 5.5)
    save_data(d1, "Fig01a_AnnualTCs_SSHS")
    cat("  Fig01a saved\n")
  }
}, error = function(e) cat("  Fig01a error:", e$message, "\n"))

tryCatch({
  d1b <- tc_year_sshs_all
  if (nrow(d1b) > 0) {
    d1b$peak_sshs <- factor(d1b$peak_sshs, levels = SSHS_LABELS)
    yr_tot <- d1b %>% group_by(tc_year) %>% summarise(total = sum(n))
    p <- ggplot(d1b, aes(tc_year, n, fill = peak_sshs)) +
      geom_col(width = 0.7) +
      geom_text(data = yr_tot, aes(tc_year, total, label = total),
                inherit.aes = FALSE, vjust = -0.5, size = 3) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = seq(min(YEARS), max(YEARS), by = 5)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Year", y = "Number of TCs",
           title = "Annual TC Count (ALL) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig01b_AnnualTCs_All"), w = 12, h = 5.5)
    save_data(d1b, "Fig01b_AnnualTCs_All")
    cat("  Fig01b saved\n")
  }
}, error = function(e) cat("  Fig01b error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS)) %>%
      arrange(desc(n_counties)) %>%
      mutate(label = factor(label, levels = rev(unique(label))))
    p <- ggplot(d, aes(n_counties, label, fill = peak_sshs_f)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = n_counties), hjust = -0.2, size = 2.5) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Affected counties", y = NULL,
           title = "Affected Counties per TC (TS+) \u2014 SSP5-8.5") +
      theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig02a_Counties_TSplus"),
              w = 9, h = max(5, nrow(d) * 0.35 + 1.5))
    save_data(d %>% as.data.frame() %>% dplyr::select(-label),
              "Fig02a_Counties_TSplus")
    cat("  Fig02a saved\n")
  }
}, error = function(e) cat("  Fig02a error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS)) %>%
      arrange(desc(n_counties)) %>%
      mutate(label = factor(label, levels = rev(unique(label))))
    p <- ggplot(d, aes(n_counties, label, fill = peak_sshs_f)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = n_counties), hjust = -0.2, size = 2.5) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Affected counties", y = NULL,
           title = "Affected Counties per TC (ALL) \u2014 SSP5-8.5") +
      theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig02b_Counties_All"),
              w = 9, h = max(5, nrow(d) * 0.35 + 1.5))
    save_data(d %>% as.data.frame() %>% dplyr::select(-label),
              "Fig02b_Counties_All")
    cat("  Fig02b saved\n")
  }
}, error = function(e) cat("  Fig02b error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS)) %>%
      arrange(desc(n_cities)) %>%
      mutate(label = factor(label, levels = rev(unique(label))))
    p <- ggplot(d, aes(n_cities, label, fill = peak_sshs_f)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = n_cities), hjust = -0.2, size = 2.5) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Affected cities", y = NULL,
           title = "Affected Cities per TC (TS+) \u2014 SSP5-8.5") +
      theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig02c_Cities_TSplus"),
              w = 9, h = max(5, nrow(d) * 0.35 + 1.5))
    save_data(d %>% as.data.frame() %>% dplyr::select(-label),
              "Fig02c_Cities_TSplus")
    cat("  Fig02c saved\n")
  }
}, error = function(e) cat("  Fig02c error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS)) %>%
      arrange(desc(n_cities)) %>%
      mutate(label = factor(label, levels = rev(unique(label))))
    p <- ggplot(d, aes(n_cities, label, fill = peak_sshs_f)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = n_cities), hjust = -0.2, size = 2.5) +
      scale_fill_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Affected cities", y = NULL,
           title = "Affected Cities per TC (ALL) \u2014 SSP5-8.5") +
      theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig02d_Cities_All"),
              w = 9, h = max(5, nrow(d) * 0.35 + 1.5))
    save_data(d %>% as.data.frame() %>% dplyr::select(-label),
              "Fig02d_Cities_All")
    cat("  Fig02d saved\n")
  }
}, error = function(e) cat("  Fig02d error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>% arrange(tc_year, tc_id) %>%
      mutate(idx = seq_len(n()), cum_pop = cumsum(pop_total),
             label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS))
    p <- ggplot(d, aes(idx, cum_pop / 1e6)) +
      geom_line(colour = "#0C2C84", linewidth = 0.8) +
      geom_point(aes(colour = peak_sshs_f), size = 2.5) +
      scale_colour_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = d$idx, labels = d$label) +
      labs(x = NULL, y = "Cumulative pop (millions)",
           title = "Cumulative Exposed Population (TS+) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
    save_both(p, file.path(DIR_FIG_EXP, "Fig03a_CumPop_TSplus"), w = 10, h = 5.5)
    save_data(d, "Fig03a_CumPop_TSplus")
    cat("  Fig03a saved\n")
  }
}, error = function(e) cat("  Fig03a error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>% arrange(tc_year, tc_id) %>%
      mutate(idx = seq_len(n()), cum_pop = cumsum(pop_total),
             label = make_tc_labels(tc_name, tc_year),
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS))
    p <- ggplot(d, aes(idx, cum_pop / 1e6)) +
      geom_line(colour = "#0C2C84", linewidth = 0.8) +
      geom_point(aes(colour = peak_sshs_f), size = 2.5) +
      scale_colour_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = d$idx, labels = d$label) +
      labs(x = NULL, y = "Cumulative pop (millions)",
           title = "Cumulative Exposed Population (ALL) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
    save_both(p, file.path(DIR_FIG_EXP, "Fig03b_CumPop_All"), w = 10, h = 5.5)
    save_data(d, "Fig03b_CumPop_All")
    cat("  Fig03b saved\n")
  }
}, error = function(e) cat("  Fig03b error:", e$message, "\n"))

tryCatch({
  if (nrow(county_daily) > 0 && any(county_daily$exposed_ts)) {
    d <- county_daily %>% filter(exposed_ts)
    d$sshs_f <- factor(d$sshs_cat, levels = SSHS_LABELS)
    present <- levels(d$sshs_f)[levels(d$sshs_f) %in% unique(as.character(d$sshs_f))]
    p <- ggplot(d, aes(ws_daily_max_pw, fill = sshs_f)) +
      geom_histogram(binwidth = 2, boundary = 0, colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = SSHS_COLORS[present], name = "SSHS", drop = TRUE) +
      geom_vline(xintercept = WS_THRESHOLD, linetype = "dashed", colour = "red") +
      labs(x = expression("Daily max pop-weighted wind (m s"^{-1}*")"),
           y = "County-day count",
           title = "Wind Speed Distribution (TS+) \u2014 SSP5-8.5") + theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig04a_WsDist_TSplus"), w = 8, h = 5)
    cat("  Fig04a saved\n")
  }
}, error = function(e) cat("  Fig04a error:", e$message, "\n"))

tryCatch({
  if (nrow(county_daily) > 0) {
    d <- county_daily
    p <- ggplot(d, aes(ws_daily_max_pw)) +
      geom_histogram(binwidth = 2, boundary = 0, fill = "#3182BD",
                     colour = "white", linewidth = 0.2) +
      labs(x = expression("Daily max pop-weighted wind (m s"^{-1}*")"),
           y = "County-day count",
           title = "Wind Speed Distribution (ALL) \u2014 SSP5-8.5") + theme_pub()
    save_both(p, file.path(DIR_FIG_EXP, "Fig04b_WsDist_All"), w = 8, h = 5)
    cat("  Fig04b saved\n")
  }
}, error = function(e) cat("  Fig04b error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)))
    age_long <- d %>%
      dplyr::select(label, pop_lt20, pop_20to44, pop_45to64, pop_ge65) %>%
      pivot_longer(-label, names_to = "age", values_to = "pop") %>%
      mutate(age = factor(age,
                          levels = c("pop_lt20","pop_20to44","pop_45to64","pop_ge65"),
                          labels = c("<20","20-44","45-64","65+")))
    age_cols <- c("<20"="#7FCDBB","20-44"="#41B6C4","45-64"="#1D91C0","65+"="#0C2C84")
    p <- ggplot(age_long, aes(label, pop / 1e6, fill = age)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = age_cols, name = "Age") +
      labs(x = NULL, y = "Exposed pop (millions)",
           title = "Age-specific Exposure per TC (TS+) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig05a_AgeExposure_TSplus"),
              w = max(8, nrow(d) * 0.7 + 2), h = 6)
    save_data(age_long, "Fig05a_AgeExposure_TSplus")
    cat("  Fig05a saved\n")
  }
}, error = function(e) cat("  Fig05a error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)))
    age_long <- d %>%
      dplyr::select(label, pop_lt20, pop_20to44, pop_45to64, pop_ge65) %>%
      pivot_longer(-label, names_to = "age", values_to = "pop") %>%
      mutate(age = factor(age,
                          levels = c("pop_lt20","pop_20to44","pop_45to64","pop_ge65"),
                          labels = c("<20","20-44","45-64","65+")))
    age_cols <- c("<20"="#7FCDBB","20-44"="#41B6C4","45-64"="#1D91C0","65+"="#0C2C84")
    p <- ggplot(age_long, aes(label, pop / 1e6, fill = age)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = age_cols, name = "Age") +
      labs(x = NULL, y = "Exposed pop (millions)",
           title = "Age-specific Exposure per TC (ALL) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig05b_AgeExposure_All"),
              w = max(18), h = 6)
    save_data(age_long, "Fig05b_AgeExposure_All")
    cat("  Fig05b saved\n")
  }
}, error = function(e) cat("  Fig05b error:", e$message, "\n"))

tryCatch({
  if (nrow(prov_summary) > 0) {
    d <- prov_summary %>% arrange(pop_total) %>%
      mutate(province = factor(province, levels = province))
    age_l <- d %>%
      dplyr::select(province, pop_lt20, pop_20to44, pop_45to64, pop_ge65) %>%
      pivot_longer(-province, names_to = "age", values_to = "pop") %>%
      mutate(age = factor(age,
                          levels = c("pop_lt20","pop_20to44","pop_45to64","pop_ge65"),
                          labels = c("<20","20-44","45-64","65+")))
    age_cols <- c("<20"="#7FCDBB","20-44"="#41B6C4","45-64"="#1D91C0","65+"="#0C2C84")
    p <- ggplot(age_l, aes(pop / 1e6, province, fill = age)) +
      geom_col(width = 0.65) +
      scale_fill_manual(values = age_cols, name = "Age") +
      labs(x = "Exposed population (millions)", y = NULL,
           title = "Provincial TC Exposure (SSP5-8.5, excl. HK/Macau/TW)") +
      theme_pub() + theme(legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig06_Provincial"),
              w = 9, h = max(5, nrow(d) * 0.35 + 1.5))
    save_data(prov_summary, "Fig06_Provincial")
    cat("  Fig06 saved\n")
  }
}, error = function(e) cat("  Fig06 error:", e$message, "\n"))

tryCatch({
  if (nrow(county_daily) > 0 && any(county_daily$exposed_ts)) {
    cal <- county_daily %>% filter(exposed_ts) %>%
      group_by(tc_id, tc_name, tc_year, date) %>%
      summarise(n_counties = n_distinct(county_code),
                mean_ws = mean(ws_daily_max_pw, na.rm = TRUE),
                .groups = "drop")
    cal <- cal %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year))
    p <- ggplot(cal, aes(date, reorder(label, tc_year), fill = mean_ws)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_gradientn(colours = pal_wind(200),
                           limits = c(WS_THRESHOLD, 70), oob = squish,
                           name = expression("Wind (m s"^{-1}*")")) +
      labs(x = "Date", y = NULL,
           title = "TC Exposure Calendar \u2014 SSP5-8.5") +
      theme_pub() + theme(axis.text.y = element_text(size = 7))
    save_both(p, file.path(DIR_FIG_EXP, "Fig07_Calendar"), w = 14, h = 8)
    save_data(cal, "Fig07_Calendar")
    cat("  Fig07 saved\n")
  }
}, error = function(e) cat("  Fig07 error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)),
             cum_cnt = cumsum(n_counties), cum_cit = cumsum(n_cities))
    cl <- d %>% dplyr::select(label, cum_cnt, cum_cit) %>%
      pivot_longer(-label, names_to = "m", values_to = "count") %>%
      mutate(m = ifelse(m == "cum_cnt", "Counties", "Cities"))
    p <- ggplot(cl, aes(label, count, colour = m, group = m)) +
      geom_line(linewidth = 0.8) + geom_point(size = 2) +
      scale_colour_manual(values = c("Counties"="#E6550D","Cities"="#3182BD"),
                          name = NULL) +
      labs(x = NULL, y = "Cumulative count",
           title = "Cumulative Affected Areas (TS+) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig08a_CumAreas_TSplus"), w = 10, h = 6)
    save_data(cl, "Fig08a_CumAreas_TSplus")
    cat("  Fig08a saved\n")
  }
}, error = function(e) cat("  Fig08a error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)),
             cum_cnt = cumsum(n_counties), cum_cit = cumsum(n_cities))
    cl <- d %>% dplyr::select(label, cum_cnt, cum_cit) %>%
      pivot_longer(-label, names_to = "m", values_to = "count") %>%
      mutate(m = ifelse(m == "cum_cnt", "Counties", "Cities"))
    p <- ggplot(cl, aes(label, count, colour = m, group = m)) +
      geom_line(linewidth = 0.8) + geom_point(size = 2) +
      scale_colour_manual(values = c("Counties"="#E6550D","Cities"="#3182BD"),
                          name = NULL) +
      labs(x = NULL, y = "Cumulative count",
           title = "Cumulative Affected Areas (ALL) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig08b_CumAreas_All"), w = 10, h = 6)
    save_data(cl, "Fig08b_CumAreas_All")
    cat("  Fig08b saved\n")
  }
}, error = function(e) cat("  Fig08b error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_exp) > 0) {
    d <- tc_summary_exp %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)),
             idx = seq_len(n()),
             c1 = cumsum(pop_lt20), c2 = cumsum(pop_20to44),
             c3 = cumsum(pop_45to64), c4 = cumsum(pop_ge65),
             cum_total = c1 + c2 + c3 + c4,
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS))
    al <- d %>% dplyr::select(idx, label, c1, c2, c3, c4) %>%
      pivot_longer(c1:c4, names_to = "age", values_to = "pop") %>%
      mutate(age = factor(age, levels = c("c4","c3","c2","c1"),
                          labels = c("65+","45-64","20-44","<20")))
    age_cols <- c("<20"="#7FCDBB","20-44"="#41B6C4","45-64"="#1D91C0","65+"="#0C2C84")
    p <- ggplot() +
      geom_area(data = al, aes(idx, pop / 1e6, fill = age), alpha = 0.85) +
      scale_fill_manual(values = age_cols, name = "Age") +
      geom_line(data = d, aes(idx, cum_total / 1e6), linewidth = 0.6,
                colour = "black") +
      geom_point(data = d, aes(idx, cum_total / 1e6, colour = peak_sshs_f),
                 size = 2.5) +
      scale_colour_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = d$idx, labels = d$label) +
      labs(x = NULL, y = "Cumulative pop (millions)",
           title = "Cumulative Exposed Population by Age (TS+) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig09a_CumPopAge_TSplus"), w = 11, h = 6)
    save_data(al, "Fig09a_CumPopAge_TSplus")
    cat("  Fig09a saved\n")
  }
}, error = function(e) cat("  Fig09a error:", e$message, "\n"))

tryCatch({
  if (nrow(tc_summary_all) > 0) {
    d <- tc_summary_all %>%
      arrange(tc_year, tc_id) %>%
      mutate(label = make_tc_labels(tc_name, tc_year),
             label = factor(label, levels = unique(label)),
             idx = seq_len(n()),
             c1 = cumsum(pop_lt20), c2 = cumsum(pop_20to44),
             c3 = cumsum(pop_45to64), c4 = cumsum(pop_ge65),
             cum_total = c1 + c2 + c3 + c4,
             peak_sshs_f = factor(peak_sshs, levels = SSHS_LABELS))
    al <- d %>% dplyr::select(idx, label, c1, c2, c3, c4) %>%
      pivot_longer(c1:c4, names_to = "age", values_to = "pop") %>%
      mutate(age = factor(age, levels = c("c4","c3","c2","c1"),
                          labels = c("65+","45-64","20-44","<20")))
    age_cols <- c("<20"="#7FCDBB","20-44"="#41B6C4","45-64"="#1D91C0","65+"="#0C2C84")
    p <- ggplot() +
      geom_area(data = al, aes(idx, pop / 1e6, fill = age), alpha = 0.85) +
      scale_fill_manual(values = age_cols, name = "Age") +
      geom_line(data = d, aes(idx, cum_total / 1e6), linewidth = 0.6,
                colour = "black") +
      geom_point(data = d, aes(idx, cum_total / 1e6, colour = peak_sshs_f),
                 size = 2.5) +
      scale_colour_manual(values = SSHS_COLORS, name = "Peak SSHS", drop = TRUE) +
      scale_x_continuous(breaks = d$idx, labels = d$label) +
      labs(x = NULL, y = "Cumulative pop (millions)",
           title = "Cumulative Exposed Population by Age (ALL) \u2014 SSP5-8.5") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "top")
    save_both(p, file.path(DIR_FIG_EXP, "Fig09b_CumPopAge_All"), w = 11, h = 6)
    save_data(al, "Fig09b_CumPopAge_All")
    cat("  Fig09b saved\n")
  }
}, error = function(e) cat("  Fig09b error:", e$message, "\n"))


# ================================================================
#  PART 19: PANEL DATASET CONSTRUCTION
# ================================================================

cat("\n=== Panel dataset construction ===\n")

if (nrow(county_daily) > 0 && any(county_daily$exposed_ts)) {
  
  exposed_counties <- county_daily %>%
    filter(exposed_ts) %>%
    dplyr::select(county_code, province, city, county,
                  province_code, city_code) %>% distinct()
  
  panel_dates <- do.call(c, lapply(YEARS, function(yr) {
    seq.Date(as.Date(paste0(yr, "-", sprintf("%02d", TC_SEASON_START), "-01")),
             as.Date(paste0(yr, "-", sprintf("%02d", TC_SEASON_END), "-31")),
             by = "day")
  }))
  
  panel_skeleton <- expand.grid(county_code = exposed_counties$county_code,
                                date = panel_dates, stringsAsFactors = FALSE)
  panel_skeleton <- panel_skeleton %>%
    left_join(exposed_counties, by = "county_code") %>%
    mutate(year = year(date), month = month(date),
           dow = wday(date), doy = yday(date))
  
  cat(sprintf("  Skeleton: %s rows\n", format(nrow(panel_skeleton), big.mark = ",")))
  
  tc_exp <- county_daily %>% filter(exposed_ts) %>%
    group_by(county_code, date) %>%
    slice_max(ws_daily_max_pw, n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(county_code, date, tc_id_exp = tc_id,
                  tc_name_exp = tc_name,
                  ws_daily_max_pw, ws_daily_mean_pw, ws_daily_max,
                  ws_daily_mean, dist_min_km, tc_vmax_max,
                  sshs_cat, landfalling)
  
  panel <- panel_skeleton %>%
    left_join(tc_exp, by = c("county_code", "date"))
  panel$tc_exposed <- !is.na(panel$tc_id_exp)
  panel$ws_daily_max_pw[!panel$tc_exposed]  <- 0
  panel$ws_daily_mean_pw[!panel$tc_exposed] <- 0
  panel$ws_daily_max[!panel$tc_exposed]     <- 0
  panel$sshs_ts_plus <- as.integer(panel$ws_daily_max_pw >= WS_THRESHOLD)
  
  pop_lookup <- county_daily %>% filter(!is.na(pop_total)) %>%
    group_by(county_code, tc_year) %>% slice(1) %>% ungroup() %>%
    dplyr::select(county_code, tc_year, pop_total, pop_lt20,
                  pop_20to44, pop_45to64, pop_ge65)
  panel <- panel %>%
    left_join(pop_lookup, by = c("county_code", "year" = "tc_year"))
  
  cat(sprintf("  Panel: %s rows | Exposed: %s\n",
              format(nrow(panel), big.mark = ","),
              format(sum(panel$tc_exposed), big.mark = ",")))
  
  
  # ================================================================
  #  PART 20: LAG STRUCTURE (VECTORIZED)
  # ================================================================
  
  cat("\n=== Lag structure (vectorized) ===\n")
  lag_start <- Sys.time()
  
  panel <- panel %>% arrange(county_code, date)
  
  cat("  Computing lag 0-21 and lead 1-3...\n")
  
  panel <- panel %>%
    group_by(county_code) %>%
    mutate(
      ws_lag0  = ws_daily_max_pw,
      ws_lag1  = lag(ws_daily_max_pw, 1),
      ws_lag2  = lag(ws_daily_max_pw, 2),
      ws_lag3  = lag(ws_daily_max_pw, 3),
      ws_lag4  = lag(ws_daily_max_pw, 4),
      ws_lag5  = lag(ws_daily_max_pw, 5),
      ws_lag6  = lag(ws_daily_max_pw, 6),
      ws_lag7  = lag(ws_daily_max_pw, 7),
      ws_lag8  = lag(ws_daily_max_pw, 8),
      ws_lag9  = lag(ws_daily_max_pw, 9),
      ws_lag10 = lag(ws_daily_max_pw, 10),
      ws_lag11 = lag(ws_daily_max_pw, 11),
      ws_lag12 = lag(ws_daily_max_pw, 12),
      ws_lag13 = lag(ws_daily_max_pw, 13),
      ws_lag14 = lag(ws_daily_max_pw, 14),
      ws_lag15 = lag(ws_daily_max_pw, 15),
      ws_lag16 = lag(ws_daily_max_pw, 16),
      ws_lag17 = lag(ws_daily_max_pw, 17),
      ws_lag18 = lag(ws_daily_max_pw, 18),
      ws_lag19 = lag(ws_daily_max_pw, 19),
      ws_lag20 = lag(ws_daily_max_pw, 20),
      ws_lag21 = lag(ws_daily_max_pw, 21),
      ws_lead1 = lead(ws_daily_max_pw, 1),
      ws_lead2 = lead(ws_daily_max_pw, 2),
      ws_lead3 = lead(ws_daily_max_pw, 3),
      tc_exp_lag0  = as.integer(tc_exposed),
      tc_exp_lag1  = lag(as.integer(tc_exposed), 1),
      tc_exp_lag2  = lag(as.integer(tc_exposed), 2),
      tc_exp_lag3  = lag(as.integer(tc_exposed), 3),
      tc_exp_lag4  = lag(as.integer(tc_exposed), 4),
      tc_exp_lag5  = lag(as.integer(tc_exposed), 5),
      tc_exp_lag6  = lag(as.integer(tc_exposed), 6),
      tc_exp_lag7  = lag(as.integer(tc_exposed), 7),
      tc_exp_lag8  = lag(as.integer(tc_exposed), 8),
      tc_exp_lag9  = lag(as.integer(tc_exposed), 9),
      tc_exp_lag10 = lag(as.integer(tc_exposed), 10),
      tc_exp_lag11 = lag(as.integer(tc_exposed), 11),
      tc_exp_lag12 = lag(as.integer(tc_exposed), 12),
      tc_exp_lag13 = lag(as.integer(tc_exposed), 13),
      tc_exp_lag14 = lag(as.integer(tc_exposed), 14),
      tc_exp_lag15 = lag(as.integer(tc_exposed), 15),
      tc_exp_lag16 = lag(as.integer(tc_exposed), 16),
      tc_exp_lag17 = lag(as.integer(tc_exposed), 17),
      tc_exp_lag18 = lag(as.integer(tc_exposed), 18),
      tc_exp_lag19 = lag(as.integer(tc_exposed), 19),
      tc_exp_lag20 = lag(as.integer(tc_exposed), 20),
      tc_exp_lag21 = lag(as.integer(tc_exposed), 21)
    ) %>%
    ungroup()
  
  panel$ws_cum_0_3  <- rowSums(panel[, paste0("ws_lag", 0:3)],  na.rm = FALSE)
  panel$ws_cum_0_7  <- rowSums(panel[, paste0("ws_lag", 0:7)],  na.rm = FALSE)
  panel$ws_cum_0_14 <- rowSums(panel[, paste0("ws_lag", 0:14)], na.rm = FALSE)
  panel$ws_cum_0_21 <- rowSums(panel[, paste0("ws_lag", 0:21)], na.rm = FALSE)
  panel$ws_ma_0_3  <- panel$ws_cum_0_3  / 4
  panel$ws_ma_0_7  <- panel$ws_cum_0_7  / 8
  panel$ws_ma_0_14 <- panel$ws_cum_0_14 / 15
  panel$ws_ma_0_21 <- panel$ws_cum_0_21 / 22
  panel$tc_window_0_3 <- as.integer(
    rowSums(panel[, paste0("tc_exp_lag", 0:3)], na.rm = TRUE) > 0)
  panel$tc_window_0_7 <- as.integer(
    rowSums(panel[, paste0("tc_exp_lag", 0:7)], na.rm = TRUE) > 0)
  
  lag_time <- as.numeric(difftime(Sys.time(), lag_start, units = "secs"))
  cat(sprintf("  Lag structure done in %.1f seconds\n", lag_time))
  
  
  # ================================================================
  #  PART 21: SAVE PANEL & QC
  # ================================================================
  
  cat("\n=== Saving panel & QC ===\n")
  
  saveRDS(panel, file.path(DIR_PANEL, "TC_Health_Panel.rds"))
  cat("  RDS saved\n")
  
  tryCatch({
    core_cols <- c("county_code","date","year","month","dow","doy",
                   "province","city","county","province_code","city_code",
                   "tc_exposed","tc_id_exp","tc_name_exp",
                   "ws_daily_max_pw","ws_daily_mean_pw","ws_daily_max",
                   "dist_min_km","tc_vmax_max","sshs_cat","landfalling",
                   "sshs_ts_plus","pop_total","pop_lt20","pop_20to44",
                   "pop_45to64","pop_ge65",
                   "ws_lag0","ws_cum_0_3","ws_cum_0_7","ws_cum_0_14",
                   "ws_cum_0_21","ws_ma_0_3","ws_ma_0_7",
                   "tc_window_0_3","tc_window_0_7")
    export_cols <- intersect(core_cols, names(panel))
    write.csv(panel[, export_cols],
              file.path(DIR_PANEL, "TC_Health_Panel.csv"),
              row.names = FALSE, na = "")
    cat("  CSV saved\n")
  }, error = function(e) cat("  CSV error:", e$message, "\n"))
  
  cat("\n  Quality checks:\n")
  n_dup <- sum(duplicated(panel[, c("county_code","date")]))
  cat(sprintf("    Duplicates: %d %s\n", n_dup,
              ifelse(n_dup == 0, "OK", "WARNING")))
  ws_rng <- range(panel$ws_daily_max_pw, na.rm = TRUE)
  cat(sprintf("    WS range: %.1f - %.1f m/s\n", ws_rng[1], ws_rng[2]))
  pop_cov <- mean(!is.na(panel$pop_total[panel$tc_exposed]))
  cat(sprintf("    Pop coverage: %.1f%%\n", pop_cov * 100))
  lag0_ok <- all(panel$ws_lag0 == panel$ws_daily_max_pw, na.rm = TRUE)
  cat(sprintf("    Lag0 == daily WS: %s\n", ifelse(lag0_ok, "OK", "MISMATCH")))
  
  qc <- data.frame(
    Check = c("Duplicates","WS_max","Pop_coverage","Lag0_ok",
              "N_panel","N_counties","N_TCs","N_LF","N_NLF_affect"),
    Value = c(n_dup, round(ws_rng[2],1), round(pop_cov,3), lag0_ok,
              nrow(panel), n_distinct(panel$county_code),
              n_tot, n_lf, n_aff),
    stringsAsFactors = FALSE)
  save_data(qc, "Quality_Check", DIR_TABLE)
  
} else {
  cat("  No exposed county-days - panel not built\n")
}


# ================================================================
#  PART 22: VARIABLE DICTIONARY
# ================================================================

cat("\n=== Variable Dictionary ===\n")

var_dict <- data.frame(
  Variable = c(
    "tc_id", "tc_name", "tc_year", "date", "county_code", "city_code",
    "province_code", "province", "city", "county",
    "tc_lon", "tc_lat", "tc_vmax", "tc_vmax_max",
    "peak_vmax", "peak_sshs",
    "landfalling", "landfall_time", "affects_china",
    "n_timesteps", "date_start", "date_end",
    "ws_daily_max_pw", "ws_daily_mean_pw", "ws_daily_max", "ws_daily_mean",
    "ws_popweighted", "ws_max", "ws_area_mean",
    "exposed_ts", "sshs_cat", "sshs_ts_plus", "tc_exposed",
    "dist_km", "dist_min_km",
    "pop_total", "pop_lt20", "pop_20to44", "pop_45to64", "pop_ge65",
    "ws_lagX", "ws_leadX", "tc_exp_lagX",
    "ws_cum_0_K", "ws_ma_0_K", "tc_window_0_K",
    "freq_all", "freq_tdts", "freq_cat12", "freq_cat35",
    "freq_td", "freq_ts", "freq_cat1", "freq_cat2",
    "freq_cat3", "freq_cat4", "freq_cat5",
    "ann_freq_all", "ann_freq_tdts", "ann_freq_cat12", "ann_freq_cat35",
    "ann_freq_td", "ann_freq_ts", "ann_freq_cat1", "ann_freq_cat2",
    "ann_freq_cat3", "ann_freq_cat4", "ann_freq_cat5",
    "max_ws", "mean_ws",
    "event_start", "event_end", "max_ws_pw", "mean_ws_pw",
    "event_sshs", "event_seq", "n_days",
    "WS_THRESHOLD", "TD", "TS", "Category 1-5"
  ),
  Description = c(
    "TC identifier: YYYY_SID_NAME_startdate_enddate",
    "TC name from HighResMIP TRACK output",
    "Year of TC occurrence",
    "Calendar date (YYYY-MM-DD)",
    "County-level administrative code (6-digit)",
    "City-level administrative code (4-digit)",
    "Province-level administrative code (2-digit)",
    "Province name (Chinese)", "City name (Chinese)", "County name (Chinese)",
    "TC centre longitude at this timestep (deg E)",
    "TC centre latitude at this timestep (deg N)",
    "TC max sustained wind at this timestep (m/s, from NC attribute)",
    "Max TC Vmax across timesteps on this day (m/s)",
    "TC lifetime peak Vmax (m/s)",
    "TC lifetime peak SSHS (TD/TS/Cat1-5)",
    "TRUE if TC centre crossed China boundary (+0.15deg buffer)",
    "Datetime of first landfall point (UTC, if applicable)",
    "TRUE if TC wind field reaches any Chinese county",
    "Number of 6-hourly timesteps processed for this TC",
    "Date of first track point", "Date of last track point",
    "Daily max pop-weighted wind speed (m/s)",
    "Daily mean pop-weighted wind speed (m/s)",
    "Daily max grid-cell wind speed (m/s)",
    "Daily mean area-averaged wind speed (m/s)",
    "Timestep-level pop-weighted wind (m/s, raw)",
    "Timestep-level max grid wind (m/s, raw)",
    "Timestep-level area-mean wind (m/s, raw)",
    "TRUE if ws_daily_max_pw >= 17.5 m/s",
    "SSHS category of ws_daily_max_pw",
    "Binary: 1 if ws_daily_max_pw >= 17.5 m/s",
    "Binary: 1 if county exposed to any TC on this day",
    "Distance from TC centre to county centroid (km)",
    "Min distance across daily timesteps (km)",
    "Total population (CMIP6 SSP5-8.5)",
    "Population aged <20 years",
    "Population aged 20-44 years",
    "Population aged 45-64 years",
    "Population aged >=65 years",
    "Lagged wind speed at lag X days (X = 0..21)",
    "Lead wind speed at lead X days (X = 1..3)",
    "Lagged TC exposure indicator at lag X days",
    "Cumulative wind sum over lags 0 to K",
    "Moving average wind over lags 0 to K",
    "Binary: any TC exposure within lag 0-K window",
    "Total TC count (2025-2050)",
    "TD+TS count", "Category 1-2 count", "Category 3-5 count",
    "TD only count", "TS only count",
    "Cat1 count", "Cat2 count", "Cat3 count", "Cat4 count", "Cat5 count",
    "Annual avg of freq_all", "Annual avg of freq_tdts",
    "Annual avg of freq_cat12", "Annual avg of freq_cat35",
    "Annual avg TD freq", "Annual avg TS freq",
    "Annual avg Cat1 freq", "Annual avg Cat2 freq",
    "Annual avg Cat3 freq", "Annual avg Cat4 freq",
    "Annual avg Cat5 freq",
    "Max wind speed (m/s)",
    "Mean wind speed (m/s)",
    "First date of exposure event", "Last date",
    "Max pop-weighted wind across event (m/s)",
    "Mean pop-weighted wind across event (m/s)",
    "SSHS of event max wind", "Sequential event number",
    "Number of days in event",
    "17.5 m/s: TS onset threshold",
    "0-17.5 m/s: Tropical Depression",
    "17.5-33.0 m/s: Tropical Storm",
    "Cat1: 33-43 | Cat2: 43-50 | Cat3: 50-58 | Cat4: 58-70 | Cat5: >=70"
  ),
  stringsAsFactors = FALSE
)

save_data(var_dict, "Variable_Dictionary_Comprehensive", DIR_TABLE)
save_data(var_dict, "Panel_Data_Dictionary", DIR_PANEL)
cat("  Variable dictionary saved\n")


# ================================================================
#  FINAL SUMMARY
# ================================================================

total_time <- as.numeric(difftime(Sys.time(), step2_start, units = "mins"))

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  STEP 2.1 v2.3 (Accelerated Future Projection, coastal-fix) COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat(sprintf("Scenario: SSP5-8.5 | Period: %d-%d\n", min(YEARS), max(YEARS)))
cat(sprintf("Wind field: Holland parametric only (no ERA5 blending)\n"))
cat(sprintf("Population: CMIP6 SSP5-8.5 gridded + county/city age data\n"))
cat(sprintf("Total time: %.1f minutes\n", total_time))
cat(sprintf("TCs: %d total | %d landfalling | %d non-LF affecting\n",
            n_tot, n_lf, n_aff))
if (exists("panel") && is.data.frame(panel)) {
  cat(sprintf("Panel: %s rows | %d counties | %s exposed days\n",
              format(nrow(panel), big.mark = ","),
              n_distinct(panel$county_code),
              format(sum(panel$tc_exposed), big.mark = ",")))
}
cat(sprintf("\nMaps:    %s\n", DIR_MAP))
cat(sprintf("Figures: %s\n", DIR_FIG_EXP))
cat(sprintf("Tables:  %s\n", DIR_TABLE))
cat(sprintf("Panel:   %s\n", DIR_PANEL))
cat("\nStep 2.1 v2.3 complete. Next: Step 3 - Health analyses.\n")




