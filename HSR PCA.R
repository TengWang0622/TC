###############################################################################
#
#                          P   R   O   J   E   C   T
#                                     of
#   Health system resilience shapes health equity and adaptation under
#                     tropical cyclones in China
#
#   Module: City-level HSR index construction (PCA, 2000-2024) and
#           Bayesian hierarchical projection under adaptation scenarios
#           (2025-2049)
#
#
#   Developed by Teng Wang, Hanxu Shi
#   Version - 20260711_V4
#
###############################################################################

############################################
#             0. Preparation
############################################

library(readxl)
library(tidyverse)
library(zoo)          
library(openxlsx)
library(brms)
library(posterior)
library(bayesplot)
library(patchwork)
library(viridis)
library(ggridges)
library(scales)
library(sf)

set.seed(2026)

# ============================ Output paths ============================

OutPath <- "/Volumes/Lenovo/HSR/Output_Structure"

dir_fig <- file.path(OutPath, "Figure")
dir_tab <- file.path(OutPath, "Table")
dir_res <- file.path(OutPath, "Result")

# New subfolder for annual HSR maps
dir_fig_maps <- file.path(dir_fig, "HSR_Annual_Maps")

for (d in c(OutPath, dir_fig, dir_tab, dir_res, dir_fig_maps)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Helper: save figure (PDF + PNG) and its underlying plotting data (xlsx)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 600)
}
save_fig_data <- function(df, name) {
  write.xlsx(as.data.frame(df), file.path(dir_tab, paste0(name, "_PlotData.xlsx")),
             rowNames = FALSE)
}

# Save map figures to the maps subfolder
save_map <- function(p, name, w = 11, h = 9) {
  tryCatch(ggsave(file.path(dir_fig_maps, paste0(name, ".pdf")), p,
                  width = w, height = h, dpi = 300, device = cairo_pdf),
           error = function(e)
             tryCatch(ggsave(file.path(dir_fig_maps, paste0(name, ".pdf")), p,
                             width = w, height = h, dpi = 300),
                      error = function(e2) cat("  PDF fail:", e2$message, "\n")))
  tryCatch(ggsave(file.path(dir_fig_maps, paste0(name, ".png")), p,
                  width = w, height = h, dpi = 600, bg = "white"),
           error = function(e) cat("  PNG fail:", e$message, "\n"))
  while (dev.cur() > 1) tryCatch(dev.off(), error = function(e) break)
}

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size, color = "grey30"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Map theme (no grid lines, clean cartographic style)
theme_map <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text          = element_text(colour = "black"),
      plot.title    = element_text(size = base_size + 2, face = "bold",
                                   hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 0.5, hjust = 0,
                                   margin = margin(b = 6), colour = "grey30"),
      plot.caption  = element_text(size = base_size - 2, hjust = 0,
                                   colour = "grey50", margin = margin(t = 8)),
      axis.title    = element_text(size = base_size, face = "bold"),
      axis.text     = element_text(size = base_size - 1, colour = "black"),
      legend.title  = element_text(size = base_size - 0.5, face = "bold"),
      legend.text   = element_text(size = base_size - 1.5),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.margin   = margin(5, 5, 5, 5, "mm"))
}

colors_group <- c("High" = "#2166AC", "Moderate" = "#F4A582", "Low" = "#B2182B")
colors_scenario <- c("No Adaptation" = "#D73027", "Autonomous" = "#FC8D59",
                     "Targeted" = "#91BFDB", "Aspirational" = "#4575B4")

# ============================ Loading files ============================

CityLevel_Pop_2016_2024 <- read_excel("/Volumes/Lenovo/Database/PopAgeSex/2016_2024_CityLevel_Pop.xlsx")
NBS_2000_2024_ShortName <- read_excel("/Volumes/Lenovo/HSR/NBS/NBS_2000_2024_ShortName.xlsx")
HSR_ShortName <- read_excel("/Volumes/Lenovo/HSR/Src/HSR_ShortName.xlsx")

# Unit harmonization for HSR (2021-2024 Medical_Assistance * 10000)
HSR_ShortName <- HSR_ShortName %>%
  mutate(Medical_Assistance = ifelse(Year >= 2021,
                                     Medical_Assistance * 10000,
                                     Medical_Assistance))

NBS_2000_2024_ShortName$CityCode <- as.character(NBS_2000_2024_ShortName$CityCode)


Year_Start <- 2000
Year_End   <- 2024
Years_Hist <- Year_Start:Year_End
Years_Proj <- 2025:2049

# ============================ Load China Shapefiles ============================

cat("=== Loading China shapefiles for mapping ===\n")

DIR_SHP <- "/Volumes/Lenovo/Database/MapExtract/Src/ChinaMap"

county_sf <- st_read(file.path(DIR_SHP, "县.shp"), quiet = TRUE)
city_sf   <- st_read(file.path(DIR_SHP, "市.shp"), quiet = TRUE)

county_sf <- st_transform(county_sf, 4326)
city_sf   <- st_transform(city_sf, 4326)

county_sf$省代码 <- as.character(county_sf$省代码)
county_sf$市代码 <- as.character(county_sf$市代码)
county_sf$县代码 <- as.character(county_sf$县代码)
city_sf$省代码   <- as.character(city_sf$省代码)
city_sf$市代码   <- as.character(city_sf$市代码)

cat(sprintf("  Counties: %d | Cities: %d\n", nrow(county_sf), nrow(city_sf)))

# Create China national boundary (union of counties)
cat("  Creating China national boundary...\n")
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

if (!is.null(china_land)) {
  cat("  China land boundary OK\n")
} else {
  cat("  WARNING: China land boundary creation FAILED\n")
}

# Map style parameters
COL_CNTY <- "grey88"
LW_CNTY  <- 0.06
COL_CITY <- "grey55"
LW_CITY  <- 0.12
COL_NATL <- "grey40"
LW_NATL  <- 0.35

CHINA_LON <- c(73, 136)
CHINA_LAT <- c(15, 55)

# Coordinate formatters
fmt_lon <- function(x) paste0(x, "\u00B0E")
fmt_lat <- function(y) ifelse(y >= 0, paste0(y, "\u00B0N"), paste0(-y, "\u00B0S"))

###############################################################################
#     PART 0. Data preparation: Population handling & Per-capita transforms
###############################################################################

# =============================================================================
# 0.1 Add total population columns to CityLevel_Pop_2016_2024
# =============================================================================

for (yr in 2016:2024) {
  col_lt20   <- paste0("<20 years_", yr)
  col_2044   <- paste0("20 to 44 years_", yr)
  col_4564   <- paste0("45 to 64 years_", yr)
  col_gt64   <- paste0(">64 years_", yr)
  col_total  <- paste0("TotalPop_", yr)
  
  CityLevel_Pop_2016_2024[[col_total]] <-
    CityLevel_Pop_2016_2024[[col_lt20]] +
    CityLevel_Pop_2016_2024[[col_2044]] +
    CityLevel_Pop_2016_2024[[col_4564]] +
    CityLevel_Pop_2016_2024[[col_gt64]]
}

# =============================================================================
# 0.2 Fill missing Province_Pop in HSR for 2021-2024
# =============================================================================

# Province-CityCode mapping from NBS
Province_City_Map <- NBS_2000_2024_ShortName %>%
  select(Province, ProvinceCode, CityCode) %>%
  distinct() %>%
  mutate(CityCode = as.character(CityCode))

# Reshape CityLevel_Pop total columns to long format
pop_total_cols <- paste0("TotalPop_", 2016:2024)

CityPop_Long <- CityLevel_Pop_2016_2024 %>%
  select(CityCode, all_of(pop_total_cols)) %>%
  mutate(CityCode = as.character(CityCode)) %>%
  pivot_longer(
    cols = all_of(pop_total_cols),
    names_to = "Year",
    values_to = "CityTotalPop"
  ) %>%
  mutate(Year = as.numeric(gsub("TotalPop_", "", Year)))

# Aggregate to province level for 2021-2024
Province_Pop_Agg <- CityPop_Long %>%
  left_join(Province_City_Map, by = "CityCode") %>%
  filter(!is.na(Province), Year >= 2021 & Year <= 2024) %>%
  group_by(Province, Year) %>%
  summarise(Province_Pop_Calc = sum(CityTotalPop, na.rm = TRUE) / 10000,
            .groups = "drop")

# Special handling for Chongqing
Chongqing_Pop <- CityPop_Long %>%
  filter(grepl("^50", CityCode), Year >= 2021 & Year <= 2024) %>%
  group_by(Year) %>%
  summarise(Province_Pop_Calc = sum(CityTotalPop, na.rm = TRUE) / 10000,
            .groups = "drop") %>%
  mutate(Province = "重庆市")

Province_Pop_Agg <- Province_Pop_Agg %>%
  filter(Province != "重庆市") %>%
  bind_rows(Chongqing_Pop)

# Fill HSR Province_Pop
HSR_ShortName <- HSR_ShortName %>%
  left_join(Province_Pop_Agg, by = c("Province", "Year")) %>%
  mutate(
    Province_Pop = ifelse(is.na(Province_Pop) & !is.na(Province_Pop_Calc),
                          Province_Pop_Calc,
                          Province_Pop)
  ) %>%
  select(-Province_Pop_Calc)

# =============================================================================
# 0.3 NBS Pop interpolation
# =============================================================================

interp_pop <- function(pop, year, max_extrap_years = 5) {
  non_na <- !is.na(pop)
  n_known <- sum(non_na)
  
  if (n_known == 0) return(pop)
  if (n_known == 1) {
    pop[is.na(pop)] <- pop[non_na][1]
    return(pop)
  }
  
  known_years <- year[non_na]
  known_pop   <- pop[non_na]
  
  # Step A: linear interpolation for internal gaps
  pop_filled <- approx(x = known_years, y = known_pop,
                       xout = year, method = "linear", rule = 1)$y
  
  # Step B: leading NAs
  leading_na <- which(is.na(pop_filled) & year < min(known_years))
  if (length(leading_na) > 0) {
    n_use <- min(n_known, 5)
    fit_years <- known_years[1:n_use]
    fit_pop   <- known_pop[1:n_use]
    lm_fit    <- lm(fit_pop ~ fit_years)
    extrap_limit <- min(known_years) - max_extrap_years
    
    for (i in leading_na) {
      if (year[i] >= extrap_limit) {
        pop_filled[i] <- predict(lm_fit, newdata = data.frame(fit_years = year[i]))
      } else {
        pop_filled[i] <- predict(lm_fit, newdata = data.frame(fit_years = extrap_limit))
      }
    }
  }
  
  # Step C: trailing NAs
  trailing_na <- which(is.na(pop_filled) & year > max(known_years))
  if (length(trailing_na) > 0) {
    n_use <- min(n_known, 5)
    fit_years <- known_years[(n_known - n_use + 1):n_known]
    fit_pop   <- known_pop[(n_known - n_use + 1):n_known]
    lm_fit    <- lm(fit_pop ~ fit_years)
    extrap_limit <- max(known_years) + max_extrap_years
    
    for (i in trailing_na) {
      if (year[i] <= extrap_limit) {
        pop_filled[i] <- predict(lm_fit, newdata = data.frame(fit_years = year[i]))
      } else {
        pop_filled[i] <- predict(lm_fit, newdata = data.frame(fit_years = extrap_limit))
      }
    }
  }
  
  # Step D: safety — population must be positive
  pop_filled <- pmax(pop_filled, min(known_pop) * 0.5)
  return(pop_filled)
}

NBS_2000_2024_ShortName <- NBS_2000_2024_ShortName %>%
  group_by(CityCode) %>%
  arrange(Year) %>%
  mutate(
    Pop_original = Pop,
    Pop = interp_pop(Pop, Year, max_extrap_years = 5)
  ) %>%
  ungroup()

# Check interpolation results
pop_check <- NBS_2000_2024_ShortName %>%
  filter(is.na(Pop_original) & !is.na(Pop)) %>%
  select(Province, CityCode, Year, Pop_original, Pop)
cat("共填充了", nrow(pop_check), "个NBS Pop缺失值\n")

pop_still_na <- NBS_2000_2024_ShortName %>% filter(is.na(Pop))
if (nrow(pop_still_na) > 0) {
  cat("WARNING:", nrow(pop_still_na),
      "city-years still have NA Pop after interpolation. These will be excluded.\n")
}

# =============================================================================
# 0.4 Per-capita transformation: HSR (provincial level)
# =============================================================================

HSR_Vars_PerCapita <- c(
  "No_Inst", "No_Hosp", "No_General_Hosp", "No_Primary", "No_CDC",
  "No_Supervision", "No_Com_Health_Station", "No_H_staff", "No_Physician",
  "No_Nurse", "No_Pharmacist", "No_Management_Staff", "No_MC_inst",
  "No_CDC_Staff", "No_Bed",
  "Tot_H_Expenditure", "Gov_H_Expenditure",
  "H_Asset", "H_Income", "H_Cost", "H_Financial_Support", "H_Staff_Expenditure",
  "No_Outpatient", "No_Check", "No_ED", "No_Hosp_Admission",
  "Tech_Policy_Consultation", "Public_Edu", "Foodborne_Disease",
  "Basic_Insurance", "Maternity_Insurance", "Medical_Assistance"
)

# Keep provincial rows only, exclude national aggregate
df_HSR <- HSR_ShortName %>%
  filter(Year %in% Years_Hist,
         !Province %in% c("全国", "中国", "合计"))

# Per-capita transform (divide by Province_Pop, unit: 万人)
df_HSR_pc <- df_HSR %>%
  mutate(across(all_of(HSR_Vars_PerCapita), ~ .x / Province_Pop))

# Identify HSR indicator columns (everything except identifiers)
HSR_ID_cols <- c("Province", "Year", "Province_Pop")
HSR_Indicator_cols <- setdiff(colnames(df_HSR_pc), HSR_ID_cols)

# =============================================================================
# 0.5 Per-capita transformation: NBS (city level)
# =============================================================================

NBS_Vars_PerCapita <- c(
  "GDP", "GDP_1", "GDP_2", "GDP_3",
  "Employment",
  "Pop_Unemployed",
  "Pop_Health_Social_Welfare",
  "Pop_Public_Management",
  "Pop_Social_Service",
  "Pop_Health_Sports_Welfare",
  "Tax",
  "Asset_Investment",
  "LocGov_Revenue", "LocGov_Expenditure",
  "Scientific_Expenditure", "Edu_Expenditure",
  "Saving",
  "No_Edu_Inst",
  "Loc_H_Inst", "Loc_Hosp_Inst",
  "Loc_Bed", "Loc_Physician",
  "TeleCom_Revenue",
  "Pop_Pension_Insurance",
  "Pop_Medical_Insurance",
  "Pop_Employment_Insurance"
)

# Filter NBS to historical years, exclude rows with NA Pop
df_NBS <- NBS_2000_2024_ShortName %>%
  filter(Year %in% Years_Hist, !is.na(Pop))

# Per-capita transform (Pop in 万人, result is "per 万人")
df_NBS_pc <- df_NBS %>%
  mutate(across(all_of(NBS_Vars_PerCapita), ~ .x / Pop))

# Identify NBS indicator columns (excluding IDs and Pop)
NBS_ID_cols <- c("Province", "ProvinceCode", "CityCode", "Year", "Pop", "Pop_original")
NBS_Indicator_cols <- setdiff(colnames(df_NBS_pc), NBS_ID_cols)

# =============================================================================
# 0.6 Merge: Assign provincial HSR per-capita indicators to each city
# =============================================================================

# Prepare HSR for merge: keep Province, Year, and all indicator columns
HSR_for_merge <- df_HSR_pc %>%
  select(Province, Year, all_of(HSR_Indicator_cols))

# Check for column name conflicts between NBS and HSR
shared_names <- intersect(NBS_Indicator_cols, HSR_Indicator_cols)
if (length(shared_names) > 0) {
  cat("WARNING: Shared indicator names between NBS and HSR:",
      paste(shared_names, collapse = ", "), "\n")
  cat("HSR columns will be suffixed with '_prov' to avoid ambiguity.\n")
  HSR_rename_map <- setNames(
    paste0(shared_names, "_prov"),
    shared_names
  )
  HSR_for_merge <- HSR_for_merge %>%
    rename(!!!HSR_rename_map)
  HSR_Indicator_cols[HSR_Indicator_cols %in% shared_names] <-
    paste0(HSR_Indicator_cols[HSR_Indicator_cols %in% shared_names], "_prov")
}

# Merge: each city inherits its province's HSR indicators
df_merged <- merge(
  df_NBS_pc,
  HSR_for_merge,
  by = c("Province", "Year"), all.x = TRUE
)

cat("Merged dataset dimensions:", nrow(df_merged), "rows x", ncol(df_merged), "cols\n")
cat("Number of unique cities:", n_distinct(df_merged$CityCode), "\n")
cat("Year range:", min(df_merged$Year), "-", max(df_merged$Year), "\n")

# Save merged dataset
saveRDS(df_merged, file.path(dir_res, "Merged_NBS_HSR_PerCapita.rds"))


###############################################################################
#                       Select indicators
###############################################################################

colnames(df_merged)

df_merged <- df_merged[, c("Province","ProvinceCode","City","CityCode","Region","Year",
                           
                           # ------------------------ Financing ------------------------
                           "GDP","H_Income","H_Cost","H_Asset","H_Staff_Expenditure","Tot_H_Expenditure",                 
                           "Outpatient_Cost","Urban_H_Expenditure","Rural_H_Expenditure","Ave_Wage","Hosp_Cost",
                           "H_Financial_Support","Gov_H_Expenditure","Saving","Asset_Investment","Personal_H_Expenditure",
                           
                           # ------------------------ Health workforce ------------------------
                           "No_Nurse","No_H_staff","No_Physician","No_Pharmacist", "Loc_Physician","Workload",                
                           "No_CDC_Staff",
                           
                           # ------------------------ Service delivery ------------------------
                           "No_Hosp","No_General_Hosp","No_Inst","Loc_H_Inst","Loc_Hosp_Inst","No_Com_Health_Station","No_Primary",
                           "No_Hosp_Admission","No_ED","No_Outpatient","No_Bed","Loc_Bed", "Bed_Rate","No_MC_inst",
                           "Prenatal_Care_Rate","Delivery_Rate",
                           "Medical_Assistance","Hosp_Day",
                           
                           # ------------------------ Information technologies ------------------------
                           "TeleCom_Revenue","Scientific_Expenditure","Tech_Policy_Consultation","No_CDC", 
                           
                           # ------------------------ Disease control and preventative interventions ------------------------
                           "longevity","Foodborne_Disease","Iodine_Prevention","Hepatitis_Incidence", "Hepatitis_A_Incidence",                  
                           "Fever_Incidence","Premarital_Check","No_Check",
                           "Infectious_AB_Mortality","Infectious_AB_Incidence","MC_Management_Rate","MC_Post_rate",
                           
                           # ------------------------ Governance and Leadership ------------------------
                           "No_Supervision","Pop_Social_Service","Pop_Health_Sports_Welfare","Pop_Health_Social_Welfare",
                           "No_Management_Staff","Pop_Public_Management",                 
                           
                           # Medical products, insurance, and technologies
                           "Basic_Insurance","Maternity_Insurance","Pop_Pension_Insurance","Pop_Medical_Insurance"                  
)]


###############################################################################
#     PART A. City-level HSR index (Hierarchical PCA composite), 2000-2024
###############################################################################

# ===================== A1. Define candidate indicators =====================

Vars_ID <- c("Region", "Province", "ProvinceCode", "City", "CityCode", "Year")
Vars_All <- setdiff(colnames(df_merged),
                    c(Vars_ID, "Pop", "Pop_original", "Province_Pop"))

cat("Total candidate indicators:", length(Vars_All), "\n")

# ===================== A2. Missingness audit & indicator screening ==========

Miss_Summary <- df_merged %>%
  summarise(across(all_of(Vars_All), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "Indicator", values_to = "Prop_Missing") %>%
  arrange(desc(Prop_Missing))

write.xlsx(Miss_Summary, file.path(dir_tab, "Indicator_Missingness_Summary.xlsx"))

Miss_Threshold <- 1   # drop indicators with >threshold missing
Vars_Keep <- Miss_Summary %>%
  filter(Prop_Missing <= Miss_Threshold) %>%
  pull(Indicator)

cat("Indicators retained after missingness screening:",
    length(Vars_Keep), "of", length(Vars_All), "\n")

# ===================== A3. Imputation =====================

df_imp <- df_merged %>%
  dplyr::select(all_of(c(Vars_ID, Vars_Keep))) %>%
  arrange(CityCode, Year) %>%
  group_by(CityCode) %>%
  mutate(across(all_of(Vars_Keep), ~ {
    if (sum(!is.na(.x)) >= 2) {
      zoo::na.approx(.x, x = Year, na.rm = FALSE, rule = 2)
    } else .x
  })) %>%
  ungroup() %>%
  group_by(Year) %>%
  mutate(across(all_of(Vars_Keep),
                ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup()

remaining_na <- sum(is.na(df_imp[, Vars_Keep]))
if (remaining_na > 0) {
  cat("INFO:", remaining_na, "NAs remain after temporal+cross-sectional imputation.",
      "Filling with global median.\n")
  df_imp <- df_imp %>%
    mutate(across(all_of(Vars_Keep),
                  ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))
}

stopifnot(!anyNA(df_imp[, Vars_Keep]))
saveRDS(df_imp, file.path(dir_res, "HSR_Indicators_Imputed_CityLevel.rds"))

# ===================== A4. Normalization (pooled across city-years) ==========

normalize_data <- function(x) {
  lower_bound <- quantile(x, 0.05, na.rm = TRUE)
  upper_bound <- quantile(x, 0.95, na.rm = TRUE)
  x[x < lower_bound] <- lower_bound
  x[x > upper_bound] <- upper_bound
  (x - lower_bound) / (upper_bound - lower_bound)
}

df_nor <- df_imp
df_nor[, Vars_Keep] <- lapply(df_nor[, Vars_Keep], normalize_data)

sd_check <- sapply(df_nor[, Vars_Keep], sd, na.rm = TRUE)
Vars_Keep <- Vars_Keep[sd_check > 1e-8]
cat("Indicators after removing zero-variance:", length(Vars_Keep), "\n")

# ===================== A5. Hierarchical PCA: Pillar-level then equal average =

# --- A5.1 Define pillar membership ---

Pillar_Definitions <- list(
  Financing = c(
    "GDP", "H_Income", "H_Cost", "H_Asset", "H_Staff_Expenditure", "Tot_H_Expenditure",
    "Outpatient_Cost", "Urban_H_Expenditure", "Rural_H_Expenditure", "Ave_Wage", "Hosp_Cost",
    "H_Financial_Support", "Gov_H_Expenditure", "Saving", "Asset_Investment", "Personal_H_Expenditure"
  ),
  Health_Workforce = c(
    "No_Nurse", "No_H_staff", "No_Physician", "No_Pharmacist", "Loc_Physician", "Workload",
    "No_CDC_Staff"
  ),
  Service_Delivery = c(
    "No_Hosp", "No_General_Hosp", "No_Inst", "Loc_H_Inst", "Loc_Hosp_Inst",
    "No_Com_Health_Station", "No_Primary",
    "No_Hosp_Admission", "No_ED", "No_Outpatient", "No_Bed", "Loc_Bed", "Bed_Rate", "No_MC_inst",
    "Prenatal_Care_Rate", "Delivery_Rate",
    "Medical_Assistance", "Hosp_Day"
  ),
  Information_Technologies = c(
    "TeleCom_Revenue", "Scientific_Expenditure", "Tech_Policy_Consultation", "No_CDC"
  ),
  Disease_Control = c(
    "longevity", "Foodborne_Disease", "Iodine_Prevention", "Hepatitis_Incidence",
    "Hepatitis_A_Incidence", "Fever_Incidence", "Premarital_Check", "No_Check",
    "Infectious_AB_Mortality", "Infectious_AB_Incidence", "MC_Management_Rate", "MC_Post_rate"
  ),
  Governance_Leadership = c(
    "No_Supervision", "Pop_Social_Service", "Pop_Health_Sports_Welfare",
    "Pop_Health_Social_Welfare", "No_Management_Staff", "Pop_Public_Management"
  ),
  Medical_Products_Insurance = c(
    "Basic_Insurance", "Maternity_Insurance", "Pop_Pension_Insurance", "Pop_Medical_Insurance"
  )
)

# Reference variables for orienting PC1 within each pillar (higher = better health system)
Pillar_Ref_Vars <- list(
  Financing             = c("GDP", "Gov_H_Expenditure", "Tot_H_Expenditure"),
  Health_Workforce      = c("No_Physician", "Loc_Physician", "No_Nurse"),
  Service_Delivery      = c("No_Bed", "Loc_Bed", "No_Hosp"),
  Information_Technologies = c("TeleCom_Revenue", "Scientific_Expenditure"),
  Disease_Control       = c("longevity", "Iodine_Prevention", "MC_Management_Rate"),
  Governance_Leadership = c("Pop_Health_Sports_Welfare", "Pop_Social_Service", "No_Supervision"),
  Medical_Products_Insurance = c("Basic_Insurance", "Pop_Medical_Insurance", "Pop_Pension_Insurance")
)

# Negative-direction indicators (higher raw value = WORSE health system)
Negative_Indicators <- c(
  "Hepatitis_Incidence", "Hepatitis_A_Incidence", "Fever_Incidence",
  "Infectious_AB_Mortality", "Infectious_AB_Incidence",
  "Foodborne_Disease", "Hosp_Day", "Personal_H_Expenditure"
)

# --- A5.2 Intersect with Vars_Keep (indicators that passed screening) ---

Pillar_Indicators <- lapply(Pillar_Definitions, function(vars) {
  intersect(vars, Vars_Keep)
})

# Report pillar composition
cat("\n========== Pillar Composition After Screening ==========\n")
for (pname in names(Pillar_Indicators)) {
  cat(sprintf("  %-30s: %d indicators\n", pname, length(Pillar_Indicators[[pname]])))
}
cat("=========================================================\n\n")

# Check for any Vars_Keep not assigned to a pillar
all_pillar_vars <- unlist(Pillar_Indicators)
unassigned <- setdiff(Vars_Keep, all_pillar_vars)
if (length(unassigned) > 0) {
  cat("WARNING: The following indicators are in Vars_Keep but not assigned to any pillar:\n")
  cat("  ", paste(unassigned, collapse = ", "), "\n")
  cat("  These will be excluded from the hierarchical HSR index.\n\n")
}

# Remove pillars with fewer than 2 indicators (PCA needs at least 2)
Pillar_Indicators <- Pillar_Indicators[sapply(Pillar_Indicators, length) >= 2]
cat("Pillars with >= 2 indicators:", length(Pillar_Indicators), "\n")

# --- A5.3 Within-pillar PCA and scoring ---

Pillar_Scores <- data.frame(row_idx = 1:nrow(df_nor))
Pillar_PCA_Results <- list()
Pillar_Summary <- data.frame()

for (pname in names(Pillar_Indicators)) {
  
  p_vars <- Pillar_Indicators[[pname]]
  cat(sprintf("\n--- Pillar: %s (%d indicators) ---\n", pname, length(p_vars)))
  
  # Extract data matrix for this pillar
  X_pillar <- df_nor[, p_vars, drop = FALSE]
  
  # Run PCA
  pca_pillar <- prcomp(X_pillar, center = TRUE, scale. = TRUE)
  
  # PC1 loadings
  pc1_loadings <- pca_pillar$rotation[, 1]
  
  # Orient PC1: use reference variable (first available)
  ref_candidates <- Pillar_Ref_Vars[[pname]]
  ref_found <- intersect(ref_candidates, p_vars)
  
  if (length(ref_found) > 0) {
    ref_var <- ref_found[1]
    if (pc1_loadings[ref_var] < 0) {
      pc1_loadings <- -pc1_loadings
      pca_pillar$rotation[, 1] <- -pca_pillar$rotation[, 1]
      cat(sprintf("  PC1 flipped based on reference variable: %s\n", ref_var))
    } else {
      cat(sprintf("  PC1 direction confirmed by reference variable: %s\n", ref_var))
    }
  } else {
    cat("  WARNING: No reference variable found. PC1 direction may need manual check.\n")
  }
  
  # Compute pillar score: X * PC1_loadings (weighted sum)
  pillar_raw_score <- as.numeric(as.matrix(X_pillar) %*% pc1_loadings)
  
  # Min-max normalize pillar score to [0, 1]
  p_min <- min(pillar_raw_score)
  p_max <- max(pillar_raw_score)
  pillar_score <- (pillar_raw_score - p_min) / (p_max - p_min)
  
  # Store
  Pillar_Scores[[pname]] <- pillar_score
  
  # Store PCA results
  eigenvalues_p <- pca_pillar$sdev^2
  Pillar_PCA_Results[[pname]] <- list(
    pca = pca_pillar,
    pc1_loadings = pc1_loadings,
    var_explained_pc1 = eigenvalues_p[1] / sum(eigenvalues_p) * 100,
    norm_params = list(min = p_min, max = p_max),
    indicators = p_vars
  )
  
  # Summary row
  Pillar_Summary <- bind_rows(Pillar_Summary, data.frame(
    Pillar = pname,
    N_Indicators = length(p_vars),
    PC1_Var_Explained_Pct = round(eigenvalues_p[1] / sum(eigenvalues_p) * 100, 2),
    PC1_Eigenvalue = round(eigenvalues_p[1], 3),
    Reference_Variable = ifelse(length(ref_found) > 0, ref_found[1], "NONE")
  ))
  
  cat(sprintf("  PC1 variance explained: %.1f%%\n",
              eigenvalues_p[1] / sum(eigenvalues_p) * 100))
}

# Remove the row_idx helper column
Pillar_Scores <- Pillar_Scores %>% select(-row_idx)

# --- A5.4 Overall HSR index: equal-weighted average across pillars ---

n_pillars <- ncol(Pillar_Scores)
HSR_overall_raw <- rowMeans(Pillar_Scores, na.rm = TRUE)

# Final min-max rescaling to [0, 1]
overall_min <- min(HSR_overall_raw)
overall_max <- max(HSR_overall_raw)
HSR_overall <- (HSR_overall_raw - overall_min) / (overall_max - overall_min)

cat(sprintf("\nOverall HSR computed as equal average of %d pillars.\n", n_pillars))

# --- A5.5 Assemble HSR City Panel ---

HSR_City_Panel <- df_nor %>%
  dplyr::select(Province, ProvinceCode, City, CityCode, Year) %>%
  bind_cols(Pillar_Scores) %>%
  mutate(HSR_raw = HSR_overall_raw,
         HSR = HSR_overall)

# Save normalization parameters for projections
HSR_norm_params <- list(
  overall_min = overall_min,
  overall_max = overall_max,
  n_pillars = n_pillars,
  pillar_results = Pillar_PCA_Results
)
saveRDS(HSR_norm_params, file.path(dir_res, "HSR_Norm_Params_Hierarchical.rds"))

saveRDS(HSR_City_Panel, file.path(dir_res, "HSR_City_Panel_2000_2024.rds"))
write.xlsx(HSR_City_Panel, file.path(dir_res, "HSR_City_Panel_2000_2024.xlsx"))

# --- A5.6 Save pillar summary and loadings ---

write.xlsx(Pillar_Summary, file.path(dir_tab, "Pillar_PCA_Summary.xlsx"))
cat("\n========== Pillar PCA Summary ==========\n")
print(Pillar_Summary)

# Save all PC1 loadings across pillars
All_Loadings <- data.frame()
for (pname in names(Pillar_PCA_Results)) {
  loadings_df <- data.frame(
    Pillar = pname,
    Indicator = names(Pillar_PCA_Results[[pname]]$pc1_loadings),
    PC1_Loading = as.numeric(Pillar_PCA_Results[[pname]]$pc1_loadings)
  )
  All_Loadings <- bind_rows(All_Loadings, loadings_df)
}
write.xlsx(All_Loadings, file.path(dir_res, "Pillar_PC1_Loadings_All.xlsx"))

# ===================== A5b. Province-level aggregation ======================

City_Pop_for_weight <- df_NBS %>%
  select(CityCode, Year, Pop) %>%
  distinct()

HSR_Prov_Panel <- HSR_City_Panel %>%
  left_join(City_Pop_for_weight, by = c("CityCode", "Year")) %>%
  group_by(Province, Year) %>%
  summarise(
    HSR = weighted.mean(HSR, w = Pop, na.rm = TRUE),
    n_cities = n(),
    .groups = "drop"
  )

saveRDS(HSR_Prov_Panel, file.path(dir_res, "HSR_Province_Panel_2000_2024.rds"))
write.xlsx(HSR_Prov_Panel, file.path(dir_res, "HSR_Province_Panel_2000_2024.xlsx"))

cat("City-level HSR panel:", n_distinct(HSR_City_Panel$CityCode), "cities,",
    n_distinct(HSR_City_Panel$Year), "years\n")

# ===================== A6. PCA diagnostics ==================================

# --- A6.1 Pillar-level diagnostics: Scree plots ---

scree_list <- list()
for (pname in names(Pillar_PCA_Results)) {
  evals <- Pillar_PCA_Results[[pname]]$pca$sdev^2
  n_pcs <- length(evals)
  scree_list[[pname]] <- data.frame(
    Pillar = pname,
    PC = 1:n_pcs,
    Eigenvalue = evals,
    CumVar = cumsum(evals) / sum(evals) * 100
  )
}
scree_all <- bind_rows(scree_list)

pA1_pillar <- ggplot(scree_all, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "grey40") +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  facet_wrap(~ Pillar, scales = "free", ncol = 3) +
  labs(x = "Principal component", y = "Eigenvalue") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_text(face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_fig(pA1_pillar, "FigS_PCA_Scree_ByPillar", 10, 8)
save_fig_data(scree_all, "FigS_PCA_Scree_ByPillar")

# --- A6.2 PC1 loading bar plot by pillar ---

All_Loadings <- All_Loadings %>%
  mutate(Direction = ifelse(PC1_Loading >= 0, "Positive", "Negative"))

pA3_pillar <- ggplot(All_Loadings, aes(x = reorder(Indicator, PC1_Loading),
                                       y = PC1_Loading, fill = Direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B")) +
  coord_flip() +
  facet_wrap(~ Pillar, scales = "free", ncol = 2) +
  labs(x = NULL, y = "PC1 loading") +
  theme_minimal(base_size = 8) +
  theme(
    axis.title = element_text(face = "bold"),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  )

save_fig(pA3_pillar, "FigS_PC1_Loadings_ByPillar", 12, 14)
save_fig_data(All_Loadings, "FigS_PC1_Loadings_ByPillar")

# --- A6.3 Pillar score correlation heatmap ---

pillar_cor <- cor(Pillar_Scores, use = "complete.obs")
pillar_cor_long <- as.data.frame(as.table(pillar_cor)) %>%
  rename(Pillar1 = Var1, Pillar2 = Var2, Correlation = Freq)

pA_cor <- ggplot(pillar_cor_long, aes(Pillar1, Pillar2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_fig(pA_cor, "FigS_Pillar_Correlation_Heatmap", 7, 6)
save_fig_data(pillar_cor_long, "FigS_Pillar_Correlation_Heatmap")

# --- A6.4 Overall PCA summary table ---

PCA_Summary <- data.frame(
  Statistic = c("Number of pillars",
                "Total indicators across all pillars",
                paste0("Pillar: ", Pillar_Summary$Pillar, " - PC1 Var%")),
  Value = c(n_pillars,
            sum(Pillar_Summary$N_Indicators),
            paste0(Pillar_Summary$PC1_Var_Explained_Pct, "%"))
)

write.xlsx(PCA_Summary, file.path(dir_tab, "PCA_Diagnostics_Summary_Hierarchical.xlsx"))

# ===================== A7. Monte Carlo weight-perturbation robustness =======

MC_years <- c(2000, 2005, 2010, 2015, 2020, 2024)
No_MC <- 1000

MC_Summary_all <- data.frame()

for (yr in MC_years) {
  
  HSR_yr <- HSR_City_Panel %>% filter(Year == yr)
  yr_idx <- which(df_nor$Year == yr)
  
  MC_hsr <- matrix(NA, nrow = length(yr_idx), ncol = No_MC)
  
  for (i in 1:No_MC) {
    pillar_scores_mc <- matrix(NA, nrow = length(yr_idx), ncol = n_pillars)
    
    for (j in seq_along(names(Pillar_PCA_Results))) {
      pname <- names(Pillar_PCA_Results)[j]
      p_vars <- Pillar_PCA_Results[[pname]]$indicators
      pc1_loads <- Pillar_PCA_Results[[pname]]$pc1_loadings
      p_min <- Pillar_PCA_Results[[pname]]$norm_params$min
      p_max <- Pillar_PCA_Results[[pname]]$norm_params$max
      
      # Perturb weights: multiply each loading by Uniform(0.5, 1.5)
      perturb <- runif(length(pc1_loads), 0.5, 1.5)
      w_perturbed <- pc1_loads * perturb
      
      X_p <- as.matrix(df_nor[yr_idx, p_vars])
      raw_p <- as.numeric(X_p %*% w_perturbed)
      pillar_scores_mc[, j] <- (raw_p - p_min) / (p_max - p_min)
    }
    
    # Equal average across pillars
    hsr_mc_raw <- rowMeans(pillar_scores_mc, na.rm = TRUE)
    MC_hsr[, i] <- (hsr_mc_raw - overall_min) / (overall_max - overall_min)
  }
  
  # Rank for each MC iteration
  MC_rank <- apply(MC_hsr, 2, function(s) rank(-s, ties.method = "min"))
  
  MC_yr <- HSR_yr %>%
    mutate(Year_label   = yr,
           Ranking_PCA  = rank(-HSR, ties.method = "min"),
           Rank_lower5  = apply(MC_rank, 1, quantile, 0.05),
           Rank_median  = apply(MC_rank, 1, quantile, 0.50),
           Rank_upper95 = apply(MC_rank, 1, quantile, 0.95)) %>%
    arrange(Ranking_PCA)
  
  MC_Summary_all <- bind_rows(MC_Summary_all, MC_yr)
}

MC_Summary_all <- MC_Summary_all %>%
  mutate(Year_label = factor(Year_label, levels = MC_years))

write.xlsx(MC_Summary_all, file.path(dir_tab, "HSR_City_Ranking_MC_Robustness_MultiYear.xlsx"))

p_mc <- ggplot(MC_Summary_all, aes(x = Ranking_PCA)) +
  geom_segment(aes(xend = Ranking_PCA, y = Rank_lower5, yend = Rank_upper95),
               color = "black", linewidth = 0.1) +
  geom_point(aes(y = Rank_lower5), color = "black", size = 0.1) +
  geom_point(aes(y = Rank_upper95), color = "black", size = 0.1) +
  geom_point(aes(y = Rank_median), color = "#2166AC", shape = 4, size = 1, stroke = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "#B2182B", linewidth = 0.8, linetype = "solid") +
  facet_wrap(~ Year_label, nrow = 2, scales = "free") +
  labs(x = "PCA-based ranking", y = "Monte Carlo ranking") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_fig(p_mc, "Fig_MC_Ranking_Robustness_City_MultiYear", 12, 6)
save_fig_data(MC_Summary_all %>% select(Year_label, Ranking_PCA, Rank_lower5, Rank_median, Rank_upper95),
              "Fig_MC_Ranking_Robustness_City_MultiYear")

# ===================== A8. Historical HSR figure ============================

# --- 8a. Year-specific tertile classification and trajectory plot -----------

HSR_City_Panel_class <- HSR_City_Panel %>%
  group_by(Year) %>%
  mutate(tertile = ntile(HSR, 3),
         Resilience_Group = factor(case_when(tertile == 3 ~ "High",
                                             tertile == 2 ~ "Moderate",
                                             tertile == 1 ~ "Low"),
                                   levels = c("High", "Moderate", "Low"))) %>%
  ungroup()

fig1_data <- HSR_City_Panel_class %>%
  group_by(Resilience_Group, Year) %>%
  summarise(mean_hsr = mean(HSR), q25 = quantile(HSR, 0.25),
            q75 = quantile(HSR, 0.75), .groups = "drop")

p1 <- ggplot(fig1_data, aes(Year, mean_hsr,
                            color = Resilience_Group, fill = Resilience_Group)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_group, name = "Resilience group") +
  scale_fill_manual(values = colors_group, name = "Resilience group") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year", y = "HSR index") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.margin = margin(4, 6, 4, 6),
    legend.key.height = unit(0.8, "lines")
  )

save_fig(p1, "Fig_Historical_HSR_Trajectories_City", 6, 4)
save_fig_data(fig1_data, "Fig_Historical_HSR_Trajectories_City")

# --- 8b. Pillar-level trajectory plot ---

Pillar_Names_Plot <- names(Pillar_Indicators)

fig1b_data <- HSR_City_Panel %>%
  select(CityCode, Year, all_of(Pillar_Names_Plot)) %>%
  pivot_longer(cols = all_of(Pillar_Names_Plot),
               names_to = "Pillar", values_to = "Score") %>%
  group_by(Pillar, Year) %>%
  summarise(mean_score = mean(Score), q25 = quantile(Score, 0.25),
            q75 = quantile(Score, 0.75), .groups = "drop")

p1b <- ggplot(fig1b_data, aes(x = Year, y = mean_score, color = Pillar, fill = Pillar)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year", y = "Pillar score") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

save_fig(p1b, "Fig_Historical_Pillar_Trajectories", 9, 5)
save_fig_data(fig1b_data, "Fig_Historical_Pillar_Trajectories")

# --- 8c. Full city-year ranking table with pillar scores and HSR level ------

# Build enriched ranking table with all pillar scores and HSR level
HSR_City_Rank_Full <- HSR_City_Panel %>%
  group_by(Year) %>%
  mutate(
    Rank = rank(-HSR, ties.method = "min"),
    # Year-specific tertile classification
    tertile = ntile(HSR, 3),
    HSR_Level = factor(case_when(
      tertile == 3 ~ "High",
      tertile == 2 ~ "Moderate",
      tertile == 1 ~ "Low"),
      levels = c("High", "Moderate", "Low"))
  ) %>%
  ungroup() %>%
  arrange(Year, Rank) %>%
  select(Year, Rank, Province, ProvinceCode, City, CityCode,
         HSR, HSR_Level,
         # All pillar dimension scores
         all_of(names(Pillar_Indicators)),
         HSR_raw)

# Verify the new columns
cat("\n========== HSR_City_Rank_Full structure ==========\n")
cat("Dimensions:", nrow(HSR_City_Rank_Full), "x", ncol(HSR_City_Rank_Full), "\n")
cat("Columns:", paste(colnames(HSR_City_Rank_Full), collapse = ", "), "\n")
cat("HSR Level distribution (all years combined):\n")
print(table(HSR_City_Rank_Full$HSR_Level))
cat("\nHSR Level distribution by year (sample):\n")
print(table(HSR_City_Rank_Full$Year[HSR_City_Rank_Full$Year %in% c(2000, 2012, 2024)],
            HSR_City_Rank_Full$HSR_Level[HSR_City_Rank_Full$Year %in% c(2000, 2012, 2024)]))

write.xlsx(HSR_City_Rank_Full, file.path(dir_tab, "HSR_City_Ranking_Full_2000_2024.xlsx"))
saveRDS(HSR_City_Rank_Full, file.path(dir_res, "HSR_City_Ranking_Full_2000_2024.rds"))

# --- 8d. Rank stability analysis ---

Rank_Stability <- HSR_City_Rank_Full %>%
  group_by(CityCode, City, Province, ProvinceCode) %>%
  summarise(
    Rank_Mean    = round(mean(Rank), 1),
    Rank_Median  = median(Rank),
    Rank_SD      = round(sd(Rank), 2),
    Rank_Min     = min(Rank),
    Rank_Max     = max(Rank),
    Rank_Range   = max(Rank) - min(Rank),
    Max_YoY_Change = max(abs(diff(Rank[order(Year)]))),
    Best_Year    = Year[which.min(Rank)],
    Worst_Year   = Year[which.max(Rank)],
    .groups = "drop"
  ) %>%
  mutate(Stability_Category = factor(case_when(
    Rank_SD <= quantile(Rank_SD, 0.33) ~ "Stable",
    Rank_SD <= quantile(Rank_SD, 0.67) ~ "Moderate",
    TRUE ~ "Volatile"),
    levels = c("Stable", "Moderate", "Volatile"))) %>%
  arrange(Rank_Mean)

write.xlsx(Rank_Stability, file.path(dir_tab, "HSR_City_Rank_Stability_Summary.xlsx"))

# Spaghetti plot
Rank_Traj <- HSR_City_Rank_Full %>%
  left_join(Rank_Stability %>% select(CityCode, Stability_Category, Rank_SD),
            by = "CityCode")

p_traj <- ggplot(Rank_Traj, aes(x = Year, y = Rank, group = CityCode,
                                color = Stability_Category)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  scale_y_reverse() +
  scale_color_manual(values = c("Stable" = "#2166AC", "Moderate" = "#F4A582", "Volatile" = "#B2182B"),
                     name = "Rank stability") +
  labs(x = "Year", y = "Rank (1 = best)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

save_fig(p_traj, "Fig_Rank_Trajectory_Spaghetti", 10, 7)
save_fig_data(Rank_Traj %>% select(Year, CityCode, City, Rank, Stability_Category),
              "Fig_Rank_Trajectory_Spaghetti")

# Scatter: Mean rank vs SD
p_scatter <- ggplot(Rank_Stability, aes(x = Rank_Mean, y = Rank_SD,
                                        color = Stability_Category)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("Stable" = "#2166AC", "Moderate" = "#F4A582", "Volatile" = "#B2182B"),
                     name = "Stability") +
  labs(x = "Mean rank (2000-2024)", y = "Rank standard deviation") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

save_fig(p_scatter, "Fig_Rank_MeanVsSD_Scatter", 9, 6)
save_fig_data(Rank_Stability %>% select(CityCode, City, Rank_Mean, Rank_SD, Stability_Category),
              "Fig_Rank_MeanVsSD_Scatter")


# ===================== A9. Annual HSR spatial distribution maps ==============
#
# Purple color scheme; one map per year (2000-2024).
# Each city polygon is filled by its HSR index value for that year.
# Maps saved to dir_fig_maps subfolder.
# ============================================================================

cat("\n=== A9: Generating annual HSR maps (purple scheme) ===\n")

# --- Purple color palette for HSR (low = light lavender, high = deep purple) ---
pal_purple <- colorRampPalette(
  c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC",
    "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D")
)

# --- Prepare city shapefile with CityCode as character for matching ---
city_map_sf <- city_sf %>%
  mutate(CityCode = as.character(市代码))

# --- Verify matching between HSR data and shapefile ---
hsr_codes <- unique(HSR_City_Panel$CityCode)
map_codes <- unique(city_map_sf$CityCode)
matched   <- intersect(hsr_codes, map_codes)
unmatched_hsr <- setdiff(hsr_codes, map_codes)

cat(sprintf("  HSR cities: %d | Map cities: %d | Matched: %d | Unmatched: %d\n",
            length(hsr_codes), length(map_codes), length(matched), length(unmatched_hsr)))

if (length(unmatched_hsr) > 0 && length(unmatched_hsr) <= 20) {
  cat("  Unmatched CityCode examples:", paste(head(unmatched_hsr, 10), collapse = ", "), "\n")
}

# --- Global HSR range for consistent color scale across all years ---
hsr_range <- range(HSR_City_Panel$HSR, na.rm = TRUE)
cat(sprintf("  HSR range: [%.4f, %.4f]\n", hsr_range[1], hsr_range[2]))

# --- Generate one map per year ---
for (yr in Years_Hist) {
  
  tryCatch({
    # Get HSR data for this year
    hsr_yr <- HSR_City_Panel %>%
      filter(Year == yr) %>%
      select(CityCode, HSR, HSR_raw) %>%
      # Also get HSR Level
      left_join(
        HSR_City_Rank_Full %>% filter(Year == yr) %>% select(CityCode, HSR_Level),
        by = "CityCode"
      )
    
    # Merge with shapefile
    city_map_yr <- city_map_sf %>%
      left_join(hsr_yr, by = "CityCode")
    
    # Create map
    p_map <- ggplot() +
      # County boundaries as faint background
      geom_sf(data = county_sf, fill = "grey97", colour = COL_CNTY,
              linewidth = LW_CNTY) +
      # City polygons filled by HSR
      geom_sf(data = city_map_yr, aes(fill = HSR),
              colour = COL_CITY, linewidth = LW_CITY) +
      # Color scale: purple gradient
      scale_fill_gradientn(
        colours = pal_purple(256),
        name = "HSR\nindex",
        limits = hsr_range,
        na.value = "grey90",
        breaks = seq(0, 1, by = 0.2),
        labels = sprintf("%.1f", seq(0, 1, by = 0.2))
      ) +
      # National boundary overlay
      {if (!is.null(china_land))
        geom_sf(data = china_land, fill = NA, colour = COL_NATL,
                linewidth = LW_NATL)
      } +
      # Map extent and coordinates
      coord_sf(xlim = CHINA_LON, ylim = CHINA_LAT, datum = st_crs(4326)) +
      scale_x_continuous(breaks = seq(80, 130, by = 10), labels = fmt_lon) +
      scale_y_continuous(breaks = seq(20, 50, by = 10), labels = fmt_lat) +
      # Labels
      labs(
        title = sprintf("Health System Resilience (HSR) Index — %d", yr),
        subtitle = sprintf("City-level composite index (hierarchical PCA, %d pillars)",
                           n_pillars),
        caption = "Source: National Bureau of Statistics & National Health Commission",
        x = "Longitude", y = "Latitude"
      ) +
      theme_map()
    
    # Save map
    save_map(p_map, sprintf("HSR_Map_%d", yr), w = 11, h = 9)
    
    if (yr %% 5 == 0 || yr == Year_Start || yr == Year_End) {
      cat(sprintf("  HSR map %d saved\n", yr))
    }
    
  }, error = function(e) {
    cat(sprintf("  Map %d error: %s\n", yr, e$message))
  })
}

cat("  All annual HSR maps complete.\n")


# --- Additional 2024 HSR Map with year-specific color scale ---
tryCatch({
  # Get HSR data for 2024
  hsr_2024 <- HSR_City_Panel %>%
    filter(Year == 2024) %>%
    select(CityCode, HSR, HSR_raw) %>%
    left_join(
      HSR_City_Rank_Full %>% filter(Year == 2024) %>% select(CityCode, HSR_Level),
      by = "CityCode"
    )
  
  # Compute 2024-specific range for color scale
  hsr_range_2024 <- range(hsr_2024$HSR, na.rm = TRUE)
  hsr_min_2024 <- floor(hsr_range_2024[1] * 10) / 10   # round down to 0.1
  hsr_max_2024 <- ceiling(hsr_range_2024[2] * 10) / 10  # round up to 0.1
  
  # Merge with shapefile
  city_map_2024 <- city_map_sf %>%
    left_join(hsr_2024, by = "CityCode")
  
  # Create map with 2024-specific color limits
  p_map_2024 <- ggplot() +
    geom_sf(data = county_sf, fill = "grey97", colour = COL_CNTY,
            linewidth = LW_CNTY) +
    geom_sf(data = city_map_2024, aes(fill = HSR),
            colour = COL_CITY, linewidth = LW_CITY) +
    scale_fill_gradientn(
      colours = pal_purple(256),
      name = "Health system\nresilience index",
      limits = c(hsr_min_2024, hsr_max_2024),
      na.value = "grey96",
      breaks = pretty(c(hsr_min_2024, hsr_max_2024), n = 5),
      labels = sprintf("%.2f", pretty(c(hsr_min_2024, hsr_max_2024), n = 5))
    ) +
    {if (!is.null(china_land))
      geom_sf(data = china_land, fill = NA, colour = COL_NATL,
              linewidth = LW_NATL)
    } +
    coord_sf(xlim = CHINA_LON, ylim = CHINA_LAT, datum = st_crs(4326)) +
    scale_x_continuous(breaks = seq(80, 130, by = 10), labels = fmt_lon) +
    scale_y_continuous(breaks = seq(20, 50, by = 10), labels = fmt_lat) +
    labs(
      #title = "Health System Resilience (HSR) Index — 2024",
      #subtitle = sprintf(
      #  "City-level composite index (hierarchical PCA, %d pillars)\nColor scaled to 2024 range [%.3f, %.3f]",
      #  n_pillars, hsr_range_2024[1], hsr_range_2024[2]
      #),
      #caption = "Source: National Bureau of Statistics & National Health Commission",
      x = "Longitude", y = "Latitude"
    ) +
    theme_map()
  
  # Save with distinct filename
  save_map(p_map_2024, "HSR_Map_2024_YearSpecificScale", w = 7, h = 5)  # 11, 9
  cat("  HSR map 2024 (year-specific scale) saved\n")
  
}, error = function(e) {
  cat(sprintf("  Map 2024 (year-specific scale) error: %s\n", e$message))
})



# --- A9b. Faceted small-multiples map (selected years) ---

cat("  Generating small-multiples overview map...\n")

selected_years <- c(2000, 2004, 2008, 2012, 2016, 2020, 2024)

hsr_selected <- HSR_City_Panel %>%
  filter(Year %in% selected_years) %>%
  select(CityCode, Year, HSR)

city_map_selected <- city_map_sf %>%
  # Cross join with selected years, then merge HSR
  crossing(Year = selected_years) %>%
  left_join(hsr_selected, by = c("CityCode", "Year"))

tryCatch({
  p_facet_map <- ggplot() +
    geom_sf(data = county_sf, fill = "grey97", colour = COL_CNTY,
            linewidth = 0.02) +
    geom_sf(data = city_map_selected, aes(fill = HSR),
            colour = COL_CITY, linewidth = 0.04) +
    scale_fill_gradientn(
      colours = pal_purple(256),
      name = "HSR\nindex",
      limits = hsr_range,
      na.value = "grey90"
    ) +
    {if (!is.null(china_land))
      geom_sf(data = china_land, fill = NA, colour = COL_NATL,
              linewidth = 0.2)
    } +
    facet_wrap(~ Year, ncol = 4) +
    coord_sf(xlim = CHINA_LON, ylim = CHINA_LAT, datum = st_crs(4326)) +
    labs(
      title = "Spatial evolution of Health System Resilience (HSR)",
      subtitle = "City-level HSR index, selected years",
      x = NULL, y = NULL
    ) +
    theme_map(base_size = 8) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  save_map(p_facet_map, "HSR_Map_Faceted_SelectedYears", w = 16, h = 10)
  cat("  Faceted overview map saved\n")
}, error = function(e) cat("  Faceted map error:", e$message, "\n"))

# --- A9c. HSR Level (categorical) map for latest year ---

cat("  Generating HSR Level categorical map (2024)...\n")

tryCatch({
  hsr_2024_level <- HSR_City_Rank_Full %>%
    filter(Year == 2024) %>%
    select(CityCode, HSR_Level)
  
  city_map_level <- city_map_sf %>%
    left_join(hsr_2024_level, by = "CityCode")
  
  level_colors <- c("High" = "#3F007D", "Moderate" = "#807DBA", "Low" = "#DADAEB")
  
  p_level_map <- ggplot() +
    geom_sf(data = county_sf, fill = "grey97", colour = COL_CNTY,
            linewidth = LW_CNTY) +
    geom_sf(data = city_map_level, aes(fill = HSR_Level),
            colour = COL_CITY, linewidth = LW_CITY) +
    scale_fill_manual(
      values = level_colors,
      name = "HSR Level",
      na.value = "grey90"
    ) +
    {if (!is.null(china_land))
      geom_sf(data = china_land, fill = NA, colour = COL_NATL,
              linewidth = LW_NATL)
    } +
    coord_sf(xlim = CHINA_LON, ylim = CHINA_LAT, datum = st_crs(4326)) +
    scale_x_continuous(breaks = seq(80, 130, by = 10), labels = fmt_lon) +
    scale_y_continuous(breaks = seq(20, 50, by = 10), labels = fmt_lat) +
    labs(
      title = "Health System Resilience Level Classification — 2024",
      subtitle = "Tertile-based classification (High / Moderate / Low)",
      caption = "Classification: within-year tertiles of overall HSR index",
      x = "Longitude", y = "Latitude"
    ) +
    theme_map()
  
  save_map(p_level_map, "HSR_Map_Level_2024", w = 11, h = 9)
  cat("  HSR Level map (2024) saved\n")
}, error = function(e) cat("  Level map error:", e$message, "\n"))

# --- A9d. Pillar-specific maps for 2024 (optional, informative) ---

cat("  Generating pillar-specific maps for 2024...\n")

for (pname in names(Pillar_Indicators)) {
  tryCatch({
    pillar_2024 <- HSR_City_Panel %>%
      filter(Year == 2024) %>%
      select(CityCode, PillarScore = all_of(pname))
    
    city_map_pillar <- city_map_sf %>%
      left_join(pillar_2024, by = "CityCode")
    
    p_pillar_map <- ggplot() +
      geom_sf(data = county_sf, fill = "grey97", colour = COL_CNTY,
              linewidth = LW_CNTY) +
      geom_sf(data = city_map_pillar, aes(fill = PillarScore),
              colour = COL_CITY, linewidth = LW_CITY) +
      scale_fill_gradientn(
        colours = pal_purple(256),
        name = "Pillar\nscore",
        limits = c(0, 1),
        na.value = "grey90"
      ) +
      {if (!is.null(china_land))
        geom_sf(data = china_land, fill = NA, colour = COL_NATL,
                linewidth = LW_NATL)
      } +
      coord_sf(xlim = CHINA_LON, ylim = CHINA_LAT, datum = st_crs(4326)) +
      scale_x_continuous(breaks = seq(80, 130, by = 10), labels = fmt_lon) +
      scale_y_continuous(breaks = seq(20, 50, by = 10), labels = fmt_lat) +
      labs(
        title = sprintf("HSR Pillar: %s — 2024", gsub("_", " ", pname)),
        subtitle = "Normalized PC1 score (0 = lowest, 1 = highest)",
        x = "Longitude", y = "Latitude"
      ) +
      theme_map()
    
    save_map(p_pillar_map, sprintf("HSR_Map_Pillar_%s_2024", pname), w = 11, h = 9)
    cat(sprintf("    Pillar map: %s saved\n", pname))
  }, error = function(e) {
    cat(sprintf("    Pillar map %s error: %s\n", pname, e$message))
  })
}


cat("\n===== Section A complete: Hierarchical PCA HSR + Maps + Enriched Ranking =====\n")
cat(sprintf("  Output folder (maps): %s\n", dir_fig_maps))
cat(sprintf("  HSR_City_Rank_Full columns: %s\n",
            paste(colnames(HSR_City_Rank_Full), collapse = ", ")))
cat(sprintf("  Total map files generated: %d annual + faceted + level + %d pillars\n",
            length(Years_Hist), length(Pillar_Indicators)))







