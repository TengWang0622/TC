####################################################################################
#                          TROPICAL CYCLONES AND HEALTH
#
#  CMIP6 SSP → 中国区县/城市人口 (年龄组) 提取
#  时间范围: 2025–2050 (台风未来预测研究)
#
####################################################################################

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(exactextractr)
library(ncdf4)
library(openxlsx)
library(ggplot2)
library(scales)
options(scipen = 6)

setwd("C:/Project/Humid")

# ==============================================================
#  0.  全局设定：输出目录 & 可视化参数 (仅在此处修改)
# ==============================================================

# ---- 0a. 输出目录 ----
out_root         <- "Output/Tropical cyclone"
out_extract      <- file.path(out_root, "Extract")
out_data         <- file.path(out_root, "Data")
out_fig          <- file.path(out_root, "Figures")
out_fig_a        <- file.path(out_fig, "Fig_a_stacked_area")
out_fig_b        <- file.path(out_fig, "Fig_b_composition")
out_fig_c        <- file.path(out_fig, "Fig_c_aging_trend")
out_fig_d        <- file.path(out_fig, "Fig_d_aging_map")
out_fig_d_county <- file.path(out_fig_d, "County")
out_fig_d_city   <- file.path(out_fig_d, "City")

for (d in c(out_root, out_extract, out_data, out_fig,
            out_fig_a, out_fig_b, out_fig_c, out_fig_d,
            out_fig_d_county, out_fig_d_city)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 0b. 可视化时间范围 ----
vis_year_min <- 2025
vis_year_max <- 2050

# ---- 0c. 地图年份 ----
map_years_general <- c(2025, 2035, 2050)
map_years_ssp585  <- c(2025, 2030, 2035, 2040, 2045, 2050)

# ---- 0d. 图片分辨率 ----
map_dpi <- 600
fig_dpi <- 300


# ==============================================================
#  A.  数 据 提 取
# ==============================================================

# ---- A1. SSP 文件 ----
ssp_files <- c(
  "SSP1-2.6" = "Src/Pop_CMIP6/grid_pop_age_gender_ssp126.nc",
  "SSP2-4.5" = "Src/Pop_CMIP6/grid_pop_age_gender_ssp245.nc",
  "SSP3-7.0" = "Src/Pop_CMIP6/grid_pop_age_gender_ssp370.nc",
  "SSP5-8.5" = "Src/Pop_CMIP6/grid_pop_age_gender_ssp585.nc"
)

# ---- A2. 坐标 ----
lon   <- seq(73.75, 135.25, by = 0.5)
lat   <- seq(16.25, 53.75,  by = 0.5)
years <- 2010:2100

# ---- A3. 区县地图 ----
map <- st_read("Src/ChinaMap/县.shp")
map <- st_transform(map, "+proj=longlat +datum=WGS84")
map <- map %>% filter(!市 %in% c("台湾省", "澳门特别行政区", "香港特别行政区"))

geo_info <- data.frame(
  Province   = map$省,
  City       = map$市,
  County     = map$县,
  CountyCode = map$县代码,
  CityCode   = map$市代码,
  stringsAsFactors = FALSE
)

# ---- A3b. 城市地图 (区县融合) ----
cat("构建城市级地图 (dissolve) ...\n")
map_city <- map %>%
  group_by(省, 市, 市代码) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  rename(Province = 省, City = 市, CityCode = 市代码) %>%
  st_make_valid()

# ---- A4. 年龄组 ----
age_groups <- list(
  "0-19"  = 1:4,
  "20-44" = 5:9,
  "45-64" = 10:13,
  "65+"   = 14:21
)

# ---- A5. 栅格函数 ----
make_raster <- function(mat, lon, lat) {
  res <- 0.5
  raster(mat,
         xmn = min(lon) - res/2, xmx = max(lon) + res/2,
         ymn = min(lat) - res/2, ymx = max(lat) + res/2,
         crs = "+proj=longlat +datum=WGS84")
}

# ---- A6. 提取主循环 ----
all_results      <- list()
all_results_city <- list()                                     # ★ NEW

for (ssp in names(ssp_files)) {
  
  cat("\n===", ssp, "===\n")
  
  nc <- nc_open(ssp_files[ssp])
  pop_all <- ncvar_get(nc, "grid_pop_age_gender")
  nc_close(nc)
  
  pop_both <- pop_all[,,,, 1] + pop_all[,,,, 2]
  rm(pop_all); gc()
  
  pop_ag <- list()
  for (ag in names(age_groups)) {
    pop_ag[[ag]] <- apply(pop_both[,, age_groups[[ag]], , drop = FALSE],
                          c(1, 2, 4), sum, na.rm = TRUE)
  }
  rm(pop_both); gc()
  
  all_list <- vector("list", length(years) * length(age_groups))
  counter  <- 1
  
  for (i in seq_along(years)) {
    cat("\r  ", years[i])
    s <- stack(lapply(names(age_groups), function(ag) {
      mat <- pop_ag[[ag]][,, i]; mat[is.na(mat)] <- 0
      make_raster(mat, lon, lat)
    }))
    ext_df <- exact_extract(s, map, fun = "sum")
    
    for (k in seq_along(names(age_groups))) {
      ag <- names(age_groups)[k]
      all_list[[counter]] <- data.frame(
        geo_info, Year = years[i], Agegroup = ag,
        Population = round(ext_df[, k], 0),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
  }
  
  result_long <- do.call(rbind, all_list)
  all_results[[ssp]] <- result_long
  
  tag <- gsub("[^A-Za-z0-9]", "", ssp)
  
  # ---- County: long & wide ----
  write.xlsx(result_long, file.path(out_extract, paste0("Pop_County_", tag, "_long.xlsx")))
  
  result_wide <- result_long %>%
    mutate(Col = paste0(Agegroup, "_", Year)) %>%
    pivot_wider(id_cols = c(Province, City, County, CountyCode, CityCode),
                names_from = Col, values_from = Population)
  write.xlsx(result_wide, file.path(out_extract, paste0("Pop_County_", tag, "_wide.xlsx")))
  
  # ---- ★ NEW: City-level aggregation + save ----
  result_city_long <- result_long %>%
    group_by(Province, City, CityCode, Year, Agegroup) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  
  all_results_city[[ssp]] <- result_city_long
  
  write.xlsx(result_city_long,
             file.path(out_extract, paste0("Pop_City_", tag, "_long.xlsx")))
  
  result_city_wide <- result_city_long %>%
    mutate(Col = paste0(Agegroup, "_", Year)) %>%
    pivot_wider(id_cols = c(Province, City, CityCode),
                names_from = Col, values_from = Population)
  write.xlsx(result_city_wide,
             file.path(out_extract, paste0("Pop_City_", tag, "_wide.xlsx")))
  # ---- ★ END NEW ----
  
  cat("\n  saved (County + City)\n")
  rm(pop_ag, all_list, result_long, result_wide,
     result_city_long, result_city_wide); gc()                 # ★ MODIFIED: 加入 city 对象
}

full_data      <- bind_rows(all_results,      .id = "SSP")
full_data_city <- bind_rows(all_results_city,  .id = "SSP")    # ★ NEW
rm(all_results, all_results_city); gc()                         # ★ MODIFIED
cat("\n== 全部提取完成 (County + City) ==\n")


# ==============================================================
#  B.  精 美 可 视 化  (2025–2050)
# ==============================================================

# ---- B0. 全局设定 ----

ssp_levels <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
age_stack  <- c("65+", "45-64", "20-44", "0-19")

full_data$SSP      <- factor(full_data$SSP,      levels = ssp_levels)
full_data$Agegroup <- factor(full_data$Agegroup,  levels = age_stack)

full_data_city$SSP      <- factor(full_data_city$SSP,      levels = ssp_levels)    # ★ NEW
full_data_city$Agegroup <- factor(full_data_city$Agegroup,  levels = age_stack)    # ★ NEW

# Lancet 配色：年龄组
age_pal <- c("0-19"  = "#0072B5",
             "20-44" = "#20854E",
             "45-64" = "#E18727",
             "65+"   = "#BC3C29")

# IPCC 配色：SSP
ssp_pal <- c("SSP1-2.6" = "#2166AC",
             "SSP2-4.5" = "#4DAF4A",
             "SSP3-7.0" = "#FF7F00",
             "SSP5-8.5" = "#D62728")

# Nature/Lancet 主题
theme_pub <- function(bs = 13) {
  theme_minimal(base_size = bs) %+replace%
    theme(
      text              = element_text(colour = "grey10"),
      plot.title        = element_text(face = "bold", size = bs + 4, hjust = 0,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(size = bs + 0.5, colour = "grey30",
                                       margin = margin(b = 12)),
      plot.caption      = element_text(size = bs - 2.5, colour = "grey50",
                                       hjust = 1, margin = margin(t = 10)),
      axis.title        = element_text(face = "bold", size = bs),
      axis.text         = element_text(size = bs - 1.5, colour = "grey20"),
      axis.line.x       = element_line(colour = "grey40", linewidth = 0.35),
      axis.ticks.x      = element_line(colour = "grey40", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.title      = element_text(face = "bold", size = bs - 0.5),
      legend.text       = element_text(size = bs - 1.5),
      legend.position   = "bottom",
      strip.text        = element_text(face = "bold", size = bs + 0.5),
      strip.background  = element_rect(fill = "grey96", colour = NA),
      plot.margin       = margin(14, 18, 10, 14)
    )
}

# ---- 全国汇总 (County → National) ----
nat <- full_data %>%
  filter(Year >= vis_year_min, Year <= vis_year_max) %>%
  group_by(SSP, Year, Agegroup) %>%
  summarise(Pop = sum(Population, na.rm = TRUE), .groups = "drop")

nat_pct <- nat %>%
  group_by(SSP, Year) %>%
  mutate(Pct = Pop / sum(Pop) * 100) %>%
  ungroup()

# ---- ★ NEW: 城市级汇总 (用于图 a/b/c 城市版) ----
nat_city <- full_data_city %>%
  filter(Year >= vis_year_min, Year <= vis_year_max) %>%
  group_by(SSP, Year, Province, City, CityCode, Agegroup) %>%
  summarise(Pop = sum(Population, na.rm = TRUE), .groups = "drop")

nat_city_pct <- nat_city %>%
  group_by(SSP, Year, Province, City, CityCode) %>%
  mutate(Total = sum(Pop),
         Pct   = Pop / Total * 100) %>%
  ungroup()
# ---- ★ END NEW ----


# ==========================================================
#  图 (a)  堆积面积图 — 每个 SSP 单独保存 (全国)
# ==========================================================

write.xlsx(nat %>% mutate(Pop_billion = Pop / 1e9),
           file.path(out_data, "Fig_a_stacked_area_data.xlsx"))

# ★ NEW: 保存城市级堆积数据
write.xlsx(nat_city %>% mutate(Pop_million = Pop / 1e6),
           file.path(out_data, "Fig_a_stacked_area_city_data.xlsx"))

for (ssp_i in ssp_levels) {
  
  ssp_tag <- gsub("[^A-Za-z0-9]", "", ssp_i)
  nat_sub <- nat %>% filter(SSP == ssp_i)
  
  pa <- ggplot(nat_sub, aes(x = Year, y = Pop / 1e9, fill = Agegroup)) +
    geom_area(alpha = 0.88, linewidth = 0.2, colour = "white") +
    scale_fill_manual(
      values = age_pal, name = "Age group",
      guide  = guide_legend(reverse = TRUE)
    ) +
    scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
                       minor_breaks = seq(vis_year_min, vis_year_max, 1),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 1.6),
                       breaks = seq(0, 1.6, 0.2),
                       labels = comma_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(title = NULL, subtitle = NULL,
         x = "Year", y = "Population (billion)") +
    theme_pub()
  
  ggsave(file.path(out_fig_a, paste0("Fig_a_", ssp_tag, ".png")), pa,
         width = 7, height = 5.5, dpi = fig_dpi, bg = "white")
}

# ★ NEW: 城市级堆积面积图 (每个 SSP × 每个城市)
#out_fig_a_city <- file.path(out_fig_a, "City")
#dir.create(out_fig_a_city, showWarnings = FALSE, recursive = TRUE)

#for (ssp_i in ssp_levels) {

#  ssp_tag  <- gsub("[^A-Za-z0-9]", "", ssp_i)
#  city_sub <- nat_city %>% filter(SSP == ssp_i)
#  cities   <- city_sub %>% distinct(Province, City, CityCode)

#  for (r in seq_len(nrow(cities))) {
#    ct_code <- cities$CityCode[r]
#    ct_name <- cities$City[r]
#    ct_tag  <- paste0(ct_code, "_", gsub("[/\\\\: ]", "", ct_name))

#    d <- city_sub %>% filter(CityCode == ct_code)

#    y_max_val <- d %>%
#      group_by(Year) %>% summarise(tot = sum(Pop), .groups = "drop") %>%
#      pull(tot) %>% max(na.rm = TRUE) / 1e6
#    y_upper <- ceiling(y_max_val / 2) * 2   # 取偶数上限
#    if (y_upper < 1) y_upper <- 1

#    pa_city <- ggplot(d, aes(x = Year, y = Pop / 1e6, fill = Agegroup)) +
#      geom_area(alpha = 0.88, linewidth = 0.2, colour = "white") +
#      scale_fill_manual(values = age_pal, name = "Age group",
#                        guide = guide_legend(reverse = TRUE)) +
#      scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
#                         minor_breaks = seq(vis_year_min, vis_year_max, 1),
#                         expand = expansion(mult = c(0.01, 0.01))) +
#      scale_y_continuous(labels = comma_format(accuracy = 0.1),
#                         expand = expansion(mult = c(0, 0.03))) +
#      labs(title = paste0(ct_name, "  |  ", ssp_i),
#           x = "Year", y = "Population (million)") +
#      theme_pub()

#    ggsave(file.path(out_fig_a_city,
#                     paste0("Fig_a_city_", ssp_tag, "_", ct_tag, ".png")),
#           pa_city, width = 7, height = 5.5, dpi = fig_dpi, bg = "white")
#  }
#}
# ★ END NEW

#cat("Fig a saved (National + City)\n")


# ==========================================================
#  图 (b)  百分比堆积面积图 — 每个 SSP 单独保存 (全国)
# ==========================================================

write.xlsx(nat_pct, file.path(out_data, "Fig_b_composition_data.xlsx"))

# ★ NEW: 保存城市级百分比数据
write.xlsx(nat_city_pct %>% select(SSP, Year, Province, City, CityCode,
                                   Agegroup, Pop, Total, Pct),
           file.path(out_data, "Fig_b_composition_city_data.xlsx"))

for (ssp_i in ssp_levels) {
  
  ssp_tag <- gsub("[^A-Za-z0-9]", "", ssp_i)
  pct_sub <- nat_pct %>% filter(SSP == ssp_i)
  
  pb <- ggplot(pct_sub, aes(x = Year, y = Pct, fill = Agegroup)) +
    geom_area(alpha = 0.88, linewidth = 0.2, colour = "white") +
    scale_fill_manual(
      values = age_pal, name = "Age group",
      guide  = guide_legend(reverse = TRUE)
    ) +
    scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
                       minor_breaks = seq(vis_year_min, vis_year_max, 1),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 10),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(title = NULL, subtitle = NULL,
         x = "Year", y = "Proportion (%)") +
    theme_pub()
  
  ggsave(file.path(out_fig_b, paste0("Fig_b_", ssp_tag, ".png")), pb,
         width = 7, height = 5.5, dpi = fig_dpi, bg = "white")
}

# ★ NEW: 城市级百分比堆积面积图 (每个 SSP × 每个城市)
#out_fig_b_city <- file.path(out_fig_b, "City")
#dir.create(out_fig_b_city, showWarnings = FALSE, recursive = TRUE)

#for (ssp_i in ssp_levels) {

#  ssp_tag  <- gsub("[^A-Za-z0-9]", "", ssp_i)
#  pct_sub  <- nat_city_pct %>% filter(SSP == ssp_i)
#  cities   <- pct_sub %>% distinct(Province, City, CityCode)

#  for (r in seq_len(nrow(cities))) {
#    ct_code <- cities$CityCode[r]
#    ct_name <- cities$City[r]
#    ct_tag  <- paste0(ct_code, "_", gsub("[/\\\\: ]", "", ct_name))

#    d <- pct_sub %>% filter(CityCode == ct_code)

#    pb_city <- ggplot(d, aes(x = Year, y = Pct, fill = Agegroup)) +
#      geom_area(alpha = 0.88, linewidth = 0.2, colour = "white") +
#      scale_fill_manual(values = age_pal, name = "Age group",
#                        guide = guide_legend(reverse = TRUE)) +
#      scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
#                         minor_breaks = seq(vis_year_min, vis_year_max, 1),
#                         expand = expansion(mult = c(0.01, 0.01))) +
#      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10),
#                         labels = function(x) paste0(x, "%"),
#                         expand = expansion(mult = c(0, 0.03))) +
#      labs(title = paste0(ct_name, "  |  ", ssp_i),
#           x = "Year", y = "Proportion (%)") +
#      theme_pub()
#    
#    ggsave(file.path(out_fig_b_city,
#                     paste0("Fig_b_city_", ssp_tag, "_", ct_tag, ".png")),
#           pb_city, width = 7, height = 5.5, dpi = fig_dpi, bg = "white")
#  }
#}
# ★ END NEW

#cat("Fig b saved (National + City)\n")


# ==========================================================
#  图 (c)  全国老龄化率折线图
# ==========================================================

aging <- full_data %>%
  filter(Year >= vis_year_min, Year <= vis_year_max) %>%
  group_by(SSP, Year) %>%
  summarise(Pop65 = sum(Population[Agegroup == "65+"], na.rm = TRUE),
            Total = sum(Population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(AgingPct = Pop65 / Total * 100)

write.xlsx(aging %>% select(SSP, Year, Pop65, Total, AgingPct),
           file.path(out_data, "Fig_c_aging_trend_data.xlsx"))

# ★ NEW: 城市级老龄化趋势数据
aging_city <- full_data_city %>%
  filter(Year >= vis_year_min, Year <= vis_year_max) %>%
  group_by(SSP, Year, Province, City, CityCode) %>%
  summarise(Pop65 = sum(Population[Agegroup == "65+"], na.rm = TRUE),
            Total = sum(Population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(AgingPct = ifelse(Total > 0, Pop65 / Total * 100, NA))

write.xlsx(aging_city,
           file.path(out_data, "Fig_c_aging_trend_city_data.xlsx"))
# ★ END NEW

key_years_all <- c(vis_year_min, 2035, vis_year_max)
key_pts <- aging %>% filter(Year %in% key_years_all)

# ---- (c-1) 全部 SSP (全国) ----
pc <- ggplot(aging, aes(Year, AgingPct, colour = SSP)) +
  geom_hline(yintercept = 14, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 21, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  annotate("text", x = vis_year_max - 1, y = 14.8,
           label = "Aged society (14 %)", size = 3.4,
           colour = "grey45", hjust = 1, fontface = "italic") +
  annotate("text", x = vis_year_max - 1, y = 21.8,
           label = "Super-aged society (21 %)", size = 3.4,
           colour = "grey45", hjust = 1, fontface = "italic") +
  geom_line(linewidth = 1.15) +
  geom_point(data = key_pts, size = 3) +
  geom_text(data = key_pts,
            aes(label = sprintf("%.1f%%", AgingPct)),
            vjust = -1.4, size = 3.5, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = ssp_pal, name = "Scenario") +
  scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.04, 0.14))) +
  labs(title    = "Population Aging Trajectory",
       subtitle = paste0("China, ", vis_year_min, "\u2013", vis_year_max,
                         "  |  Proportion of population aged 65+"),
       caption  = "Dashed lines: UN aging-society thresholds  |  Data: CMIP6 SSP",
       x = NULL, y = "Aging ratio (%)") +
  theme_pub(bs = 14) +
  theme(
    legend.position   = c(0.15, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.88),
                                     colour = "grey75", linewidth = 0.3),
    legend.key.width  = unit(1.5, "cm"),
    panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.25)
  )

ggsave(file.path(out_fig_c, "Fig_c_aging_trend_allSSP.png"), pc,
       width = 11, height = 7.5, dpi = fig_dpi, bg = "white")
ggsave(file.path(out_fig_c, "Fig_c_aging_trend_allSSP.pdf"), pc,
       width = 11, height = 7.5, device = cairo_pdf)

# ---- (c-2) 仅 SSP5-8.5 (全国) ----
aging_585 <- aging %>% filter(SSP == "SSP5-8.5")

key_years_585 <- c(2025, 2030, 2035, 2040, 2045, 2050)
key_585 <- aging_585 %>% filter(Year %in% key_years_585)

pc585 <- ggplot(aging_585, aes(Year, AgingPct)) +
  geom_hline(yintercept = 14, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 21, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  annotate("text", x = vis_year_max - 1, y = 14.8,
           label = "Aged society (14 %)", size = 3.4,
           colour = "grey45", hjust = 1, fontface = "italic") +
  annotate("text", x = vis_year_max - 1, y = 21.8,
           label = "Super-aged society (21 %)", size = 3.4,
           colour = "grey45", hjust = 1, fontface = "italic") +
  geom_line(linewidth = 1.15, colour = ssp_pal["SSP5-8.5"]) +
  geom_point(data = key_585, size = 3, colour = ssp_pal["SSP5-8.5"]) +
  geom_text(data = key_585,
            aes(label = sprintf("%.1f%%", AgingPct)),
            vjust = -1.4, size = 3.5, fontface = "bold",
            colour = ssp_pal["SSP5-8.5"]) +
  scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.04, 0.14))) +
  labs(title    = "Population Aging Trajectory — SSP5-8.5",
       subtitle = paste0("China, ", vis_year_min, "\u2013", vis_year_max,
                         "  |  Proportion of population aged 65+"),
       caption  = "Dashed lines: UN aging-society thresholds  |  Data: CMIP6 SSP",
       x = NULL, y = "Aging ratio (%)") +
  theme_pub(bs = 14) +
  theme(panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.25))

ggsave(file.path(out_fig_c, "Fig_c_aging_trend_SSP585.png"), pc585,
       width = 11, height = 7.5, dpi = fig_dpi, bg = "white")
ggsave(file.path(out_fig_c, "Fig_c_aging_trend_SSP585.pdf"), pc585,
       width = 11, height = 7.5, device = cairo_pdf)

# ★ NEW: (c-3) 城市级老龄化折线 — 每个 SSP × 每个城市
#out_fig_c_city <- file.path(out_fig_c, "City")
#dir.create(out_fig_c_city, showWarnings = FALSE, recursive = TRUE)

#for (ssp_i in ssp_levels) {

#  ssp_tag   <- gsub("[^A-Za-z0-9]", "", ssp_i)
#  aging_sub <- aging_city %>% filter(SSP == ssp_i)
#  cities    <- aging_sub %>% distinct(Province, City, CityCode)

#  for (r in seq_len(nrow(cities))) {
#    ct_code <- cities$CityCode[r]
#    ct_name <- cities$City[r]
#    ct_tag  <- paste0(ct_code, "_", gsub("[/\\\\: ]", "", ct_name))

#    d <- aging_sub %>% filter(CityCode == ct_code)
#    key_d <- d %>% filter(Year %in% key_years_585)

#    y_lo <- floor(min(d$AgingPct, na.rm = TRUE) / 5) * 5
#    y_hi <- ceiling(max(d$AgingPct, na.rm = TRUE) / 5) * 5
#    if (y_hi < 21) y_hi <- max(y_hi, 25)

#    pc_city <- ggplot(d, aes(Year, AgingPct)) +
#      geom_hline(yintercept = 14, linetype = "dashed",
#                 colour = "grey60", linewidth = 0.3) +
#      geom_hline(yintercept = 21, linetype = "dashed",
#                 colour = "grey60", linewidth = 0.3) +
#      geom_line(linewidth = 1.1, colour = ssp_pal[ssp_i]) +
#      geom_point(data = key_d, size = 2.5, colour = ssp_pal[ssp_i]) +
#      geom_text(data = key_d,
#                aes(label = sprintf("%.1f%%", AgingPct)),
#                vjust = -1.3, size = 3.2, fontface = "bold",
#                colour = ssp_pal[ssp_i]) +
#      scale_x_continuous(breaks = seq(vis_year_min, vis_year_max, 5),
#                         expand = expansion(mult = c(0.02, 0.04))) +
#      scale_y_continuous(labels = function(x) paste0(x, "%"),
#                         expand = expansion(mult = c(0.04, 0.14))) +
#      labs(title = paste0(ct_name, "  |  ", ssp_i),
#           subtitle = "Proportion of population aged 65+",
#           x = NULL, y = "Aging ratio (%)") +
#      theme_pub(bs = 13)

#    ggsave(file.path(out_fig_c_city,
#                     paste0("Fig_c_city_", ssp_tag, "_", ct_tag, ".png")),
#           pc_city, width = 8, height = 6, dpi = fig_dpi, bg = "white")
#  }
#}
# ★ END NEW

#cat("Fig c saved (National + City)\n")


# ==========================================================
#  图 (d)  老龄化率空间分布 — 逐 SSP × Year, County + City
# ==========================================================

# ---- 简化地图 ----
map_county_slim <- st_simplify(map, dTolerance = 0.002, preserveTopology = TRUE) %>%
  select(CountyCode = 县代码)

map_city_slim <- st_simplify(map_city, dTolerance = 0.005, preserveTopology = TRUE)

# ---- 红色连续色阶 ----
red_pal <- colorRampPalette(
  c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
    "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
)(256)

# ---- 合并所有需要的年份 ----
all_map_years <- sort(unique(c(map_years_general, map_years_ssp585)))

# ---- 区县级老龄化率 ----
county_ag <- full_data %>%
  filter(Year %in% all_map_years) %>%
  group_by(SSP, Year, Province, City, County, CountyCode, CityCode) %>%
  summarise(Pop65 = sum(Population[Agegroup == "65+"], na.rm = TRUE),
            Total = sum(Population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(AgingPct = ifelse(Total > 0, Pop65 / Total * 100, NA))

write.xlsx(county_ag, file.path(out_data, "Fig_d_aging_county_data.xlsx"))

# ---- 城市级老龄化率 ----
city_ag <- full_data_city %>%                                    # ★ MODIFIED: 改用 full_data_city
  filter(Year %in% all_map_years) %>%
  group_by(SSP, Year, Province, City, CityCode) %>%
  summarise(Pop65 = sum(Population[Agegroup == "65+"], na.rm = TRUE),
            Total = sum(Population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(AgingPct = ifelse(Total > 0, Pop65 / Total * 100, NA))

write.xlsx(city_ag, file.path(out_data, "Fig_d_aging_city_data.xlsx"))

# ---- 绘图函数 (优化版) ----
plot_aging_map <- function(map_data, ssp_name, yr, level_tag,
                           fill_pal, fill_lim = c(0, 50),
                           border_lw = 0.05) {
  
  ggplot(map_data) +
    geom_sf(aes(fill = AgingPct),
            colour = "grey80", linewidth = border_lw) +
    scale_fill_gradientn(
      colours  = fill_pal,
      limits   = fill_lim,
      oob      = squish,
      na.value = "grey85",
      name     = "Aging ratio (%)",
      guide    = guide_colorbar(
        barwidth        = unit(10, "cm"),
        barheight       = unit(0.45, "cm"),
        ticks           = FALSE,
        frame.colour    = "grey40",
        frame.linewidth = 0.3,
        title.position  = "top",
        direction       = "horizontal"
      )
    ) +
    coord_sf(xlim = c(73, 136), ylim = c(16, 55), expand = FALSE) +
    labs(title    = paste0(ssp_name, "  |  ", yr),
         subtitle = NULL) +
    theme_pub(bs = 12) +
    theme(
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(face = "bold", size = 11),
      legend.text      = element_text(size = 10),
      legend.margin    = margin(t = 4, b = 2),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
}

# ---- 逐 SSP × Year 绘图保存 ----
for (ssp_i in ssp_levels) {
  
  ssp_tag <- gsub("[^A-Za-z0-9]", "", ssp_i)
  
  if (ssp_i == "SSP5-8.5") {
    yrs <- map_years_ssp585
  } else {
    yrs <- map_years_general
  }
  
  for (yr in yrs) {
    
    cat("  Map:", ssp_i, yr, "\n")
    
    # --- 区县级 ---
    dat_c <- county_ag %>% filter(SSP == ssp_i, Year == yr)
    map_c <- left_join(map_county_slim, dat_c, by = "CountyCode")
    
    p_county <- plot_aging_map(map_c, ssp_i, yr, "County",
                               fill_pal = red_pal, fill_lim = c(0, 50),
                               border_lw = 0.08)
    
    ggsave(file.path(out_fig_d_county,
                     paste0("Fig_d_county_", ssp_tag, "_", yr, ".png")),
           p_county, width = 8, height = 7, dpi = map_dpi, bg = "white")
    
    # --- 城市级 ---
    dat_ct <- city_ag %>% filter(SSP == ssp_i, Year == yr)
    map_ct <- left_join(map_city_slim, dat_ct, by = "CityCode")
    
    p_city <- plot_aging_map(map_ct, ssp_i, yr, "City",
                             fill_pal = red_pal, fill_lim = c(0, 50),
                             border_lw = 0.15)
    
    ggsave(file.path(out_fig_d_city,
                     paste0("Fig_d_city_", ssp_tag, "_", yr, ".png")),
           p_city, width = 8, height = 7, dpi = map_dpi, bg = "white")
  }
}

cat("Fig d saved\n")
cat("\n== 全部完成 ==\n")


