###############################################################################
#     PART B (V5). Bayesian hierarchical projection (2025-2049) — City level
#
###############################################################################

# ===================== B0. Setup: classification, themes, linetypes =========

HSR_2024_class <- HSR_City_Panel %>%
  filter(Year == Year_End) %>%
  mutate(tertile = ntile(HSR, 3),
         Resilience_Group = factor(case_when(tertile == 3 ~ "High",
                                             tertile == 2 ~ "Moderate",
                                             TRUE ~ "Low"),
                                   levels = c("High", "Moderate", "Low"))) %>%
  dplyr::select(Province, CityCode, HSR_2024 = HSR, Resilience_Group)

write.xlsx(HSR_2024_class,
           file.path(dir_res, "HSR_2024_Tertile_Classification_City.xlsx"))

theme_pub2 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size, color = "grey30"),
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Linetypes: dotted / dot-dash / dashed / solid (No Adaptation -> Aspirational)
scenario_linetypes <- c("No Adaptation" = "dotted",
                        "Autonomous"    = "dotdash",
                        "Targeted"      = "dashed",
                        "Aspirational"  = "solid")

# ===================== B1. Model data =====================

City_Map <- HSR_City_Panel %>%
  distinct(Province, CityCode) %>%
  arrange(CityCode) %>%
  mutate(City_ID = sprintf("C%04d", row_number()))
write.xlsx(City_Map, file.path(dir_res, "City_ID_Mapping.xlsx"))

Prov_Map <- HSR_City_Panel %>%
  distinct(Province) %>%
  arrange(Province) %>%
  mutate(Province_ID = sprintf("P%02d", row_number()))
write.xlsx(Prov_Map, file.path(dir_res, "Province_ID_Mapping.xlsx"))

model_df <- HSR_City_Panel %>%
  left_join(City_Map, by = c("Province", "CityCode")) %>%
  left_join(Prov_Map, by = "Province") %>%
  arrange(City_ID, Year) %>%
  group_by(City_ID) %>%
  mutate(hsr_lag = lag(HSR, 1),
         delta_hsr = HSR - hsr_lag,
         time_index = Year - Year_Start) %>%
  ungroup() %>%
  filter(!is.na(delta_hsr))

time_mean <- mean(model_df$time_index)
time_sd   <- sd(model_df$time_index)
model_df  <- model_df %>%
  mutate(time_std = (time_index - time_mean) / time_sd)

# ===================== B2. Bayesian hierarchical beta-convergence model =====

model_priors <- c(
  prior(normal(0, 0.1), class = "b"),
  prior(normal(0.005, 0.01), class = "Intercept"),
  prior(exponential(20), class = "sd"),
  prior(exponential(50), class = "sigma")
)

fit_convergence <- brm(
  delta_hsr ~ hsr_lag + time_std + (1 + time_std | Province_ID) + (1 | City_ID),
  data = model_df, family = gaussian(), prior = model_priors,
  chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 2026,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  backend = "rstan"
)

saveRDS(fit_convergence, file.path(dir_res, "brms_fit_convergence_city.rds"))

sink(file.path(dir_res, "brms_model_summary_city.txt"))
print(summary(fit_convergence))
rhat_vals <- brms::rhat(fit_convergence)
cat("\nMax Rhat:", max(rhat_vals, na.rm = TRUE), "\n")
cat("All Rhat < 1.01:", all(rhat_vals < 1.01, na.rm = TRUE), "\n")
sink()

post <- as_draws_df(fit_convergence)
beta_tab <- data.frame(
  Parameter = "beta (hsr_lag, convergence)",
  Mean = mean(post$b_hsr_lag),
  CI_low = quantile(post$b_hsr_lag, 0.025),
  CI_high = quantile(post$b_hsr_lag, 0.975),
  P_beta_lt_0 = mean(post$b_hsr_lag < 0)
)
write.xlsx(beta_tab, file.path(dir_tab, "Beta_Convergence_Posterior_City.xlsx"))

pp <- pp_check(fit_convergence, ndraws = 100) +
  labs(title = "Posterior predictive check (annual \u0394HSR, city-level)",
       subtitle = "Black: observed density (y); blue: 100 replicated datasets (y_rep)") +
  theme_pub2()
save_fig(pp, "FigS_PP_Check_City", 7, 5)

# ===================== B2. Posterior Predictive Check Plot ====================

pp <- pp_check(fit_convergence, ndraws = 100, type = "dens_overlay") +
  scale_x_continuous(
    name = expression(Delta * "HSR (annual change in HSR index)"),
    breaks = pretty
  ) +
  scale_y_continuous(
    name = "Probability density"
  ) +
  #labs(
  #  title = "Posterior Predictive Check: Bayesian Hierarchical Beta-Convergence Model",
  #  subtitle = "Black: observed density of annual \u0394HSR; Blue: 100 simulated datasets from posterior predictive distribution"
  #) +
  theme_pub2() +
  theme(
    legend.position = "right"
  )

save_fig(pp, "FigS_PP_Check_City", 5, 3.5)

# ===================== B3. Forward projection engine ========================
# Trend damping: standardized time index frozen at end-of-training value.
# HSR uncapped above; lower bound 0.

project_hsr_city <- function(post, n_draws, hsr_start, city_ids,
                             city_to_prov, years_proj,
                             time_mean, time_sd, year0 = 2000,
                             t_cap_year = 2024, damp_trend = TRUE,
                             lb = 0, seed = 2026) {
  set.seed(seed)
  draw_ids <- sample(seq_len(nrow(post)), n_draws)
  n_cities <- length(city_ids); n_years <- length(years_proj)
  
  unique_provs <- unique(unname(city_to_prov))
  re_prov_int  <- as.matrix(post[draw_ids,
                                 paste0("r_Province_ID[", unique_provs, ",Intercept]")])
  re_prov_time <- as.matrix(post[draw_ids,
                                 paste0("r_Province_ID[", unique_provs, ",time_std]")])
  colnames(re_prov_int)  <- unique_provs
  colnames(re_prov_time) <- unique_provs
  
  re_city_int <- as.matrix(post[draw_ids,
                                paste0("r_City_ID[", city_ids, ",Intercept]")])
  
  b0    <- post$b_Intercept[draw_ids]
  b_lag <- post$b_hsr_lag[draw_ids]
  b_t   <- post$b_time_std[draw_ids]
  sig   <- post$sigma[draw_ids]
  
  prov_for_city <- city_to_prov[city_ids]
  ts_cap <- (t_cap_year - year0 - time_mean) / time_sd
  
  arr <- array(NA_real_, dim = c(n_draws, n_cities, n_years),
               dimnames = list(NULL, city_ids, as.character(years_proj)))
  
  for (d in seq_len(n_draws)) {
    h <- hsr_start
    for (t in seq_len(n_years)) {
      ts_raw <- (years_proj[t] - year0 - time_mean) / time_sd
      ts <- if (damp_trend) min(ts_raw, ts_cap) else ts_raw
      mu <- b0[d] +
        re_prov_int[d, prov_for_city] +
        re_city_int[d, ] +
        b_lag[d] * h +
        (b_t[d] + re_prov_time[d, prov_for_city]) * ts
      h <- pmax(h + rnorm(n_cities, mu, sig[d]), lb)
      arr[d, , t] <- h
    }
  }
  arr
}

hsr_start_2024 <- HSR_City_Panel %>%
  filter(Year == Year_End) %>%
  left_join(City_Map, by = c("Province", "CityCode")) %>%
  left_join(Prov_Map, by = "Province") %>%
  arrange(City_ID)

city_to_prov <- setNames(hsr_start_2024$Province_ID, hsr_start_2024$City_ID)

N_Draws <- 1000

hsr_proj_auto <- project_hsr_city(
  post, N_Draws,
  hsr_start = hsr_start_2024$HSR,
  city_ids = hsr_start_2024$City_ID,
  city_to_prov = city_to_prov,
  years_proj = Years_Proj,
  time_mean = time_mean, time_sd = time_sd,
  t_cap_year = Year_End, damp_trend = TRUE
)
saveRDS(hsr_proj_auto, file.path(dir_res, "HSR_Projection_Autonomous_CityDraws.rds"))

hsr_proj_auto_undamped <- project_hsr_city(
  post, N_Draws,
  hsr_start = hsr_start_2024$HSR,
  city_ids = hsr_start_2024$City_ID,
  city_to_prov = city_to_prov,
  years_proj = Years_Proj,
  time_mean = time_mean, time_sd = time_sd,
  t_cap_year = Year_End, damp_trend = FALSE
)

qs <- c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975)
summarize_proj <- function(arr, base_df, years_proj) {
  map_dfr(seq_along(years_proj), function(t) {
    map_dfr(seq_len(nrow(base_df)), function(i) {
      v <- arr[, i, t]
      tibble(City_ID = base_df$City_ID[i],
             Year = years_proj[t],
             hsr_mean = mean(v),
             !!!setNames(as.list(quantile(v, qs)),
                         c("q025","q10","q25","hsr_median","q75","q90","q975")))
    })
  })
}

proj_auto_summary <- summarize_proj(hsr_proj_auto, hsr_start_2024, Years_Proj) %>%
  left_join(City_Map, by = "City_ID") %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode")

proj_auto_summary_undamped <-
  summarize_proj(hsr_proj_auto_undamped, hsr_start_2024, Years_Proj) %>%
  left_join(City_Map, by = "City_ID") %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode")

write.xlsx(proj_auto_summary,
           file.path(dir_res, "HSR_Projection_Autonomous_City_Summary.xlsx"))
write.xlsx(proj_auto_summary_undamped,
           file.path(dir_res, "HSR_Projection_Autonomous_City_Summary_UndampedTrend.xlsx"))

# ===================== B4. Adaptation scenarios ==============================
# DISTRIBUTIONAL targets via rank-preserving quantile mapping, delivered
# through a smooth additive wedge on top of the autonomous trajectory.
#
#   Targeted:     Low-group cities mapped onto the Moderate group's
#                 autonomous-2049 distribution.
#   Aspirational: Low + Moderate cities mapped onto the High group's
#                 autonomous-2049 distribution.
#
#   HSR_scen(c,t) = HSR_auto(c,t) + w(t) * gap_c
#     w(t)  = (t - 2024) / 25            (0 at 2024, 1 at 2049)
#     gap_c = max(QM_target_c - HSR_auto(c,2049), 0)
#
# Properties: immediate separation from Autonomous in 2025; smooth monotone
# paths (no pmax kink); 2049 cross-city dispersion of intervened cities
# matches the reference group's distribution (no interval collapse); no city
# is made worse off than under autonomous adaptation.

auto_2049 <- proj_auto_summary %>% filter(Year == 2049)

ref_mod_dist  <- auto_2049 %>%
  filter(Resilience_Group == "Moderate") %>% pull(hsr_median)
ref_high_dist <- auto_2049 %>%
  filter(Resilience_Group == "High") %>% pull(hsr_median)

# --- Quantile mapping: Targeted (Low -> Moderate distribution) ---
targ_map <- auto_2049 %>%
  filter(Resilience_Group == "Low") %>%
  mutate(p_rank = (rank(hsr_median, ties.method = "average") - 0.5) / n(),
         hsr_target_2049 = as.numeric(quantile(ref_mod_dist, p_rank, type = 7)),
         gap_targ = pmax(hsr_target_2049 - hsr_median, 0)) %>%
  dplyr::select(CityCode, hsr_target_2049, gap_targ)

# --- Quantile mapping: Aspirational (Low + Moderate -> High distribution) ---
asp_map <- auto_2049 %>%
  filter(Resilience_Group %in% c("Low", "Moderate")) %>%
  mutate(p_rank = (rank(hsr_median, ties.method = "average") - 0.5) / n(),
         hsr_target_2049 = as.numeric(quantile(ref_high_dist, p_rank, type = 7)),
         gap_asp = pmax(hsr_target_2049 - hsr_median, 0)) %>%
  dplyr::select(CityCode, hsr_target_2049, gap_asp)

# Reference distribution summary (for reporting)
Ref_Dist_Summary <- tibble(
  Reference = c("Moderate group, autonomous 2049 (Targeted reference)",
                "High group, autonomous 2049 (Aspirational reference)"),
  N_cities = c(length(ref_mod_dist), length(ref_high_dist)),
  Mean   = c(mean(ref_mod_dist), mean(ref_high_dist)),
  P25    = c(quantile(ref_mod_dist, 0.25), quantile(ref_high_dist, 0.25)),
  Median = c(median(ref_mod_dist), median(ref_high_dist)),
  P75    = c(quantile(ref_mod_dist, 0.75), quantile(ref_high_dist, 0.75))
)
write.xlsx(Ref_Dist_Summary, file.path(dir_tab, "Scenario_Reference_Distributions.xlsx"))
write.xlsx(targ_map, file.path(dir_tab, "QuantileMapping_Targeted.xlsx"))
write.xlsx(asp_map, file.path(dir_tab, "QuantileMapping_Aspirational.xlsx"))

# (1) No adaptation: HSR fixed at 2024
scen_none <- expand_grid(CityCode = hsr_start_2024$CityCode, Year = Years_Proj) %>%
  left_join(hsr_start_2024 %>% dplyr::select(CityCode, HSR_2024 = HSR),
            by = "CityCode") %>%
  transmute(CityCode, Year, HSR_proj = HSR_2024, Scenario = "No Adaptation")

# (2) Autonomous: posterior median trajectory (trend-damped)
scen_auto <- proj_auto_summary %>%
  transmute(CityCode, Year, HSR_proj = hsr_median, Scenario = "Autonomous")

# (3) Targeted: autonomous + linear wedge to quantile-mapped Moderate target
scen_targ <- proj_auto_summary %>%
  dplyr::select(CityCode, Year, hsr_auto = hsr_median) %>%
  left_join(targ_map %>% dplyr::select(CityCode, gap_targ), by = "CityCode") %>%
  mutate(gap_targ = replace_na(gap_targ, 0),
         HSR_proj = hsr_auto + (Year - Year_End) / 25 * gap_targ,
         Scenario = "Targeted") %>%
  dplyr::select(CityCode, Year, HSR_proj, Scenario)

# (4) Aspirational: autonomous + linear wedge to quantile-mapped High target
scen_asp <- proj_auto_summary %>%
  dplyr::select(CityCode, Year, hsr_auto = hsr_median) %>%
  left_join(asp_map %>% dplyr::select(CityCode, gap_asp), by = "CityCode") %>%
  mutate(gap_asp = replace_na(gap_asp, 0),
         HSR_proj = hsr_auto + (Year - Year_End) / 25 * gap_asp,
         Scenario = "Aspirational") %>%
  dplyr::select(CityCode, Year, HSR_proj, Scenario)

scenarios_df <- bind_rows(scen_none, scen_auto, scen_targ, scen_asp) %>%
  left_join(City_Map %>% dplyr::select(CityCode, Province), by = "CityCode") %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  mutate(Scenario = factor(Scenario, levels = names(colors_scenario)))

saveRDS(scenarios_df, file.path(dir_res, "HSR_Projection_AllScenarios_City.rds"))
write.xlsx(scenarios_df, file.path(dir_res, "HSR_Projection_AllScenarios_City.xlsx"))

# --- 2024 anchor rows (connect history and projection in figures) ---
anchor_2024 <- expand_grid(CityCode = hsr_start_2024$CityCode,
                           Scenario = names(colors_scenario)) %>%
  left_join(hsr_start_2024 %>% dplyr::select(CityCode, Province, HSR_2024 = HSR),
            by = "CityCode") %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  transmute(CityCode, Year = Year_End, HSR_proj = HSR_2024,
            Scenario = factor(Scenario, levels = names(colors_scenario)),
            Province, Resilience_Group)

scenarios_plot <- bind_rows(scenarios_df, anchor_2024) %>%
  arrange(Scenario, CityCode, Year)

HSR_Full_Panel <- bind_rows(
  HSR_City_Panel %>%
    mutate(Scenario = "Historical") %>%
    dplyr::select(Province, CityCode, Year, HSR, Scenario),
  scenarios_df %>%
    dplyr::select(Province, CityCode, Year, HSR = HSR_proj, Scenario))
write.xlsx(HSR_Full_Panel, file.path(dir_res, "HSR_Full_Panel_City_2000_2049.xlsx"))
saveRDS(HSR_Full_Panel, file.path(dir_res, "HSR_Full_Panel_City_2000_2049.rds"))

# ===================== B5. Figures ==========================================

# --- p2: Scenario trajectories (mean +/- IQR), with historical band ---------

fig2_data <- scenarios_plot %>%
  group_by(Scenario, Year) %>%
  summarise(mean_hsr = mean(HSR_proj), q25 = quantile(HSR_proj, 0.25),
            q75 = quantile(HSR_proj, 0.75), .groups = "drop")

hist_band <- HSR_City_Panel %>%
  group_by(Year) %>%
  summarise(mean_hsr = mean(HSR), q25 = quantile(HSR, 0.25),
            q75 = quantile(HSR, 0.75), .groups = "drop")

p2 <- ggplot() +
  geom_ribbon(data = hist_band, aes(Year, ymin = q25, ymax = q75),
              fill = "grey60", alpha = 0.25) +
  geom_line(data = hist_band, aes(Year, mean_hsr),
            color = "grey25", linewidth = 1) +
  geom_ribbon(data = fig2_data,
              aes(Year, ymin = q25, ymax = q75, fill = Scenario),
              alpha = 0.13) +
  geom_line(data = fig2_data,
            aes(Year, mean_hsr, color = Scenario), linewidth = 1.1) +
  geom_vline(xintercept = Year_End, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(#subtitle = "Lines: cross-city mean; ribbons: interquartile range (middle 50% of cities)",
    x = "Year", y = "HSR index") +
  theme_pub2() + theme(legend.position = "right")
save_fig(p2, "Fig_Scenario_Trajectories_City", 7, 4.5)
save_fig_data(fig2_data, "Fig_Scenario_Trajectories_City")
save_fig_data(hist_band, "Fig_Scenario_Trajectories_City_Historical")

# --- p3: Fan chart by resilience group, all four scenarios ------------------
# Historical lines are cross-city group MEANS (fixed 2024 classification).
# Bands are Bayesian credible intervals of the autonomous projection.

fig3_band_anchor <- HSR_City_Panel %>%
  filter(Year == Year_End) %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  group_by(Resilience_Group) %>%
  summarise(m = mean(HSR), .groups = "drop") %>%
  transmute(Resilience_Group, Year = Year_End,
            hsr_median = m, q025 = m, q10 = m, q25 = m,
            q75 = m, q90 = m, q975 = m)

fig3_band <- proj_auto_summary %>%
  group_by(Resilience_Group, Year) %>%
  summarise(across(c(hsr_median, q025, q10, q25, q75, q90, q975), mean),
            .groups = "drop") %>%
  bind_rows(fig3_band_anchor) %>%
  arrange(Resilience_Group, Year)

fig3_lines <- scenarios_plot %>%
  group_by(Scenario, Resilience_Group, Year) %>%
  summarise(mean_hsr = mean(HSR_proj), .groups = "drop")

hist_group_fixed <- HSR_City_Panel %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  group_by(Resilience_Group, Year) %>%
  summarise(mean_hsr = mean(HSR), .groups = "drop")

make_fan_plot <- function(band, lines, hist, groups, show_legend = TRUE) {
  ggplot(band %>% filter(Resilience_Group %in% groups), aes(Year)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = Resilience_Group), alpha = 0.10) +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = Resilience_Group), alpha = 0.18) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = Resilience_Group), alpha = 0.28) +
    geom_line(data = hist %>% filter(Resilience_Group %in% groups),
              aes(Year, mean_hsr, color = Resilience_Group),
              linewidth = 0.9, show.legend = FALSE) +
    geom_line(data = lines %>% filter(Resilience_Group %in% groups),
              aes(Year, mean_hsr, color = Resilience_Group,
                  linetype = Scenario), linewidth = 0.8) +
    geom_vline(xintercept = Year_End, linetype = "dashed", color = "grey40") +
    facet_wrap(~ Resilience_Group) +
    scale_color_manual(values = colors_group, guide = "none") +
    scale_fill_manual(values = colors_group, guide = "none") +
    scale_linetype_manual(values = scenario_linetypes, name = "Scenario") +
    labs(#subtitle = "Bands: 50% / 80% / 95% credible intervals (autonomous adaptation); historical lines: cross-city group means",
      x = "Year", y = "HSR index") +
    theme_pub2() +
    theme(legend.position = if (show_legend) "bottom" else "none",
          strip.text = element_text(face = "bold"))
}

p3 <- make_fan_plot(fig3_band, fig3_lines, hist_group_fixed,
                    c("High", "Moderate", "Low"))
save_fig(p3, "Fig_Fan_Chart_AllScenarios_City", 9, 4)
save_fig_data(fig3_band, "Fig_Fan_Chart_AllScenarios_City_Bands")
save_fig_data(fig3_lines, "Fig_Fan_Chart_AllScenarios_City_Lines")
save_fig_data(hist_group_fixed, "Fig_Fan_Chart_AllScenarios_City_Historical")

for (g in c("High", "Moderate", "Low")) {
  pg <- make_fan_plot(fig3_band, fig3_lines, hist_group_fixed, g)
  save_fig(pg, paste0("Fig_Fan_Chart_AllScenarios_City_", g), 4, 4)
}

# --- p4: Single-panel overlay, COMMON bandwidth + GLOBAL height scale -------
# Previous version used two geom_density_ridges layers, each estimating its
# own joint bandwidth and normalizing heights to the tallest density WITHIN
# the layer -> identical No-Adaptation data appeared vertically compressed.
# Fix: compute all densities manually with ONE common bandwidth, scale by ONE
# global maximum, and draw with geom_ridgeline (absolute heights).

fig4_base <- HSR_2024_class %>%
  dplyr::select(CityCode, HSR = HSR_2024) %>%
  crossing(Scenario = names(colors_scenario)) %>%
  mutate(Distribution = "2024 (baseline)")

fig4_2049 <- scenarios_df %>%
  filter(Year == 2049) %>%
  transmute(CityCode, HSR = HSR_proj, Scenario = as.character(Scenario),
            Distribution = "2049 (projected)")

fig4_all <- bind_rows(fig4_base, fig4_2049)

bw_common <- bw.nrd0(fig4_all$HSR)   # single bandwidth for ALL densities
x_lo <- min(fig4_all$HSR) - 3 * bw_common
x_hi <- max(fig4_all$HSR) + 3 * bw_common

dens_df <- fig4_all %>%
  group_by(Distribution, Scenario) %>%
  group_modify(~ {
    d <- density(.x$HSR, bw = bw_common, from = x_lo, to = x_hi, n = 512)
    tibble(HSR = d$x, dens = d$y)
  }) %>%
  ungroup()

global_max <- max(dens_df$dens)
scen_levels_p4 <- rev(names(colors_scenario))   # Aspirational at bottom row

dens_df <- dens_df %>%
  mutate(Scenario = factor(Scenario, levels = scen_levels_p4),
         y_pos = as.numeric(Scenario),
         h = dens / global_max * 0.95)          # heights on one absolute scale

p4 <- ggplot() +
  geom_ridgeline(data = dens_df %>% filter(Distribution == "2024 (baseline)"),
                 aes(x = HSR, y = y_pos, height = h,
                     group = Scenario, fill = Scenario),
                 alpha = 0.55, color = "white") +
  geom_ridgeline(data = dens_df %>% filter(Distribution == "2049 (projected)"),
                 aes(x = HSR, y = y_pos, height = h,
                     group = Scenario, color = Scenario),
                 fill = NA, linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(breaks = seq_along(scen_levels_p4),
                     labels = scen_levels_p4) +
  scale_fill_manual(values = colors_scenario, guide = "none") +
  scale_color_manual(values = colors_scenario, guide = "none") +
  labs(#subtitle = "Filled: 2024 baseline; dashed outline: 2049 projection.\nCommon kernel bandwidth and global height scale (densities directly comparable)",
    x = "HSR index", y = NULL) +
  theme_pub2()
save_fig(p4, "Fig_Scenario_Distributions_2024_vs_2049_City", 3, 3)
save_fig_data(dens_df, "Fig_Scenario_Distributions_2024_vs_2049_City")

# --- Summary table by period ---
summary_table <- scenarios_df %>%
  mutate(Period = case_when(Year <= 2030 ~ "2025-2030",
                            Year <= 2040 ~ "2031-2040",
                            TRUE ~ "2041-2049")) %>%
  group_by(Scenario, Period) %>%
  summarise(mean_HSR = mean(HSR_proj), sd_HSR = sd(HSR_proj),
            min_HSR = min(HSR_proj), max_HSR = max(HSR_proj), .groups = "drop")
write.xlsx(summary_table, file.path(dir_tab, "HSR_Summary_By_Scenario_Period_City.xlsx"))

# ===================== B6. Dual-window leave-future-out cross-validation ====

run_lfo_cv <- function(train_end, test_years, label) {
  
  train_df <- model_df %>% filter(Year <= train_end)
  
  fit_cv <- update(fit_convergence, newdata = train_df,
                   chains = 4, iter = 3000, warmup = 1500,
                   cores = 4, seed = 2026)
  post_cv <- as_draws_df(fit_cv)
  
  hsr_base <- HSR_City_Panel %>%
    filter(Year == train_end) %>%
    left_join(City_Map, by = c("Province", "CityCode")) %>%
    left_join(Prov_Map, by = "Province") %>%
    arrange(City_ID)
  
  c2p <- setNames(hsr_base$Province_ID, hsr_base$City_ID)
  
  proj <- project_hsr_city(
    post_cv, 1000,
    hsr_start = hsr_base$HSR,
    city_ids = hsr_base$City_ID,
    city_to_prov = c2p,
    years_proj = test_years,
    time_mean = mean(train_df$time_index),
    time_sd = sd(train_df$time_index),
    t_cap_year = train_end, damp_trend = TRUE
  )
  
  cv_summary <- map_dfr(seq_along(test_years), function(t) {
    map_dfr(seq_len(nrow(hsr_base)), function(i) {
      v <- proj[, i, t]
      tibble(City_ID = hsr_base$City_ID[i], Year = test_years[t],
             pred_median = median(v),
             pred_q025 = quantile(v, 0.025), pred_q975 = quantile(v, 0.975))
    })
  }) %>%
    left_join(City_Map, by = "City_ID") %>%
    left_join(HSR_City_Panel %>%
                dplyr::select(CityCode, Year, obs_hsr = HSR),
              by = c("CityCode", "Year")) %>%
    mutate(CV_Window = label)
  
  metrics <- cv_summary %>%
    filter(!is.na(obs_hsr)) %>%
    summarise(CV_Window = label,
              Train_End = train_end,
              Test_Years = paste(range(test_years), collapse = "-"),
              RMSE = sqrt(mean((obs_hsr - pred_median)^2)),
              MAE = mean(abs(obs_hsr - pred_median)),
              Correlation = cor(obs_hsr, pred_median),
              Coverage_95 = mean(obs_hsr >= pred_q025 & obs_hsr <= pred_q975),
              Mean_Bias = mean(pred_median - obs_hsr))
  
  list(summary = cv_summary, metrics = metrics)
}

cv1 <- run_lfo_cv(2019, 2020:2024, "CV1: train 2000-2019 (pandemic stress test)")
cv2 <- run_lfo_cv(2022, 2023:2024, "CV2: train 2000-2022 (operational accuracy)")

val_metrics_all <- bind_rows(cv1$metrics, cv2$metrics)
cv_summary_all  <- bind_rows(cv1$summary, cv2$summary)

write.xlsx(val_metrics_all, file.path(dir_tab, "LFO_CV_Validation_Metrics_City.xlsx"))
write.xlsx(cv_summary_all, file.path(dir_tab, "LFO_CV_Predictions_City.xlsx"))

make_cv_plot <- function(cvres) {
  m <- cvres$metrics
  ggplot(cvres$summary %>% filter(!is.na(obs_hsr)),
         aes(obs_hsr, pred_median)) +
    geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975),
                  color = "grey70", width = 0, alpha = 0.3) +
    geom_point(alpha = 0.4, color = "#2C3E50", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    annotate("text",
             x = min(cvres$summary$obs_hsr, na.rm = TRUE),
             y = max(cvres$summary$pred_q975, na.rm = TRUE),
             hjust = 0, vjust = 1, size = 3.3,
             label = paste0("r = ", round(m$Correlation, 3),
                            "\nRMSE = ", round(m$RMSE, 3),
                            "\nMean bias = ", round(m$Mean_Bias, 3),
                            "\n95% coverage = ",
                            round(m$Coverage_95 * 100, 1), "%")) +
    labs(#title = m$CV_Window,
      x = "Observed HSR", y = "Predicted HSR (posterior median)") +
    coord_equal() + theme_pub2()
}

p5a <- make_cv_plot(cv1)
p5b <- make_cv_plot(cv2)

# Separate figures
save_fig(p5a, "Fig_LFO_CV_Window1_PandemicStressTest", 6, 4)
save_fig(p5b, "Fig_LFO_CV_Window2_OperationalAccuracy", 6, 4)

# Combined figure
p5 <- p5a + p5b #+ plot_annotation(tag_levels = "a")
save_fig(p5, "Fig_LFO_CrossValidation_DualWindow_City", 6, 3.5)
save_fig_data(cv_summary_all, "Fig_LFO_CrossValidation_DualWindow_City")

# ===================== B7. Sensitivity: alternative projection methods ======
# Lines: cross-city group MEANS of per-city point projections (posterior
# medians for Bayesian variants; deterministic fits for linear/logistic).
# Ribbons: cross-city IQR. Historical 2000-2024 shown with IQR, connected
# at the 2024 anchor; dashed reference line at 2024.

proj_linear <- HSR_City_Panel %>%
  group_by(CityCode) %>%
  group_modify(~ {
    fit <- lm(HSR ~ Year, data = .x)
    tibble(Year = Years_Proj,
           HSR_linear = pmax(predict(fit,
                                     newdata = data.frame(Year = Years_Proj)), 0))
  }) %>% ungroup()

proj_logistic <- HSR_City_Panel %>%
  group_by(CityCode) %>%
  group_modify(~ {
    out <- tryCatch({
      nf <- nls(HSR ~ K / (1 + exp(-r * (Year - t0))), data = .x,
                start = list(K = max(.x$HSR) * 1.2, r = 0.1, t0 = 2012),
                control = nls.control(maxiter = 500))
      predict(nf, newdata = data.frame(Year = Years_Proj))
    }, error = function(e) {
      predict(lm(HSR ~ Year, data = .x),
              newdata = data.frame(Year = Years_Proj))
    })
    tibble(Year = Years_Proj, HSR_logistic = pmax(out, 0))
  }) %>% ungroup()

sens_methods <- proj_auto_summary %>%
  dplyr::select(CityCode, Year, HSR_bayes_damped = hsr_median, Resilience_Group) %>%
  left_join(proj_auto_summary_undamped %>%
              dplyr::select(CityCode, Year, HSR_bayes_extrap = hsr_median),
            by = c("CityCode", "Year")) %>%
  left_join(proj_linear, by = c("CityCode", "Year")) %>%
  left_join(proj_logistic, by = c("CityCode", "Year"))
write.xlsx(sens_methods, file.path(dir_tab, "Sensitivity_ProjectionMethods_City.xlsx"))

method_levels <- c("Bayesian, trend-damped (primary)",
                   "Bayesian, extrapolated trend",
                   "Linear extrapolation",
                   "Logistic growth")

sens_long <- sens_methods %>%
  pivot_longer(c(HSR_bayes_damped, HSR_bayes_extrap, HSR_linear, HSR_logistic),
               names_to = "Method", values_to = "HSR") %>%
  mutate(Method = recode(Method,
                         HSR_bayes_damped = "Bayesian, trend-damped (primary)",
                         HSR_bayes_extrap = "Bayesian, extrapolated trend",
                         HSR_linear = "Linear extrapolation",
                         HSR_logistic = "Logistic growth"),
         Method = factor(Method, levels = method_levels)) %>%
  dplyr::select(CityCode, Resilience_Group, Year, Method, HSR)

# 2024 anchor: observed values, replicated for every method
anchor_methods <- HSR_City_Panel %>%
  filter(Year == Year_End) %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  dplyr::select(CityCode, Resilience_Group, HSR) %>%
  crossing(Method = factor(method_levels, levels = method_levels)) %>%
  mutate(Year = Year_End)

fig6_data <- bind_rows(sens_long, anchor_methods) %>%
  group_by(Method, Resilience_Group, Year) %>%
  summarise(mean_hsr = mean(HSR, na.rm = TRUE),
            q25 = quantile(HSR, 0.25, na.rm = TRUE),
            q75 = quantile(HSR, 0.75, na.rm = TRUE), .groups = "drop")

hist6_band <- HSR_City_Panel %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  group_by(Resilience_Group, Year) %>%
  summarise(mean_hsr = mean(HSR),
            q25 = quantile(HSR, 0.25),
            q75 = quantile(HSR, 0.75), .groups = "drop")

make_method_plot <- function(groups, show_legend = TRUE) {
  ggplot() +
    geom_ribbon(data = hist6_band %>% filter(Resilience_Group %in% groups),
                aes(Year, ymin = q25, ymax = q75),
                fill = "grey60", alpha = 0.25) +
    geom_line(data = hist6_band %>% filter(Resilience_Group %in% groups),
              aes(Year, mean_hsr), color = "grey25", linewidth = 0.9) +
    geom_ribbon(data = fig6_data %>% filter(Resilience_Group %in% groups),
                aes(Year, ymin = q25, ymax = q75, fill = Method),
                alpha = 0.12) +
    geom_line(data = fig6_data %>% filter(Resilience_Group %in% groups),
              aes(Year, mean_hsr, color = Method , linetype = Method
              ),
              linewidth = 0.9) +
    geom_vline(xintercept = Year_End, linetype = "dashed", color = "grey40") +
    facet_wrap(~ Resilience_Group) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(#subtitle = "Lines: cross-city group means; ribbons: interquartile range across cities;\ngrey: historical (2000-2024)",
      x = "Year", y = "HSR index"
    ) +
    theme_pub2() +
    theme(legend.position = if (show_legend) "bottom" else "none",
          strip.text = element_text(face = "bold"))
}

p6 <- make_method_plot(c("High", "Moderate", "Low"))
save_fig(p6, "Fig_Sensitivity_Methods_City", 9, 4)
save_fig_data(fig6_data, "Fig_Sensitivity_Methods_City")
save_fig_data(hist6_band, "Fig_Sensitivity_Methods_City_Historical")

for (g in c("High", "Moderate", "Low")) {
  pg <- make_method_plot(g)
  save_fig(pg, paste0("Fig_Sensitivity_Methods_City_", g), 4, 4)
}

# ===================== B8. Sensitivity: intervention ambition ================
# Ambition redefined as GAP-CLOSURE FRACTION gamma of the quantile-mapped
# Targeted gap:  HSR = auto + w(t) * gamma * gap_c.
# gamma = 1 reproduces the main Targeted scenario; trajectories separate from
# Autonomous immediately and smoothly for all gamma > 0.

ambition_grid <- tibble(
  Ambition = c("25% gap closure", "50% gap closure", "75% gap closure",
               "100% gap closure (main Targeted)", "125% gap closure"),
  gamma = c(0.25, 0.50, 0.75, 1.00, 1.25)
)
write.xlsx(ambition_grid, file.path(dir_tab, "Sensitivity_AmbitionGrid.xlsx"))

sens_target <- map_dfr(seq_len(nrow(ambition_grid)), function(k) {
  g <- ambition_grid$gamma[k]
  proj_auto_summary %>%
    dplyr::select(CityCode, Year, hsr_auto = hsr_median, Resilience_Group) %>%
    left_join(targ_map %>% dplyr::select(CityCode, gap_targ), by = "CityCode") %>%
    mutate(gap_targ = replace_na(gap_targ, 0),
           HSR_target_scen = hsr_auto + (Year - Year_End) / 25 * g * gap_targ,
           Ambition = ambition_grid$Ambition[k])
}) %>%
  left_join(City_Map %>% dplyr::select(CityCode, Province), by = "CityCode") %>%
  dplyr::select(CityCode, Province, Year, Ambition, HSR_target_scen,
                Resilience_Group)
write.xlsx(sens_target, file.path(dir_tab, "Sensitivity_InterventionTargets_City.xlsx"))

# 2024 anchor for the Low group
anchor_low_2024 <- HSR_City_Panel %>%
  filter(Year == Year_End) %>%
  left_join(HSR_2024_class %>% dplyr::select(CityCode, Resilience_Group),
            by = "CityCode") %>%
  filter(Resilience_Group == "Low") %>%
  summarise(mean_hsr = mean(HSR)) %>%
  crossing(Ambition = ambition_grid$Ambition) %>%
  mutate(Year = Year_End)

fig7_data <- sens_target %>%
  filter(Resilience_Group == "Low") %>%
  group_by(Ambition, Year) %>%
  summarise(mean_hsr = mean(HSR_target_scen), .groups = "drop") %>%
  bind_rows(anchor_low_2024) %>%
  arrange(Ambition, Year)

fig7_data <- sens_target %>%
  filter(Resilience_Group == "Low") %>%
  group_by(Ambition, Year) %>%
  summarise(mean_hsr = mean(HSR_target_scen), .groups = "drop") %>%
  bind_rows(anchor_low_2024) %>%
  arrange(Ambition, Year) %>%
  mutate(Ambition = factor(Ambition, levels = ambition_grid$Ambition))

fig7_data <- sens_target %>%
  filter(Resilience_Group == "Low") %>%
  group_by(Ambition, Year) %>%
  summarise(mean_hsr = mean(HSR_target_scen), .groups = "drop") %>%
  bind_rows(anchor_low_2024) %>%
  arrange(Ambition, Year) %>%
  mutate(Ambition = factor(Ambition, levels = rev(ambition_grid$Ambition)))


p7 <- ggplot(fig7_data, aes(Year, mean_hsr, color = Ambition)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = Year_End, linetype = "dashed", color = "grey40") +
  scale_color_viridis_d(option = "plasma", end = 0.9, name = "Ambition") +
  labs(# subtitle = "Fraction of the quantile-mapped gap to the Moderate reference distribution closed by 2049",
    x = "Year", y = "Mean HSR (low-resilience group)"
  ) +
  theme_pub2() + theme(legend.position = "right", legend.justification = "top")
save_fig(p7, "Fig_Sensitivity_Targets_City", 10, 4)
save_fig_data(fig7_data, "Fig_Sensitivity_Targets_City")

# ===================== B9. Province-level summary (for downstream) ==========

Prov_Proj_Summary <- scenarios_df %>%
  left_join(City_Pop_for_weight %>% filter(Year == Year_End) %>%
              dplyr::select(CityCode, Pop), by = "CityCode") %>%
  group_by(Province, Year, Scenario) %>%
  summarise(HSR_proj = weighted.mean(HSR_proj, w = Pop, na.rm = TRUE),
            .groups = "drop")
write.xlsx(Prov_Proj_Summary,
           file.path(dir_res, "HSR_Province_Projection_AllScenarios.xlsx"))
saveRDS(Prov_Proj_Summary,
        file.path(dir_res, "HSR_Province_Projection_AllScenarios.rds"))

cat("\n==============================================\n")
cat("ANALYSIS COMPLETE (V5).\n")
cat("  - Distributional targets via rank-preserving quantile mapping\n")
cat("  - Smooth additive-wedge intervention trajectories (no pmax kink)\n")
cat("  - p4 rebuilt: common bandwidth + global height scale\n")
cat("  - Ambition sensitivity as gap-closure fractions\n")
cat("==============================================\n")





