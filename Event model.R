#############################################################################################
#
#                          TROPICAL CYCLONES AND HEALTH
#                  Event Model: Bayesian INLA Analysis (City Level, v7.0)
#
#############################################################################################

# ================================================================
#  PART 0: CONFIGURATION & UTILITY FUNCTIONS
# ================================================================

library(dplyr); library(tidyr); library(lubridate); library(ggplot2)
library(openxlsx); library(stringr); library(INLA)

# ===================== USER INPUT =====================
Outcome_var   <- "All_cause"
Analysis_type <- "Disease_highHSR"

# "All_cause","Infectious","Circulatory","Respiratory",
# "Injuries","Mental","Nervous","Genitourinary",
# "Neoplasms","Endocrine","Nutritional","Metabolic",
# "Age_0_19","Age_20_44","Age_45_64","Age_65plus",
# "Male","Female"

# "Disease","Age_Stratified","Sex_Stratified"

# ======================================================

setwd("C:/Project/Tropical cyclone")

# --- Subgroup analysis: population column mapping ---
POP_COL_MAP <- c(
  "All_cause"    = "pop_total",  "Infectious"  = "pop_total",
  
  "Circulatory"  = "pop_total",  "Respiratory" = "pop_total",
  "Injuries"     = "pop_total",  "Mental"      = "pop_total",
  "Nervous"      = "pop_total",  "Genitourinary" = "pop_total",
  "Neoplasms"    = "pop_total",  "Endocrine"   = "pop_total",
  "Nutritional"  = "pop_total",  "Metabolic"   = "pop_total",
  "Age_0_19"     = "pop_lt20",   "Age_20_44"   = "pop_20to44",
  "Age_45_64"    = "pop_45to64", "Age_65plus"  = "pop_ge65",
  "Male"         = "pop_total",  "Female"      = "pop_total"
)

if (Analysis_type %in% c("Age_Stratified", "Sex_Stratified")) {
  PATH_EVENT_PANEL <- "Result TC-health CityLevel/Event_Panel_Filtered_City_AgeSex.rds"
} else {
  PATH_EVENT_PANEL <- "Result TC-health CityLevel/Event_Panel_Filtered_City.rds"
}
cat(sprintf("  Event panel path: %s\n", PATH_EVENT_PANEL))

DIR_OUTPUT    <- file.path("Result TC-health CityLevel", Analysis_type, Outcome_var)
DIR_FIG       <- file.path(DIR_OUTPUT, "Figures")
DIR_TABLE_OUT <- file.path(DIR_OUTPUT, "Tables")
for (d in c(DIR_OUTPUT, DIR_FIG, DIR_TABLE_OUT))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Parameters ----
PRE_DAYS <- 90; MAX_LAG_DAYS <- 180; WINDOW_SIZE <- 30
N_PRE_WINDOWS  <- PRE_DAYS %/% WINDOW_SIZE
N_POST_WINDOWS <- MAX_LAG_DAYS %/% WINDOW_SIZE
window_labels  <- c(paste0("-", N_PRE_WINDOWS:1, "mo"), "Active",
                    paste0("+", 1:N_POST_WINDOWS, "mo"))
window_display <- c(paste0("-", N_PRE_WINDOWS:1), "Active",
                    paste0("+", 1:N_POST_WINDOWS))
POST_WINDOWS   <- c("Active", paste0("+", 1:N_POST_WINDOWS, "mo"))

REF_MAIN  <- "-3mo"
REF_SENS1 <- "-2mo"
N_SAMPLES <- 1000
HSR_STEP  <- 0.1
MIN_COVERAGE <- 0.5

if (Analysis_type %in% c("Age_Stratified", "Sex_Stratified")) {
  MIN_TOTAL_COUNT <- 100
  cat(sprintf("  Subgroup analysis: MIN_TOTAL_COUNT lowered to %d\n", MIN_TOTAL_COUNT))
} else {
  MIN_TOTAL_COUNT <- 500
}

# ---- Priors ----
PRIOR_FIXED <- list(mean = 0, prec = 0.001,
                    mean.intercept = 0, prec.intercept = 0.001)
CTRL_FAMILY <- list(hyper = list(size = list(prior = "pc.mgamma", param = c(7))))

cat("Prior specification:\n")
cat("  Fixed effects  : Gaussian N(0, 1000)\n")
cat("  Random effects : PC priors, see formula definitions\n")
cat("  NB dispersion  : pc.mgamma(7)\n")

# ---- Plotting utilities ----
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(text = element_text(colour = "black"),
          axis.title = element_text(size = base_size, face = "bold"),
          axis.text = element_text(size = base_size - 1, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size = base_size, face = "bold"),
          legend.background = element_blank())
}

save_both <- function(p, path_base, w = 7, h = 5.5) {
  p_trans <- p +
    theme(plot.background  = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  tryCatch(ggsave(paste0(path_base, ".pdf"), p_trans, width = w, height = h,
                  dpi = 300, device = cairo_pdf, bg = "transparent"), error = function(e) NULL)
  tryCatch(ggsave(paste0(path_base, ".png"), p_trans, width = w, height = h,
                  dpi = 300, bg = "transparent"), error = function(e) NULL)
  invisible(NULL)
}

save_zoom_transparent <- function(p, path_base, w = 4.4, h = 3) {
  ggsave(paste0(path_base, "_Zoom.png"), plot = p +
           theme(plot.background = element_rect(fill = "transparent", colour = NA),
                 panel.background = element_rect(fill = "transparent", colour = NA),
                 legend.background = element_rect(fill = "transparent", colour = NA),
                 legend.box.background = element_rect(fill = "transparent", colour = NA)),
         width = w, height = h, dpi = 300, bg = "transparent")
}

# ================================================================
#  UTILITY FUNCTIONS
# ================================================================

extract_window_effects <- function(fit, ref_window, window_labels) {
  if (is.null(fit)) return(NULL)
  fs <- fit$summary.fixed
  out <- lapply(window_labels, function(wl) {
    if (wl == ref_window)
      return(data.frame(window = wl, ERR_pct = 0, lower = 0, upper = 0,
                        significant = "ref"))
    i <- which(rownames(fs) == paste0("window_f", wl))
    if (length(i) != 1) return(NULL)
    data.frame(window = wl,
               ERR_pct = (exp(fs$mean[i]) - 1) * 100,
               lower   = (exp(fs$`0.025quant`[i]) - 1) * 100,
               upper   = (exp(fs$`0.975quant`[i]) - 1) * 100,
               significant = ifelse(fs$`0.025quant`[i] > 0 |
                                      fs$`0.975quant`[i] < 0, "*", ""))
  })
  out <- bind_rows(out)
  out$window <- factor(out$window, levels = window_labels)
  out[order(out$window), ]
}

extract_interaction_effects <- function(fit, modifier, window_labels,
                                        ref_window, scale_by = 1) {
  if (is.null(fit)) return(NULL)
  fs <- fit$summary.fixed
  out <- lapply(setdiff(window_labels, ref_window), function(wl) {
    i <- which(rownames(fs) == paste0("window_f", wl, ":", modifier))
    if (length(i) != 1) return(NULL)
    data.frame(window = wl,
               theta      = fs$mean[i],
               theta_lo   = fs$`0.025quant`[i],
               theta_hi   = fs$`0.975quant`[i],
               pct_change    = (exp(fs$mean[i] * scale_by) - 1) * 100,
               pct_change_lo = (exp(fs$`0.025quant`[i] * scale_by) - 1) * 100,
               pct_change_hi = (exp(fs$`0.975quant`[i] * scale_by) - 1) * 100,
               significant = ifelse(fs$`0.025quant`[i] > 0 |
                                      fs$`0.975quant`[i] < 0, "*", ""))
  })
  out <- bind_rows(out)
  if (nrow(out) == 0) return(NULL)
  out$window <- factor(out$window, levels = window_labels)
  out[order(out$window), ]
}

get_coef_draws <- function(fit, coef_names, n = N_SAMPLES, seed = 123) {
  ps <- inla.posterior.sample(n, fit, seed = seed)
  ln <- rownames(ps[[1]]$latent)
  idx <- sapply(coef_names, function(nm) {
    hit <- which(ln == paste0(nm, ":1"))
    if (length(hit) == 0) hit <- which(ln == nm)
    if (length(hit) == 0) NA_integer_ else hit[1]
  })
  if (any(is.na(idx))) {
    cat("  WARNING: not found:", paste(coef_names[is.na(idx)], collapse = ", "), "\n")
    return(NULL)
  }
  m <- sapply(ps, function(s) s$latent[idx, 1])
  if (is.null(dim(m))) m <- matrix(m, nrow = 1)
  rownames(m) <- coef_names
  m
}

summarise_draws <- function(x) data.frame(
  median = median(x), lower_2.5 = quantile(x, 0.025, names = FALSE),
  upper_97.5 = quantile(x, 0.975, names = FALSE), prob_positive = mean(x > 0))

fmt_p <- function(p) ifelse(p >= 0.999, ">0.999",
                            ifelse(p <= 0.001, "<0.001", sprintf("%.3f", p)))

ptxt <- function(x, dir = c("gt0", "lt0")) {
  dir <- match.arg(dir)
  if (is.null(x) || length(x) == 0) return(NA_character_)
  p <- if (dir == "gt0") mean(x > 0) else mean(x < 0)
  sprintf("P(%s0) = %s", ifelse(dir == "gt0", ">", "<"), fmt_p(p))
}

# --- NEW: Compute window-specific ΔERR from beta and theta draw matrices ---
compute_window_derr <- function(M_beta, M_theta, scale_by = 1) {
  # M_beta: [n_windows x N_SAMPLES] matrix of beta draws
  
  # M_theta: [n_windows x N_SAMPLES] matrix of theta draws
  # Returns: data.frame with window-level ΔERR summary
  n_w <- nrow(M_beta)
  out <- lapply(seq_len(n_w), function(w) {
    rr_w <- exp(M_beta[w, ])
    derr_w <- rr_w * (exp(M_theta[w, ] * scale_by) - 1) * 100
    data.frame(
      window_idx  = w,
      derr_median = median(derr_w),
      derr_lo     = quantile(derr_w, 0.025, names = FALSE),
      derr_hi     = quantile(derr_w, 0.975, names = FALSE),
      significant = ifelse(quantile(derr_w, 0.025) > 0 |
                             quantile(derr_w, 0.975) < 0, "*", "")
    )
  })
  bind_rows(out)
}

# --- NEW: Compute cumulative ΔERR from cumulative beta and theta draws ---
compute_cum_derr <- function(beta_cum_draws, theta_cum_draws, scale_by = 1) {
  rr_base <- exp(beta_cum_draws)
  derr <- rr_base * (exp(theta_cum_draws * scale_by) - 1) * 100
  derr
}


# ================================================================
#  PART 1: DATA PREPARATION
# ================================================================

cat("\n========== PART 1: DATA PREPARATION ==========\n")
event_panel <- readRDS(PATH_EVENT_PANEL)

df <- event_panel %>% 
  filter(!is.na(.data[[Outcome_var]]), .data[[Outcome_var]] >= 0,
         obs_days > 0, data_coverage >= MIN_COVERAGE)

pop_col_name <- POP_COL_MAP[Outcome_var]
if (is.na(pop_col_name)) pop_col_name <- "pop_total"

if (pop_col_name != "pop_total" && pop_col_name %in% names(df)) {
  df$pop_denom <- df[[pop_col_name]]
  n_miss_denom <- sum(is.na(df$pop_denom) | df$pop_denom <= 0, na.rm = TRUE)
  if (n_miss_denom > 0) {
    cat(sprintf("  %d rows missing '%s'; falling back to pop_total\n",
                n_miss_denom, pop_col_name))
    idx_miss <- is.na(df$pop_denom) | df$pop_denom <= 0
    df$pop_denom[idx_miss] <- df$pop_total[idx_miss]
  }
  cat(sprintf("  Using population denominator: %s\n", pop_col_name))
} else {
  df$pop_denom <- df$pop_total
  cat(sprintf("  Using population denominator: pop_total\n"))
}
df <- df %>% filter(!is.na(pop_denom), pop_denom > 0)

event_total <- df %>% group_by(event_id) %>%
  summarise(tot = sum(.data[[Outcome_var]], na.rm = TRUE), .groups = "drop")
low_events <- event_total$event_id[event_total$tot < MIN_TOTAL_COUNT]
cat(sprintf("  Min-count filter (>=%d): removing %d of %d events\n",
            MIN_TOTAL_COUNT, length(low_events), n_distinct(df$event_id)))
df <- df %>% filter(!(event_id %in% low_events))

df$Y             <- df[[Outcome_var]]
df$log_offset    <- log(df$obs_days) + log(df$pop_denom / 1e5)
df$city_id       <- as.integer(factor(df$CityCode))
df$year_idx      <- as.integer(factor(df$tc_year))
df$obs_cal_month <- as.integer(month(df$window_start))
df$window_f      <- relevel(factor(as.character(df$window_f),
                                   levels = window_labels), ref = REF_MAIN)

df$temp_round <- as.integer(round(df$temp_mean))
temp_min_val  <- min(df$temp_round, na.rm = TRUE)
df$temp_idx   <- df$temp_round - temp_min_val + 1L
df$temp_idx[is.na(df$temp_idx)] <- as.integer(median(df$temp_idx, na.rm = TRUE))

df$HSR_raw <- ifelse(!is.na(df$HSR) & df$HSR > 0, df$HSR, NA_real_)
HSR_center <- mean(df$HSR_raw, na.rm = TRUE)
df$HSR_c   <- df$HSR_raw - HSR_center

EXP_center <- mean(df$H_log, na.rm = TRUE)
EXP_sd     <- sd(df$H_log, na.rm = TRUE)
df$EXP_std <- (df$H_log - EXP_center) / EXP_sd

df$period_f <- factor(ifelse(df$tc_year <= 2020, "P2016_2020", "P2021_2024"),
                      levels = c("P2016_2020", "P2021_2024"))

cat(sprintf("  Final: N=%s | events=%d | cities=%d\n",
            format(nrow(df), big.mark = ","), n_distinct(df$event_id),
            n_distinct(df$CityCode)))
cat(sprintf("  Population denominator: %s | range=[%.0f, %.0f]\n",
            pop_col_name, min(df$pop_denom), max(df$pop_denom)))
cat(sprintf("  HSR: center=%.3f | observed range=[%.3f, %.3f] | available=%.1f%%\n",
            HSR_center, min(df$HSR_raw, na.rm = TRUE), max(df$HSR_raw, na.rm = TRUE),
            100 * mean(!is.na(df$HSR_raw))))
cat(sprintf("  EXP: center(H_log)=%.3f | SD=%.3f | range=[%.3f, %.3f]\n",
            EXP_center, EXP_sd, min(df$H_log, na.rm = TRUE), max(df$H_log, na.rm = TRUE)))
cat(sprintf("  cor(HSR, Experience) = %.2f\n",
            cor(df$HSR_c, df$EXP_std, use = "complete.obs")))

RE_TERMS <- paste(
  'f(city_id, model = "iid",
     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
  'f(obs_cal_month, model = "rw1", cyclic = TRUE, constr = TRUE,
     hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))',
  'f(year_idx, model = "rw1",
     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
  'f(temp_idx, model = "rw2", constr = TRUE,
     hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))',
  sep = " + ")
RE_TERMS_NOYEAR <- paste(
  'f(city_id, model = "iid",
     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
  'f(obs_cal_month, model = "rw1", cyclic = TRUE, constr = TRUE,
     hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))',
  'f(temp_idx, model = "rw2", constr = TRUE,
     hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))',
  sep = " + ")

F_EVENT   <- as.formula(paste("Y ~ window_f +", RE_TERMS))
F_HSR_INT <- as.formula(paste("Y ~ window_f * HSR_c +", RE_TERMS))
F_EXP_INT <- as.formula(paste("Y ~ window_f * EXP_std +", RE_TERMS))
F_PERIOD  <- as.formula(paste("Y ~ window_f * period_f +", RE_TERMS_NOYEAR))


# ================================================================
#  PART 2: MAIN MODEL
# ================================================================

cat("\n========== PART 2: MAIN MODEL ==========\n")

t0 <- Sys.time()
fit_main <- inla(F_EVENT, family = "nbinomial", data = df, offset = log_offset,
                 control.compute   = list(config = TRUE, waic = TRUE, dic = TRUE),
                 control.predictor = list(link = 1),
                 control.inla      = list(strategy = "simplified.laplace"),
                 control.fixed     = PRIOR_FIXED,
                 control.family    = CTRL_FAMILY,
                 num.threads = "4:1", verbose = FALSE)
cat(sprintf("  Time: %.1f min\n", difftime(Sys.time(), t0, units = "mins")))

res_main <- extract_window_effects(fit_main, REF_MAIN, window_labels)
print(res_main)
print(fit_main$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")])

# ---------------------------------------------------------------
# Fig2.1-01: Window-specific excess relative risk (ERR, %) of
#   post-tropical cyclone hospitalisation relative to the pre-TC
#   reference period (-3 months). The dashed vertical line indicates
#   the active TC exposure period. Error bars denote 95% credible
#   intervals from a Bayesian negative binomial event-study model.
# ---------------------------------------------------------------
pd <- res_main; pd$x_pos <- as.numeric(pd$window)
active_x <- which(window_labels == "Active")
p_main <- ggplot(pd, aes(x = x_pos, y = ERR_pct)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9B59B6", alpha = 0.20) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15,
                colour = "#9B59B6", linewidth = 0.7) +
  geom_line(colour = "#6C3483", linewidth = 0.8) +
  geom_point(colour = "#6C3483", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)", y = "Excess relative risk (ERR, %)") + theme_pub()
save_both(p_main, file.path(DIR_FIG,
                            sprintf("Fig2.1-01_TC_%s_EventStudy_Main", Outcome_var)), w = 6, h = 4)
write.xlsx(pd, file.path(DIR_TABLE_OUT,
                         sprintf("Fig2.1-01_TC_%s_EventStudy_Main_PlotData.xlsx", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-02: Posterior marginal densities of the window-specific
#   excess relative risk (ERR, %) for each post-TC period (Active
#   through +6 months). Annotations show posterior median, 95%
#   credible interval, and the probability of ERR > 0.
# ---------------------------------------------------------------
marg_list <- lapply(POST_WINDOWS, function(wl) {
  m <- fit_main$marginals.fixed[[paste0("window_f", wl)]]
  if (is.null(m)) return(NULL)
  data.frame(window = wl, ERR = (exp(m[, "x"]) - 1) * 100, density = m[, "y"])
})
marg_df <- bind_rows(marg_list)
marg_df$window <- factor(marg_df$window, levels = POST_WINDOWS)

marg_stats <- bind_rows(lapply(POST_WINDOWS, function(wl) {
  m <- fit_main$marginals.fixed[[paste0("window_f", wl)]]
  if (is.null(m)) return(NULL)
  data.frame(
    window       = wl,
    ERR_median   = (exp(inla.qmarginal(0.5,   m)) - 1) * 100,
    ERR_mean     = (inla.emarginal(exp, m) - 1) * 100,
    ERR_lo_2.5   = (exp(inla.qmarginal(0.025, m)) - 1) * 100,
    ERR_hi_97.5  = (exp(inla.qmarginal(0.975, m)) - 1) * 100,
    prob_ERR_pos = 1 - inla.pmarginal(0, m))
}))
marg_stats$window <- factor(marg_stats$window, levels = POST_WINDOWS)
cat("\n  Posterior summaries for reporting:\n"); print(marg_stats)

marg_stats$anno <- sprintf("%.2f%% (%.2f, %.2f)\nP(ERR > 0) = %s",
                           marg_stats$ERR_median, marg_stats$ERR_lo_2.5,
                           marg_stats$ERR_hi_97.5, fmt_p(marg_stats$prob_ERR_pos))

p_marg <- ggplot(marg_df, aes(x = ERR, y = density)) +
  geom_area(fill = "#9B59B6", alpha = 0.35) +
  geom_vline(data = marg_stats, aes(xintercept = ERR_median),
             colour = "#6C3483", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = marg_stats, aes(x = Inf, y = Inf, label = anno),
            inherit.aes = FALSE, hjust = 1.02, vjust = 1.25,
            size = 2.5, lineheight = 0.95) +
  facet_wrap(~ window, scales = "free") +
  labs(x = "ERR (%)", y = "Posterior density") + theme_pub()
save_both(p_marg, file.path(DIR_FIG,
                            sprintf("Fig2.1-02_TC_%s_PosteriorDensity_Windows", Outcome_var)), w = 9, h = 5)
write.xlsx(marg_df, file.path(DIR_TABLE_OUT,
                              sprintf("Fig2.1-02_TC_%s_PosteriorDensity_Windows_PlotData.xlsx", Outcome_var)))
write.xlsx(marg_stats, file.path(DIR_TABLE_OUT,
                                 sprintf("Tab2.1-01_TC_%s_Posterior_Window_Stats.xlsx", Outcome_var)))


# ================================================================
#  PART 3: SENSITIVITY ANALYSES + WAIC
# ================================================================

cat("\n========== PART 3: SENSITIVITY MODELS ==========\n")

df_s1 <- df
df_s1$window_f <- relevel(factor(as.character(df_s1$window_f),
                                 levels = window_labels), ref = REF_SENS1)
fit_sens1 <- inla(F_EVENT, family = "nbinomial", data = df_s1, offset = log_offset,
                  control.compute   = list(config = TRUE, waic = TRUE, dic = TRUE),
                  control.predictor = list(link = 1),
                  control.inla      = list(strategy = "simplified.laplace"),
                  control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                  num.threads = "4:1", verbose = FALSE)
res_sens1 <- extract_window_effects(fit_sens1, REF_SENS1, window_labels)
cat("\n  --- Sens1 (ref = -2mo) ---\n"); print(res_sens1)

fit_sens2 <- inla(F_EVENT, family = "nbinomial", data = df, offset = log_offset,
                  control.compute   = list(config = TRUE, waic = TRUE, dic = TRUE),
                  control.predictor = list(link = 1),
                  control.inla      = list(strategy = "gaussian"),
                  control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                  num.threads = "4:1", verbose = FALSE)
res_sens2 <- extract_window_effects(fit_sens2, REF_MAIN, window_labels)
cat("\n  --- Sens2 (gaussian) ---\n"); print(res_sens2)

df_s3 <- df %>% filter(!(tc_year %in% 2020:2022))
df_s3$city_id  <- as.integer(factor(df_s3$CityCode))
df_s3$year_idx <- as.integer(factor(df_s3$tc_year))
fit_sens3 <- inla(F_EVENT, family = "nbinomial", data = df_s3, offset = log_offset,
                  control.compute   = list(config = TRUE, waic = TRUE, dic = TRUE),
                  control.predictor = list(link = 1),
                  control.inla      = list(strategy = "simplified.laplace"),
                  control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                  num.threads = "4:1", verbose = FALSE)
res_sens3 <- extract_window_effects(fit_sens3, REF_MAIN, window_labels)
cat("\n  --- Sens3 (excl. 2020-2022) ---\n"); print(res_sens3)

get_waic <- function(fit) if (!is.null(fit$waic)) round(fit$waic$waic, 1) else NA
waic_df <- data.frame(
  Model = c("Main", "Sens1 (ref=-2mo)", "Sens2 (gaussian)", "Sens3 (excl.2020-2022)"),
  WAIC  = c(get_waic(fit_main), get_waic(fit_sens1),
            get_waic(fit_sens2), get_waic(fit_sens3)))
print(waic_df)


# ================================================================
#  PART 4: CUMULATIVE ERR (JOINT POSTERIOR SAMPLING)
# ================================================================

cat("\n========== PART 4: CUMULATIVE ERR ==========\n")

last_sig <- "Active"
for (wl in POST_WINDOWS) {
  r <- res_main[as.character(res_main$window) == wl, ]
  if (nrow(r) == 1 && r$significant == "*") last_sig <- wl
}
sig_windows <- POST_WINDOWS[1:which(POST_WINDOWS == last_sig)]
cumulative_periods <- list(Active_to_PostFull = POST_WINDOWS,
                           Active_to_Significant = sig_windows)
cat(sprintf("  Last significant window: %s\n", last_sig))

compute_cum <- function(fit, ref, seed) {
  out <- list(results = NULL, samples = list())
  for (nm in names(cumulative_periods)) {
    coefs <- paste0("window_f", setdiff(cumulative_periods[[nm]], ref))
    M <- get_coef_draws(fit, coefs, N_SAMPLES, seed)
    if (is.null(M)) next
    err <- (exp(colSums(M)) - 1) * 100
    out$results <- rbind(out$results,
                         cbind(data.frame(period = nm,
                                          n_windows = length(cumulative_periods[[nm]])),
                               summarise_draws(err)))
    out$samples[[nm]] <- err
  }
  out
}
cum_main  <- compute_cum(fit_main,  REF_MAIN,  123)
cum_sens1 <- compute_cum(fit_sens1, REF_SENS1, 124)
cum_sens2 <- compute_cum(fit_sens2, REF_MAIN,  125)
cum_sens3 <- compute_cum(fit_sens3, REF_MAIN,  126)
cat("\n  Main model cumulative ERR:\n"); print(cum_main$results)

# ---------------------------------------------------------------
# Fig2.1-03: Posterior density of the cumulative post-TC excess
#   relative risk (ERR, %), defined as (exp(Σβ_w) − 1) × 100%
#   summed over Active to +6 months. The solid line marks the
#   posterior median; the dashed line marks ERR = 0.
# ---------------------------------------------------------------
cum_draw_df <- data.frame(ERR = cum_main$samples[["Active_to_PostFull"]])
cum_x <- cum_draw_df$ERR
anno_cum <- sprintf("%.2f%% (%.2f, %.2f)\n%s",
                    median(cum_x), quantile(cum_x, .025, names = FALSE),
                    quantile(cum_x, .975, names = FALSE), ptxt(cum_x, "gt0"))
p_cum <- ggplot(cum_draw_df, aes(x = ERR)) +
  geom_density(fill = "#9B59B6", alpha = 0.35, colour = "#9B59B6") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = median(cum_x), colour = "#9B59B6") +
  annotate("text", x = Inf, y = Inf, label = anno_cum,
           hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
  labs(x = sprintf("Cumulative Post-TC ERR (%%)"),
       y = "Posterior density") + theme_pub()
save_both(p_cum, file.path(DIR_FIG,
                           sprintf("Fig2.1-03_TC_%s_CumERR_Posterior", Outcome_var)), w = 2.5, h = 2)
save_zoom_transparent(p_cum, file.path(DIR_FIG,
                                       sprintf("Fig2.1-03_TC_%s_CumERR_Posterior", Outcome_var)))
write.xlsx(cum_draw_df, file.path(DIR_TABLE_OUT,
                                  sprintf("Fig2.1-03_TC_%s_CumERR_PosteriorDraws.xlsx", Outcome_var)))


# ================================================================
#  PART 5: FORMAL MODEL COMPARISON
# ================================================================

cat("\n========== PART 5: FORMAL MODEL COMPARISON ==========\n")

compare_two <- function(sA, sB, labA, labB, period) {
  if (!period %in% names(sA) || !period %in% names(sB)) return(NULL)
  d <- sA[[period]] - sB[[period]]; p <- mean(sA[[period]] > sB[[period]])
  cat(sprintf("  %s vs %s [%s]: P(A>B)=%.3f | diff=%+.2f%% [%+.2f, %+.2f]\n",
              labA, labB, period, p, median(d), quantile(d, .025), quantile(d, .975)))
  data.frame(comparison = paste(labA, "vs", labB), period = period,
             prob_A_gt_B = p, diff_median = median(d),
             diff_lo = quantile(d, .025, names = FALSE),
             diff_hi = quantile(d, .975, names = FALSE),
             agreement = ifelse(p > 0.2 & p < 0.8, "Agree", "Differ"))
}
comp_all <- bind_rows(lapply(names(cumulative_periods), function(cp) bind_rows(
  compare_two(cum_main$samples, cum_sens1$samples, "Main", "Sens1", cp),
  compare_two(cum_main$samples, cum_sens2$samples, "Main", "Sens2", cp),
  compare_two(cum_main$samples, cum_sens3$samples, "Main", "Sens3", cp))))


# ================================================================
#  PART 6a: HSR × WINDOW INTERACTION (ALL CITIES)
# ================================================================

cat("\n========== PART 6a: HSR x WINDOW INTERACTION ==========\n")

df_hsr <- df %>% filter(!is.na(HSR_c))
df_hsr$city_id  <- as.integer(factor(df_hsr$CityCode))
df_hsr$year_idx <- as.integer(factor(df_hsr$tc_year))

fit_hsr <- inla(F_HSR_INT, family = "nbinomial", data = df_hsr, offset = log_offset,
                control.compute   = list(config = TRUE, waic = TRUE),
                control.predictor = list(link = 1),
                control.inla      = list(strategy = "simplified.laplace"),
                control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                num.threads = "4:1", verbose = FALSE)

hsr_int <- extract_interaction_effects(fit_hsr, "HSR_c", window_labels,
                                       REF_MAIN, scale_by = HSR_STEP)
cat("\n  HSR interaction (per +0.1 HSR):\n"); print(hsr_int)

# ---------------------------------------------------------------
# Fig2.1-04: Window-specific percentage change in the rate ratio
#   (% RR change) per +0.1 unit increase in the health system
#   resilience (HSR) index. Each point represents the posterior
#   mean of (exp(θ_w × 0.1) − 1) × 100% for window w, with 95%
#   credible intervals. Negative values indicate HSR is protective.
# ---------------------------------------------------------------
pd <- hsr_int; pd$x_pos <- as.numeric(pd$window)
p_theta <- ggplot(pd, aes(x = x_pos, y = pct_change)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = pct_change_lo, ymax = pct_change_hi),
                width = 0.15, colour = "#2166AC", linewidth = 0.7) +
  geom_point(colour = "#2166AC", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = "Change in window-specific RR per +0.1 HSR (%)") + theme_pub()
save_both(p_theta, file.path(DIR_FIG,
                             sprintf("Fig2.1-04_TC_%s_HSR_Interaction_RRpct", Outcome_var)), w = 6, h = 4)
write.xlsx(pd, file.path(DIR_TABLE_OUT,
                         sprintf("Fig2.1-04_TC_%s_HSR_Interaction_RRpct_PlotData.xlsx", Outcome_var)))

# --- Extract beta and theta draws for full HSR model ---
beta_names  <- paste0("window_f", POST_WINDOWS)
theta_names <- paste0("window_f", POST_WINDOWS, ":HSR_c")
M_beta  <- get_coef_draws(fit_hsr, beta_names,  N_SAMPLES, 321)
M_theta <- get_coef_draws(fit_hsr, theta_names, N_SAMPLES, 321)

beta_cum_draws  <- colSums(M_beta)
theta_cum_draws <- colSums(M_theta)

# --- Cumulative % RR change per +0.1 HSR ---
hsr_cum_rr_draws <- (exp(theta_cum_draws * HSR_STEP) - 1) * 100
hsr_cum_mod <- summarise_draws(hsr_cum_rr_draws)
cat(sprintf("\n  Cumulative %%RR change per +0.1 HSR: %+.2f%% [%+.2f, %+.2f], P(protective)=%.3f\n",
            hsr_cum_mod$median, hsr_cum_mod$lower_2.5, hsr_cum_mod$upper_97.5,
            1 - hsr_cum_mod$prob_positive))

# --- [NEW] Cumulative ΔERR (pp) per +0.1 HSR ---
hsr_cum_derr_draws <- compute_cum_derr(beta_cum_draws, theta_cum_draws, scale_by = HSR_STEP)
hsr_cum_derr_summ  <- summarise_draws(hsr_cum_derr_draws)
cat(sprintf("  Cumulative ΔERR per +0.1 HSR: %+.2f pp [%+.2f, %+.2f], P(protective)=%.3f\n",
            hsr_cum_derr_summ$median, hsr_cum_derr_summ$lower_2.5,
            hsr_cum_derr_summ$upper_97.5, 1 - hsr_cum_derr_summ$prob_positive))

# --- Implied effect of observed 0.2 HSR growth ---
d02 <- (exp(theta_cum_draws * 0.2) - 1) * 100
cat(sprintf("  Implied effect of +0.2 HSR (%%RR): %+.2f%% [%+.2f, %+.2f]\n",
            median(d02), quantile(d02, .025), quantile(d02, .975)))
d02_derr <- compute_cum_derr(beta_cum_draws, theta_cum_draws, scale_by = 0.2)
cat(sprintf("  Implied effect of +0.2 HSR (ΔERR): %+.2f pp [%+.2f, %+.2f]\n",
            median(d02_derr), quantile(d02_derr, .025), quantile(d02_derr, .975)))

# ---------------------------------------------------------------
# Fig2.1-05: Posterior density of the cumulative relative change
#   in the post-TC rate ratio (% RR change) per +0.1 unit HSR,
#   estimated from the full-sample HSR × window interaction model.
#   Negative values indicate that higher HSR attenuates the
#   cumulative post-TC hospitalisation risk.
# ---------------------------------------------------------------
tc_draws <- hsr_cum_rr_draws
anno_tc <- sprintf("%.2f%% (%.2f, %.2f)\n%s",
                   median(tc_draws), quantile(tc_draws, .025, names = FALSE),
                   quantile(tc_draws, .975, names = FALSE), ptxt(tc_draws, "lt0"))
p_tc <- ggplot(data.frame(x = tc_draws), aes(x = x)) +
  geom_density(fill = "#9B59B6", alpha = 0.35, colour = "#9B59B6") +
  geom_vline(xintercept = median(tc_draws), colour = "#9B59B6", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = Inf, y = Inf, label = anno_tc,
           hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
  labs(x = "Relative change in cumulative post-TC rate ratio per +0.1 HSR (%)",
       y = "Posterior density") + theme_pub()
save_both(p_tc, file.path(DIR_FIG,
                          sprintf("Fig2.1-05_TC_%s_HSR_CumRRpct_Posterior", Outcome_var)), w = 5, h = 4)
save_zoom_transparent(p_tc, file.path(DIR_FIG,
                                      sprintf("Fig2.1-05_TC_%s_HSR_CumRRpct_Posterior", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-05b: Posterior density of the cumulative absolute change
#   in post-TC ERR (ΔERR, percentage points) per +0.1 unit HSR.
#   ΔERR = RR_base × (exp(Θ_cum × 0.1) − 1) × 100, where RR_base
#   is the cumulative rate ratio at the mean HSR. Negative values
#   indicate that higher HSR reduces the absolute excess risk.
# ---------------------------------------------------------------
anno_derr_hsr <- sprintf("%.2f %% (%.2f, %.2f)\n%s",
                         median(hsr_cum_derr_draws), quantile(hsr_cum_derr_draws, .025, names = FALSE),
                         quantile(hsr_cum_derr_draws, .975, names = FALSE),
                         ptxt(hsr_cum_derr_draws, "lt0"))
p_derr_hsr <- ggplot(data.frame(x = hsr_cum_derr_draws), aes(x = x)) +
  geom_density(fill = "#9B59B6", alpha = 0.35, colour = "#9B59B6") +
  geom_vline(xintercept = median(hsr_cum_derr_draws), colour = "#9B59B6", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = Inf, y = Inf, label = anno_derr_hsr,
           hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
  labs(x = "ΔCumulative ERR per +0.1 HSR (%)",
       y = "Posterior density") + theme_pub()
save_both(p_derr_hsr, file.path(DIR_FIG,
                                sprintf("Fig2.1-05b_TC_%s_HSR_CumDERR_Posterior", Outcome_var)), w = 3, h = 2)
save_zoom_transparent(p_derr_hsr, file.path(DIR_FIG,
                                            sprintf("Fig2.1-05b_TC_%s_HSR_CumDERR_Posterior", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-05c: Window-specific absolute change in ERR (ΔERR, pp)
#   per +0.1 unit HSR across the post-TC timeline. ΔERR_w =
#   exp(β_w) × (exp(θ_w × 0.1) − 1) × 100. This quantifies the
#   absolute attenuation of window-specific ERR by HSR.
# ---------------------------------------------------------------
hsr_win_derr <- compute_window_derr(M_beta, M_theta, scale_by = HSR_STEP)
hsr_win_derr$window <- POST_WINDOWS
hsr_win_derr$x_pos  <- as.numeric(factor(POST_WINDOWS, levels = POST_WINDOWS)) +
  N_PRE_WINDOWS + 1

p_derr_hsr_win <- ggplot(hsr_win_derr, aes(x = x_pos, y = derr_median)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = derr_lo, ymax = derr_hi),
                width = 0.15, colour = "#2166AC", linewidth = 0.7) +
  geom_point(colour = "#2166AC", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = "Window-specific ΔERR per +0.1 HSR (pp)") + theme_pub()
save_both(p_derr_hsr_win, file.path(DIR_FIG,
                                    sprintf("Fig2.1-05c_TC_%s_HSR_WindowDERR", Outcome_var)), w = 6, h = 4)
write.xlsx(hsr_win_derr, file.path(DIR_TABLE_OUT,
                                   sprintf("Fig2.1-05c_TC_%s_HSR_WindowDERR_PlotData.xlsx", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-06: Predicted cumulative post-TC ERR (%) as a function
#   of the HSR index, with 95% credible band. The rug plot shows
#   the observed distribution of city-year HSR values.
# ---------------------------------------------------------------
hsr_grid <- seq(min(df_hsr$HSR_raw, na.rm = TRUE),
                max(df_hsr$HSR_raw, na.rm = TRUE), length.out = 80)

hsr_curve <- bind_rows(lapply(hsr_grid, function(h) {
  err <- (exp(beta_cum_draws + (h - HSR_center) * theta_cum_draws) - 1) * 100
  data.frame(HSR = h, ERR = median(err),
             lo = quantile(err, .025, names = FALSE),
             hi = quantile(err, .975, names = FALSE))
}))
p_curve <- ggplot(hsr_curve, aes(x = HSR, y = ERR)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#9B59B6", alpha = 0.2) +
  geom_line(colour = "#9B59B6", linewidth = 0.9) +
  geom_rug(data = df_hsr %>% distinct(CityCode, tc_year, HSR_raw),
           aes(x = HSR_raw), inherit.aes = FALSE, alpha = 0.15,
           length = unit(0.02, "npc")) +
  labs(x = "Health system resilience index",
       y = "Cumulative post-TC ERR (%)") + theme_pub()
save_both(p_curve, file.path(DIR_FIG,
                             sprintf("Fig2.1-06_TC_%s_CumERR_vs_HSR", Outcome_var)), w = 5, h = 4)
write.xlsx(hsr_curve, file.path(DIR_TABLE_OUT,
                                sprintf("Fig2.1-06_TC_%s_CumERR_vs_HSR_PlotData.xlsx", Outcome_var)))

saveRDS(list(beta_cum_draws = beta_cum_draws, theta_cum_draws = theta_cum_draws,
             HSR_center = HSR_center,
             HSR_observed_range = range(df_hsr$HSR_raw, na.rm = TRUE),
             cum_rr_pct_draws = hsr_cum_rr_draws,
             cum_derr_draws   = hsr_cum_derr_draws,
             note = "cumERR(h) = (exp(beta_cum + (h - HSR_center)*theta_cum) - 1)*100"),
        file.path(DIR_OUTPUT, sprintf("TC_%s_HSR_Projection_Draws.rds", Outcome_var)))


# ================================================================
#  PART 6b: HSR-STRATIFIED MODELS
# ================================================================

cat("\n========== PART 6b: HSR-STRATIFIED MODELS ==========\n")

hsr_strata <- list()
for (lvl in c("High", "Moderate", "Low")) {
  d <- df_hsr %>% filter(HSR_Level == lvl)
  if (nrow(d) < 200 || n_distinct(d$event_id) < 10) {
    cat(sprintf("  HSR=%s: insufficient data, skipped\n", lvl)); next }
  d$city_id  <- as.integer(factor(d$CityCode))
  d$year_idx <- as.integer(factor(d$tc_year))
  cat(sprintf("\n  --- HSR = %s | N=%s | events=%d ---\n", lvl,
              format(nrow(d), big.mark = ","), n_distinct(d$event_id)))
  fit_s <- inla(F_EVENT, family = "nbinomial", data = d, offset = log_offset,
                control.compute   = list(config = TRUE, waic = TRUE),
                control.predictor = list(link = 1),
                control.inla      = list(strategy = "simplified.laplace"),
                control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                num.threads = "4:1", verbose = FALSE)
  we <- extract_window_effects(fit_s, REF_MAIN, window_labels)
  print(we)
  M <- get_coef_draws(fit_s, paste0("window_f", POST_WINDOWS), N_SAMPLES,
                      200 + which(c("High","Moderate","Low") == lvl))
  hsr_strata[[lvl]] <- list(window_effects = we,
                            cum_draws = (exp(colSums(M)) - 1) * 100,
                            beta_cum_draws = colSums(M),
                            weight = sum(d$Y, na.rm = TRUE))
}

hsr_comp <- data.frame()
for (pair in list(c("High","Low"), c("High","Moderate"), c("Moderate","Low"))) {
  a <- hsr_strata[[pair[1]]]$cum_draws; b <- hsr_strata[[pair[2]]]$cum_draws
  if (is.null(a) || is.null(b)) next
  d <- a - b
  cat(sprintf("  HSR %s vs %s: P(%s < %s)=%.3f | diff=%+.2f%% [%+.2f, %+.2f]\n",
              pair[1], pair[2], pair[1], pair[2], mean(a < b),
              median(d), quantile(d, .025), quantile(d, .975)))
  hsr_comp <- rbind(hsr_comp, data.frame(
    comparison = paste("HSR", pair[1], "vs", pair[2]),
    prob_first_lower = mean(a < b), diff_median = median(d),
    diff_lo = quantile(d, .025, names = FALSE),
    diff_hi = quantile(d, .975, names = FALSE)))
}

if (length(hsr_strata) == 3) {
  wts <- sapply(hsr_strata, `[[`, "weight"); wts <- wts / sum(wts)
  pooled <- Reduce(`+`, Map(function(s, w) s$cum_draws * w, hsr_strata, wts))
  cat(sprintf("\n  Weighted pooling of strata: %+.2f%% [%+.2f, %+.2f]\n",
              median(pooled), quantile(pooled, .025), quantile(pooled, .975)))
}

# ---------------------------------------------------------------
# Fig2.1-07: Window-specific ERR (%) stratified by HSR level
#   (High, Moderate, Low). Each line represents the posterior mean
#   ERR with 95% CrI from stratum-specific Bayesian event-study
#   models. Comparison across strata illustrates how health system
#   resilience modulates the post-TC hospitalisation trajectory.
# ---------------------------------------------------------------
if (length(hsr_strata) >= 2) {
  pd <- bind_rows(lapply(names(hsr_strata), function(l) {
    e <- hsr_strata[[l]]$window_effects; e$HSR_Level <- l
    e$x_pos <- as.numeric(e$window); e }))
  p_strat <- ggplot(pd, aes(x = x_pos, y = ERR_pct, colour = HSR_Level)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15,
                  position = position_dodge(0.4), linewidth = 0.6) +
    geom_point(position = position_dodge(0.4), size = 2.2) +
    geom_line(position = position_dodge(0.4), linewidth = 0.7) +
    scale_colour_manual(values = c(High = "#2166AC", Moderate = "#F4A582",
                                   Low = "#B2182B"), name = "HSR level") +
    scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
    labs(x = "Lag period (month)", y = "ERR (%)") + theme_pub()
  save_both(p_strat, file.path(DIR_FIG,
                               sprintf("Fig2.1-07_TC_%s_EventStudy_HSR_Strata", Outcome_var)), w = 9, h = 5)
  write.xlsx(pd, file.path(DIR_TABLE_OUT,
                           sprintf("Fig2.1-07_TC_%s_EventStudy_HSR_Strata_PlotData.xlsx", Outcome_var)))
}


# ================================================================
#  PART 7a: EXPERIENCE × WINDOW INTERACTION (STANDARDISED SCALE)
# ================================================================

cat("\n========== PART 7a: EXPERIENCE x WINDOW INTERACTION ==========\n")

fit_exp <- inla(F_EXP_INT, family = "nbinomial", data = df, offset = log_offset,
                control.compute   = list(config = TRUE, waic = TRUE),
                control.predictor = list(link = 1),
                control.inla      = list(strategy = "simplified.laplace"),
                control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                num.threads = "4:1", verbose = FALSE)
exp_int <- extract_interaction_effects(fit_exp, "EXP_std", window_labels, REF_MAIN)
cat("\n  Experience interaction (per +1 SD experience):\n"); print(exp_int)

M_exp <- get_coef_draws(fit_exp, paste0("window_f", POST_WINDOWS, ":EXP_std"),
                        N_SAMPLES, 331)
exp_cum_draws <- (exp(colSums(M_exp)) - 1) * 100
exp_cum_mod <- summarise_draws(exp_cum_draws)
cat(sprintf("  Cumulative modification per +1 SD experience: %+.2f%% [%+.2f, %+.2f], P(adaptation)=%.3f\n",
            exp_cum_mod$median, exp_cum_mod$lower_2.5, exp_cum_mod$upper_97.5,
            1 - exp_cum_mod$prob_positive))


# ================================================================
#  PART 7a-2: EXPERIENCE INTERACTION — RAW H_log SCALE
# ================================================================

cat("\n========== PART 7a-2: EXPERIENCE (RAW H_log SCALE) ==========\n")

EXPC_center <- mean(df$H_log, na.rm = TRUE)
df$EXP_c    <- df$H_log - EXPC_center

F_EXP_INT_RAW <- as.formula(paste("Y ~ window_f * EXP_c +", RE_TERMS))

fit_exp_raw <- inla(F_EXP_INT_RAW, family = "nbinomial", data = df, offset = log_offset,
                    control.compute   = list(config = TRUE, waic = TRUE),
                    control.predictor = list(link = 1),
                    control.inla      = list(strategy = "simplified.laplace"),
                    control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                    num.threads = "4:1", verbose = FALSE)

exp_int_raw <- extract_interaction_effects(fit_exp_raw, "EXP_c", window_labels,
                                           REF_MAIN, scale_by = 1)
cat("\n  Experience interaction (per +1 unit H_log):\n"); print(exp_int_raw)

beta_names_e  <- paste0("window_f", POST_WINDOWS)
theta_names_e <- paste0("window_f", POST_WINDOWS, ":EXP_c")
M_beta_e  <- get_coef_draws(fit_exp_raw, beta_names_e,  N_SAMPLES, 351)
M_theta_e <- get_coef_draws(fit_exp_raw, theta_names_e, N_SAMPLES, 351)
beta_cum_e  <- colSums(M_beta_e)
theta_cum_e <- colSums(M_theta_e)

# --- Cumulative % RR change per +1 H_log ---
exp_cum_raw_draws <- (exp(theta_cum_e) - 1) * 100
exp_cum_raw <- summarise_draws(exp_cum_raw_draws)
cat(sprintf("  Cumulative %%RR change per +1 H_log: %+.2f%% [%+.2f, %+.2f], P(adaptation)=%.3f\n",
            exp_cum_raw$median, exp_cum_raw$lower_2.5, exp_cum_raw$upper_97.5,
            1 - exp_cum_raw$prob_positive))

# --- [NEW] Cumulative ΔERR (pp) per +1 H_log ---
exp_cum_derr_draws <- compute_cum_derr(beta_cum_e, theta_cum_e, scale_by = 1)
exp_cum_derr_summ  <- summarise_draws(exp_cum_derr_draws)
cat(sprintf("  Cumulative ΔERR per +1 H_log: %+.2f pp [%+.2f, %+.2f], P(adaptation)=%.3f\n",
            exp_cum_derr_summ$median, exp_cum_derr_summ$lower_2.5,
            exp_cum_derr_summ$upper_97.5, 1 - exp_cum_derr_summ$prob_positive))

# ---------------------------------------------------------------
# Fig2.1-08a: Window-specific percentage change in the rate ratio
#   (% RR change) per +1 unit increase in the log-transformed
#   historical TC experience index (H_log). Negative values
#   indicate that prior TC exposure is protective for that window.
# ---------------------------------------------------------------
pd_exp_win <- exp_int_raw; pd_exp_win$x_pos <- as.numeric(pd_exp_win$window)
p_exp_win_rr <- ggplot(pd_exp_win, aes(x = x_pos, y = pct_change)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = pct_change_lo, ymax = pct_change_hi),
                width = 0.15, colour = "#E18727", linewidth = 0.7) +
  geom_point(colour = "#E18727", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = "Change in window-specific RR per +1 H_log (%)") + theme_pub()
save_both(p_exp_win_rr, file.path(DIR_FIG,
                                  sprintf("Fig2.1-08a_TC_%s_EXP_Interaction_RRpct", Outcome_var)), w = 6, h = 4)
write.xlsx(pd_exp_win, file.path(DIR_TABLE_OUT,
                                 sprintf("Fig2.1-08a_TC_%s_EXP_Interaction_RRpct_PlotData.xlsx", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-08b: Window-specific absolute change in ERR (ΔERR, pp)
#   per +1 unit H_log. ΔERR_w = exp(β_w) × (exp(θ_w) − 1) × 100.
#   This captures the absolute ERR reduction attributable to prior
#   TC experience for each post-TC window.
# ---------------------------------------------------------------
exp_win_derr <- compute_window_derr(M_beta_e, M_theta_e, scale_by = 1)
exp_win_derr$window <- POST_WINDOWS
exp_win_derr$x_pos  <- as.numeric(factor(POST_WINDOWS, levels = POST_WINDOWS)) +
  N_PRE_WINDOWS + 1

p_exp_win_derr <- ggplot(exp_win_derr, aes(x = x_pos, y = derr_median)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = derr_lo, ymax = derr_hi),
                width = 0.15, colour = "#E18727", linewidth = 0.7) +
  geom_point(colour = "#E18727", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = "Window-specific ΔERR per +1 H_log (pp)") + theme_pub()
save_both(p_exp_win_derr, file.path(DIR_FIG,
                                    sprintf("Fig2.1-08b_TC_%s_EXP_WindowDERR", Outcome_var)), w = 6, h = 4)
write.xlsx(exp_win_derr, file.path(DIR_TABLE_OUT,
                                   sprintf("Fig2.1-08b_TC_%s_EXP_WindowDERR_PlotData.xlsx", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-08c: Posterior density of the cumulative % RR change per
#   +1 unit H_log (full-sample experience interaction model).
#   Negative values indicate experience-based adaptation.
# ---------------------------------------------------------------
anno_exp_rr <- sprintf("%.2f%% (%.2f, %.2f)\n%s",
                       median(exp_cum_raw_draws), quantile(exp_cum_raw_draws, .025, names = FALSE),
                       quantile(exp_cum_raw_draws, .975, names = FALSE),
                       ptxt(exp_cum_raw_draws, "lt0"))
p_exp_cum_rr <- ggplot(data.frame(x = exp_cum_raw_draws), aes(x = x)) +
  geom_density(fill = "#E18727", alpha = 0.35, colour = "#E18727") +
  geom_vline(xintercept = median(exp_cum_raw_draws), colour = "#E18727", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = Inf, y = Inf, label = anno_exp_rr,
           hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
  labs(x = "Cumulative % RR change per +1 H_log (%)",
       y = "Posterior density") + theme_pub()
save_both(p_exp_cum_rr, file.path(DIR_FIG,
                                  sprintf("Fig2.1-08c_TC_%s_EXP_CumRRpct_Posterior", Outcome_var)), w = 5, h = 4)
save_zoom_transparent(p_exp_cum_rr, file.path(DIR_FIG,
                                              sprintf("Fig2.1-08c_TC_%s_EXP_CumRRpct_Posterior", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-08d: Posterior density of the cumulative ΔERR (pp) per
#   +1 unit H_log. ΔERR = RR_base × (exp(Θ_cum) − 1) × 100.
#   Negative values indicate that higher historical TC experience
#   reduces the absolute cumulative excess hospitalisation risk.
# ---------------------------------------------------------------
anno_exp_derr <- sprintf("%.2f pp (%.2f, %.2f)\n%s",
                         median(exp_cum_derr_draws), quantile(exp_cum_derr_draws, .025, names = FALSE),
                         quantile(exp_cum_derr_draws, .975, names = FALSE),
                         ptxt(exp_cum_derr_draws, "lt0"))
p_exp_cum_derr <- ggplot(data.frame(x = exp_cum_derr_draws), aes(x = x)) +
  geom_density(fill = "#E18727", alpha = 0.35, colour = "#E18727") +
  geom_vline(xintercept = median(exp_cum_derr_draws), colour = "#E18727", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = Inf, y = Inf, label = anno_exp_derr,
           hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
  labs(x = "ΔCumulative ERR per +1 unit \nlog-scaled TC experience index (%)",
       y = "Posterior density") + theme_pub()
save_both(p_exp_cum_derr, file.path(DIR_FIG,
                                    sprintf("Fig2.1-08d_TC_%s_EXP_CumDERR_Posterior", Outcome_var)), w = 3, h = 2)
save_zoom_transparent(p_exp_cum_derr, file.path(DIR_FIG,
                                                sprintf("Fig2.1-08d_TC_%s_EXP_CumDERR_Posterior", Outcome_var)))

# ---------------------------------------------------------------
# Fig2.1-08e: Predicted cumulative post-TC ERR (%) as a function
#   of the historical TC experience index (H_log), with 95%
#   credible band and rug plot of observed city-year values.
# ---------------------------------------------------------------
exp_grid <- seq(min(df$H_log, na.rm = TRUE),
                max(df$H_log, na.rm = TRUE), length.out = 80)

exp_curve <- bind_rows(lapply(exp_grid, function(h) {
  err <- (exp(beta_cum_e + (h - EXPC_center) * theta_cum_e) - 1) * 100
  data.frame(H_log = h, H_raw = exp(h) - 1, ERR = median(err),
             lo = quantile(err, .025, names = FALSE),
             hi = quantile(err, .975, names = FALSE))
}))

p_expcurve <- ggplot(exp_curve, aes(x = H_log, y = ERR)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#E18727", alpha = 0.2) +
  geom_line(colour = "#E18727", linewidth = 0.9) +
  geom_rug(data = df %>% distinct(CityCode, tc_year, H_log),
           aes(x = H_log), inherit.aes = FALSE, alpha = 0.15,
           length = unit(0.02, "npc")) +
  scale_x_continuous(
    sec.axis = sec_axis(~ exp(.) - 1,
                        name = "Historical TC experience index (raw H)",
                        breaks = scales::pretty_breaks(5))) +
  labs(x = expression("Historical TC experience index (log"[e]*"(H+1))"),
       y = "Cumulative post-TC ERR (%)") + theme_pub()
save_both(p_expcurve, file.path(DIR_FIG,
                                sprintf("Fig2.1-08e_TC_%s_CumERR_vs_Experience", Outcome_var)), w = 5, h = 4)
write.xlsx(exp_curve, file.path(DIR_TABLE_OUT,
                                sprintf("Fig2.1-08e_TC_%s_CumERR_vs_Experience_PlotData.xlsx", Outcome_var)))



# ---------------------------------------------------------------
# Fig2.1-08e2: Predicted cumulative post-TC ERR (%) as a function
#   of the log-transformed historical TC experience index
#   (log_e(H+1)), with 95% credible band. This is a simplified
#   version of Fig2.1-08e without the secondary raw-scale axis,
#   emphasising the near-linear relationship on the log scale
#   that arises from the linear interaction specification. The rug
#   plot shows the observed distribution of city-year H_log values.
# ---------------------------------------------------------------
p_expcurve_logonly <- ggplot(exp_curve, aes(x = H_log, y = ERR)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#E18727", alpha = 0.2) +
  geom_line(colour = "#E18727", linewidth = 0.9) +
  geom_rug(data = df %>% distinct(CityCode, tc_year, H_log),
           aes(x = H_log), inherit.aes = FALSE, alpha = 0.15,
           length = unit(0.02, "npc")) +
  labs(x = expression("Historical TC experience index (log"[e]*"(H+1))"),
       y = "Cumulative post-TC ERR (%)") + theme_pub()
save_both(p_expcurve_logonly, file.path(DIR_FIG,
                                        sprintf("Fig2.1-08e2_TC_%s_CumERR_vs_Experience_LogOnly", Outcome_var)),
          w = 5, h = 4)



# ---------------------------------------------------------------
# Fig2.1-08f: Predicted cumulative post-TC ERR (%) as a function
#   of the historical TC experience index on the raw (untransformed)
#   scale (H), with 95% credible band. Unlike Fig2.1-08e which
#   plots against log(H+1), this figure reveals the diminishing
#   marginal returns of TC experience on the natural scale — the
#   relationship is logarithmic because the model is linear in
#   log(H+1). The rug plot shows observed city-year values of H.
# ---------------------------------------------------------------
exp_grid_raw <- seq(min(exp(df$H_log) - 1, na.rm = TRUE),
                    max(exp(df$H_log) - 1, na.rm = TRUE), length.out = 100)

exp_curve_raw <- bind_rows(lapply(exp_grid_raw, function(h_raw) {
  h_log <- log(h_raw + 1)
  err <- (exp(beta_cum_e + (h_log - EXPC_center) * theta_cum_e) - 1) * 100
  data.frame(H_raw = h_raw, H_log = h_log, ERR = median(err),
             lo = quantile(err, .025, names = FALSE),
             hi = quantile(err, .975, names = FALSE))
}))

p_expcurve_raw <- ggplot(exp_curve_raw, aes(x = H_raw, y = ERR)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#E18727", alpha = 0.2) +
  geom_line(colour = "#E18727", linewidth = 0.9) +
  geom_rug(data = df %>% distinct(CityCode, tc_year, H_log) %>%
             mutate(H_raw = exp(H_log) - 1),
           aes(x = H_raw), inherit.aes = FALSE, alpha = 0.15,
           length = unit(0.02, "npc")) +
  labs(x = "Historical TC experience index (H, raw scale)",
       y = "Cumulative post-TC ERR (%)") + theme_pub()
save_both(p_expcurve_raw, file.path(DIR_FIG,
                                    sprintf("Fig2.1-08f_TC_%s_CumERR_vs_Experience_RawScale", Outcome_var)),
          w = 5, h = 4)
write.xlsx(exp_curve_raw, file.path(DIR_TABLE_OUT,
                                    sprintf("Fig2.1-08f_TC_%s_CumERR_vs_Experience_RawScale_PlotData.xlsx", Outcome_var)))





saveRDS(list(beta_cum_draws = beta_cum_e, theta_cum_draws = theta_cum_e,
             EXP_center = EXPC_center, EXP_sd = EXP_sd,
             EXP_observed_range = range(df$H_log, na.rm = TRUE),
             cum_rr_pct_draws = exp_cum_raw_draws,
             cum_derr_draws   = exp_cum_derr_draws,
             note = "cumERR(h) = (exp(beta_cum + (h - EXP_center)*theta_cum) - 1)*100"),
        file.path(DIR_OUTPUT, sprintf("TC_%s_EXP_Projection_Draws.rds", Outcome_var)))


# ================================================================
#  PART 7b: EXPERIENCE & HSR INTERACTION WITHIN HSR STRATA [v7.0]
# ================================================================

cat("\n========== PART 7b: EXPERIENCE & HSR INTERACTION WITHIN HSR STRATA ==========\n")

adapt_strata     <- list()
adapt_strata_hsr <- list()

for (lvl in c("High", "Moderate", "Low")) {
  d <- df_hsr %>% filter(HSR_Level == lvl)
  if (nrow(d) < 200 || n_distinct(d$event_id) < 10) {
    cat(sprintf("  HSR=%s: insufficient data, skipped\n", lvl)); next
  }
  d$city_id  <- as.integer(factor(d$CityCode))
  d$year_idx <- as.integer(factor(d$tc_year))
  
  cat(sprintf("\n  --- HSR = %s | N=%s | events=%d | HSR range=[%.3f, %.3f] ---\n",
              lvl, format(nrow(d), big.mark = ","), n_distinct(d$event_id),
              min(d$HSR_raw, na.rm = TRUE), max(d$HSR_raw, na.rm = TRUE)))
  
  # ---- 7b-1: Experience × window interaction (standardised) ----
  fit_a <- inla(F_EXP_INT, family = "nbinomial", data = d, offset = log_offset,
                control.compute   = list(config = TRUE),
                control.predictor = list(link = 1),
                control.inla      = list(strategy = "simplified.laplace"),
                control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                num.threads = "4:1", verbose = FALSE)
  
  int_eff_std <- extract_interaction_effects(fit_a, "EXP_std", window_labels, REF_MAIN)
  
  int_eff_hlog <- NULL
  if (!is.null(int_eff_std) && nrow(int_eff_std) > 0) {
    int_eff_hlog <- int_eff_std
    int_eff_hlog$theta          <- int_eff_std$theta       / EXP_sd
    int_eff_hlog$theta_lo       <- int_eff_std$theta_lo    / EXP_sd
    int_eff_hlog$theta_hi       <- int_eff_std$theta_hi    / EXP_sd
    int_eff_hlog$pct_change     <- (exp(int_eff_std$theta       / EXP_sd) - 1) * 100
    int_eff_hlog$pct_change_lo  <- (exp(int_eff_std$theta_lo    / EXP_sd) - 1) * 100
    int_eff_hlog$pct_change_hi  <- (exp(int_eff_std$theta_hi    / EXP_sd) - 1) * 100
    int_eff_hlog$significant    <- ifelse(int_eff_hlog$theta_lo > 0 |
                                            int_eff_hlog$theta_hi < 0, "*", "")
  }
  
  # --- [NEW v7.0] Extract BOTH beta and theta draws for ΔERR computation ---
  beta_names_s  <- paste0("window_f", POST_WINDOWS)
  theta_names_s <- paste0("window_f", POST_WINDOWS, ":EXP_std")
  M_beta_s  <- get_coef_draws(fit_a, beta_names_s,  N_SAMPLES,
                              400 + which(c("High","Moderate","Low") == lvl))
  M_theta_s <- get_coef_draws(fit_a, theta_names_s, N_SAMPLES,
                              400 + which(c("High","Moderate","Low") == lvl))
  
  if (is.null(M_theta_s) || is.null(M_beta_s)) next
  
  Theta_std  <- colSums(M_theta_s)
  Theta_hlog <- Theta_std / EXP_sd
  beta_cum_s <- colSums(M_beta_s)
  
  # Compute %RR change draws
  pct_std  <- (exp(Theta_std)  - 1) * 100
  pct_hlog <- (exp(Theta_hlog) - 1) * 100
  
  # --- [NEW v7.0] Compute ΔERR draws ---
  derr_std  <- compute_cum_derr(beta_cum_s, Theta_std,  scale_by = 1)
  derr_hlog <- compute_cum_derr(beta_cum_s, Theta_hlog, scale_by = 1)
  
  adapt_strata[[lvl]] <- list(
    draws_std    = Theta_std,
    draws_hlog   = Theta_hlog,
    beta_cum_draws = beta_cum_s,
    cum_std_pct  = pct_std,
    cum_hlog_pct = pct_hlog,
    cum_std_derr  = derr_std,
    cum_hlog_derr = derr_hlog,
    int_eff_std  = int_eff_std,
    int_eff_hlog = int_eff_hlog
  )
  
  cat(sprintf("  HSR=%s [EXP per 1SD]:     %%RR = %+.2f [%+.2f, %+.2f] | ΔERR = %+.2f pp [%+.2f, %+.2f]\n",
              lvl, median(pct_std), quantile(pct_std, .025), quantile(pct_std, .975),
              median(derr_std), quantile(derr_std, .025), quantile(derr_std, .975)))
  cat(sprintf("  HSR=%s [EXP per 1 H_log]: %%RR = %+.2f [%+.2f, %+.2f] | ΔERR = %+.2f pp [%+.2f, %+.2f]\n",
              lvl, median(pct_hlog), quantile(pct_hlog, .025), quantile(pct_hlog, .975),
              median(derr_hlog), quantile(derr_hlog, .025), quantile(derr_hlog, .975)))
  
  # ---- 7b-2: HSR × window interaction within this stratum ----
  HSR_center_local <- mean(d$HSR_raw, na.rm = TRUE)
  d$HSR_c_local    <- d$HSR_raw - HSR_center_local
  F_HSR_LOCAL <- as.formula(paste("Y ~ window_f * HSR_c_local +", RE_TERMS))
  
  fit_h <- tryCatch({
    inla(F_HSR_LOCAL, family = "nbinomial", data = d, offset = log_offset,
         control.compute   = list(config = TRUE),
         control.predictor = list(link = 1),
         control.inla      = list(strategy = "simplified.laplace"),
         control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
         num.threads = "4:1", verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("    HSR=%s: within-stratum HSR model failed: %s\n",
                lvl, conditionMessage(e))); NULL
  })
  
  if (!is.null(fit_h)) {
    hsr_int_local <- extract_interaction_effects(fit_h, "HSR_c_local",
                                                 window_labels, REF_MAIN,
                                                 scale_by = HSR_STEP)
    
    # --- [NEW v7.0] Extract beta AND theta draws for within-stratum HSR ---
    beta_names_h  <- paste0("window_f", POST_WINDOWS)
    theta_names_h <- paste0("window_f", POST_WINDOWS, ":HSR_c_local")
    M_beta_h <- tryCatch(
      get_coef_draws(fit_h, beta_names_h, N_SAMPLES,
                     500 + which(c("High","Moderate","Low") == lvl)),
      error = function(e) NULL)
    M_theta_h <- tryCatch(
      get_coef_draws(fit_h, theta_names_h, N_SAMPLES,
                     500 + which(c("High","Moderate","Low") == lvl)),
      error = function(e) NULL)
    
    if (!is.null(M_theta_h) && !is.null(M_beta_h)) {
      hsr_cum_local_rr_draws   <- (exp(colSums(M_theta_h) * HSR_STEP) - 1) * 100
      beta_cum_h_local         <- colSums(M_beta_h)
      hsr_cum_local_derr_draws <- compute_cum_derr(beta_cum_h_local,
                                                   colSums(M_theta_h),
                                                   scale_by = HSR_STEP)
      
      adapt_strata_hsr[[lvl]] <- list(
        int_eff      = hsr_int_local,
        cum_draws    = hsr_cum_local_rr_draws,
        cum_derr_draws = hsr_cum_local_derr_draws,
        beta_cum_draws = beta_cum_h_local,
        HSR_center   = HSR_center_local,
        HSR_range    = range(d$HSR_raw, na.rm = TRUE)
      )
      
      cat(sprintf("    HSR=%s [HSR per +0.1, within]: %%RR = %+.2f [%+.2f, %+.2f] | ΔERR = %+.2f pp [%+.2f, %+.2f]\n",
                  lvl,
                  median(hsr_cum_local_rr_draws), quantile(hsr_cum_local_rr_draws, .025),
                  quantile(hsr_cum_local_rr_draws, .975),
                  median(hsr_cum_local_derr_draws), quantile(hsr_cum_local_derr_draws, .025),
                  quantile(hsr_cum_local_derr_draws, .975)))
    }
  }
}

# --- Comparison: experience adaptation between strata ---
adapt_comp <- data.frame()
if (!is.null(adapt_strata[["High"]]) && !is.null(adapt_strata[["Low"]])) {
  dHL <- adapt_strata[["High"]]$draws_std - adapt_strata[["Low"]]$draws_std
  pHL <- mean(adapt_strata[["High"]]$draws_std < adapt_strata[["Low"]]$draws_std)
  cat(sprintf("\n  P(Theta_High < Theta_Low) = %.3f\n", pHL))
  adapt_comp <- data.frame(comparison = "Theta_High vs Theta_Low",
                           prob_High_more_protective = pHL,
                           diff_median = median(dHL),
                           diff_lo = quantile(dHL, .025, names = FALSE),
                           diff_hi = quantile(dHL, .975, names = FALSE))
}

# ---------------------------------------------------------------
# Fig2.1-09: Posterior densities of the cumulative experience-
#   adaptation coefficient (Θ_cum, per +1 SD) stratified by HSR
#   level. More negative Θ_cum indicates stronger experience-based
#   adaptation (prior TC exposure more effectively reduces future
#   post-TC hospitalisation within that HSR stratum).
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  theta_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, Theta = adapt_strata[[l]]$draws_std)))
  theta_df$HSR_Level <- factor(theta_df$HSR_Level,
                               levels = c("High", "Moderate", "Low"))
  p_adapt <- ggplot(theta_df, aes(x = Theta, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High = "#2166AC", Moderate = "#F4A582",
                                 Low = "#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High = "#2166AC", Moderate = "#F4A582",
                                   Low = "#B2182B"), name = "HSR level") +
    labs(x = expression(Theta[cum]~"(per +1 SD experience, log scale)"),
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.85, 0.8))
  save_both(p_adapt, file.path(DIR_FIG,
                               sprintf("Fig2.1-09_TC_%s_Adaptation_Theta_SD", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_adapt, file.path(DIR_FIG,
                                           sprintf("Fig2.1-09_TC_%s_Adaptation_Theta_SD", Outcome_var)))
  write.xlsx(theta_df, file.path(DIR_TABLE_OUT,
                                 sprintf("Fig2.1-09_TC_%s_Adaptation_Theta_SD_Draws.xlsx", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09b: Posterior densities of the cumulative experience-
#   adaptation coefficient (Θ_cum, per +1 unit H_log) stratified
#   by HSR level.
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  theta_hlog_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, Theta = adapt_strata[[l]]$draws_hlog)))
  theta_hlog_df$HSR_Level <- factor(theta_hlog_df$HSR_Level,
                                    levels = c("High", "Moderate", "Low"))
  p_adapt_hlog <- ggplot(theta_hlog_df, aes(x = Theta, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High = "#2166AC", Moderate = "#F4A582",
                                 Low = "#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High = "#2166AC", Moderate = "#F4A582",
                                   Low = "#B2182B"), name = "HSR level") +
    labs(x = expression(Theta[cum]~"(per +1 H_log, log scale)"),
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.85, 0.8))
  save_both(p_adapt_hlog, file.path(DIR_FIG,
                                    sprintf("Fig2.1-09b_TC_%s_Adaptation_Theta_Hlog", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_adapt_hlog, file.path(DIR_FIG,
                                                sprintf("Fig2.1-09b_TC_%s_Adaptation_Theta_Hlog", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09c: Posterior densities of the cumulative % RR change
#   per +1 SD experience, stratified by HSR level.
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  rr_sd_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, pct = adapt_strata[[l]]$cum_std_pct)))
  rr_sd_df$HSR_Level <- factor(rr_sd_df$HSR_Level, levels = c("High","Moderate","Low"))
  p_adapt_rr_sd <- ggplot(rr_sd_df, aes(x = pct, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    labs(x = "Cumulative % RR change per +1 SD experience (%)",
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.15, 0.8))
  save_both(p_adapt_rr_sd, file.path(DIR_FIG,
                                     sprintf("Fig2.1-09c_TC_%s_Adaptation_RRpct_SD", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_adapt_rr_sd, file.path(DIR_FIG,
                                                 sprintf("Fig2.1-09c_TC_%s_Adaptation_RRpct_SD", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09d: Posterior densities of the cumulative % RR change
#   per +1 unit H_log, stratified by HSR level.
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  rr_hlog_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, pct = adapt_strata[[l]]$cum_hlog_pct)))
  rr_hlog_df$HSR_Level <- factor(rr_hlog_df$HSR_Level, levels = c("High","Moderate","Low"))
  p_adapt_rr_hlog <- ggplot(rr_hlog_df, aes(x = pct, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    labs(x = "Cumulative % RR change per +1 H_log (%)",
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.15, 0.8))
  save_both(p_adapt_rr_hlog, file.path(DIR_FIG,
                                       sprintf("Fig2.1-09d_TC_%s_Adaptation_RRpct_Hlog", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_adapt_rr_hlog, file.path(DIR_FIG,
                                                   sprintf("Fig2.1-09d_TC_%s_Adaptation_RRpct_Hlog", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09e: Posterior densities of the cumulative ΔERR (pp) per
#   +1 SD experience, stratified by HSR level. ΔERR = RR_base ×
#   (exp(Θ_cum) − 1) × 100. More negative values indicate greater
#   absolute reduction in post-TC hospitalisation burden.
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  derr_sd_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, derr = adapt_strata[[l]]$cum_std_derr)))
  derr_sd_df$HSR_Level <- factor(derr_sd_df$HSR_Level, levels = c("High","Moderate","Low"))
  p_adapt_derr_sd <- ggplot(derr_sd_df, aes(x = derr, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    labs(x = "Cumulative ΔERR per +1 SD experience (pp)",
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.15, 0.8))
  save_both(p_adapt_derr_sd, file.path(DIR_FIG,
                                       sprintf("Fig2.1-09e_TC_%s_Adaptation_DERR_SD", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_adapt_derr_sd, file.path(DIR_FIG,
                                                   sprintf("Fig2.1-09e_TC_%s_Adaptation_DERR_SD", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09f: Posterior densities of the cumulative ΔERR (pp) per
#   +1 unit H_log, stratified by HSR level.
# ---------------------------------------------------------------
if (length(adapt_strata) >= 2) {
  derr_hlog_df <- bind_rows(lapply(names(adapt_strata), function(l)
    data.frame(HSR_Level = l, derr = adapt_strata[[l]]$cum_hlog_derr)))
  derr_hlog_df$HSR_Level <- factor(derr_hlog_df$HSR_Level, levels = c("High","Moderate","Low"))
  p_adapt_derr_hlog <- ggplot(derr_hlog_df, aes(x = derr, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    labs(x = "ΔCumulative ERR per +1 unit \nlog-scaled TC experience index (%)",
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.82, 0.7))
  save_both(p_adapt_derr_hlog, file.path(DIR_FIG,
                                         sprintf("Fig2.1-09f_TC_%s_Adaptation_DERR_Hlog", Outcome_var)), w = 3.5, h = 2.5)
  save_zoom_transparent(p_adapt_derr_hlog, file.path(DIR_FIG,
                                                     sprintf("Fig2.1-09f_TC_%s_Adaptation_DERR_Hlog", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09g: Posterior densities of the within-stratum HSR
#   modification (% RR change per +0.1 HSR), stratified by HSR
#   level. This shows residual HSR variation effects within each
#   stratum.
# ---------------------------------------------------------------
if (length(adapt_strata_hsr) >= 2) {
  hsr_within_df <- bind_rows(lapply(names(adapt_strata_hsr), function(l)
    data.frame(HSR_Level = l, pct = adapt_strata_hsr[[l]]$cum_draws)))
  hsr_within_df$HSR_Level <- factor(hsr_within_df$HSR_Level, levels = c("High","Moderate","Low"))
  p_hsr_within <- ggplot(hsr_within_df, aes(x = pct, fill = HSR_Level, colour = HSR_Level)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
    labs(x = "Within-stratum: % RR change per +0.1 HSR (%)",
         y = "Posterior density") + theme_pub() +
    theme(legend.position = c(0.15, 0.8))
  save_both(p_hsr_within, file.path(DIR_FIG,
                                    sprintf("Fig2.1-09g_TC_%s_WithinHSR_RRpct", Outcome_var)), w = 5, h = 4)
  save_zoom_transparent(p_hsr_within, file.path(DIR_FIG,
                                                sprintf("Fig2.1-09g_TC_%s_WithinHSR_RRpct", Outcome_var)))
}

# ---------------------------------------------------------------
# Fig2.1-09h: Posterior densities of the within-stratum HSR ΔERR
#   (pp per +0.1 HSR), stratified by HSR level.
# ---------------------------------------------------------------
if (length(adapt_strata_hsr) >= 2) {
  hsr_derr_within_df <- bind_rows(lapply(names(adapt_strata_hsr), function(l) {
    if (!is.null(adapt_strata_hsr[[l]]$cum_derr_draws))
      data.frame(HSR_Level = l, derr = adapt_strata_hsr[[l]]$cum_derr_draws)
    else NULL
  }))
  if (nrow(hsr_derr_within_df) > 0) {
    hsr_derr_within_df$HSR_Level <- factor(hsr_derr_within_df$HSR_Level,
                                           levels = c("High","Moderate","Low"))
    p_hsr_derr_within <- ggplot(hsr_derr_within_df, aes(x = derr, fill = HSR_Level, colour = HSR_Level)) +
      geom_density(alpha = 0.3, linewidth = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_fill_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
      scale_colour_manual(values = c(High="#2166AC", Moderate="#F4A582", Low="#B2182B"), name = "HSR level") +
      labs(x = "Within-stratum: ΔERR per +0.1 HSR (pp)",
           y = "Posterior density") + theme_pub() +
      theme(legend.position = c(0.15, 0.8))
    save_both(p_hsr_derr_within, file.path(DIR_FIG,
                                           sprintf("Fig2.1-09h_TC_%s_WithinHSR_DERR", Outcome_var)), w = 5, h = 4)
    save_zoom_transparent(p_hsr_derr_within, file.path(DIR_FIG,
                                                       sprintf("Fig2.1-09h_TC_%s_WithinHSR_DERR", Outcome_var)))
  }
}

# ---------------------------------------------------------------
# Fig2.1-10: Posterior density of the difference in Θ_cum between
#   High-HSR and Low-HSR strata (Θ_High − Θ_Low). Negative values
#   indicate that cities with high HSR exhibit stronger experience-
#   based adaptation (greater protective effect per SD experience).
# ---------------------------------------------------------------
if (nrow(adapt_comp) > 0) {
  d_hl <- adapt_strata[["High"]]$draws_std - adapt_strata[["Low"]]$draws_std
  anno_hl <- sprintf("Median diff = %.3f (%.3f, %.3f)\n%s",
                     median(d_hl), quantile(d_hl, .025, names = FALSE),
                     quantile(d_hl, .975, names = FALSE), ptxt(d_hl, "lt0"))
  p_diff <- ggplot(data.frame(d = d_hl), aes(x = d)) +
    geom_density(fill = "#6C3483", alpha = 0.35, colour = "#6C3483") +
    geom_vline(xintercept = median(d_hl), colour = "#6C3483", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    annotate("text", x = Inf, y = Inf, label = anno_hl,
             hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
    labs(x = expression(Theta[High] - Theta[Low]~"(negative = HSR facilitates adaptation)"),
         y = "Posterior density") + theme_pub()
  save_both(p_diff, file.path(DIR_FIG,
                              sprintf("Fig2.1-10_TC_%s_Adaptation_ThetaDiff_Posterior", Outcome_var)),
            w = 5, h = 4)
  save_zoom_transparent(p_diff, file.path(DIR_FIG,
                                          sprintf("Fig2.1-10_TC_%s_Adaptation_ThetaDiff_Posterior", Outcome_var)))
}


# ================================================================
#  PART 7c: PERIOD SPLIT 2016-2020 vs 2021-2024
# ================================================================

cat("\n========== PART 7c: PERIOD SPLIT ==========\n")

fit_per <- inla(F_PERIOD, family = "nbinomial", data = df, offset = log_offset,
                control.compute   = list(config = TRUE),
                control.predictor = list(link = 1),
                control.inla      = list(strategy = "simplified.laplace"),
                control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                num.threads = "4:1", verbose = FALSE)
per_int <- extract_interaction_effects(fit_per, "period_fP2021_2024",
                                       window_labels, REF_MAIN)
print(per_int)
M_per <- get_coef_draws(fit_per,
                        paste0("window_f", POST_WINDOWS, ":period_fP2021_2024"),
                        N_SAMPLES, 341)
if (!is.null(M_per)) {
  per_draws <- (exp(colSums(M_per)) - 1) * 100
  cat(sprintf("  Cumulative post-TC RR, 2021-24 vs 2016-20: %+.2f%% [%+.2f, %+.2f], P(decline)=%.3f\n",
              median(per_draws), quantile(per_draws, .025),
              quantile(per_draws, .975), mean(per_draws < 0)))
}


# ================================================================
#  PART 8: SAVE ALL RESULTS [v7.0]
# ================================================================

cat("\n========== PART 8: SAVE ==========\n")
wb <- createWorkbook()
add_sheet <- function(nm, x) if (!is.null(x) && NROW(x) > 0) {
  addWorksheet(wb, nm); writeData(wb, nm, x) }
add_sheet("Main_Window_ERR",   res_main)
add_sheet("Main_Posterior_Stats", marg_stats)
add_sheet("Sens1_Window_ERR",  res_sens1)
add_sheet("Sens2_Window_ERR",  res_sens2)
add_sheet("Sens3_Window_ERR",  res_sens3)
add_sheet("WAIC",              waic_df)
add_sheet("Cumulative_ERR",    bind_rows(
  mutate(cum_main$results, model = "Main"),  mutate(cum_sens1$results, model = "Sens1"),
  mutate(cum_sens2$results, model = "Sens2"), mutate(cum_sens3$results, model = "Sens3")))
add_sheet("Model_Comparison",  comp_all)
add_sheet("HSR_Interaction",   hsr_int)
add_sheet("HSR_CumRRpct",     as.data.frame(hsr_cum_mod))
add_sheet("HSR_CumDERR",      as.data.frame(hsr_cum_derr_summ))
add_sheet("HSR_Window_DERR",   hsr_win_derr)
add_sheet("HSR_ERR_Curve",     hsr_curve)
add_sheet("HSR_Strata_Comparison", hsr_comp)
add_sheet("EXP_Interaction",   exp_int)
add_sheet("EXP_Interaction_Raw", exp_int_raw)
add_sheet("EXP_CumRRpct",     as.data.frame(exp_cum_raw))
add_sheet("EXP_CumDERR",      as.data.frame(exp_cum_derr_summ))
add_sheet("EXP_Window_DERR",   exp_win_derr)
add_sheet("EXP_ERR_Curve",     exp_curve)
add_sheet("Adaptation_HSRxEXP", adapt_comp)
add_sheet("Period_Interaction", per_int)

# Within-stratum summaries
for (lv in names(adapt_strata)) {
  s <- adapt_strata[[lv]]
  row_s <- data.frame(
    HSR_Level = lv,
    EXP_SD_RRpct_median = median(s$cum_std_pct),
    EXP_SD_RRpct_lo = quantile(s$cum_std_pct, .025, names = FALSE),
    EXP_SD_RRpct_hi = quantile(s$cum_std_pct, .975, names = FALSE),
    EXP_Hlog_RRpct_median = median(s$cum_hlog_pct),
    EXP_Hlog_RRpct_lo = quantile(s$cum_hlog_pct, .025, names = FALSE),
    EXP_Hlog_RRpct_hi = quantile(s$cum_hlog_pct, .975, names = FALSE),
    EXP_SD_DERR_median = median(s$cum_std_derr),
    EXP_SD_DERR_lo = quantile(s$cum_std_derr, .025, names = FALSE),
    EXP_SD_DERR_hi = quantile(s$cum_std_derr, .975, names = FALSE),
    EXP_Hlog_DERR_median = median(s$cum_hlog_derr),
    EXP_Hlog_DERR_lo = quantile(s$cum_hlog_derr, .025, names = FALSE),
    EXP_Hlog_DERR_hi = quantile(s$cum_hlog_derr, .975, names = FALSE))
  if (lv == names(adapt_strata)[1]) strata_exp_summ <- row_s
  else strata_exp_summ <- rbind(strata_exp_summ, row_s)
}
add_sheet("Strata_EXP_Summary", strata_exp_summ)

for (lv in names(adapt_strata_hsr)) {
  s <- adapt_strata_hsr[[lv]]
  row_h <- data.frame(
    HSR_Level = lv,
    HSR_RRpct_median = median(s$cum_draws),
    HSR_RRpct_lo = quantile(s$cum_draws, .025, names = FALSE),
    HSR_RRpct_hi = quantile(s$cum_draws, .975, names = FALSE),
    HSR_DERR_median = if (!is.null(s$cum_derr_draws)) median(s$cum_derr_draws) else NA,
    HSR_DERR_lo = if (!is.null(s$cum_derr_draws)) quantile(s$cum_derr_draws, .025, names = FALSE) else NA,
    HSR_DERR_hi = if (!is.null(s$cum_derr_draws)) quantile(s$cum_derr_draws, .975, names = FALSE) else NA)
  if (lv == names(adapt_strata_hsr)[1]) strata_hsr_summ <- row_h
  else strata_hsr_summ <- rbind(strata_hsr_summ, row_h)
}
add_sheet("Strata_HSR_Summary", strata_hsr_summ)

saveWorkbook(wb, file.path(DIR_TABLE_OUT,
                           sprintf("Tab2.1-00_TC_%s_Step2.1_AllResults.xlsx", Outcome_var)),
             overwrite = TRUE)

# --- Save posterior RDS (comprehensive, for Step 2.1S) ---
hsr_strata_cum <- if (length(hsr_strata) > 0) {
  lapply(hsr_strata, function(s) s$cum_draws)
} else { list() }

saveRDS(list(outcome        = Outcome_var,
             analysis_type  = Analysis_type,
             pop_col        = pop_col_name,
             window_results = res_main,
             cum_samples    = cum_main$samples,
             hsr_strata_cum_draws = hsr_strata_cum,
             hsr_strata_beta_cum  = lapply(hsr_strata, function(s) s$beta_cum_draws),
             hsr_curve      = hsr_curve,
             hsr_cum_mod    = hsr_cum_mod,
             hsr_cum_derr_summ = hsr_cum_derr_summ,
             adapt_strata_draws = adapt_strata,
             adapt_strata_hsr   = adapt_strata_hsr,
             adapt_comp     = adapt_comp,
             HSR_center     = HSR_center,
             EXP_center     = EXPC_center,
             EXP_sd         = EXP_sd,
             overall_hsr_cum_rr_draws  = hsr_cum_rr_draws,
             overall_hsr_cum_derr_draws = hsr_cum_derr_draws,
             overall_exp_cum_rr_draws  = exp_cum_raw_draws,
             overall_exp_cum_derr_draws = exp_cum_derr_draws,
             created        = Sys.time()),
        file.path(DIR_OUTPUT, sprintf("TC_%s_Step2.1_Posterior.rds", Outcome_var)))
cat("  Step 2.1 posterior data saved.\n")


# ================================================================
#  PART 9: SUMMARY TABLE [v7.0 — with Scale & Interpretation]
# ================================================================

cat("\n========== PART 9: SUMMARY TABLE ==========\n")

SUMMARY_COLS <- c("Analysis", "Fig", "Scale", "-3", "-2", "-1", "Active",
                  as.character(1:6), "Cumulative", "P_Cumulative", "Interpretation")

fmt_est_sum <- function(est, lo, hi, digits = 2)
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
          est, lo, hi)

WIN_COL_MAP <- setNames(c("-3", "-2", "-1", "Active", as.character(1:6)),
                        window_labels)

make_summary_row <- function(analysis, fig, scale = "", df_win = NULL,
                             win_col = "window", est_col = "ERR_pct",
                             lo_col = "lower", hi_col = "upper",
                             ref_label = NULL, cum_triple = NULL,
                             p_text = NULL, interpretation = "",
                             digits = 2) {
  row <- as.data.frame(as.list(setNames(rep(NA_character_, length(SUMMARY_COLS)),
                                        SUMMARY_COLS)),
                       check.names = FALSE, stringsAsFactors = FALSE)
  row[["Analysis"]]       <- analysis
  row[["Fig"]]            <- fig
  row[["Scale"]]          <- scale
  row[["Interpretation"]] <- interpretation
  if (!is.null(df_win)) {
    need <- c(win_col, est_col, lo_col, hi_col)
    miss <- setdiff(need, names(df_win))
    if (length(miss) > 0) {
      cat(sprintf("  [Summary] '%s' missing cols: %s\n", analysis, paste(miss, collapse = ", ")))
      return(row)
    }
    for (k in seq_len(nrow(df_win))) {
      wl <- as.character(df_win[[win_col]][k])
      if (!wl %in% names(WIN_COL_MAP)) next
      if (!is.null(ref_label) && wl == ref_label) {
        row[[WIN_COL_MAP[[wl]]]] <- "0 (ref)"
      } else {
        row[[WIN_COL_MAP[[wl]]]] <- fmt_est_sum(df_win[[est_col]][k],
                                                df_win[[lo_col]][k],
                                                df_win[[hi_col]][k], digits)
      }
    }
    if (!is.null(ref_label) && is.na(row[[WIN_COL_MAP[[ref_label]]]]))
      row[[WIN_COL_MAP[[ref_label]]]] <- "0 (ref)"
  }
  if (!is.null(cum_triple) && length(cum_triple) == 3 && all(is.finite(cum_triple)))
    row[["Cumulative"]] <- fmt_est_sum(cum_triple[1], cum_triple[2],
                                       cum_triple[3], digits)
  if (!is.null(p_text) && !is.na(p_text))
    row[["P_Cumulative"]] <- p_text
  row
}

cum_from_results <- function(cum_obj, period = "Active_to_PostFull") {
  r <- cum_obj$results[cum_obj$results$period == period, ]
  if (nrow(r) != 1) return(NULL)
  as.numeric(r[1, c("median", "lower_2.5", "upper_97.5")])
}
cum_from_err_draws <- function(x) {
  if (is.null(x)) return(NULL)
  c(median(x), quantile(x, .025, names = FALSE), quantile(x, .975, names = FALSE))
}

rows_21 <- list()
add_row <- function(...) {
  r <- tryCatch(make_summary_row(...),
                error = function(e) { cat("  [Summary ERROR]:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) rows_21[[length(rows_21) + 1]] <<- r
}

## ---------- Main model ----------
main_cum_draws <- cum_main$samples[["Active_to_PostFull"]]
add_row(analysis = "Main", fig = "Fig2.1-01/02/03",
        scale = "ERR (%)",
        df_win = res_main, ref_label = REF_MAIN,
        cum_triple = cum_from_results(cum_main),
        p_text = ptxt(main_cum_draws, "gt0"),
        interpretation = "Cumulative post-TC ERR (Active to +6mo). Positive = elevated risk.")

## ---------- Sensitivity ----------
add_row(analysis = "Sen_AltRef_-2mo", fig = "Tab",
        scale = "ERR (%)",
        df_win = res_sens1, ref_label = REF_SENS1,
        cum_triple = cum_from_results(cum_sens1),
        p_text = ptxt(cum_sens1$samples[["Active_to_PostFull"]], "gt0"),
        interpretation = "Robustness check: alternative reference window (-2mo).")
add_row(analysis = "Sen_Gaussian", fig = "Tab",
        scale = "ERR (%)",
        df_win = res_sens2, ref_label = REF_MAIN,
        cum_triple = cum_from_results(cum_sens2),
        p_text = ptxt(cum_sens2$samples[["Active_to_PostFull"]], "gt0"),
        interpretation = "Robustness check: Gaussian approximation strategy.")
add_row(analysis = "Sen_PandemicExcl", fig = "Tab",
        scale = "ERR (%)",
        df_win = res_sens3, ref_label = REF_MAIN,
        cum_triple = cum_from_results(cum_sens3),
        p_text = ptxt(cum_sens3$samples[["Active_to_PostFull"]], "gt0"),
        interpretation = "Robustness check: COVID-19 period (2020-2022) excluded.")

## ---------- HSR strata (cumulative ERR) ----------
for (lv in c("High", "Moderate", "Low")) {
  s <- hsr_strata[[lv]]
  if (is.null(s)) next
  add_row(analysis = paste0(lv, "-HSR_CumERR"), fig = "Fig2.1-07",
          scale = "ERR (%)",
          df_win = s$window_effects, ref_label = REF_MAIN,
          cum_triple = cum_from_err_draws(s$cum_draws),
          p_text = ptxt(s$cum_draws, "gt0"),
          interpretation = sprintf("Cumulative ERR within %s-HSR stratum.", lv))
}

## ---------- Overall HSR interaction (% RR) ----------
add_row(analysis = "Overall_HSR_per0.1_RRpct", fig = "Fig2.1-04/05",
        scale = "% RR change",
        df_win = hsr_int, est_col = "pct_change",
        lo_col = "pct_change_lo", hi_col = "pct_change_hi",
        ref_label = REF_MAIN,
        cum_triple = as.numeric(hsr_cum_mod[1, c("median", "lower_2.5", "upper_97.5")]),
        p_text = ptxt(hsr_cum_rr_draws, "lt0"),
        interpretation = "Relative change in cumulative RR per +0.1 HSR. Negative = HSR protective.")

## ---------- Overall HSR interaction (ΔERR) ----------
add_row(analysis = "Overall_HSR_per0.1_DERR", fig = "Fig2.1-05b/05c",
        scale = "ΔERR (pp)",
        cum_triple = as.numeric(hsr_cum_derr_summ[1, c("median", "lower_2.5", "upper_97.5")]),
        p_text = ptxt(hsr_cum_derr_draws, "lt0"),
        interpretation = "Absolute cumulative ERR change (pp) per +0.1 HSR. Negative = HSR reduces absolute excess risk.")

## ---------- Overall EXP interaction (% RR, per 1 SD) ----------
add_row(analysis = "Overall_EXP_per1SD_RRpct", fig = "Tab",
        scale = "% RR change",
        df_win = exp_int, est_col = "pct_change",
        lo_col = "pct_change_lo", hi_col = "pct_change_hi",
        ref_label = REF_MAIN,
        cum_triple = as.numeric(exp_cum_mod[1, c("median", "lower_2.5", "upper_97.5")]),
        p_text = ptxt(exp_cum_draws, "lt0"),
        interpretation = "Relative change in cumulative RR per +1 SD experience. Negative = adaptation.")

## ---------- Overall EXP interaction (% RR, per 1 H_log) ----------
add_row(analysis = "Overall_EXP_per1Hlog_RRpct", fig = "Fig2.1-08a/08c",
        scale = "% RR change",
        df_win = exp_int_raw, est_col = "pct_change",
        lo_col = "pct_change_lo", hi_col = "pct_change_hi",
        ref_label = REF_MAIN,
        cum_triple = as.numeric(exp_cum_raw[1, c("median", "lower_2.5", "upper_97.5")]),
        p_text = ptxt(exp_cum_raw_draws, "lt0"),
        interpretation = "Relative change in cumulative RR per +1 unit log(H+1). Negative = adaptation.")

## ---------- Overall EXP interaction (ΔERR, per 1 H_log) ----------
add_row(analysis = "Overall_EXP_per1Hlog_DERR", fig = "Fig2.1-08b/08d",
        scale = "ΔERR (pp)",
        cum_triple = as.numeric(exp_cum_derr_summ[1, c("median", "lower_2.5", "upper_97.5")]),
        p_text = ptxt(exp_cum_derr_draws, "lt0"),
        interpretation = "Absolute cumulative ERR change (pp) per +1 H_log. Negative = experience reduces excess risk.")

## ---------- Within-stratum Experience (all metrics, per stratum) ----------
for (lv in c("High", "Moderate", "Low")) {
  s <- adapt_strata[[lv]]
  if (is.null(s)) next
  
  add_row(analysis = paste0(lv, "-HSR_EXP_per1SD_RRpct"), fig = "Fig2.1-09c",
          scale = "% RR change",
          df_win = s$int_eff_std, est_col = "pct_change",
          lo_col = "pct_change_lo", hi_col = "pct_change_hi", ref_label = REF_MAIN,
          cum_triple = cum_from_err_draws(s$cum_std_pct),
          p_text = ptxt(s$draws_std, "lt0"),
          interpretation = sprintf("Within %s-HSR: %%RR change per +1 SD experience.", lv))
  
  add_row(analysis = paste0(lv, "-HSR_EXP_per1SD_DERR"), fig = "Fig2.1-09e",
          scale = "ΔERR (pp)",
          cum_triple = cum_from_err_draws(s$cum_std_derr),
          p_text = ptxt(s$cum_std_derr, "lt0"),
          interpretation = sprintf("Within %s-HSR: absolute ΔERR (pp) per +1 SD experience.", lv))
  
  add_row(analysis = paste0(lv, "-HSR_EXP_per1Hlog_RRpct"), fig = "Fig2.1-09d",
          scale = "% RR change",
          df_win = s$int_eff_hlog, est_col = "pct_change",
          lo_col = "pct_change_lo", hi_col = "pct_change_hi", ref_label = REF_MAIN,
          cum_triple = cum_from_err_draws(s$cum_hlog_pct),
          p_text = ptxt(s$draws_hlog, "lt0"),
          interpretation = sprintf("Within %s-HSR: %%RR change per +1 H_log.", lv))
  
  add_row(analysis = paste0(lv, "-HSR_EXP_per1Hlog_DERR"), fig = "Fig2.1-09f",
          scale = "ΔERR (pp)",
          cum_triple = cum_from_err_draws(s$cum_hlog_derr),
          p_text = ptxt(s$cum_hlog_derr, "lt0"),
          interpretation = sprintf("Within %s-HSR: absolute ΔERR (pp) per +1 H_log.", lv))
}

## ---------- Within-stratum HSR (per stratum) ----------
for (lv in c("High", "Moderate", "Low")) {
  s <- adapt_strata_hsr[[lv]]
  if (is.null(s)) next
  
  add_row(analysis = paste0(lv, "-HSR_HSR_per0.1_within_RRpct"), fig = "Fig2.1-09g",
          scale = "% RR change",
          df_win = s$int_eff, est_col = "pct_change",
          lo_col = "pct_change_lo", hi_col = "pct_change_hi", ref_label = REF_MAIN,
          cum_triple = cum_from_err_draws(s$cum_draws),
          p_text = ptxt(s$cum_draws, "lt0"),
          interpretation = sprintf("Within %s-HSR: %%RR change per +0.1 HSR variation.", lv))
  
  if (!is.null(s$cum_derr_draws)) {
    add_row(analysis = paste0(lv, "-HSR_HSR_per0.1_within_DERR"), fig = "Fig2.1-09h",
            scale = "ΔERR (pp)",
            cum_triple = cum_from_err_draws(s$cum_derr_draws),
            p_text = ptxt(s$cum_derr_draws, "lt0"),
            interpretation = sprintf("Within %s-HSR: absolute ΔERR (pp) per +0.1 HSR variation.", lv))
  }
}

## ---------- Period ----------
if (exists("per_int") && !is.null(per_int))
  add_row(analysis = "Period_2021-24_vs_2016-20", fig = "Tab",
          scale = "% RR change",
          df_win = per_int, est_col = "pct_change",
          lo_col = "pct_change_lo", hi_col = "pct_change_hi", ref_label = REF_MAIN,
          cum_triple = if (exists("per_draws")) cum_from_err_draws(per_draws) else NULL,
          p_text = if (exists("per_draws")) ptxt(per_draws, "lt0") else NA_character_,
          interpretation = "Temporal change: cumulative %RR for 2021-24 relative to 2016-20.")

## ---------- Compile and save ----------
summary_21 <- do.call(rbind, rows_21)
print(summary_21[, c("Analysis", "Scale", "Cumulative", "P_Cumulative", "Interpretation")],
      right = FALSE)
write.xlsx(summary_21,
           file.path(DIR_TABLE_OUT,
                     sprintf("Tab2.1-Summary_TC_%s_AllAnalyses.xlsx", Outcome_var)),
           overwrite = TRUE)


# ================================================================
#  PART 9b: POSTERIOR PROBABILITY TABLE [v7.0]
# ================================================================

cat("\n========== PART 9b: POSTERIOR PROBABILITY TABLE ==========\n")

prow <- function(Analysis, Quantity, Statement, P, Interpretation)
  data.frame(Analysis = Analysis, Quantity = Quantity, Statement = Statement,
             Probability = round(P, 4), Interpretation = Interpretation,
             stringsAsFactors = FALSE)

pv <- list()

for (k in seq_len(nrow(marg_stats)))
  pv[[length(pv) + 1]] <- prow("Main",
                               paste0("ERR, window ", marg_stats$window[k]), "P(ERR > 0)",
                               marg_stats$prob_ERR_pos[k],
                               "Posterior probability that hospitalisation is elevated in this window")

pv[[length(pv) + 1]] <- prow("Main", "Cumulative ERR (Active to +6mo)",
                             "P(ERR > 0)", mean(main_cum_draws > 0),
                             "Posterior probability of a net post-TC increase in hospitalisation")

for (nm in c("Sens1", "Sens2", "Sens3")) {
  s <- get(paste0("cum_", tolower(nm)))$samples[["Active_to_PostFull"]]
  if (!is.null(s))
    pv[[length(pv) + 1]] <- prow(nm, "Cumulative ERR (Active to +6mo)",
                                 "P(ERR > 0)", mean(s > 0), "Robustness check")
}

for (lv in names(hsr_strata))
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR"), "Cumulative ERR",
                               "P(ERR > 0)", mean(hsr_strata[[lv]]$cum_draws > 0),
                               "Within-stratum cumulative effect")
if (exists("hsr_comp") && nrow(hsr_comp) > 0)
  for (k in seq_len(nrow(hsr_comp)))
    pv[[length(pv) + 1]] <- prow("HSR strata", hsr_comp$comparison[k],
                                 "P(first < second)", hsr_comp$prob_first_lower[k],
                                 "Higher HSR stratum suffers smaller ERR")

# Overall HSR
pv[[length(pv) + 1]] <- prow("Overall_HSR", "Cum %RR per +0.1",
                             "P(< 0)", mean(hsr_cum_rr_draws < 0),
                             "Higher HSR attenuates cumulative RR (relative)")
pv[[length(pv) + 1]] <- prow("Overall_HSR", "Cum ΔERR per +0.1",
                             "P(< 0)", mean(hsr_cum_derr_draws < 0),
                             "Higher HSR attenuates cumulative ERR (absolute)")

# Overall EXP
pv[[length(pv) + 1]] <- prow("Overall_EXP_SD", "Cum %RR per +1SD",
                             "P(< 0)", mean(exp_cum_draws < 0),
                             "Experience-based adaptation (standardised)")
pv[[length(pv) + 1]] <- prow("Overall_EXP_Hlog", "Cum %RR per +1 H_log",
                             "P(< 0)", mean(exp_cum_raw_draws < 0),
                             "Experience-based adaptation (raw log scale)")
pv[[length(pv) + 1]] <- prow("Overall_EXP_Hlog", "Cum ΔERR per +1 H_log",
                             "P(< 0)", mean(exp_cum_derr_draws < 0),
                             "Experience reduces absolute cumulative ERR")

# Within-stratum EXP adaptation
for (lv in names(adapt_strata)) {
  s <- adapt_strata[[lv]]
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_EXP_SD"),
                               "Theta_cum (per 1 SD)", "P(< 0)",
                               mean(s$draws_std < 0),
                               sprintf("Experience adaptation within %s-HSR (per SD)", lv))
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_EXP_Hlog"),
                               "Theta_cum (per 1 H_log)", "P(< 0)",
                               mean(s$draws_hlog < 0),
                               sprintf("Experience adaptation within %s-HSR (per H_log)", lv))
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_EXP_SD_DERR"),
                               "ΔERR per +1 SD", "P(< 0)",
                               mean(s$cum_std_derr < 0),
                               sprintf("Absolute ERR reduction per +1 SD within %s-HSR", lv))
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_EXP_Hlog_DERR"),
                               "ΔERR per +1 H_log", "P(< 0)",
                               mean(s$cum_hlog_derr < 0),
                               sprintf("Absolute ERR reduction per +1 H_log within %s-HSR", lv))
}

# Within-stratum HSR
for (lv in names(adapt_strata_hsr)) {
  s <- adapt_strata_hsr[[lv]]
  pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_within"),
                               "%RR per +0.1 HSR (within)", "P(< 0)",
                               mean(s$cum_draws < 0),
                               sprintf("Within %s: +0.1 HSR attenuates RR", lv))
  if (!is.null(s$cum_derr_draws))
    pv[[length(pv) + 1]] <- prow(paste0(lv, "-HSR_within_DERR"),
                                 "ΔERR per +0.1 HSR (within)", "P(< 0)",
                                 mean(s$cum_derr_draws < 0),
                                 sprintf("Within %s: +0.1 HSR reduces absolute ERR", lv))
}

# Theta comparison across strata
if (exists("adapt_comp") && nrow(adapt_comp) > 0)
  pv[[length(pv) + 1]] <- prow("Adaptation_HighvsLow", "Theta_High vs Theta_Low",
                               "P(Theta_High < Theta_Low)", adapt_comp$prob_High_more_protective[1],
                               "High HSR facilitates stronger experience-based adaptation")

# Period
if (exists("per_draws"))
  pv[[length(pv) + 1]] <- prow("Period_2021-24",
                               "Cumulative RR ratio", "P(< 0)", mean(per_draws < 0),
                               "Post-TC risk declined over time")

pvals_21 <- bind_rows(pv)
print(pvals_21, right = FALSE)
write.xlsx(pvals_21,
           file.path(DIR_TABLE_OUT,
                     sprintf("Tab2.1-P_TC_%s_PosteriorProbabilities.xlsx", Outcome_var)),
           overwrite = TRUE)

cat(sprintf("  Summary: %d rows | Posterior probability table: %d rows\n",
            nrow(summary_21), nrow(pvals_21)))
cat(sprintf("\n  Step 2.1 (v7.0) FINISHED for Outcome=%s, Analysis=%s\n",
            Outcome_var, Analysis_type))











