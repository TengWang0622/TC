#############################################################################################
#
#                          TROPICAL CYCLONES AND HEALTH
#                Pairwise Subgroup Posterior Comparison (City Level, v2.3)
#
#  Developed by Teng Wang  @ The University of Hong Kong
#              Hanxu Shi  @ Peking University
#  Version: 2.3
#  Last updated: 2026-08-18
#
#############################################################################################


# ================================================================
#  PART 0: CONFIGURATION
# ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(stringr)

setwd("C:/Project/Tropical cyclone")
DIR_BASE <- "Result TC-health CityLevel"

# ===================== USER INPUT: SPECIFY TWO GROUPS =====================
#
#  REF  = Reference group  (Group A)
#  COMP = Comparison group (Group B)
#
#  HSR:       "All" | "High" | "Moderate" | "Low"
#             "All" = use the overall cumulative ERR (all cities pooled)
#             "High"/"Moderate"/"Low" = use the HSR-stratified cumulative ERR
#
#  ANALYSIS:  "Disease" | "Age_Stratified" | "Sex_Stratified"
#
#  OUTCOME:   any outcome variable from Step 2.1
#             Disease: "All_cause","Infectious","Circulatory","Respiratory",
#                      "Injuries","Mental","Nervous","Genitourinary",
#                      "Neoplasms","Endocrine","Nutritional","Metabolic"
#             Age:     "Age_0_19","Age_20_44","Age_45_64","Age_65plus"
#             Sex:     "Male","Female"
#
# --- Reference group ---
REF_HSR       <- "All"
REF_ANALYSIS  <- "Sex_Stratified"
REF_OUTCOME   <- "Female"

# --- Comparison group ---
COMP_HSR      <- "All"
COMP_ANALYSIS <- "Sex_Stratified"
COMP_OUTCOME  <- "Male"

# =========================================================================

# ---- Parameters (must match Step 2.1) ----
PRE_DAYS       <- 90
MAX_LAG_DAYS   <- 180
WINDOW_SIZE    <- 30
N_PRE_WINDOWS  <- PRE_DAYS %/% WINDOW_SIZE          # 3
N_POST_WINDOWS <- MAX_LAG_DAYS %/% WINDOW_SIZE      # 6
HSR_STEP       <- 0.1
N_SAMPLES      <- 1000

window_labels  <- c(paste0("-", N_PRE_WINDOWS:1, "mo"), "Active",
                    paste0("+", 1:N_POST_WINDOWS, "mo"))
window_display <- c(paste0("-", N_PRE_WINDOWS:1), "Active",
                    paste0("+", 1:N_POST_WINDOWS))
POST_WINDOWS   <- c("Active", paste0("+", 1:N_POST_WINDOWS, "mo"))

# ---- Derived labels ----
REF_LABEL  <- paste0(REF_HSR, "_", REF_ANALYSIS, "_", REF_OUTCOME)
COMP_LABEL <- paste0(COMP_HSR, "_", COMP_ANALYSIS, "_", COMP_OUTCOME)

make_display_name <- function(hsr, analysis, outcome) {
  out_clean <- gsub("_", " ", outcome)
  out_clean <- gsub("Age ", "", out_clean)
  if (hsr == "All") return(out_clean)
  paste0(hsr, "-HSR ", out_clean)
}
REF_DISPLAY  <- make_display_name(REF_HSR, REF_ANALYSIS, REF_OUTCOME)
COMP_DISPLAY <- make_display_name(COMP_HSR, COMP_ANALYSIS, COMP_OUTCOME)

# ---- Output paths ----
DIR_AGESEX_COMP <- file.path(DIR_BASE, "AgeSex_Comparison")
DIR_COMP <- file.path(DIR_AGESEX_COMP,
                      sprintf("Comparison_%s_vs_%s", REF_LABEL, COMP_LABEL))
DIR_FIG  <- file.path(DIR_COMP, "Figures")
DIR_TAB  <- file.path(DIR_COMP, "Tables")
for (d in c(DIR_AGESEX_COMP, DIR_COMP, DIR_FIG, DIR_TAB))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Colours ----
COL_REF  <- "#2166AC"   # Blue
COL_COMP <- "#B2182B"   # Red
COL_DIFF <- "#6C3483"   # Purple

# ---- Print header ----
cat(sprintf("\n%s\n", paste(rep("=", 70), collapse = "")))
cat("  Step 2.1S: Pairwise Subgroup Posterior Comparison (v2.3)\n")
cat(sprintf("%s\n", paste(rep("=", 70), collapse = "")))
cat(sprintf("  Reference group:  %s  [Display: '%s']\n", REF_LABEL, REF_DISPLAY))
cat(sprintf("  Comparison group: %s  [Display: '%s']\n", COMP_LABEL, COMP_DISPLAY))
cat(sprintf("  Output:           %s\n", DIR_COMP))
cat(sprintf("%s\n\n", paste(rep("=", 70), collapse = "")))

# ---- Utility functions ----

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(text = element_text(colour = "black"),
          plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
          axis.title = element_text(size = base_size, face = "bold"),
          axis.text = element_text(size = base_size - 1, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = base_size, face = "bold"),
          legend.background = element_blank())
}

save_both <- function(p, path_base, w = 7, h = 5.5) {
  tryCatch(ggsave(paste0(path_base, ".pdf"), p, width = w, height = h,
                  dpi = 300, device = cairo_pdf), error = function(e) NULL)
  tryCatch(ggsave(paste0(path_base, ".png"), p, width = w, height = h,
                  dpi = 300, bg = "white"), error = function(e) NULL)
  invisible(NULL)
}

fmt_p <- function(p) {
  ifelse(p >= 0.999, ">0.999",
         ifelse(p <= 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_est <- function(m, lo, hi, d = 2) {
  
  sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), m, lo, hi)
}


# ================================================================
#  PART 1: LOAD POSTERIOR DATA
# ================================================================

cat("=== PART 1: Loading posterior data ===\n")

load_posterior <- function(analysis, outcome) {
  path <- file.path(DIR_BASE, analysis, outcome,
                    sprintf("TC_%s_Step2.1_Posterior.rds", outcome))
  if (!file.exists(path)) {
    stop(sprintf("  ERROR: Posterior RDS not found:\n    %s\n  Run Step 2.1 for %s/%s first.",
                 path, analysis, outcome))
  }
  cat(sprintf("  Loaded: %s\n", path))
  readRDS(path)
}

load_projection <- function(analysis, outcome, type) {
  path <- file.path(DIR_BASE, analysis, outcome,
                    sprintf("TC_%s_%s_Projection_Draws.rds", outcome, type))
  if (file.exists(path)) {
    cat(sprintf("  Loaded: %s\n", path))
    return(readRDS(path))
  }
  cat(sprintf("  Optional file not found: TC_%s_%s_Projection_Draws.rds\n", outcome, type))
  NULL
}

ref_post  <- load_posterior(REF_ANALYSIS, REF_OUTCOME)
comp_post <- load_posterior(COMP_ANALYSIS, COMP_OUTCOME)

ref_hsr_proj  <- load_projection(REF_ANALYSIS, REF_OUTCOME, "HSR")
comp_hsr_proj <- load_projection(COMP_ANALYSIS, COMP_OUTCOME, "HSR")
ref_exp_proj  <- load_projection(REF_ANALYSIS, REF_OUTCOME, "EXP")
comp_exp_proj <- load_projection(COMP_ANALYSIS, COMP_OUTCOME, "EXP")


# ================================================================
#  PART 2: EXTRACT ALL DRAWS
# ================================================================

cat("\n=== PART 2: Extracting posterior draws ===\n")

# ---- 2a. Cumulative ERR draws (in % scale) ----
extract_cum <- function(post, hsr) {
  if (hsr == "All") {
    d <- post$cum_samples[["Active_to_PostFull"]]
    if (!is.null(d)) return(d)
    if (length(post$cum_samples) > 0) return(post$cum_samples[[1]])
    return(NULL)
  }
  if (!is.null(post$hsr_strata_cum_draws)) return(post$hsr_strata_cum_draws[[hsr]])
  NULL
}

ref_cum  <- extract_cum(ref_post, REF_HSR)
comp_cum <- extract_cum(comp_post, COMP_HSR)
has_cum  <- !is.null(ref_cum) && !is.null(comp_cum)

# ---- 2b. Within-stratum experience adaptation draws ----
# [v2.3 FIX] Now also extracts pre-computed ΔERR draws directly from Step 2.1
extract_adapt <- function(post, hsr) {
  if (hsr == "All") return(NULL)
  if (is.null(post$adapt_strata_draws)) return(NULL)
  ad <- post$adapt_strata_draws[[hsr]]
  if (is.null(ad)) return(NULL)
  if (is.list(ad) && "draws_std" %in% names(ad)) {
    return(list(
      theta_std  = ad$draws_std,
      theta_hlog = if ("draws_hlog" %in% names(ad)) ad$draws_hlog else NULL,
      beta_cum_draws = if ("beta_cum_draws" %in% names(ad)) ad$beta_cum_draws else NULL,
      pct_std    = if ("cum_std_pct" %in% names(ad)) ad$cum_std_pct
      else (exp(ad$draws_std) - 1) * 100,
      pct_hlog   = if ("cum_hlog_pct" %in% names(ad)) ad$cum_hlog_pct
      else if (!is.null(ad$draws_hlog)) (exp(ad$draws_hlog) - 1) * 100 else NULL,
      # [v2.3] Direct extraction of pre-computed ΔERR from Step 2.1
      derr_std   = if ("cum_std_derr" %in% names(ad)) ad$cum_std_derr else NULL,
      derr_hlog  = if ("cum_hlog_derr" %in% names(ad)) ad$cum_hlog_derr else NULL
    ))
  }
  if (is.numeric(ad)) {
    return(list(theta_std = ad, theta_hlog = NULL, beta_cum_draws = NULL,
                pct_std = (exp(ad) - 1) * 100, pct_hlog = NULL,
                derr_std = NULL, derr_hlog = NULL))
  }
  NULL
}

ref_adapt  <- extract_adapt(ref_post, REF_HSR)
comp_adapt <- extract_adapt(comp_post, COMP_HSR)
has_adapt_std  <- !is.null(ref_adapt$theta_std) && !is.null(comp_adapt$theta_std)
has_adapt_hlog <- !is.null(ref_adapt$theta_hlog) && !is.null(comp_adapt$theta_hlog)

# ---- 2c. Within-stratum HSR modification draws ----
# [v2.3 FIX] Now extracts both %RR and ΔERR draws
extract_within_hsr <- function(post, hsr) {
  if (hsr == "All") return(NULL)
  if (is.null(post$adapt_strata_hsr)) return(NULL)
  s <- post$adapt_strata_hsr[[hsr]]
  if (is.null(s)) return(NULL)
  list(
    cum_rr_pct   = if (!is.null(s$cum_draws)) s$cum_draws else NULL,
    cum_derr     = if (!is.null(s$cum_derr_draws)) s$cum_derr_draws else NULL
  )
}

ref_within_hsr_obj  <- extract_within_hsr(ref_post, REF_HSR)
comp_within_hsr_obj <- extract_within_hsr(comp_post, COMP_HSR)

ref_within_hsr  <- ref_within_hsr_obj$cum_rr_pct
comp_within_hsr <- comp_within_hsr_obj$cum_rr_pct
has_within_hsr  <- !is.null(ref_within_hsr) && !is.null(comp_within_hsr)

# [v2.3] Direct extraction of within-stratum HSR ΔERR
ref_within_hsr_derr  <- ref_within_hsr_obj$cum_derr
comp_within_hsr_derr <- comp_within_hsr_obj$cum_derr
has_within_hsr_derr  <- !is.null(ref_within_hsr_derr) && !is.null(comp_within_hsr_derr)

# ---- 2d. Overall HSR modification draws ----
ref_hsr_theta  <- if (!is.null(ref_hsr_proj)) ref_hsr_proj$theta_cum_draws else NULL
comp_hsr_theta <- if (!is.null(comp_hsr_proj)) comp_hsr_proj$theta_cum_draws else NULL
has_hsr_mod    <- !is.null(ref_hsr_theta) && !is.null(comp_hsr_theta)

# ---- 2e. Overall Experience modification draws ----
ref_exp_theta  <- if (!is.null(ref_exp_proj)) ref_exp_proj$theta_cum_draws else NULL
comp_exp_theta <- if (!is.null(comp_exp_proj)) comp_exp_proj$theta_cum_draws else NULL
has_exp_mod    <- !is.null(ref_exp_theta) && !is.null(comp_exp_theta)

# ---- 2f. Window-specific ERR summaries ----
ref_win  <- if (REF_HSR == "All") ref_post$window_results else NULL
comp_win <- if (COMP_HSR == "All") comp_post$window_results else NULL
has_win  <- !is.null(ref_win) && !is.null(comp_win)

# ---- 2g. Compute % RR change scale ----
if (has_hsr_mod) {
  ref_hsr_pct  <- (exp(ref_hsr_theta * HSR_STEP) - 1) * 100
  comp_hsr_pct <- (exp(comp_hsr_theta * HSR_STEP) - 1) * 100
}

if (has_exp_mod) {
  ref_exp_pct  <- (exp(ref_exp_theta) - 1) * 100   # per +1 unit H_log
  comp_exp_pct <- (exp(comp_exp_theta) - 1) * 100
}

# ---- 2h. ΔERR draws — DIRECT EXTRACTION from Step 2.1 [v2.3 FIX] ----
# Within-stratum adaptation ΔERR (per +1 SD) — directly from Step 2.1
has_adapt_derr_sd <- !is.null(ref_adapt$derr_std) && !is.null(comp_adapt$derr_std)
if (has_adapt_derr_sd) {
  ref_adapt_derr_sd  <- ref_adapt$derr_std
  comp_adapt_derr_sd <- comp_adapt$derr_std
  cat("  [v2.3] Within-stratum adaptation ΔERR (SD): extracted directly from Step 2.1\n")
}

# Within-stratum adaptation ΔERR (per +1 H_log) — directly from Step 2.1
has_adapt_derr_hlog <- !is.null(ref_adapt$derr_hlog) && !is.null(comp_adapt$derr_hlog)
if (has_adapt_derr_hlog) {
  ref_adapt_derr_hlog  <- ref_adapt$derr_hlog
  comp_adapt_derr_hlog <- comp_adapt$derr_hlog
  cat("  [v2.3] Within-stratum adaptation ΔERR (H_log): extracted directly from Step 2.1\n")
}

# Within-stratum HSR ΔERR — already extracted above (from adapt_strata_hsr$cum_derr_draws)
if (has_within_hsr_derr) {
  cat("  [v2.3] Within-stratum HSR ΔERR: extracted directly from Step 2.1\n")
}

# Overall HSR ΔERR — directly from posterior RDS
ref_hsr_derr  <- ref_post$overall_hsr_cum_derr_draws
comp_hsr_derr <- comp_post$overall_hsr_cum_derr_draws
has_hsr_derr  <- !is.null(ref_hsr_derr) && !is.null(comp_hsr_derr)
if (has_hsr_derr) {
  cat("  [v2.3] Overall HSR ΔERR: extracted directly from Step 2.1\n")
}

# Overall EXP ΔERR — directly from posterior RDS
ref_exp_derr  <- ref_post$overall_exp_cum_derr_draws
comp_exp_derr <- comp_post$overall_exp_cum_derr_draws
has_exp_derr  <- !is.null(ref_exp_derr) && !is.null(comp_exp_derr)
if (has_exp_derr) {
  cat("  [v2.3] Overall EXP ΔERR: extracted directly from Step 2.1\n")
}

# ---- 2i. Fallback: re-compute ΔERR only if pre-computed not available ----
# This uses beta_cum_draws from the SAME interaction model (consistent baseline)
compute_derr_consistent <- function(beta_cum_draws, theta_cum_draws, scale_by = 1) {
  # Both beta_cum and theta_cum MUST be from the SAME model
  rr_base <- exp(beta_cum_draws)
  derr <- rr_base * (exp(theta_cum_draws * scale_by) - 1) * 100
  return(derr)
}

# Fallback for within-stratum adaptation ΔERR (SD) if not pre-computed
if (!has_adapt_derr_sd && has_adapt_std) {
  if (!is.null(ref_adapt$beta_cum_draws) && !is.null(comp_adapt$beta_cum_draws)) {
    ref_adapt_derr_sd  <- compute_derr_consistent(ref_adapt$beta_cum_draws,
                                                  ref_adapt$theta_std, scale_by = 1)
    comp_adapt_derr_sd <- compute_derr_consistent(comp_adapt$beta_cum_draws,
                                                  comp_adapt$theta_std, scale_by = 1)
    has_adapt_derr_sd  <- TRUE
    cat("  [v2.3] Within-stratum ΔERR (SD): re-computed from consistent beta+theta\n")
  }
}

# Fallback for within-stratum adaptation ΔERR (H_log) if not pre-computed
if (!has_adapt_derr_hlog && has_adapt_hlog) {
  if (!is.null(ref_adapt$beta_cum_draws) && !is.null(comp_adapt$beta_cum_draws)) {
    ref_adapt_derr_hlog  <- compute_derr_consistent(ref_adapt$beta_cum_draws,
                                                    ref_adapt$theta_hlog, scale_by = 1)
    comp_adapt_derr_hlog <- compute_derr_consistent(comp_adapt$beta_cum_draws,
                                                    comp_adapt$theta_hlog, scale_by = 1)
    has_adapt_derr_hlog  <- TRUE
    cat("  [v2.3] Within-stratum ΔERR (H_log): re-computed from consistent beta+theta\n")
  }
}

# Fallback for overall HSR ΔERR if not pre-computed in posterior RDS
if (!has_hsr_derr && has_hsr_mod) {
  if (!is.null(ref_hsr_proj$beta_cum_draws) && !is.null(comp_hsr_proj$beta_cum_draws)) {
    ref_hsr_derr  <- compute_derr_consistent(ref_hsr_proj$beta_cum_draws,
                                             ref_hsr_theta, scale_by = HSR_STEP)
    comp_hsr_derr <- compute_derr_consistent(comp_hsr_proj$beta_cum_draws,
                                             comp_hsr_theta, scale_by = HSR_STEP)
    has_hsr_derr  <- TRUE
    cat("  [v2.3] Overall HSR ΔERR: re-computed from projection draws (consistent)\n")
  }
}

# Fallback for overall EXP ΔERR if not pre-computed in posterior RDS
if (!has_exp_derr && has_exp_mod) {
  if (!is.null(ref_exp_proj$beta_cum_draws) && !is.null(comp_exp_proj$beta_cum_draws)) {
    ref_exp_derr  <- compute_derr_consistent(ref_exp_proj$beta_cum_draws,
                                             ref_exp_theta, scale_by = 1)
    comp_exp_derr <- compute_derr_consistent(comp_exp_proj$beta_cum_draws,
                                             comp_exp_theta, scale_by = 1)
    has_exp_derr  <- TRUE
    cat("  [v2.3] Overall EXP ΔERR: re-computed from projection draws (consistent)\n")
  }
}

# ---- 2j. Print availability summary ----
avail <- data.frame(
  Code = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"),
  Comparison = c(
    "Window-specific ERR (visual comparison)",
    "Cumulative ERR",
    "Within-stratum adaptation Theta_cum (per +1 SD, log scale)",
    "Within-stratum adaptation Theta_cum (per +1 H_log, log scale)",
    "Within-stratum adaptation: % RR change per +1 SD",
    "Within-stratum adaptation: % RR change per +1 H_log",
    "Within-stratum adaptation: ΔERR per +1 SD",
    "Within-stratum adaptation: ΔERR per +1 H_log",
    "Within-stratum HSR: % RR change per +0.1 HSR",
    "Within-stratum HSR: ΔERR per +0.1 HSR",
    "Overall HSR: % RR change per +0.1 HSR",
    "Overall HSR: ΔERR per +0.1 HSR",
    "Overall EXP: % RR change per +1 H_log",
    "Overall EXP: ΔERR per +1 H_log"
  ),
  Available = c(has_win, has_cum, has_adapt_std, has_adapt_hlog,
                has_adapt_std, has_adapt_hlog,
                has_adapt_derr_sd, has_adapt_derr_hlog,
                has_within_hsr, has_within_hsr_derr,
                has_hsr_mod, has_hsr_derr,
                has_exp_mod, has_exp_derr),
  stringsAsFactors = FALSE
)

cat("\n  Availability summary:\n")
print(avail, row.names = FALSE, right = FALSE)
cat("\n")


# ================================================================
#  PART 3: PAIRWISE COMPARISONS
# ================================================================

cat("\n=== PART 3: Pairwise posterior comparisons ===\n")

compare_draws <- function(draws_A, draws_B, label_A, label_B, quantity, fig_ref = NA,
                          scale_name = "value") {
  if (is.null(draws_A) || is.null(draws_B)) return(NULL)
  n <- min(length(draws_A), length(draws_B))
  dA <- draws_A[1:n]; dB <- draws_B[1:n]
  diff_d <- dA - dB
  
  data.frame(
    Quantity       = quantity,
    Fig            = fig_ref,
    Scale          = scale_name,
    Group_Ref      = label_A,
    Group_Comp     = label_B,
    Ref_median     = median(dA),
    Ref_lo         = quantile(dA, .025, names = FALSE),
    Ref_hi         = quantile(dA, .975, names = FALSE),
    Ref_formatted  = fmt_est(median(dA), quantile(dA, .025, names = FALSE),
                             quantile(dA, .975, names = FALSE)),
    Comp_median    = median(dB),
    Comp_lo        = quantile(dB, .025, names = FALSE),
    Comp_hi        = quantile(dB, .975, names = FALSE),
    Comp_formatted = fmt_est(median(dB), quantile(dB, .025, names = FALSE),
                             quantile(dB, .975, names = FALSE)),
    Diff_median    = median(diff_d),
    Diff_lo        = quantile(diff_d, .025, names = FALSE),
    Diff_hi        = quantile(diff_d, .975, names = FALSE),
    Diff_formatted = fmt_est(median(diff_d), quantile(diff_d, .025, names = FALSE),
                             quantile(diff_d, .975, names = FALSE)),
    P_Ref_gt_Comp  = mean(diff_d > 0),
    P_Comp_gt_Ref  = mean(diff_d < 0),
    N_draws        = n,
    stringsAsFactors = FALSE
  )
}

comp_results <- list()

# --- B. Cumulative ERR ---
if (has_cum) {
  comp_results[["CumERR"]] <- compare_draws(
    ref_cum, comp_cum, REF_DISPLAY, COMP_DISPLAY,
    "Cumulative post-TC ERR (%)", fig_ref = "FigS-02/03", scale_name = "ERR (%)")
  cat(sprintf("\n  --- B. Cumulative ERR: Ref=%s, Comp=%s, Diff=%s ---\n",
              comp_results[["CumERR"]]$Ref_formatted,
              comp_results[["CumERR"]]$Comp_formatted,
              comp_results[["CumERR"]]$Diff_formatted))
}

# --- C. Within-stratum adaptation Theta per +1 SD (log scale) ---
if (has_adapt_std) {
  comp_results[["Adapt_Theta_SD"]] <- compare_draws(
    ref_adapt$theta_std, comp_adapt$theta_std, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum: Theta_cum (per +1 SD, log scale)",
    fig_ref = "FigS-04/05", scale_name = "log-scale coefficient")
  cat(sprintf("\n  --- C. Adapt Theta (SD): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_Theta_SD"]]$Ref_formatted,
              comp_results[["Adapt_Theta_SD"]]$Comp_formatted))
}

# --- D. Within-stratum adaptation Theta per +1 H_log (log scale) ---
if (has_adapt_hlog) {
  comp_results[["Adapt_Theta_Hlog"]] <- compare_draws(
    ref_adapt$theta_hlog, comp_adapt$theta_hlog, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum: Theta_cum (per +1 H_log, log scale)",
    fig_ref = "FigS-06/07", scale_name = "log-scale coefficient")
  cat(sprintf("\n  --- D. Adapt Theta (H_log): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_Theta_Hlog"]]$Ref_formatted,
              comp_results[["Adapt_Theta_Hlog"]]$Comp_formatted))
}

# --- E. Within-stratum adaptation: % RR change per +1 SD ---
if (has_adapt_std) {
  comp_results[["Adapt_RRpct_SD"]] <- compare_draws(
    ref_adapt$pct_std, comp_adapt$pct_std, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum adaptation: % RR change per +1 SD experience",
    fig_ref = "FigS-08/09", scale_name = "% RR change")
  cat(sprintf("\n  --- E. Adapt %%RR (SD): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_RRpct_SD"]]$Ref_formatted,
              comp_results[["Adapt_RRpct_SD"]]$Comp_formatted))
}

# --- F. Within-stratum adaptation: % RR change per +1 H_log ---
if (has_adapt_hlog) {
  comp_results[["Adapt_RRpct_Hlog"]] <- compare_draws(
    ref_adapt$pct_hlog, comp_adapt$pct_hlog, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum adaptation: % RR change per +1 H_log",
    fig_ref = "FigS-10/11", scale_name = "% RR change")
  cat(sprintf("\n  --- F. Adapt %%RR (H_log): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_RRpct_Hlog"]]$Ref_formatted,
              comp_results[["Adapt_RRpct_Hlog"]]$Comp_formatted))
}

# --- G. Within-stratum adaptation: ΔERR per +1 SD ---
if (has_adapt_derr_sd) {
  comp_results[["Adapt_dERR_SD"]] <- compare_draws(
    ref_adapt_derr_sd, comp_adapt_derr_sd, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum adaptation: ΔERR (pp) per +1 SD experience",
    fig_ref = "FigS-12/13", scale_name = "ΔERR (pp)")
  cat(sprintf("\n  --- G. Adapt ΔERR (SD): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_dERR_SD"]]$Ref_formatted,
              comp_results[["Adapt_dERR_SD"]]$Comp_formatted))
}

# --- H. Within-stratum adaptation: ΔERR per +1 H_log ---
if (has_adapt_derr_hlog) {
  comp_results[["Adapt_dERR_Hlog"]] <- compare_draws(
    ref_adapt_derr_hlog, comp_adapt_derr_hlog, REF_DISPLAY, COMP_DISPLAY,
    "Within-stratum adaptation: ΔERR (pp) per +1 H_log",
    fig_ref = "FigS-14/15", scale_name = "ΔERR (pp)")
  cat(sprintf("\n  --- H. Adapt ΔERR (H_log): Ref=%s, Comp=%s ---\n",
              comp_results[["Adapt_dERR_Hlog"]]$Ref_formatted,
              comp_results[["Adapt_dERR_Hlog"]]$Comp_formatted))
}

# --- I. Within-stratum HSR: % RR change per +0.1 HSR ---
if (has_within_hsr) {
  comp_results[["WithinHSR_RRpct"]] <- compare_draws(
    ref_within_hsr, comp_within_hsr, REF_DISPLAY, COMP_DISPLAY,
    sprintf("Within-stratum HSR: %% RR change per +%.1f HSR", HSR_STEP),
    fig_ref = "FigS-16/17", scale_name = "% RR change")
  cat(sprintf("\n  --- I. Within HSR %%RR: Ref=%s, Comp=%s ---\n",
              comp_results[["WithinHSR_RRpct"]]$Ref_formatted,
              comp_results[["WithinHSR_RRpct"]]$Comp_formatted))
}

# --- J. Within-stratum HSR: ΔERR per +0.1 HSR ---
if (has_within_hsr_derr) {
  comp_results[["WithinHSR_dERR"]] <- compare_draws(
    ref_within_hsr_derr, comp_within_hsr_derr, REF_DISPLAY, COMP_DISPLAY,
    sprintf("Within-stratum HSR: ΔERR (pp) per +%.1f HSR", HSR_STEP),
    fig_ref = "FigS-18/19", scale_name = "ΔERR (pp)")
  cat(sprintf("\n  --- J. Within HSR ΔERR: Ref=%s, Comp=%s ---\n",
              comp_results[["WithinHSR_dERR"]]$Ref_formatted,
              comp_results[["WithinHSR_dERR"]]$Comp_formatted))
}

# --- K. Overall HSR: % RR change per +0.1 HSR ---
if (has_hsr_mod) {
  comp_results[["HSR_Overall_RRpct"]] <- compare_draws(
    ref_hsr_pct, comp_hsr_pct, REF_DISPLAY, COMP_DISPLAY,
    sprintf("Overall HSR: %% RR change per +%.1f HSR (full-sample)", HSR_STEP),
    fig_ref = "FigS-20/21", scale_name = "% RR change")
  cat(sprintf("\n  --- K. Overall HSR %%RR: Ref=%s, Comp=%s ---\n",
              comp_results[["HSR_Overall_RRpct"]]$Ref_formatted,
              comp_results[["HSR_Overall_RRpct"]]$Comp_formatted))
}

# --- L. Overall HSR: ΔERR per +0.1 HSR ---
if (has_hsr_derr) {
  comp_results[["HSR_Overall_dERR"]] <- compare_draws(
    ref_hsr_derr, comp_hsr_derr, REF_DISPLAY, COMP_DISPLAY,
    sprintf("Overall HSR: ΔERR (pp) per +%.1f HSR (full-sample)", HSR_STEP),
    fig_ref = "FigS-22/23", scale_name = "ΔERR (pp)")
  cat(sprintf("\n  --- L. Overall HSR ΔERR: Ref=%s, Comp=%s ---\n",
              comp_results[["HSR_Overall_dERR"]]$Ref_formatted,
              comp_results[["HSR_Overall_dERR"]]$Comp_formatted))
}

# --- M. Overall EXP: % RR change per +1 H_log ---
if (has_exp_mod) {
  comp_results[["EXP_Overall_RRpct"]] <- compare_draws(
    ref_exp_pct, comp_exp_pct, REF_DISPLAY, COMP_DISPLAY,
    "Overall EXP: % RR change per +1 H_log (full-sample)",
    fig_ref = "FigS-24/25", scale_name = "% RR change")
  cat(sprintf("\n  --- M. Overall EXP %%RR: Ref=%s, Comp=%s ---\n",
              comp_results[["EXP_Overall_RRpct"]]$Ref_formatted,
              comp_results[["EXP_Overall_RRpct"]]$Comp_formatted))
}

# --- N. Overall EXP: ΔERR per +1 H_log ---
if (has_exp_derr) {
  comp_results[["EXP_Overall_dERR"]] <- compare_draws(
    ref_exp_derr, comp_exp_derr, REF_DISPLAY, COMP_DISPLAY,
    "Overall EXP: ΔERR (pp) per +1 H_log (full-sample)",
    fig_ref = "FigS-26/27", scale_name = "ΔERR (pp)")
  cat(sprintf("\n  --- N. Overall EXP ΔERR: Ref=%s, Comp=%s ---\n",
              comp_results[["EXP_Overall_dERR"]]$Ref_formatted,
              comp_results[["EXP_Overall_dERR"]]$Comp_formatted))
}


# ================================================================
#  PART 4: FIGURES
# ================================================================

cat("\n=== PART 4: Generating figures ===\n")

fig_counter <- 0
active_x <- which(window_labels == "Active")

# ---- Helper: comparison density plot ----
plot_comparison <- function(df_long, x_var, fill_var, xlab, figname,
                            vline_zero = TRUE, legend_pos = c(0.2, 0.85),
                            anno_text = NULL, w = 6, h = 4) {
  p <- ggplot(df_long, aes(x = .data[[x_var]], fill = .data[[fill_var]],
                           colour = .data[[fill_var]])) +
    geom_density(alpha = 0.3, linewidth = 0.7)
  
  if (vline_zero) p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
  
  p <- p +
    scale_fill_manual(values = setNames(c(COL_REF, COL_COMP),
                                        c(REF_DISPLAY, COMP_DISPLAY)),
                      name = "Group") +
    scale_colour_manual(values = setNames(c(COL_REF, COL_COMP),
                                          c(REF_DISPLAY, COMP_DISPLAY)),
                        name = "Group") +
    labs(x = xlab, y = "Posterior density") +
    theme_pub() +
    theme(legend.position = legend_pos)
  
  if (!is.null(anno_text))
    p <- p + annotate("text", x = Inf, y = Inf, label = anno_text,
                      hjust = 1.05, vjust = 1.3, size = 2.6, lineheight = 1.1)
  
  save_both(p, file.path(DIR_FIG, figname), w = w, h = h)
  cat(sprintf("    %s saved\n", figname))
  p
}

# ---- Helper: difference density plot ----
plot_difference <- function(diff_vec, xlab, figname, anno = NULL, w = 6, h = 4) {
  med <- median(diff_vec)
  if (is.null(anno)) {
    anno <- sprintf("Median diff: %s\nP(%s > %s) = %s",
                    fmt_est(med, quantile(diff_vec, .025),
                            quantile(diff_vec, .975)),
                    REF_DISPLAY, COMP_DISPLAY,
                    fmt_p(mean(diff_vec > 0)))
  }
  
  p <- ggplot(data.frame(x = diff_vec), aes(x = x)) +
    geom_density(fill = COL_DIFF, alpha = 0.35, colour = COL_DIFF, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = med, colour = COL_DIFF, linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf, label = anno,
             hjust = 1.05, vjust = 1.3, size = 2.8, lineheight = 1.0) +
    labs(x = xlab, y = "Posterior density") +
    theme_pub()
  
  save_both(p, file.path(DIR_FIG, figname), w = w, h = h)
  cat(sprintf("    %s saved\n", figname))
  p
}

# ---- Helper: make comparison annotation ----
make_anno <- function(draws_ref, draws_comp, d = 2) {
  sprintf("%s: %s\n%s: %s",
          REF_DISPLAY,
          fmt_est(median(draws_ref), quantile(draws_ref, .025),
                  quantile(draws_ref, .975), d = d),
          COMP_DISPLAY,
          fmt_est(median(draws_comp), quantile(draws_comp, .025),
                  quantile(draws_comp, .975), d = d))
}

# ---- Helper: make difference annotation with adaptation interpretation ----
make_diff_anno_adapt <- function(diff_vec) {
  sprintf("Median diff: %s\nP(%s more adaptive) = %s\nP(%s more adaptive) = %s",
          fmt_est(median(diff_vec), quantile(diff_vec, .025),
                  quantile(diff_vec, .975), d = 3),
          REF_DISPLAY, fmt_p(mean(diff_vec < 0)),
          COMP_DISPLAY, fmt_p(mean(diff_vec > 0)))
}

# ---- Helper: generate paired comparison + difference figures ----
generate_pair <- function(ref_draws, comp_draws, x_var_name, xlab_comp, xlab_diff,
                          figname_comp, figname_diff, d = 2,
                          diff_anno_fn = NULL) {
  n <- min(length(ref_draws), length(comp_draws))
  df <- bind_rows(
    data.frame(Group = REF_DISPLAY,  val = ref_draws[1:n]),
    data.frame(Group = COMP_DISPLAY, val = comp_draws[1:n])
  )
  df$Group <- factor(df$Group, levels = c(REF_DISPLAY, COMP_DISPLAY))
  names(df)[2] <- x_var_name
  
  plot_comparison(df, x_var_name, "Group", xlab_comp, figname_comp,
                  anno_text = make_anno(ref_draws, comp_draws, d = d))
  
  diff_v <- ref_draws[1:n] - comp_draws[1:n]
  
  if (!is.null(diff_anno_fn)) {
    da <- diff_anno_fn(diff_v)
  } else {
    da <- NULL
  }
  
  plot_difference(diff_v, xlab_diff, figname_diff, anno = da)
  
  write.xlsx(df, file.path(DIR_TAB, paste0(figname_comp, "_PlotData.xlsx")))
  write.xlsx(data.frame(difference = diff_v),
             file.path(DIR_TAB, paste0(figname_diff, "_Draws.xlsx")))
}


# ---------------------------------------------------------------
# FigS-01: Window-specific excess relative risk (ERR, %) of post-
#   tropical cyclone hospitalisation compared between subgroups
#   across the full event timeline (−3 to +6 months relative to TC
#   landfall). Each line represents the posterior mean ERR with 95%
#   credible intervals from stratum-specific Bayesian negative
#   binomial event-study models. The dashed vertical line marks the
#   active TC exposure period. Pre-TC windows provide a visual
#   assessment of the parallel trends assumption.
# ---------------------------------------------------------------
if (has_win) {
  fig_counter <- fig_counter + 1
  figname <- sprintf("FigS-%02d_WindowERR_Comparison", fig_counter)
  
  pd_ref  <- ref_win  %>% mutate(Group = REF_DISPLAY, x_pos = as.numeric(window))
  pd_comp <- comp_win %>% mutate(Group = COMP_DISPLAY, x_pos = as.numeric(window))
  pd_both <- bind_rows(pd_ref, pd_comp)
  pd_both$Group <- factor(pd_both$Group, levels = c(REF_DISPLAY, COMP_DISPLAY))
  
  p_win <- ggplot(pd_both, aes(x = x_pos, y = ERR_pct, colour = Group)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15,
                  position = position_dodge(0.4), linewidth = 0.6) +
    geom_point(position = position_dodge(0.4), size = 2.2) +
    geom_line(position = position_dodge(0.4), linewidth = 0.7) +
    scale_colour_manual(values = setNames(c(COL_REF, COL_COMP),
                                          c(REF_DISPLAY, COMP_DISPLAY)),
                        name = "Subgroup") +
    scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
    labs(x = "Lag period (month)", y = "ERR (%)") +
    theme_pub() +
    theme(legend.position = c(0.15, 0.85))
  
  save_both(p_win, file.path(DIR_FIG, figname), w = 7, h = 4.5)
  write.xlsx(pd_both, file.path(DIR_TAB, paste0(figname, "_PlotData.xlsx")))
  cat(sprintf("    %s saved\n", figname))
} else {
  cat("  FigS-01 (Window ERR): SKIPPED\n")
}


# ---------------------------------------------------------------
# FigS-02: Posterior density comparison of the cumulative post-TC
#   excess relative risk (ERR, %) between subgroups, defined as
#   (exp(Σβ_w) − 1) × 100% summed over all post-TC windows (Active
#   through +6 months). Each density represents the full marginal
#   posterior from a separate stratum-specific Bayesian event-study
#   model. Positive values indicate net elevated hospitalisation
#   risk following TC exposure. The dashed line marks ERR = 0.
#
# FigS-03: Posterior density of the pairwise difference in
#   cumulative post-TC ERR (percentage points) between subgroups
#   (Reference minus Comparison), derived from paired posterior
#   draws. The solid vertical line marks the posterior median; the
#   dashed line marks zero (null hypothesis of no difference).
#   Annotations report the median difference, 95% credible interval,
#   and the posterior probability that one group exceeds the other.
# ---------------------------------------------------------------
if (has_cum) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_CumERR_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_CumERR_Difference", fig_counter)
  
  generate_pair(ref_cum, comp_cum, "ERR", "Cumulative post-TC ERR (%)",
                sprintf("Diff in cumulative ERR (pp): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-04: Posterior density comparison of the cumulative
#   experience-adaptation coefficient (Θ_cum, per +1 SD of
#   historical TC experience) between subgroups, estimated from
#   within-stratum experience × window interaction models. Θ_cum =
#   Σ θ_w is the sum of window-specific interaction terms on the
#   log-rate-ratio scale. More negative values indicate stronger
#   experience-based adaptation, whereby prior TC exposure more
#   effectively attenuates future post-TC hospitalisation risk.
#
# FigS-05: Posterior density of the pairwise difference in Θ_cum
#   (per +1 SD) between subgroups (Reference minus Comparison).
#   Negative difference indicates the reference group exhibits
#   stronger per-SD experiential adaptation. Annotations report
#   the posterior probability that each group is more adaptive.
# ---------------------------------------------------------------
if (has_adapt_std) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptTheta_SD_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptTheta_SD_Difference", fig_counter)
  
  generate_pair(ref_adapt$theta_std, comp_adapt$theta_std, "Theta",
                expression(Theta[cum]~"(per +1 SD experience, log scale)"),
                sprintf("Diff in Theta_cum (per SD): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff, d = 3,
                diff_anno_fn = make_diff_anno_adapt)
}


# ---------------------------------------------------------------
# FigS-06: Posterior density comparison of the cumulative
#   experience-adaptation coefficient (Θ_cum, per +1 unit of
#   log-transformed historical TC experience, log_e(H+1)) between
#   subgroups. Unlike the standardised (per-SD) metric in FigS-04,
#   this preserves the original scale: per +1 unit H_log corresponds
#   to an approximately e-fold (2.72×) increase in cumulative
#   historical TC energy exposure. More negative Θ_cum indicates
#   stronger adaptation per unit of log-experience.
#
# FigS-07: Posterior density of the pairwise difference in Θ_cum
#   (per +1 H_log) between subgroups (Reference minus Comparison).
#   Negative difference indicates the reference group derives
#   stronger protection per unit of log-scaled historical TC
#   experience.
# ---------------------------------------------------------------
if (has_adapt_hlog) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptTheta_Hlog_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptTheta_Hlog_Difference", fig_counter)
  
  generate_pair(ref_adapt$theta_hlog, comp_adapt$theta_hlog, "Theta",
                expression(Theta[cum]~"(per +1 H_log, log scale)"),
                sprintf("Diff in Theta_cum (per H_log): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff, d = 3,
                diff_anno_fn = make_diff_anno_adapt)
}


# ---------------------------------------------------------------
# FigS-08: Posterior density comparison of the within-stratum
#   experience-adaptation effect expressed as the relative change
#   (%) in the cumulative post-TC rate ratio per +1 SD of
#   historical TC experience: (exp(Θ_cum) − 1) × 100%. This is
#   the ratio-of-rate-ratios (RRR) on the percentage scale,
#   quantifying how much the cumulative RR is multiplicatively
#   reduced per SD of prior experience. Negative values indicate
#   adaptation. This metric is independent of the baseline ERR
#   magnitude.
#
# FigS-09: Posterior density of the pairwise difference in % RR
#   change (per +1 SD experience) between subgroups (Reference
#   minus Comparison). More negative difference indicates stronger
#   relative adaptation in the reference group.
# ---------------------------------------------------------------
if (has_adapt_std) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptRRpct_SD_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptRRpct_SD_Difference", fig_counter)
  
  generate_pair(ref_adapt$pct_std, comp_adapt$pct_std, "pct",
                "% change in cumulative RR per +1 SD experience",
                sprintf("Diff in %%RR change (per SD): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-10: Posterior density comparison of the within-stratum
#   experience-adaptation effect expressed as the relative change
#   (%) in the cumulative rate ratio per +1 unit H_log (log-
#   transformed historical TC experience index). Interpretation
#   parallels FigS-08 but on the original log-transformed
#   experience scale rather than the standardised (per-SD) scale.
#   Negative values indicate experience-based adaptation.
#
# FigS-11: Posterior density of the pairwise difference in % RR
#   change (per +1 H_log) between subgroups (Reference minus
#   Comparison). Negative difference indicates the reference group
#   exhibits stronger per-unit adaptation on the log-experience
#   scale.
# ---------------------------------------------------------------
if (has_adapt_hlog) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptRRpct_Hlog_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptRRpct_Hlog_Difference", fig_counter)
  
  generate_pair(ref_adapt$pct_hlog, comp_adapt$pct_hlog, "pct",
                "% change in cumulative RR per +1 H_log",
                sprintf("Diff in %%RR change (per H_log): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-12: Posterior density comparison of the absolute change in
#   cumulative ERR (ΔERR, percentage points) per +1 SD of
#   historical TC experience between subgroups. ΔERR = exp(β_cum) ×
#   (exp(Θ_cum) − 1) × 100%, where β_cum and Θ_cum are jointly
#   estimated from the same within-stratum interaction model. Unlike
#   the % RR change (FigS-08), ΔERR depends on the baseline post-TC
#   effect: for subgroups with higher baseline cumulative ERR, the
#   same Θ_cum translates to a larger absolute ERR reduction.
#   Negative values indicate that experience reduces the absolute
#   post-TC excess hospitalisation rate.
#
# FigS-13: Posterior density of the pairwise difference in ΔERR
#   (per +1 SD experience) between subgroups (Reference minus
#   Comparison). More negative difference indicates the reference
#   group derives a larger absolute reduction in cumulative ERR
#   from each additional SD of historical TC experience.
# ---------------------------------------------------------------
if (has_adapt_derr_sd) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptDERR_SD_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptDERR_SD_Difference", fig_counter)
  
  generate_pair(ref_adapt_derr_sd, comp_adapt_derr_sd, "DERR",
                "ΔERR (pp): Absolute change in cumulative ERR per +1 SD experience",
                sprintf("Diff in ΔERR (per SD): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-14: Posterior density comparison of the absolute change in
#   cumulative ERR (ΔERR, percentage points) per +1 unit H_log
#   between subgroups. ΔERR is evaluated on the original log-
#   transformed experience scale, capturing how much the absolute
#   post-TC excess risk decreases per unit increase in log(H+1).
#   Both β_cum and Θ_cum are extracted from the same within-stratum
#   interaction model to ensure consistency.
#
# FigS-15: Posterior density of the pairwise difference in ΔERR
#   (per +1 H_log) between subgroups (Reference minus Comparison).
#   More negative values indicate that the reference group benefits
#   from a larger absolute ERR reduction per unit of log-scaled
#   historical TC experience.
# ---------------------------------------------------------------
if (has_adapt_derr_hlog) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_AdaptDERR_Hlog_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_AdaptDERR_Hlog_Difference", fig_counter)
  
  generate_pair(ref_adapt_derr_hlog, comp_adapt_derr_hlog, "DERR",
                "ΔERR (pp): Absolute change in cumulative ERR per +1 H_log",
                sprintf("Diff in ΔERR (per H_log): %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-16: Posterior density comparison of the within-stratum HSR
#   modification effect, expressed as the relative change (%) in
#   the cumulative rate ratio per +0.1 unit increase in the health
#   system resilience (HSR) index. This quantifies how residual HSR
#   variation among cities already classified within the same HSR
#   stratum modifies the cumulative post-TC hospitalisation risk.
#   Negative values indicate that higher within-stratum HSR is
#   protective. This metric is a ratio-of-rate-ratios and does not
#   depend on the baseline ERR magnitude.
#
# FigS-17: Posterior density of the pairwise difference in within-
#   stratum HSR % RR change between subgroups (Reference minus
#   Comparison). Positive difference indicates that within-stratum
#   HSR variation is less protective (or more deleterious) for the
#   reference group relative to the comparison group.
# ---------------------------------------------------------------
if (has_within_hsr) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_WithinHSR_RRpct_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_WithinHSR_RRpct_Difference", fig_counter)
  
  generate_pair(ref_within_hsr, comp_within_hsr, "pct",
                sprintf("Within-stratum: %% change in cumulative RR per +%.1f HSR", HSR_STEP),
                sprintf("Diff in within-HSR %%RR change: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-18: Posterior density comparison of the within-stratum HSR
#   modification effect expressed as the absolute change in
#   cumulative ERR (ΔERR, percentage points) per +0.1 unit HSR.
#   ΔERR = exp(β_cum) × (exp(Θ_cum × 0.1) − 1) × 100%, with β_cum
#   and Θ_cum from the same within-stratum HSR interaction model.
#   Negative values indicate that higher within-stratum HSR reduces
#   the absolute post-TC excess hospitalisation rate. Unlike % RR
#   change, ΔERR is sensitive to the baseline effect size.
#
# FigS-19: Posterior density of the pairwise difference in within-
#   stratum HSR ΔERR between subgroups (Reference minus Comparison).
#   More positive difference indicates that HSR is less effective at
#   reducing the absolute post-TC ERR in the reference group.
# ---------------------------------------------------------------
if (has_within_hsr_derr) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_WithinHSR_DERR_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_WithinHSR_DERR_Difference", fig_counter)
  
  generate_pair(ref_within_hsr_derr, comp_within_hsr_derr, "DERR",
                sprintf("Within-stratum: ΔERR (pp) per +%.1f HSR", HSR_STEP),
                sprintf("Diff in within-HSR ΔERR: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-20: Posterior density comparison of the overall HSR
#   modification effect (% change in cumulative rate ratio per +0.1
#   unit HSR) between subgroups, estimated from full-sample
#   HSR × window interaction models fitted to all cities regardless
#   of HSR stratum. This captures the population-average gradient
#   of cumulative post-TC risk with respect to health system
#   resilience. Negative values indicate that HSR is protective
#   against post-TC hospitalisation.
#
# FigS-21: Posterior density of the pairwise difference in overall
#   HSR % RR change between subgroups (Reference minus Comparison).
#   Positive values indicate that the cumulative rate ratio for the
#   reference group is less attenuated by HSR than that of the
#   comparison group.
# ---------------------------------------------------------------
if (has_hsr_mod) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_OverallHSR_RRpct_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_OverallHSR_RRpct_Difference", fig_counter)
  
  generate_pair(ref_hsr_pct, comp_hsr_pct, "pct",
                sprintf("Overall: %% change in cumulative RR per +%.1f HSR", HSR_STEP),
                sprintf("Diff in overall HSR %%RR change: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-22: Posterior density comparison of the overall HSR ΔERR
#   (percentage points per +0.1 HSR) between subgroups. ΔERR
#   translates the relative rate-ratio change into the absolute ERR
#   scale, computed as exp(β_cum) × (exp(Θ_cum × 0.1) − 1) × 100%
#   with both β_cum and Θ_cum from the same full-sample HSR
#   interaction model. Differences between subgroups reflect both
#   the modifier strength (Θ) and the baseline post-TC risk (β).
#
# FigS-23: Posterior density of the pairwise difference in overall
#   HSR ΔERR between subgroups (Reference minus Comparison). More
#   negative difference indicates the comparison group benefits
#   more (in absolute ERR terms) from health system resilience.
# ---------------------------------------------------------------
if (has_hsr_derr) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_OverallHSR_DERR_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_OverallHSR_DERR_Difference", fig_counter)
  
  generate_pair(ref_hsr_derr, comp_hsr_derr, "DERR",
                sprintf("Overall: ΔERR (pp) per +%.1f HSR", HSR_STEP),
                sprintf("Diff in overall HSR ΔERR: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-24: Posterior density comparison of the overall experience
#   modification effect (% change in cumulative rate ratio per +1
#   unit H_log) between subgroups, from full-sample experience ×
#   window interaction models. This estimates how the cumulative
#   post-TC hospitalisation rate ratio changes per unit increase in
#   log-transformed historical TC experience across the entire
#   sample. Negative values indicate experience-based adaptation
#   (prior TC exposure is protective).
#
# FigS-25: Posterior density of the pairwise difference in overall
#   experience % RR change between subgroups (Reference minus
#   Comparison). Negative difference indicates the reference group
#   derives greater relative adaptation benefit per unit of log-
#   scaled historical TC experience.
# ---------------------------------------------------------------
if (has_exp_mod) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_OverallEXP_RRpct_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_OverallEXP_RRpct_Difference", fig_counter)
  
  generate_pair(ref_exp_pct, comp_exp_pct, "pct",
                "Overall: % change in cumulative RR per +1 H_log",
                sprintf("Diff in overall EXP %%RR change: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}


# ---------------------------------------------------------------
# FigS-26: Posterior density comparison of the overall experience
#   ΔERR (percentage points per +1 unit H_log) between subgroups.
#   This converts the relative % RR change into the absolute ERR
#   scale using exp(β_cum) × (exp(Θ_cum) − 1) × 100%, where both
#   parameters are jointly estimated from the same full-sample
#   experience interaction model. This quantifies how much the
#   absolute cumulative excess hospitalisation risk changes per
#   unit of log-experience, conditional on each subgroup's baseline
#   post-TC effect.
#
# FigS-27: Posterior density of the pairwise difference in overall
#   experience ΔERR between subgroups (Reference minus Comparison).
#   More negative difference indicates the reference group has a
#   larger absolute ERR reduction per unit of log-scaled historical
#   TC experience.
# ---------------------------------------------------------------
if (has_exp_derr) {
  fig_counter <- fig_counter + 1
  figname_comp <- sprintf("FigS-%02d_OverallEXP_DERR_Comparison", fig_counter)
  fig_counter <- fig_counter + 1
  figname_diff <- sprintf("FigS-%02d_OverallEXP_DERR_Difference", fig_counter)
  
  generate_pair(ref_exp_derr, comp_exp_derr, "DERR",
                "Overall: ΔERR (pp) per +1 H_log",
                sprintf("Diff in overall EXP ΔERR: %s minus %s",
                        REF_DISPLAY, COMP_DISPLAY),
                figname_comp, figname_diff)
}

cat(sprintf("\n  Total figures generated: %d\n", fig_counter))


# ================================================================
#  PART 5: SUMMARY TABLE
# ================================================================

cat("\n=== PART 5: Summary comparison table ===\n")

summary_all <- bind_rows(comp_results)

# Add interpretation column
summary_all$Interpretation <- sapply(seq_len(nrow(summary_all)), function(i) {
  r <- summary_all[i, ]
  p_ab <- r$P_Ref_gt_Comp
  if (p_ab > 0.975 || p_ab < 0.025)
    return("Strong evidence of difference (P > 0.975 or < 0.025)")
  if (p_ab > 0.95  || p_ab < 0.05)
    return("Moderate evidence of difference (P in 0.95-0.975)")
  if (p_ab > 0.90  || p_ab < 0.10)
    return("Weak evidence of difference (P in 0.90-0.95)")
  "No clear difference"
})

cat("\n")
print(summary_all %>%
        dplyr::select(Quantity, Fig, Scale, Ref_formatted, Comp_formatted,
                      Diff_formatted, P_Ref_gt_Comp, Interpretation),
      right = FALSE, row.names = FALSE)

write.xlsx(summary_all,
           file.path(DIR_TAB, "Summary_Comparison_AllMetrics.xlsx"),
           overwrite = TRUE)
cat("  Summary table saved\n")


# ================================================================
#  PART 6: POSTERIOR PROBABILITY TABLE
# ================================================================

cat("\n=== PART 6: Posterior probability table ===\n")

prow <- function(Analysis, Quantity, Statement, P, Interpretation)
  data.frame(Analysis = Analysis, Quantity = Quantity, Statement = Statement,
             Probability = round(P, 4), Interpretation = Interpretation,
             stringsAsFactors = FALSE)

pv <- list()

# --- Group-level P(ERR > 0) ---
if (!is.null(ref_cum))
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Cumulative ERR", "P(ERR > 0)",
    mean(ref_cum > 0),
    "Post-TC hospitalisation elevated in reference group")
if (!is.null(comp_cum))
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Cumulative ERR", "P(ERR > 0)",
    mean(comp_cum > 0),
    "Post-TC hospitalisation elevated in comparison group")

# --- Cumulative ERR comparison ---
if (has_cum) {
  n_c <- min(length(ref_cum), length(comp_cum))
  d_c <- ref_cum[1:n_c] - comp_cum[1:n_c]
  pv[[length(pv) + 1]] <- prow(
    "Cumulative ERR comparison",
    sprintf("%s minus %s", REF_DISPLAY, COMP_DISPLAY),
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_c > 0),
    "Reference group has higher cumulative post-TC ERR")
}

# --- Within-stratum adaptation per SD ---
if (has_adapt_std) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Adaptation Theta_cum (per +1 SD)",
    "P(Theta < 0)", mean(ref_adapt$theta_std < 0),
    "Evidence of experience-based adaptation in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Adaptation Theta_cum (per +1 SD)",
    "P(Theta < 0)", mean(comp_adapt$theta_std < 0),
    "Evidence of experience-based adaptation in comparison group")
  n_a <- min(length(ref_adapt$theta_std), length(comp_adapt$theta_std))
  d_a <- ref_adapt$theta_std[1:n_a] - comp_adapt$theta_std[1:n_a]
  pv[[length(pv) + 1]] <- prow(
    "Adaptation comparison (per +1 SD)",
    "Theta difference",
    sprintf("P(%s more adaptive)", REF_DISPLAY),
    mean(d_a < 0),
    "Reference group has stronger experience-based adaptation (per SD)")
}

# --- Within-stratum adaptation per H_log ---
if (has_adapt_hlog) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Adaptation Theta_cum (per +1 H_log)",
    "P(Theta < 0)", mean(ref_adapt$theta_hlog < 0),
    "Evidence of adaptation (H_log scale) in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Adaptation Theta_cum (per +1 H_log)",
    "P(Theta < 0)", mean(comp_adapt$theta_hlog < 0),
    "Evidence of adaptation (H_log scale) in comparison group")
  n_a <- min(length(ref_adapt$theta_hlog), length(comp_adapt$theta_hlog))
  d_ah <- ref_adapt$theta_hlog[1:n_a] - comp_adapt$theta_hlog[1:n_a]
  pv[[length(pv) + 1]] <- prow(
    "Adaptation comparison (per +1 H_log)",
    "Theta difference",
    sprintf("P(%s more adaptive)", REF_DISPLAY),
    mean(d_ah < 0),
    "Reference group has stronger adaptation (per H_log)")
}

# --- Within-stratum adaptation ΔERR per SD ---
if (has_adapt_derr_sd) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Adaptation ΔERR (per +1 SD)",
    "P(ΔERR < 0)", mean(ref_adapt_derr_sd < 0),
    "Experience reduces absolute ERR in reference group (per SD)")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Adaptation ΔERR (per +1 SD)",
    "P(ΔERR < 0)", mean(comp_adapt_derr_sd < 0),
    "Experience reduces absolute ERR in comparison group (per SD)")
  n_a <- min(length(ref_adapt_derr_sd), length(comp_adapt_derr_sd))
  d_a <- ref_adapt_derr_sd[1:n_a] - comp_adapt_derr_sd[1:n_a]
  pv[[length(pv) + 1]] <- prow(
    "Adaptation ΔERR comparison (per +1 SD)",
    "ΔERR difference",
    sprintf("P(%s more negative ΔERR)", REF_DISPLAY),
    mean(d_a < 0),
    "Reference group has larger absolute ERR reduction per SD experience")
}

# --- Within-stratum adaptation ΔERR per H_log ---
if (has_adapt_derr_hlog) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Adaptation ΔERR (per +1 H_log)",
    "P(ΔERR < 0)", mean(ref_adapt_derr_hlog < 0),
    "Experience reduces absolute ERR in reference group (per H_log)")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Adaptation ΔERR (per +1 H_log)",
    "P(ΔERR < 0)", mean(comp_adapt_derr_hlog < 0),
    "Experience reduces absolute ERR in comparison group (per H_log)")
  n_a <- min(length(ref_adapt_derr_hlog), length(comp_adapt_derr_hlog))
  d_a <- ref_adapt_derr_hlog[1:n_a] - comp_adapt_derr_hlog[1:n_a]
  pv[[length(pv) + 1]] <- prow(
    "Adaptation ΔERR comparison (per +1 H_log)",
    "ΔERR difference",
    sprintf("P(%s more negative ΔERR)", REF_DISPLAY),
    mean(d_a < 0),
    "Reference group has larger absolute ERR reduction per H_log")
}

# --- Within-stratum HSR ---
if (has_within_hsr) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Within-stratum HSR (per +0.1)",
    "P(HSR protective: %RR < 0)", mean(ref_within_hsr < 0),
    "Within-stratum HSR reduces post-TC RR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Within-stratum HSR (per +0.1)",
    "P(HSR protective: %RR < 0)", mean(comp_within_hsr < 0),
    "Within-stratum HSR reduces post-TC RR in comparison group")
  n_w <- min(length(ref_within_hsr), length(comp_within_hsr))
  d_w <- ref_within_hsr[1:n_w] - comp_within_hsr[1:n_w]
  pv[[length(pv) + 1]] <- prow(
    "Within-stratum HSR comparison",
    "%RR change difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_w > 0),
    "HSR less protective for reference group (within-stratum)")
}

# --- Within-stratum HSR ΔERR ---
if (has_within_hsr_derr) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Within-stratum HSR ΔERR (per +0.1)",
    "P(ΔERR < 0)", mean(ref_within_hsr_derr < 0),
    "Within-stratum HSR reduces absolute ERR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Within-stratum HSR ΔERR (per +0.1)",
    "P(ΔERR < 0)", mean(comp_within_hsr_derr < 0),
    "Within-stratum HSR reduces absolute ERR in comparison group")
  n_w <- min(length(ref_within_hsr_derr), length(comp_within_hsr_derr))
  d_w <- ref_within_hsr_derr[1:n_w] - comp_within_hsr_derr[1:n_w]
  pv[[length(pv) + 1]] <- prow(
    "Within-stratum HSR ΔERR comparison",
    "ΔERR difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_w > 0),
    "HSR less effective at reducing absolute ERR for reference group")
}

# --- Overall HSR ---
if (has_hsr_mod) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Overall HSR (per +0.1)",
    "P(HSR protective: %RR < 0)", mean(ref_hsr_pct < 0),
    "Overall HSR reduces post-TC RR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Overall HSR (per +0.1)",
    "P(HSR protective: %RR < 0)", mean(comp_hsr_pct < 0),
    "Overall HSR reduces post-TC RR in comparison group")
  n_h <- min(length(ref_hsr_pct), length(comp_hsr_pct))
  d_h <- ref_hsr_pct[1:n_h] - comp_hsr_pct[1:n_h]
  pv[[length(pv) + 1]] <- prow(
    "Overall HSR comparison",
    "%RR change difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_h > 0),
    "HSR less protective for reference group (overall)")
}

# --- Overall HSR ΔERR ---
if (has_hsr_derr) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Overall HSR ΔERR (per +0.1)",
    "P(ΔERR < 0)", mean(ref_hsr_derr < 0),
    "Overall HSR reduces absolute ERR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Overall HSR ΔERR (per +0.1)",
    "P(ΔERR < 0)", mean(comp_hsr_derr < 0),
    "Overall HSR reduces absolute ERR in comparison group")
  n_h <- min(length(ref_hsr_derr), length(comp_hsr_derr))
  d_h <- ref_hsr_derr[1:n_h] - comp_hsr_derr[1:n_h]
  pv[[length(pv) + 1]] <- prow(
    "Overall HSR ΔERR comparison",
    "ΔERR difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_h > 0),
    "HSR less effective at reducing absolute ERR for reference group (overall)")
}

# --- Overall EXP ---
if (has_exp_mod) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Overall EXP (per +1 H_log)",
    "P(EXP protective: %RR < 0)", mean(ref_exp_pct < 0),
    "Overall experience reduces post-TC RR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Overall EXP (per +1 H_log)",
    "P(EXP protective: %RR < 0)", mean(comp_exp_pct < 0),
    "Overall experience reduces post-TC RR in comparison group")
  n_e <- min(length(ref_exp_pct), length(comp_exp_pct))
  d_e <- ref_exp_pct[1:n_e] - comp_exp_pct[1:n_e]
  pv[[length(pv) + 1]] <- prow(
    "Overall EXP comparison",
    "%RR change difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_e > 0),
    "Experience less protective for reference group (overall)")
}

# --- Overall EXP ΔERR ---
if (has_exp_derr) {
  pv[[length(pv) + 1]] <- prow(
    sprintf("REF: %s", REF_DISPLAY), "Overall EXP ΔERR (per +1 H_log)",
    "P(ΔERR < 0)", mean(ref_exp_derr < 0),
    "Overall experience reduces absolute ERR in reference group")
  pv[[length(pv) + 1]] <- prow(
    sprintf("COMP: %s", COMP_DISPLAY), "Overall EXP ΔERR (per +1 H_log)",
    "P(ΔERR < 0)", mean(comp_exp_derr < 0),
    "Overall experience reduces absolute ERR in comparison group")
  n_e <- min(length(ref_exp_derr), length(comp_exp_derr))
  d_e <- ref_exp_derr[1:n_e] - comp_exp_derr[1:n_e]
  pv[[length(pv) + 1]] <- prow(
    "Overall EXP ΔERR comparison",
    "ΔERR difference",
    sprintf("P(%s > %s)", REF_DISPLAY, COMP_DISPLAY),
    mean(d_e > 0),
    "Experience less effective at reducing absolute ERR for reference group (overall)")
}

# Theta comparison across strata (adaptation strength)
if (has_adapt_std) {
  n_a <- min(length(ref_adapt$theta_std), length(comp_adapt$theta_std))
  d_a <- ref_adapt$theta_std[1:n_a] - comp_adapt$theta_std[1:n_a]
  pv[[length(pv) + 1]] <- prow(
    "Adaptation_HighvsLow", "Theta_Ref vs Theta_Comp",
    sprintf("P(%s more protective per SD)", REF_DISPLAY),
    mean(d_a < 0),
    "Reference HSR stratum facilitates stronger experience-based adaptation")
}

pvals <- bind_rows(pv)
cat("\n")
print(pvals, right = FALSE, row.names = FALSE)

write.xlsx(pvals,
           file.path(DIR_TAB, "PosteriorProbabilities.xlsx"),
           overwrite = TRUE)
cat("  Posterior probability table saved\n")


# ================================================================
#  PART 7: WINDOW-SPECIFIC ERR COMPARISON TABLE
# ================================================================

cat("\n=== PART 7: Window-specific ERR comparison table ===\n")

if (has_win) {
  win_merged <- ref_win %>%
    dplyr::select(window, ERR_pct, lower, upper) %>%
    rename(Ref_ERR = ERR_pct, Ref_lo = lower, Ref_hi = upper) %>%
    left_join(
      comp_win %>%
        dplyr::select(window, ERR_pct, lower, upper) %>%
        rename(Comp_ERR = ERR_pct, Comp_lo = lower, Comp_hi = upper),
      by = "window"
    ) %>%
    mutate(
      Ref_formatted  = sprintf("%.2f%% (%.2f, %.2f)", Ref_ERR, Ref_lo, Ref_hi),
      Comp_formatted = sprintf("%.2f%% (%.2f, %.2f)", Comp_ERR, Comp_lo, Comp_hi),
      Diff_point     = Ref_ERR - Comp_ERR
    )
  
  names(win_merged)[names(win_merged) == "Ref_formatted"]  <- paste0(REF_DISPLAY, "_ERR")
  names(win_merged)[names(win_merged) == "Comp_formatted"] <- paste0(COMP_DISPLAY, "_ERR")
  
  cat("  Window-specific differences are point-estimate-based.\n")
  print(win_merged %>% dplyr::select(window, ends_with("_ERR"), Diff_point),
        right = FALSE, row.names = FALSE)
  
  write.xlsx(win_merged,
             file.path(DIR_TAB, "WindowERR_Comparison_Summary.xlsx"),
             overwrite = TRUE)
  cat("  Window-specific comparison table saved\n")
} else {
  cat("  Window-specific ERR comparison: SKIPPED (requires HSR = 'All' for both groups)\n")
}


# ================================================================
#  PART 8: MASTER EXCEL WORKBOOK
# ================================================================

cat("\n=== PART 8: Master Excel workbook ===\n")

wb <- createWorkbook()
add_sheet <- function(nm, x) {
  if (!is.null(x) && NROW(x) > 0) {
    addWorksheet(wb, nm)
    writeData(wb, nm, x)
  }
}

meta <- data.frame(
  Field = c("Reference group (REF)", "Comparison group (COMP)",
            "REF_HSR", "REF_ANALYSIS", "REF_OUTCOME",
            "COMP_HSR", "COMP_ANALYSIS", "COMP_OUTCOME",
            "REF display name", "COMP display name",
            "N_SAMPLES", "HSR_STEP", "Created",
            "", "--- KEY TERMINOLOGY ---",
            "Overall",
            "Within-stratum",
            "ERR (%)",
            "% RR change",
            "ΔERR (pp)",
            "Theta_cum (log scale)",
            "",
            "--- CONVERSION FORMULA ---",
            "Formula",
            "",
            "--- v2.3 BUG FIX NOTE ---",
            "Fix description"),
  Value = c(REF_LABEL, COMP_LABEL,
            REF_HSR, REF_ANALYSIS, REF_OUTCOME,
            COMP_HSR, COMP_ANALYSIS, COMP_OUTCOME,
            REF_DISPLAY, COMP_DISPLAY,
            as.character(N_SAMPLES), as.character(HSR_STEP),
            as.character(Sys.time()),
            "",
            "",
            "Interaction model fitted to ALL cities (full-sample, ignoring HSR strata)",
            "Interaction model fitted ONLY within a specific HSR stratum (High/Moderate/Low)",
            "Excess Relative Risk = (exp(beta_cum) - 1) x 100%. The cumulative post-TC risk.",
            "(exp(theta_cum x Delta) - 1) x 100%. Ratio-of-rate-ratios. Baseline-INDEPENDENT.",
            "exp(beta_cum) x (exp(theta_cum x Delta) - 1) x 100. Both beta and theta from SAME model.",
            "Sum of window-specific interaction coefficients. More negative = more protective.",
            "",
            "",
            "ΔERR (pp) = exp(beta_cum) x (exp(theta_cum x Delta) - 1) x 100. beta_cum and theta_cum MUST be from the same model.",
            "",
            "",
            "v2.3 fixes critical bug: ΔERR draws now extracted directly from Step 2.1 posterior (same-model baseline). Prior v2.2 incorrectly used base-model cumulative ERR as baseline for interaction-model theta, causing ΔERR inconsistency with Step 2.1 results."),
  stringsAsFactors = FALSE
)
add_sheet("Metadata", meta)
add_sheet("Summary_AllMetrics", summary_all)
add_sheet("PosteriorProbabilities", pvals)
add_sheet("Availability", avail)
if (has_win && exists("win_merged")) add_sheet("WindowERR_Comparison", win_merged)

for (nm in names(comp_results)) {
  sheet_name <- substr(nm, 1, 31)
  add_sheet(sheet_name, comp_results[[nm]])
}

wb_path <- file.path(DIR_TAB,
                     sprintf("MasterSummary_%s_vs_%s.xlsx",
                             REF_LABEL, COMP_LABEL))
saveWorkbook(wb, wb_path, overwrite = TRUE)
cat(sprintf("  Master workbook saved: %s\n", basename(wb_path)))


# ================================================================
#  PART 9: SAVE RDS
# ================================================================

cat("\n=== PART 9: Save RDS ===\n")

saveRDS(list(
  ref_label    = REF_LABEL,
  comp_label   = COMP_LABEL,
  ref_display  = REF_DISPLAY,
  comp_display = COMP_DISPLAY,
  ref_config   = list(HSR = REF_HSR, Analysis = REF_ANALYSIS, Outcome = REF_OUTCOME),
  comp_config  = list(HSR = COMP_HSR, Analysis = COMP_ANALYSIS, Outcome = COMP_OUTCOME),
  comp_results = comp_results,
  pvals        = pvals,
  draws = list(
    ref_cum              = ref_cum,
    comp_cum             = comp_cum,
    ref_adapt            = ref_adapt,
    comp_adapt           = comp_adapt,
    ref_within_hsr       = ref_within_hsr,
    comp_within_hsr      = comp_within_hsr,
    ref_within_hsr_derr  = if (has_within_hsr_derr) ref_within_hsr_derr else NULL,
    comp_within_hsr_derr = if (has_within_hsr_derr) comp_within_hsr_derr else NULL,
    ref_adapt_derr_sd    = if (has_adapt_derr_sd) ref_adapt_derr_sd else NULL,
    comp_adapt_derr_sd   = if (has_adapt_derr_sd) comp_adapt_derr_sd else NULL,
    ref_adapt_derr_hlog  = if (has_adapt_derr_hlog) ref_adapt_derr_hlog else NULL,
    comp_adapt_derr_hlog = if (has_adapt_derr_hlog) comp_adapt_derr_hlog else NULL,
    ref_hsr_pct          = if (has_hsr_mod) ref_hsr_pct else NULL,
    comp_hsr_pct         = if (has_hsr_mod) comp_hsr_pct else NULL,
    ref_hsr_derr         = if (has_hsr_derr) ref_hsr_derr else NULL,
    comp_hsr_derr        = if (has_hsr_derr) comp_hsr_derr else NULL,
    ref_exp_pct          = if (has_exp_mod) ref_exp_pct else NULL,
    comp_exp_pct         = if (has_exp_mod) comp_exp_pct else NULL,
    ref_exp_derr         = if (has_exp_derr) ref_exp_derr else NULL,
    comp_exp_derr        = if (has_exp_derr) comp_exp_derr else NULL
  ),
  availability = avail,
  version      = "2.3",
  created      = Sys.time()
), file.path(DIR_COMP, sprintf("Step2.1S_%s_vs_%s.rds", REF_LABEL, COMP_LABEL)))
cat("  RDS saved\n")


# ================================================================
#  PART 10: AUTO-GENERATED RESULTS TEXT
# ================================================================

cat("\n=== PART 10: Auto-generated results text ===\n")

txt_parts <- list()

# Cumulative ERR
if (has_cum) {
  n_c <- min(length(ref_cum), length(comp_cum))
  d_c <- ref_cum[1:n_c] - comp_cum[1:n_c]
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("The cumulative post-TC ERR was %s%% (95%% CrI: %s, %s) for %s and ",
           "%s%% (%s, %s) for %s. The posterior difference (%s minus %s) was ",
           "%s percentage points (%s, %s); P(%s > %s) = %s."),
    sprintf("%.2f", median(ref_cum)),
    sprintf("%.2f", quantile(ref_cum, .025)),
    sprintf("%.2f", quantile(ref_cum, .975)),
    REF_DISPLAY,
    sprintf("%.2f", median(comp_cum)),
    sprintf("%.2f", quantile(comp_cum, .025)),
    sprintf("%.2f", quantile(comp_cum, .975)),
    COMP_DISPLAY,
    REF_DISPLAY, COMP_DISPLAY,
    sprintf("%.2f", median(d_c)),
    sprintf("%.2f", quantile(d_c, .025)),
    sprintf("%.2f", quantile(d_c, .975)),
    REF_DISPLAY, COMP_DISPLAY,
    fmt_p(mean(d_c > 0)))
}

# HSR modification (both scales)
if (has_hsr_mod) {
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("The overall HSR modification (per +0.1 HSR) was %s%% RR change for %s ",
           "and %s%% for %s."),
    fmt_est(median(ref_hsr_pct), quantile(ref_hsr_pct, .025),
            quantile(ref_hsr_pct, .975)),
    REF_DISPLAY,
    fmt_est(median(comp_hsr_pct), quantile(comp_hsr_pct, .025),
            quantile(comp_hsr_pct, .975)),
    COMP_DISPLAY)
}

if (has_hsr_derr) {
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("Translating to the ERR scale, each +0.1 HSR changed the cumulative ERR by ",
           "%s pp for %s and %s pp for %s."),
    fmt_est(median(ref_hsr_derr), quantile(ref_hsr_derr, .025),
            quantile(ref_hsr_derr, .975)),
    REF_DISPLAY,
    fmt_est(median(comp_hsr_derr), quantile(comp_hsr_derr, .025),
            quantile(comp_hsr_derr, .975)),
    COMP_DISPLAY)
}

# Experience modification
if (has_exp_mod) {
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("The overall experience modification (per +1 H_log) was %s%% RR change for %s ",
           "and %s%% for %s."),
    fmt_est(median(ref_exp_pct), quantile(ref_exp_pct, .025),
            quantile(ref_exp_pct, .975)),
    REF_DISPLAY,
    fmt_est(median(comp_exp_pct), quantile(comp_exp_pct, .025),
            quantile(comp_exp_pct, .975)),
    COMP_DISPLAY)
}

if (has_exp_derr) {
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("On the ERR scale, each +1 H_log changed the cumulative ERR by ",
           "%s pp for %s and %s pp for %s."),
    fmt_est(median(ref_exp_derr), quantile(ref_exp_derr, .025),
            quantile(ref_exp_derr, .975)),
    REF_DISPLAY,
    fmt_est(median(comp_exp_derr), quantile(comp_exp_derr, .025),
            quantile(comp_exp_derr, .975)),
    COMP_DISPLAY)
}

# Adaptation
if (has_adapt_std) {
  n_a <- min(length(ref_adapt$theta_std), length(comp_adapt$theta_std))
  d_a <- ref_adapt$theta_std[1:n_a] - comp_adapt$theta_std[1:n_a]
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("The within-stratum adaptation coefficient (Theta_cum per +1 SD) was ",
           "%s for %s and %s for %s ",
           "(P(%s more adaptive) = %s)."),
    fmt_est(median(ref_adapt$theta_std), quantile(ref_adapt$theta_std, .025),
            quantile(ref_adapt$theta_std, .975), d = 3),
    REF_DISPLAY,
    fmt_est(median(comp_adapt$theta_std), quantile(comp_adapt$theta_std, .025),
            quantile(comp_adapt$theta_std, .975), d = 3),
    COMP_DISPLAY,
    REF_DISPLAY, fmt_p(mean(d_a < 0)))
}

# Adaptation ΔERR (per H_log)
if (has_adapt_derr_hlog) {
  n_a <- min(length(ref_adapt_derr_hlog), length(comp_adapt_derr_hlog))
  d_a <- ref_adapt_derr_hlog[1:n_a] - comp_adapt_derr_hlog[1:n_a]
  txt_parts[[length(txt_parts) + 1]] <- sprintf(
    paste0("The within-stratum ΔERR per +1 H_log was ",
           "%s pp for %s and %s pp for %s; ",
           "difference (Ref minus Comp) = %s pp, ",
           "P(%s more negative) = %s."),
    fmt_est(median(ref_adapt_derr_hlog), quantile(ref_adapt_derr_hlog, .025),
            quantile(ref_adapt_derr_hlog, .975)),
    REF_DISPLAY,
    fmt_est(median(comp_adapt_derr_hlog), quantile(comp_adapt_derr_hlog, .025),
            quantile(comp_adapt_derr_hlog, .975)),
    COMP_DISPLAY,
    fmt_est(median(d_a), quantile(d_a, .025), quantile(d_a, .975)),
    REF_DISPLAY, fmt_p(mean(d_a < 0)))
}

results_text <- paste(txt_parts, collapse = " ")
if (nchar(results_text) > 0) {
  writeLines(results_text, file.path(DIR_TAB, "Results_Text.txt"))
  cat("\n  ", results_text, "\n")
}


# ================================================================
#  FINISH
# ================================================================

cat(sprintf("\n%s\n", paste(rep("=", 70), collapse = "")))
cat("  Step 2.1S (v2.3) complete.\n")
cat(sprintf("  Reference:  %s [%s]\n", REF_LABEL, REF_DISPLAY))
cat(sprintf("  Comparison: %s [%s]\n", COMP_LABEL, COMP_DISPLAY))
cat(sprintf("  Output:     %s\n", DIR_COMP))
cat(sprintf("  Figures: %d  |  Excel tables: %d\n",
            length(list.files(DIR_FIG, pattern = "\\.(pdf|png)$")),
            length(list.files(DIR_TAB, pattern = "\\.xlsx$"))))
cat(sprintf("%s\n", paste(rep("=", 70), collapse = "")))












