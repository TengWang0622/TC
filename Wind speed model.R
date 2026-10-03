#############################################################################################
#
#                          TROPICAL CYCLONES AND HEALTH
#                    Wind Speed Dose-Response Model (City Level, v3.0)
#
#############################################################################################

cat("\n========== STEP 2.2: WIND SPEED DOSE-RESPONSE (v3.0) ==========\n")

DIR_OUT_WS <- file.path("Result TC-health CityLevel", Analysis_type, Outcome_var, "WindSpeed")
DIR_FIG_WS <- file.path(DIR_OUT_WS, "Figures")
DIR_TAB_WS <- file.path(DIR_OUT_WS, "Tables")
for (d in c(DIR_OUT_WS, DIR_FIG_WS, DIR_TAB_WS))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ================================================================
#  2.2 PART 1: 构建窗口特异的风速变量（参考窗口 = -3mo）
# ================================================================

REF_WS <- REF_MAIN   # "-3mo", 与 event model 一致

ws_suffix_all <- c(paste0("m", N_PRE_WINDOWS:1), "act", paste0("p", 1:N_POST_WINDOWS))
names(ws_suffix_all) <- window_labels
nonref_labels <- setdiff(window_labels, REF_WS)   # 除参考期外的 9 个窗口

df_ws <- df %>% filter(!is.na(max_ws_pw), max_ws_pw >= 0)
df_ws$city_id  <- as.integer(factor(df_ws$CityCode))
df_ws$year_idx <- as.integer(factor(df_ws$tc_year))
# window_f 已在 PART 1 relevel 到 REF_MAIN

V_SCALE <- 10   # 三次项用 u = V/10, 避免数值过大

# --- 风速主效应 (事件内常数, 所有行都取 V): 吸收事件间水平混杂 ---
df_ws$ws_main <- df_ws$max_ws_pw
df_ws$v1_main <- df_ws$max_ws_pw / V_SCALE
df_ws$v2_main <- df_ws$v1_main^2
df_ws$v3_main <- df_ws$v1_main^3

# --- 窗口特异风速项 (仅非参考窗口): DID 交互 ---
for (wl in nonref_labels) {
  sfx <- ws_suffix_all[wl]
  ind <- as.character(df_ws$window) == wl
  df_ws[[paste0("ws_", sfx)]] <- ifelse(ind, df_ws$max_ws_pw, 0)
  u <- ifelse(ind, df_ws$max_ws_pw / V_SCALE, 0)
  df_ws[[paste0("v1_", sfx)]] <- u
  df_ws[[paste0("v2_", sfx)]] <- u^2
  df_ws[[paste0("v3_", sfx)]] <- u^3
}
lin_terms <- paste0("ws_", ws_suffix_all[nonref_labels])
cub_terms <- as.vector(t(outer(c("v1_", "v2_", "v3_"),
                               ws_suffix_all[nonref_labels], paste0)))

F_WS_LIN <- as.formula(paste("Y ~ window_f + ws_main +",
                             paste(lin_terms, collapse = " + "), "+", RE_TERMS))
F_WS_CUB <- as.formula(paste("Y ~ window_f + v1_main + v2_main + v3_main +",
                             paste(cub_terms, collapse = " + "), "+", RE_TERMS))

# ================================================================
#  2.2 PART 2: 主模型（线性 DID）
# ================================================================

t0 <- Sys.time()
fit_ws_lin <- inla(F_WS_LIN, family = "nbinomial", data = df_ws, offset = log_offset,
                   control.compute   = list(config = TRUE, waic = TRUE, dic = TRUE),
                   control.predictor = list(link = 1),
                   control.inla      = list(strategy = "simplified.laplace"),
                   control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                   num.threads = "4:1", verbose = FALSE)
cat(sprintf("  Linear DID model time: %.1f min\n", difftime(Sys.time(), t0, units = "mins")))

# --- 提取 delta_w (参考窗口 = 0), per 1 m/s 与 per 5 m/s ---
fs <- fit_ws_lin$summary.fixed
ws_lin_eff <- bind_rows(lapply(window_labels, function(wl) {
  if (wl == REF_WS)
    return(data.frame(window = wl, ERR_1ms = 0, ERR_1ms_lo = 0, ERR_1ms_hi = 0,
                      ERR_5ms = 0, ERR_5ms_lo = 0, ERR_5ms_hi = 0,
                      significant = "ref"))
  i <- which(rownames(fs) == paste0("ws_", ws_suffix_all[wl]))
  if (length(i) != 1) return(NULL)
  data.frame(window = wl,
             ERR_1ms    = (exp(fs$mean[i]) - 1) * 100,
             ERR_1ms_lo = (exp(fs$`0.025quant`[i]) - 1) * 100,
             ERR_1ms_hi = (exp(fs$`0.975quant`[i]) - 1) * 100,
             ERR_5ms    = (exp(5 * fs$mean[i]) - 1) * 100,
             ERR_5ms_lo = (exp(5 * fs$`0.025quant`[i]) - 1) * 100,
             ERR_5ms_hi = (exp(5 * fs$`0.975quant`[i]) - 1) * 100,
             significant = ifelse(fs$`0.025quant`[i] > 0 |
                                    fs$`0.975quant`[i] < 0, "*", "")) }))
ws_lin_eff$window <- factor(ws_lin_eff$window, levels = window_labels)
cat("\n  Per-window delta (per 1 m/s, relative to -3mo of same event).\n")
cat("  Falsification: -2mo & -1mo should now be ~0 (CrI crossing zero):\n")
print(ws_lin_eff)

# eta (风速主效应): 只打印诊断, 不解释
i_eta <- which(rownames(fs) == "ws_main")
cat(sprintf("\n  [Diagnostic] eta (ws_main, absorbs between-event confounding): %.4f [%.4f, %.4f]\n",
            fs$mean[i_eta], fs$`0.025quant`[i_eta], fs$`0.975quant`[i_eta]))

# --- Fig2.2-01: 每窗口 ERR per 1 m/s ---
pd <- ws_lin_eff; pd$x_pos <- as.numeric(pd$window)
p_1ms <- ggplot(pd, aes(x = x_pos, y = ERR_1ms)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ERR_1ms_lo, ymax = ERR_1ms_hi), width = 0.15,
                colour = "#0072B5", linewidth = 0.7) +
  geom_point(colour = "#0072B5", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = expression("ERR per 1 m s"^{-1}*" relative to reference window (%)")) +
  theme_pub()
save_both(p_1ms, file.path(DIR_FIG_WS,
                           sprintf("Fig2.2-01_TC_%s_WS_ERR_per1ms", Outcome_var)), w = 6, h = 4)

# --- Fig2.2-02: 每窗口 ERR per 5 m/s ---
p_5ms <- ggplot(pd, aes(x = x_pos, y = ERR_5ms)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = ERR_5ms_lo, ymax = ERR_5ms_hi), fill = "#0072B5", alpha = 0.18) +
  geom_errorbar(aes(ymin = ERR_5ms_lo, ymax = ERR_5ms_hi), width = 0.15,
                colour = "#0072B5", linewidth = 0.7) +
  geom_line(colour = "#00509B", linewidth = 0.8) +
  geom_point(colour = "#00509B", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = expression("ERR per 5 m s"^{-1}*" relative to reference window (%)")) +
  theme_pub()
save_both(p_5ms, file.path(DIR_FIG_WS,
                           sprintf("Fig2.2-02_TC_%s_WS_ERR_per5ms", Outcome_var)), w = 6, h = 4)
write.xlsx(pd, file.path(DIR_TAB_WS,
                         sprintf("Fig2.2-01-02_TC_%s_WS_PerWindow_PlotData.xlsx", Outcome_var)))

# ================================================================
#  2.2 PART 3: 累积剂量-反应（联合后验抽样）
#  cumERR(V) = exp( sum_post beta_w + V * sum_post delta_w ) - 1
#  同一 seed 抽 beta 与 delta → 保持联合后验相关性
# ================================================================

post_beta_names  <- paste0("window_f", POST_WINDOWS)
post_delta_names <- paste0("ws_", ws_suffix_all[POST_WINDOWS])
M_B <- get_coef_draws(fit_ws_lin, post_beta_names,  N_SAMPLES, 501)
M_D <- get_coef_draws(fit_ws_lin, post_delta_names, N_SAMPLES, 501)
B_draws <- colSums(M_B)   # 累积窗口主效应 (V=0 时的纯时间效应)
S_draws <- colSums(M_D)   # 累积风速斜率 (log 尺度, per 1 m/s)

cat(sprintf("\n  Cumulative post-TC slope: per 1 m/s = %+.3f%% [%+.3f, %+.3f]; per 5 m/s = %+.2f%% [%+.2f, %+.2f]\n",
            median((exp(S_draws) - 1) * 100),
            quantile((exp(S_draws) - 1) * 100, .025),
            quantile((exp(S_draws) - 1) * 100, .975),
            median((exp(5 * S_draws) - 1) * 100),
            quantile((exp(5 * S_draws) - 1) * 100, .025),
            quantile((exp(5 * S_draws) - 1) * 100, .975)))

# --- Fig2.2-03: 累积斜率后验密度（实线 = 后验中位数） ---
sd_draws <- (exp(S_draws) - 1) * 100
p_sd <- ggplot(data.frame(x = sd_draws), aes(x = x)) +
  geom_density(fill = "#0072B5", alpha = 0.35, colour = "#0072B5") +
  geom_vline(xintercept = median(sd_draws), colour = "#0072B5", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(x = expression("Cumulative post-TC ERR per 1 m s"^{-1}*" (%)"),
       y = "Posterior density") + theme_pub()
save_both(p_sd, file.path(DIR_FIG_WS,
                          sprintf("Fig2.2-03_TC_%s_WS_CumSlope_Posterior", Outcome_var)), w = 5, h = 4)

save_both(p_sd, file.path(DIR_FIG_WS,
                          sprintf("Fig2.2-03_TC_%s_WS_CumSlope_Posterior_Zoom", Outcome_var)), w = 3.1, h = 2.5)

# --- 剂量-反应曲线网格 (截断到观测风速 99 分位数, 避免外推区) ---
wind_grid <- seq(0, quantile(df_ws$max_ws_pw, 0.998, na.rm = TRUE), by = 0.5)

curve_lin <- bind_rows(lapply(wind_grid, function(v) {
  err <- (exp(B_draws + v * S_draws) - 1) * 100
  data.frame(wind = v, model = "Linear", ERR = median(err),
             lo = quantile(err, .025, names = FALSE),
             hi = quantile(err, .975, names = FALSE)) }))

# ================================================================
#  2.2 PART 4: 敏感性 1 — 三次多项式 DID 模型
# ================================================================

t0 <- Sys.time()
fit_ws_cub <- inla(F_WS_CUB, family = "nbinomial", data = df_ws, offset = log_offset,
                   control.compute   = list(config = TRUE, waic = TRUE),
                   control.predictor = list(link = 1),
                   control.inla      = list(strategy = "simplified.laplace"),
                   control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                   num.threads = "4:1", verbose = FALSE)
cat(sprintf("  Cubic DID model time: %.1f min\n", difftime(Sys.time(), t0, units = "mins")))

post_sfx <- ws_suffix_all[POST_WINDOWS]
c1 <- paste0("v1_", post_sfx); c2 <- paste0("v2_", post_sfx); c3 <- paste0("v3_", post_sfx)
M_Bc  <- get_coef_draws(fit_ws_cub, paste0("window_f", POST_WINDOWS), N_SAMPLES, 502)
M_cub <- get_coef_draws(fit_ws_cub, c(c1, c2, c3), N_SAMPLES, 502)
Bc <- colSums(M_Bc)
S1 <- colSums(M_cub[c1, , drop = FALSE])
S2 <- colSums(M_cub[c2, , drop = FALSE])
S3 <- colSums(M_cub[c3, , drop = FALSE])

curve_cub <- bind_rows(lapply(wind_grid, function(v) {
  u <- v / V_SCALE
  err <- (exp(Bc + S1*u + S2*u^2 + S3*u^3) - 1) * 100
  data.frame(wind = v, model = "Cubic", ERR = median(err),
             lo = quantile(err, .025, names = FALSE),
             hi = quantile(err, .975, names = FALSE)) }))

# --- Fig2.2-04: 剂量-反应曲线叠加（线性 vs 三次） ---
#  注: 指数链接 + 高风速处数据稀疏 → cubic 的后验在高 V 处右偏,
#      中位数天然贴近下界, 属正常现象, caption 中说明即可。
curve_all <- bind_rows(curve_lin, curve_cub)
p_dr <- ggplot(curve_all, aes(x = wind, y = ERR, colour = model, fill = model)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = 17.5, linetype = "dashed", colour = "grey50") +
  #annotate("text", x = 17.5, y = max(curve_all$hi) * 0.95, label = "TS threshold = 17.5",
  #         hjust = -0.1, size = 3, colour = "black") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = df_ws %>% distinct(event_id, max_ws_pw),
           aes(x = max_ws_pw), inherit.aes = FALSE, alpha = 0.1,
           length = unit(0.02, "npc")) +
  scale_colour_manual(values = c(Linear = "#0072B5", Cubic = "#BC3C29"), name = NULL) +
  scale_fill_manual(values = c(Linear = "#0072B5", Cubic = "#BC3C29"), name = NULL) +
  labs(x = expression("Maximum sustained wind speed (m s"^{-1}*")"),
       y = "Cumulative post-TC ERR (%)") +
  theme_pub() + theme(legend.position = c(0.15, 0.85))
save_both(p_dr, file.path(DIR_FIG_WS,
                          sprintf("Fig2.2-04_TC_%s_WS_DoseResponse_LinVsCubic", Outcome_var)), w = 6, h = 4.5)
write.xlsx(curve_all, file.path(DIR_TAB_WS,
                                sprintf("Fig2.2-04_TC_%s_WS_DoseResponse_PlotData.xlsx", Outcome_var)))

# ================================================================
#  2.2 PART 5: 敏感性 2 — 线性 DID, 剔除 2020-2022
# ================================================================

df_ws_noP <- df_ws %>% filter(!(tc_year %in% 2020:2022))
df_ws_noP$city_id  <- as.integer(factor(df_ws_noP$CityCode))
df_ws_noP$year_idx <- as.integer(factor(df_ws_noP$tc_year))
fit_ws_noP <- inla(F_WS_LIN, family = "nbinomial", data = df_ws_noP,
                   offset = log_offset,
                   control.compute   = list(config = TRUE, waic = TRUE),
                   control.predictor = list(link = 1),
                   control.inla      = list(strategy = "simplified.laplace"),
                   control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                   num.threads = "4:1", verbose = FALSE)
fs2 <- fit_ws_noP$summary.fixed
ws_noP_eff <- bind_rows(lapply(nonref_labels, function(wl) {
  i <- which(rownames(fs2) == paste0("ws_", ws_suffix_all[wl]))
  if (length(i) != 1) return(NULL)
  data.frame(window = wl,
             ERR_1ms    = (exp(fs2$mean[i]) - 1) * 100,
             ERR_1ms_lo = (exp(fs2$`0.025quant`[i]) - 1) * 100,
             ERR_1ms_hi = (exp(fs2$`0.975quant`[i]) - 1) * 100) }))
cat("\n  Excl. 2020-2022 (per 1 m/s, relative to reference window):\n")
print(ws_noP_eff)

# ================================================================
#  2.2 PART 6: 诊断 — 旧版无参考模型（展示事件间混杂）
#  预期: pre 期斜率整体下移且显著为负 = 事件层面混杂的证据,
#  与主模型对比可写入 Supplementary
# ================================================================

for (wl in window_labels) {   # 旧版需要全部 10 个窗口的 ws 项 (含参考)
  sfx <- ws_suffix_all[wl]
  ind <- as.character(df_ws$window) == wl
  df_ws[[paste0("wsN_", sfx)]] <- ifelse(ind, df_ws$max_ws_pw, 0)
}
naive_terms <- paste0("wsN_", ws_suffix_all)
F_WS_NAIVE <- as.formula(paste("Y ~", paste(naive_terms, collapse = " + "),
                               "+", RE_TERMS))
fit_ws_naive <- inla(F_WS_NAIVE, family = "nbinomial", data = df_ws,
                     offset = log_offset,
                     control.compute   = list(config = TRUE, waic = TRUE),
                     control.predictor = list(link = 1),
                     control.inla      = list(strategy = "simplified.laplace"),
                     control.fixed     = PRIOR_FIXED, control.family = CTRL_FAMILY,
                     num.threads = "4:1", verbose = FALSE)
fsN <- fit_ws_naive$summary.fixed
ws_naive_eff <- bind_rows(lapply(window_labels, function(wl) {
  i <- which(rownames(fsN) == paste0("wsN_", ws_suffix_all[wl]))
  if (length(i) != 1) return(NULL)
  data.frame(window = wl,
             ERR_1ms    = (exp(fsN$mean[i]) - 1) * 100,
             ERR_1ms_lo = (exp(fsN$`0.025quant`[i]) - 1) * 100,
             ERR_1ms_hi = (exp(fsN$`0.975quant`[i]) - 1) * 100) }))
ws_naive_eff$window <- factor(ws_naive_eff$window, levels = window_labels)
cat("\n  [Diagnostic] Naive model (no within-event anchoring);\n")
cat("  negative pre-window slopes indicate between-event confounding:\n")
print(ws_naive_eff)

pdN <- ws_naive_eff; pdN$x_pos <- as.numeric(pdN$window)
p_naive <- ggplot(pdN, aes(x = x_pos, y = ERR_1ms)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ERR_1ms_lo, ymax = ERR_1ms_hi), width = 0.15,
                colour = "grey40", linewidth = 0.7) +
  geom_point(colour = "grey30", size = 2.5) +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = expression("ERR per 1 m s"^{-1}*" vs V = 0 baseline (%)")) +
  theme_pub()
save_both(p_naive, file.path(DIR_FIG_WS,
                             sprintf("Fig2.2-05_TC_%s_WS_NaiveModel_Diagnostic", Outcome_var)), w = 6, h = 4)
write.xlsx(ws_naive_eff, file.path(DIR_TAB_WS,
                                   sprintf("Fig2.2-05_TC_%s_WS_NaiveModel_PlotData.xlsx", Outcome_var)))

# ================================================================
#  2.2 PART 7: 保存
# ================================================================

wb2 <- createWorkbook()
addWorksheet(wb2, "Linear_PerWindow");  writeData(wb2, "Linear_PerWindow", ws_lin_eff)
addWorksheet(wb2, "Linear_Excl20_22");  writeData(wb2, "Linear_Excl20_22", ws_noP_eff)
addWorksheet(wb2, "Naive_Diagnostic");  writeData(wb2, "Naive_Diagnostic", ws_naive_eff)
addWorksheet(wb2, "DoseResponse");      writeData(wb2, "DoseResponse", curve_all)
addWorksheet(wb2, "WAIC")
writeData(wb2, "WAIC", data.frame(
  Model = c("Linear DID (main)", "Cubic DID", "Linear DID excl.2020-2022",
            "Naive (diagnostic)"),
  WAIC  = c(get_waic(fit_ws_lin), get_waic(fit_ws_cub),
            get_waic(fit_ws_noP), get_waic(fit_ws_naive))))
saveWorkbook(wb2, file.path(DIR_TAB_WS,
                            sprintf("Tab2.2-00_TC_%s_Step2.2_WindSpeed_Results.xlsx", Outcome_var)),
             overwrite = TRUE)
saveRDS(list(linear = ws_lin_eff, excl2020_2022 = ws_noP_eff,
             naive_diagnostic = ws_naive_eff, curves = curve_all,
             cum_slope_draws = S_draws, cum_beta_draws = B_draws,
             created = Sys.time()),
        file.path(DIR_OUT_WS, sprintf("TC_%s_Step2.2_WindSpeed.rds", Outcome_var)))
cat("  Step 2.2 complete.\n")













# ================================================================
#  2.2 PART 7b: 代表性风速轨迹 + 与 Step 2.1 的一致性检查
#  模型中 ws_main 与 ws_w 均未中心化 → V 的原点即 0 m/s, V_CENTER = 0
#  ERR(w | V) = exp(beta_w + delta_w * V) - 1  (相对同一事件的参考窗口)
# ================================================================

if (!exists("active_x")) active_x <- which(window_labels == "Active")
if (!exists("fmt_p")) {   # 若未运行 2.1 块 A
  fmt_p <- function(p) ifelse(p >= 0.999, ">0.999",
                              ifelse(p <= 0.001, "<0.001", sprintf("%.3f", p)))
  ptxt <- function(x, dir = c("gt0", "lt0")) {
    dir <- match.arg(dir)
    if (is.null(x) || length(x) == 0) return(NA_character_)
    p <- if (dir == "gt0") mean(x > 0) else mean(x < 0)
    sprintf("P(%s0) = %s", ifelse(dir == "gt0", ">", "<"), fmt_p(p))
  }
}

V_CENTER <- 0

# 代表性风速: 事件级 max_ws_pw 的 P25/P50/P75/P95
ev_ws  <- df_ws %>% distinct(event_id, max_ws_pw)
V_REPS <- data.frame(level = c("P25", "P50", "P75", "P95"),
                     V = as.numeric(quantile(ev_ws$max_ws_pw,
                                             c(.25, .50, .75, .95), na.rm = TRUE)))
V_MEAN <- mean(ev_ws$max_ws_pw, na.rm = TRUE)
cat(sprintf("  Representative winds: %s | mean = %.1f m/s\n",
            paste(sprintf("%s=%.1f", V_REPS$level, V_REPS$V), collapse = ", "), V_MEAN))

# 全部非参考窗口 beta 与 delta 的联合后验抽样 (同一 seed 保持相关性)
all_beta_names  <- paste0("window_f", nonref_labels)
all_delta_names <- paste0("ws_", ws_suffix_all[nonref_labels])
M_B_all <- get_coef_draws(fit_ws_lin, all_beta_names,  N_SAMPLES, 505)
M_D_all <- get_coef_draws(fit_ws_lin, all_delta_names, N_SAMPLES, 505)

# 各代表性风速下的窗口轨迹
traj_df <- bind_rows(lapply(seq_len(nrow(V_REPS)), function(j) {
  v <- V_REPS$V[j]
  bind_rows(lapply(nonref_labels, function(wl) {
    dr  <- M_B_all[paste0("window_f", wl), ] +
      v * M_D_all[paste0("ws_", ws_suffix_all[wl]), ]
    err <- (exp(dr) - 1) * 100
    data.frame(wind_level = V_REPS$level[j], wind = v, window = wl,
               ERR = median(err), lo = quantile(err, .025, names = FALSE),
               hi = quantile(err, .975, names = FALSE))
  }))
}))
traj_df <- bind_rows(traj_df,
                     data.frame(wind_level = V_REPS$level, wind = V_REPS$V,
                                window = REF_WS, ERR = 0, lo = 0, hi = 0))
traj_df$window     <- factor(traj_df$window, levels = window_labels)
traj_df$wind_level <- factor(traj_df$wind_level, levels = V_REPS$level)

# --- Fig2.2-06: 代表性风速下的 event-study 轨迹 ---
pdT <- traj_df; pdT$x_pos <- as.numeric(pdT$window)
p_traj <- ggplot(pdT, aes(x = x_pos, y = ERR, colour = wind_level)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = active_x, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                position = position_dodge(0.5), linewidth = 0.5) +
  geom_line(position = position_dodge(0.5), linewidth = 0.7) +
  geom_point(position = position_dodge(0.5), size = 2) +
  scale_colour_manual(values = c(P25 = "#0072B5", P50 = "#20854E",
                                 P75 = "#E18727", P95 = "#BC3C29"),
                      labels = sprintf("%s (%.1f m/s)", V_REPS$level, V_REPS$V),
                      name = "Wind speed") +
  scale_x_continuous(breaks = seq_along(window_labels), labels = window_display) +
  labs(x = "Lag period (month)",
       y = "ERR relative to reference window (%)") + theme_pub()
save_both(p_traj, file.path(DIR_FIG_WS,
                            sprintf("Fig2.2-06_TC_%s_WS_Trajectories_RepWinds", Outcome_var)), w = 8, h = 5)
write.xlsx(traj_df, file.path(DIR_TAB_WS,
                              sprintf("Fig2.2-06_TC_%s_WS_Trajectories_PlotData.xlsx", Outcome_var)))

# 一致性检查: 平均风速处各窗口 ERR (应与 Step 2.1 主模型的 res_main 量级接近)
beta_check <- bind_rows(lapply(nonref_labels, function(wl) {
  dr  <- M_B_all[paste0("window_f", wl), ] +
    V_MEAN * M_D_all[paste0("ws_", ws_suffix_all[wl]), ]
  err <- (exp(dr) - 1) * 100
  data.frame(window = wl,
             ERR_at_meanV    = median(err),
             ERR_at_meanV_lo = quantile(err, .025, names = FALSE),
             ERR_at_meanV_hi = quantile(err, .975, names = FALSE))
}))
cat("\n  Consistency check (ERR at mean wind, cf. Step 2.1 res_main):\n")
print(beta_check)

# ================================================================
#  2.2 PART 8: 汇总输出表 + 后验概率表
# ================================================================

cat("\n========== 2.2 PART 8: SUMMARY TABLE + POSTERIOR PROBABILITIES ==========\n")

# ---- helper (与 Step 2.1 PART 9 同版本, 此处无条件重定义以保证一致) ----
SUMMARY_COLS <- c("Analysis", "Fig", "-3", "-2", "-1", "Active",
                  as.character(1:6), "Cumulative", "P_Cumulative")
fmt_est <- function(est, lo, hi, digits = 2)
  sprintf(paste0("%.", digits, "f%% (%.", digits, "f, %.", digits, "f)"),
          est, lo, hi)
WIN_COL_MAP <- setNames(c("-3", "-2", "-1", "Active", as.character(1:6)),
                        window_labels)
make_summary_row <- function(analysis, fig, df_win = NULL,
                             win_col = "window", est_col = "ERR_pct",
                             lo_col = "lower", hi_col = "upper",
                             ref_label = NULL, cum_triple = NULL,
                             p_text = NULL, digits = 2) {
  row <- as.data.frame(as.list(setNames(rep(NA_character_, length(SUMMARY_COLS)),
                                        SUMMARY_COLS)),
                       check.names = FALSE, stringsAsFactors = FALSE)
  row[["Analysis"]] <- analysis
  row[["Fig"]]      <- fig
  if (!is.null(df_win)) {
    need <- c(win_col, est_col, lo_col, hi_col)
    miss <- setdiff(need, names(df_win))
    if (length(miss) > 0)
      stop(sprintf("[Summary] '%s' 缺少列: %s | 实际列名: %s",
                   analysis, paste(miss, collapse = ", "),
                   paste(names(df_win), collapse = ", ")), call. = FALSE)
    for (k in seq_len(nrow(df_win))) {
      wl <- as.character(df_win[[win_col]][k])
      if (!wl %in% names(WIN_COL_MAP)) next
      if (!is.null(ref_label) && wl == ref_label) {
        row[[WIN_COL_MAP[[wl]]]] <- "0 (ref)"
      } else {
        row[[WIN_COL_MAP[[wl]]]] <- fmt_est(df_win[[est_col]][k],
                                            df_win[[lo_col]][k],
                                            df_win[[hi_col]][k], digits)
      }
    }
    if (!is.null(ref_label) && is.na(row[[WIN_COL_MAP[[ref_label]]]]))
      row[[WIN_COL_MAP[[ref_label]]]] <- "0 (ref)"
  }
  if (!is.null(cum_triple) && length(cum_triple) == 3 && all(is.finite(cum_triple)))
    row[["Cumulative"]] <- fmt_est(cum_triple[1], cum_triple[2],
                                   cum_triple[3], digits)
  if (!is.null(p_text) && !is.na(p_text))
    row[["P_Cumulative"]] <- p_text
  row
}
# log 尺度后验样本 → ERR% 三元组
cum_from_log_draws <- function(x) {
  err <- (exp(x) - 1) * 100
  c(median(err), quantile(err, .025, names = FALSE),
    quantile(err, .975, names = FALSE))
}

rows_22 <- list()
add_row22 <- function(...) {
  r <- tryCatch(make_summary_row(...),
                error = function(e) { cat(" ", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) rows_22[[length(rows_22) + 1]] <<- r
}

## 1) 剂量梯度 per 1 m/s
add_row22(analysis = "WS_Gradient_per1ms", fig = "Fig2.2-01",
          df_win = ws_lin_eff, est_col = "ERR_1ms",
          lo_col = "ERR_1ms_lo", hi_col = "ERR_1ms_hi",
          ref_label = REF_WS,
          cum_triple = cum_from_log_draws(S_draws),
          p_text = ptxt(S_draws, "gt0"))

## 2) 剂量梯度 per 5 m/s
add_row22(analysis = "WS_Gradient_per5ms", fig = "Fig2.2-02",
          df_win = ws_lin_eff, est_col = "ERR_5ms",
          lo_col = "ERR_5ms_lo", hi_col = "ERR_5ms_hi",
          ref_label = REF_WS,
          cum_triple = cum_from_log_draws(5 * S_draws),
          p_text = ptxt(5 * S_draws, "gt0"))

## 3) 平均风速处的窗口 ERR (与 Step 2.1 一致性检查)
cum_meanV <- B_draws + V_MEAN * S_draws
add_row22(analysis = sprintf("ERR_at_MeanWind_%.1fms_ConsistencyCheck", V_MEAN),
          fig = "Table only",
          df_win = beta_check, est_col = "ERR_at_meanV",
          lo_col = "ERR_at_meanV_lo", hi_col = "ERR_at_meanV_hi",
          ref_label = REF_WS,
          cum_triple = cum_from_log_draws(cum_meanV),
          p_text = ptxt(cum_meanV, "gt0"))

## 4) 各代表性风速下的窗口轨迹 + 该风速下的总累积
for (j in seq_len(nrow(V_REPS))) {
  lv <- as.character(V_REPS$level[j]); v <- V_REPS$V[j]
  cum_v <- B_draws + (v - V_CENTER) * S_draws
  add_row22(analysis = sprintf("WS_ERR_at_%s_%.1fms", lv, v),
            fig = "Fig2.2-06",
            df_win = traj_df[traj_df$wind_level == lv, ],
            est_col = "ERR", lo_col = "lo", hi_col = "hi",
            ref_label = REF_WS,
            cum_triple = cum_from_log_draws(cum_v),
            p_text = ptxt(cum_v, "gt0"))
}

## 5) 敏感性: 剔除 2020-2022 (per 1 m/s)
S_noP_draws <- colSums(get_coef_draws(fit_ws_noP, post_delta_names,
                                      N_SAMPLES, 601))
add_row22(analysis = "Sen_PandemicExcluded_per1ms", fig = "Table only",
          df_win = ws_noP_eff, est_col = "ERR_1ms",
          lo_col = "ERR_1ms_lo", hi_col = "ERR_1ms_hi",
          ref_label = REF_WS,
          cum_triple = cum_from_log_draws(S_noP_draws),
          p_text = ptxt(S_noP_draws, "gt0"))

## 6) 诊断: 旧版无参考模型 (系数非事件内对比, 不给 Cumulative)
add_row22(analysis = "Diag_NaiveModel_per1ms", fig = "Fig2.2-05",
          df_win = ws_naive_eff, est_col = "ERR_1ms",
          lo_col = "ERR_1ms_lo", hi_col = "ERR_1ms_hi",
          ref_label = NULL)

## 7) 敏感性: 三次多项式 — 报 P50/P75/P95 风速下的总累积
for (qq in c("P50", "P75", "P95")) {
  j <- which(as.character(V_REPS$level) == qq)
  if (length(j) != 1) next
  v  <- V_REPS$V[j]
  uc <- (v - V_CENTER) / V_SCALE
  cub_draws <- Bc + S1 * uc + S2 * uc^2 + S3 * uc^3
  add_row22(analysis = sprintf("Sen_Cubic_CumERR_at_%s_%.1fms", qq, v),
            fig = "Fig2.2-04",
            cum_triple = cum_from_log_draws(cub_draws),
            p_text = ptxt(cub_draws, "gt0"))
}

## ---------- 汇总并输出 ----------
summary_22 <- do.call(rbind, rows_22)
print(summary_22, right = FALSE)
write.xlsx(summary_22,
           file.path(DIR_TAB_WS,
                     sprintf("Tab2.2-Summary_TC_%s_AllAnalyses.xlsx", Outcome_var)),
           overwrite = TRUE)

# ================================================================
#  2.2 PART 8b: 后验概率总表 (pvals_22)
# ================================================================

prow <- function(Analysis, Quantity, Statement, P, Interpretation)
  data.frame(Analysis = Analysis, Quantity = Quantity, Statement = Statement,
             Probability = round(P, 4), Interpretation = Interpretation,
             stringsAsFactors = FALSE)

pv2 <- list()

# Falsification: pre 期 delta (来自边际后验, 精确值; 期望 ~0.5)
for (wl in setdiff(window_labels[1:N_PRE_WINDOWS], REF_WS)) {
  m <- fit_ws_lin$marginals.fixed[[paste0("ws_", ws_suffix_all[wl])]]
  if (!is.null(m))
    pv2[[length(pv2) + 1]] <- prow("Main_LinearDID",
                                   paste0("delta, window ", wl), "P(delta > 0)", 1 - inla.pmarginal(0, m),
                                   "Falsification test: should be near 0.5 (no pre-TC wind gradient)")
}
# Post 期各窗口 delta
for (wl in POST_WINDOWS) {
  m <- fit_ws_lin$marginals.fixed[[paste0("ws_", ws_suffix_all[wl])]]
  if (!is.null(m))
    pv2[[length(pv2) + 1]] <- prow("Main_LinearDID",
                                   paste0("delta, window ", wl), "P(delta > 0)", 1 - inla.pmarginal(0, m),
                                   "Posterior probability of a positive wind-speed gradient in this window")
}
# 累积斜率与代表性风速处累积 ERR
pv2[[length(pv2) + 1]] <- prow("Main_LinearDID", "Cumulative slope (per 1 m/s)",
                               "P(> 0)", mean(S_draws > 0),
                               "Posterior probability of a positive cumulative dose-response")
for (j in seq_len(nrow(V_REPS))) {
  v <- V_REPS$V[j]
  pv2[[length(pv2) + 1]] <- prow("Main_LinearDID",
                                 sprintf("Cumulative ERR at %s (%.1f m/s)", V_REPS$level[j], v),
                                 "P(> 0)", mean((B_draws + (v - V_CENTER) * S_draws) > 0),
                                 "Posterior probability of net post-TC increase at this wind speed")
}
# Cubic
for (qq in c("P50", "P75", "P95")) {
  j <- which(as.character(V_REPS$level) == qq)
  if (length(j) != 1) next
  uc <- (V_REPS$V[j] - V_CENTER) / V_SCALE
  pv2[[length(pv2) + 1]] <- prow("Sen_Cubic",
                                 sprintf("Cumulative ERR at %s", qq), "P(> 0)",
                                 mean((Bc + S1 * uc + S2 * uc^2 + S3 * uc^3) > 0),
                                 "Robustness of dose-response to nonlinear specification")
}
# Sens noP
pv2[[length(pv2) + 1]] <- prow("Sen_PandemicExcluded",
                               "Cumulative slope (per 1 m/s)", "P(> 0)", mean(S_noP_draws > 0),
                               "Robustness excluding 2020-2022")
# Naive 诊断: pre 期斜率显著为负 = 事件间混杂证据
for (wl in window_labels[1:N_PRE_WINDOWS]) {
  m <- fit_ws_naive$marginals.fixed[[paste0("wsN_", ws_suffix_all[wl])]]
  if (!is.null(m))
    pv2[[length(pv2) + 1]] <- prow("Diag_Naive",
                                   paste0("slope, window ", wl), "P(< 0)", inla.pmarginal(0, m),
                                   "High value indicates between-event confounding in the naive model")
}

pvals_22 <- bind_rows(pv2)
print(pvals_22, right = FALSE)
write.xlsx(pvals_22,
           file.path(DIR_TAB_WS,
                     sprintf("Tab2.2-P_TC_%s_PosteriorProbabilities.xlsx", Outcome_var)),
           overwrite = TRUE)
cat(sprintf("  Step 2.2 summary: %d rows | probability table: %d rows\n",
            nrow(summary_22), nrow(pvals_22)))


