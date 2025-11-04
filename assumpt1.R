# ===========================
# mediation_assumption_check.R
# Checks key assumptions for linear mediation (continuous X, M, Y)
# Saves CSV + TXT report + diagnostic plots to /Users/jamesjackson/Desktop/dir
# ===========================

# --------- 0) Packages ----------
required_packages <- c(
  "dplyr","readr","tidyr","stringr","purrr","broom","ggplot2",
  "car","lmtest"
)
missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# --------- 1) User config ----------
in_csv  <- "/Users/jamesjackson/Desktop/dir/working.csv"
out_dir <- "/Users/jamesjackson/Desktop/dir"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# >>> EDIT THESE NAMES AS NEEDED <<<
X_var <- "t_sacs"           # predictor
M_var <- "atq_f_total"      # mediator
Y_var <- "dass_dep"         # outcome
covariates <- c()           # e.g., c("Age")

# Filenames
out_summary_csv <- file.path(out_dir, "med_assumption_summary.csv")
out_report_txt  <- file.path(out_dir, "med_assumption_report.txt")

plot_histqq_png <- function(name) file.path(out_dir, paste0("assump_", name, "_hist_qq.png"))
plot_scatter_png <- function(name) file.path(out_dir, paste0("assump_scatter_", name, ".png"))
plot_resid_png  <- function(name) file.path(out_dir, paste0("assump_resid_", name, ".png"))
plot_cook_png   <- function(name) file.path(out_dir, paste0("assump_cooks_", name, ".png"))

# --------- 2) Load data ----------
cat("Reading:", in_csv, "\n")
df <- readr::read_csv(in_csv, show_col_types = FALSE)

all_needed <- unique(c(X_var, M_var, Y_var, covariates))
missing_cols <- setdiff(all_needed, names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in working.csv: ", paste(missing_cols, collapse = ", "))
}

# Keep only variables of interest
d <- df %>% dplyr::select(dplyr::all_of(all_needed))

# Ensure X, M, Y are numeric
to_num <- function(x) { if (is.numeric(x)) x else suppressWarnings(readr::parse_number(as.character(x))) }
d <- d %>%
  dplyr::mutate(
    !!X_var := to_num(.data[[X_var]]),
    !!M_var := to_num(.data[[M_var]]),
    !!Y_var := to_num(.data[[Y_var]])
  )

# Coerce numeric covariates where possible; leave others as-is (they’ll be ignored in VIF if non-numeric factors)
if (length(covariates) > 0) {
  for (cv in covariates) {
    if (!is.numeric(d[[cv]]) && !is.factor(d[[cv]])) {
      d[[cv]] <- suppressWarnings(readr::parse_number(as.character(d[[cv]])))
    }
  }
}

# --------- 3) Missingness / usable Ns ----------
pair_n <- function(a,b) sum(stats::complete.cases(d[,c(a,b)]))
n_xm <- pair_n(X_var, M_var)
n_xy <- pair_n(X_var, Y_var)
n_my <- pair_n(M_var, Y_var)
n_all_a <- sum(stats::complete.cases(d[, c(X_var, M_var, covariates)]))
n_all_b <- sum(stats::complete.cases(d[, c(X_var, M_var, Y_var, covariates)]))

# --------- 4) Bivariate links & significance (Pearson) ----------
pearson_test <- function(a,b) {
  tmp <- d[, c(a,b)]
  tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]
  if (nrow(tmp) < 3) return(list(r=NA_real_, p=NA_real_, n=nrow(tmp)))
  res <- suppressWarnings(cor.test(tmp[[1]], tmp[[2]], method = "pearson"))
  list(r=unname(res$estimate), p=res$p.value, n=nrow(tmp))
}
res_xm <- pearson_test(X_var, M_var)
res_xy <- pearson_test(X_var, Y_var)
res_my <- pearson_test(M_var, Y_var)

# --------- 5) Fit mediation component models ----------
# a-path: M ~ X + C
fmla_a <- as.formula(paste(M_var, "~", paste(c(X_var, covariates), collapse = " + ")))
# b/c'-paths: Y ~ X + M + C
fmla_b <- as.formula(paste(Y_var, "~", paste(c(X_var, M_var, covariates), collapse = " + ")))
# c-path (optional): Y ~ X + C
fmla_c <- as.formula(paste(Y_var, "~", paste(c(X_var, covariates), collapse = " + ")))

d_a <- d %>% dplyr::select(all_of(c(M_var, X_var, covariates))) %>% tidyr::drop_na()
d_b <- d %>% dplyr::select(all_of(c(Y_var, X_var, M_var, covariates))) %>% tidyr::drop_na()
d_c <- d %>% dplyr::select(all_of(c(Y_var, X_var, covariates))) %>% tidyr::drop_na()

fit_a <- lm(fmla_a, data = d_a)
fit_b <- lm(fmla_b, data = d_b)
fit_c <- lm(fmla_c, data = d_c)

k_a <- length(coef(fit_a))                   # parameters incl. intercept
k_b <- length(coef(fit_b))
n_a <- nrow(d_a); n_b <- nrow(d_b)

# --------- 6) Assumption diagnostics ----------
# 6.1 Linearity visuals: scatter with loess
make_scatter <- function(x, y, name) {
  dd <- d %>% dplyr::select(all_of(c(x,y))) %>% tidyr::drop_na()
  if (nrow(dd) < 3) return(invisible(NULL))
  p <- ggplot(dd, aes_string(x=x, y=y)) +
    geom_point(alpha=.7) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = paste0("Scatter: ", y, " ~ ", x)) +
    theme_minimal()
  ggplot2::ggsave(plot_scatter_png(name), p, width = 6, height = 5, dpi = 300)
}
make_scatter(X_var, M_var, "M_on_X")
make_scatter(X_var, Y_var, "Y_on_X")
make_scatter(M_var, Y_var, "Y_on_M")

# 6.2 Hist + QQ for X, M, Y
make_histqq <- function(v) {
  vec <- d[[v]]
  vec <- vec[!is.na(vec)]
  if (length(vec) < 3) return(invisible(NULL))
  dfp <- data.frame(val = vec)
  p1 <- ggplot(dfp, aes(val)) + geom_histogram(bins = 30) + theme_minimal() + labs(title=paste("Histogram:", v))
  p2 <- ggplot(dfp, aes(sample = val)) + stat_qq() + stat_qq_line() + theme_minimal() + labs(title=paste("QQ plot:", v))
  # arrange side-by-side
  g <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1,1))
  ggplot2::ggsave(plot_histqq_png(v), g, width = 9, height = 4.5, dpi = 300)
}
# cowplot for arranging
if (!"cowplot" %in% installed.packages()[,"Package"]) install.packages("cowplot", dependencies = TRUE)
library(cowplot)
purrr::walk(c(X_var, M_var, Y_var), make_histqq)

# 6.3 Residual diagnostics (normality & homoscedasticity) for fit_a, fit_b
diag_model <- function(fit, name) {
  res <- rstandard(fit)
  fitvals <- fitted(fit)
  # Resid vs Fitted
  p <- ggplot(data.frame(fit=fitvals, res=res), aes(fit, res)) +
    geom_point(alpha=.7) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste("Residuals vs Fitted:", name), x="Fitted", y="Std. residual") +
    theme_minimal()
  ggplot2::ggsave(plot_resid_png(name), p, width = 6, height = 5, dpi = 300)
  
  # Breusch–Pagan for heteroskedasticity
  bp <- tryCatch(lmtest::bptest(fit), error = function(e) NULL)
  
  # Shapiro–Wilk on residuals (n<=5000)
  sw <- if (length(res) >= 3 && length(res) <= 5000) tryCatch(shapiro.test(res), error=function(e) NULL) else NULL
  
  list(bp_p = if (!is.null(bp)) bp$p.value else NA_real_,
       sw_p = if (!is.null(sw)) sw$p.value else NA_real_)
}
d_a_diag <- diag_model(fit_a, "M_on_X")
d_b_diag <- diag_model(fit_b, "Y_on_XM")

# 6.4 Multicollinearity (VIF) for outcome model
vif_vals <- tryCatch(car::vif(fit_b), error = function(e) NA)
max_vif <- if (all(is.na(vif_vals))) NA_real_ else max(vif_vals, na.rm = TRUE)

# 6.5 Influence & outliers (fit_b)
cooks <- cooks.distance(fit_b)
lev   <- hatvalues(fit_b)
stud  <- rstudent(fit_b)
thr_cook <- 4 / n_b
thr_lev  <- 2 * k_b / n_b
thr_stud <- 3

flags <- (cooks > thr_cook) | (lev > thr_lev) | (abs(stud) > thr_stud)
n_flag <- sum(flags, na.rm = TRUE)

# Plot Cook's distances
cook_df <- data.frame(i = seq_along(cooks), cook = cooks)
p_cook <- ggplot(cook_df, aes(i, cook)) +
  geom_col() +
  geom_hline(yintercept = thr_cook, color = "red") +
  labs(title = "Cook's Distance (Outcome Model)", y = "Cook's D", x = "Observation") +
  theme_minimal()
ggplot2::ggsave(plot_cook_png("Y_on_XM"), p_cook, width = 7, height = 4.5, dpi = 300)

# --------- 7) Summaries & flags ----------
# Simple heuristics:
flag_nonlin <- NA  # visual; user checks scatter & residual smooth
flag_hetero_a <- !is.na(d_a_diag$bp_p) && d_a_diag$bp_p < .05
flag_hetero_b <- !is.na(d_b_diag$bp_p) && d_b_diag$bp_p < .05
flag_nonNorm_a <- !is.na(d_a_diag$sw_p) && d_a_diag$sw_p < .05
flag_nonNorm_b <- !is.na(d_b_diag$sw_p) && d_b_diag$sw_p < .05
flag_vif <- !is.na(max_vif) && max_vif >= 5  # >5 often concerning; >10 severe

coefs_a <- broom::tidy(fit_a)
coefs_b <- broom::tidy(fit_b)
r2_a <- broom::glance(fit_a)$r.squared
r2_b <- broom::glance(fit_b)$r.squared

summary_tbl <- dplyr::tibble(
  check = c(
    "Pairwise N: X-M","Pairwise N: X-Y","Pairwise N: M-Y",
    "Model a N (M~X+C)","Model b N (Y~X+M+C)",
    "Linearity (visual)","Homoscedasticity a (BP p)","Homoscedasticity b (BP p)",
    "Normality residuals a (SW p)","Normality residuals b (SW p)",
    "Max VIF (Y~X+M+C)","Influential cases (count)",
    "Bivariate r: X-M (p)","Bivariate r: X-Y (p)","Bivariate r: M-Y (p)",
    "a-path: X→M (coef, p)","b-path: M→Y|X (coef, p)","c'-path: X→Y|M (coef, p)"
  ),
  value = c(
    n_xm, n_xy, n_my,
    n_a, n_b,
    "See scatter & residual plots",
    round(d_a_diag$bp_p, 4), round(d_b_diag$bp_p, 4),
    round(d_a_diag$sw_p, 4), round(d_b_diag$sw_p, 4),
    round(max_vif, 3), n_flag,
    paste0(round(res_xm$r,3), " (p=", signif(res_xm$p,3), ")"),
    paste0(round(res_xy$r,3), " (p=", signif(res_xy$p,3), ")"),
    paste0(round(res_my$r,3), " (p=", signif(res_my$p,3), ")"),
    {
      ca <- coefs_a %>% dplyr::filter(term == X_var)
      if (nrow(ca)==0) NA_character_ else paste0(round(ca$estimate,3), " (p=", signif(ca$p.value,3), ")")
    },
    {
      cb <- coefs_b %>% dplyr::filter(term == M_var)
      if (nrow(cb)==0) NA_character_ else paste0(round(cb$estimate,3), " (p=", signif(cb$p.value,3), ")")
    },
    {
      cc <- coefs_b %>% dplyr::filter(term == X_var)
      if (nrow(cc)==0) NA_character_ else paste0(round(cc$estimate,3), " (p=", signif(cc$p.value,3), ")")
    }
  ),
  flag = c(
    NA, NA, NA,
    NA, NA,
    NA,
    ifelse(flag_hetero_a, "P<.05 (heteroskedastic)", "OK"),
    ifelse(flag_hetero_b, "P<.05 (heteroskedastic)", "OK"),
    ifelse(flag_nonNorm_a, "P<.05 (non-normal)", "OK"),
    ifelse(flag_nonNorm_b, "P<.05 (non-normal)", "OK"),
    ifelse(flag_vif, "High VIF", "OK"),
    ifelse(n_flag>0, paste0(n_flag, " flagged"), "OK"),
    NA, NA, NA,
    NA, NA, NA
  )
)

readr::write_csv(summary_tbl, out_summary_csv)

# --------- 8) Human-readable report ----------
sink(out_report_txt)
cat("Mediation Assumption Check\n")
cat("====================================\n")
cat("X (predictor):   ", X_var, "\n")
cat("M (mediator):    ", M_var, "\n")
cat("Y (outcome):     ", Y_var, "\n")
cat("Covariates:      ", ifelse(length(covariates)==0, "(none)", paste(covariates, collapse=", ")), "\n\n")

cat("Usable Ns:\n")
cat("- Pairwise X-M:", n_xm, "  X-Y:", n_xy, "  M-Y:", n_my, "\n")
cat("- Model a N (M~X+C):", n_a, "   Model b N (Y~X+M+C):", n_b, "\n\n")

cat("Bivariate Pearson (two-sided):\n")
cat(sprintf("  X–M: r=%.3f (p=%s)\n", res_xm$r, signif(res_xm$p,3)))
cat(sprintf("  X–Y: r=%.3f (p=%s)\n", res_xy$r, signif(res_xy$p,3)))
cat(sprintf("  M–Y: r=%.3f (p=%s)\n\n", res_my$r, signif(res_my$p,3)))

cat("Model a: M ~ X + covariates\n")
print(broom::tidy(fit_a))
cat(sprintf("R^2 = %.3f\n\n", r2_a))

cat("Model b: Y ~ X + M + covariates\n")
print(broom::tidy(fit_b))
cat(sprintf("R^2 = %.3f\n\n", r2_b))

cat("Homoscedasticity (Breusch–Pagan):\n")
cat("  Model a p =", signif(d_a_diag$bp_p, 3), "\n")
cat("  Model b p =", signif(d_b_diag$bp_p, 3), "\n\n")

cat("Normality of residuals (Shapiro–Wilk):\n")
cat("  Model a p =", ifelse(is.na(d_a_diag$sw_p), "NA (n>5000 or error)", signif(d_a_diag$sw_p, 3)), "\n")
cat("  Model b p =", ifelse(is.na(d_b_diag$sw_p), "NA (n>5000 or error)", signif(d_b_diag$sw_p, 3)), "\n\n")

cat("Multicollinearity (VIF) — Outcome model:\n")
if (all(is.na(vif_vals))) {
  cat("  VIF not available (non-numeric predictors only or error)\n\n")
} else {
  print(vif_vals)
  cat("  Max VIF =", round(max_vif,3), "\n\n")
}

cat("Influential observations (Outcome model):\n")
cat(sprintf("  Cook's D threshold: > %.3f; leverage threshold: > %.3f; |studentized| > %d\n",
            thr_cook, thr_lev, thr_stud))
cat("  Count flagged:", n_flag, "\n")
if (n_flag > 0) {
  idx <- which(flags)
  out_infl <- data.frame(
    row = idx,
    cooksD = cooks[idx],
    leverage = lev[idx],
    rstudent = stud[idx]
  )
  print(head(out_infl, 20))
  # also save a CSV of all flagged
  readr::write_csv(out_infl, file.path(out_dir, "med_influential_cases.csv"))
  cat("\n  Full list saved to med_influential_cases.csv\n")
}
cat("\nPlots saved:\n")
cat(" - ", plot_histqq_png(X_var), "\n", sep="")
cat(" - ", plot_histqq_png(M_var), "\n", sep="")
cat(" - ", plot_histqq_png(Y_var), "\n", sep="")
cat(" - ", plot_scatter_png("M_on_X"), "\n", sep="")
cat(" - ", plot_scatter_png("Y_on_X"), "\n", sep="")
cat(" - ", plot_scatter_png("Y_on_M"), "\n", sep="")
cat(" - ", plot_resid_png("M_on_X"), "\n", sep="")
cat(" - ", plot_resid_png("Y_on_XM"), "\n", sep="")
cat(" - ", plot_cook_png("Y_on_XM"), "\n", sep="")

cat("\nGuidance:\n")
cat("- Linearity: Inspect scatter plots & residual smooths. Curvature suggests transform or interaction.\n")
cat("- Homoscedasticity: BP p<.05 -> consider robust SEs (e.g., HC3) or transform.\n")
cat("- Residual normality: SW p<.05 -> large samples OK; else consider robust/bootstrapped mediation.\n")
cat("- Multicollinearity: VIF >~5 is concerning; consider centering or removing redundant predictors.\n")
cat("- Influential cases: Review Cook's/leverage/rstudent; decide on data quality vs. robustness.\n")
cat("- Causal assumptions (unmeasured confounding, correct temporal order) are *substantive* and not testable statistically.\n")
sink()

# --------- 9) Done ----------
readr::write_csv(summary_tbl, out_summary_csv)
cat("Saved summary CSV:", out_summary_csv, "\n")
cat("Saved report TXT:", out_report_txt, "\n")
cat("✓ Assumption checks completed.\n")
