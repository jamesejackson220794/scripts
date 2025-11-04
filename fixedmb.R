# file: scripts/final_analysis_pipeline.R

# --------------------------- Setup ---------------------------

setup_environment <- function(packages) {
  # Why: Ensure reproducibility and non-interactive installs.
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  suppressPackageStartupMessages(invisible(lapply(packages, require, character.only = TRUE)))
}

required_packages <- c(
  "tidyverse", "lavaan", "blavaan", "psych", "semTools",
  "MVN", "naniar", "openxlsx", "BaylorEdPsych", "ggplot2"
)

setup_environment(required_packages)

# --------------------- Paths and I/O helpers ------------------

root_dir <- "/Users/jamesjackson/Desktop/dir"
output_dir <- file.path(root_dir, "final_analysis_pipeline_outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Output directory created at:", output_dir, "\n")
} else {
  cat("Output directory already exists at:", output_dir, "\n")
}

setwd(root_dir)
cat("Working directory set to:", getwd(), "\n")

save_output <- function(df, filename, path = output_dir) {
  # Why: Dual-format outputs ease audit and sharing.
  stopifnot(is.data.frame(df))
  csv_path <- file.path(path, paste0(filename, ".csv"))
  xlsx_path <- file.path(path, paste0(filename, ".xlsx"))
  write.csv(df, csv_path, row.names = FALSE)
  openxlsx::write.xlsx(df, xlsx_path, rowNames = FALSE)
  cat("Saved:", basename(csv_path), "and", basename(xlsx_path), "\n")
}

# ------------------ 1) Data intake & prep --------------------

if (!file.exists("raw.csv")) stop("raw.csv not found in ", getwd())

# Keep header from row 1. Drop the second physical row post-read.
raw_data <- read.csv("raw.csv", header = TRUE, stringsAsFactors = FALSE, check.names = TRUE)
if (nrow(raw_data) < 1) stop("raw.csv has no data rows.")
raw_data <- raw_data[-1, , drop = FALSE]  # drop import ID row
cat("Raw data loaded with", nrow(raw_data), "rows and", ncol(raw_data), "columns.\n")

# Item names
sacs_items      <- paste0("SACS_", 1:10)
wemwbs_items    <- paste0("WEMWBS_", 1:14)
dass_items      <- paste0("DASS.21_", 1:21)   # header already standardized by read.csv
cfq_items       <- paste0("CFQ_", 1:7)
atq_b_items     <- paste0("ATQ_", 1:15, "b")
demographic_vars <- c("Age", "COB", "COR", "Rship", "Ethn", "Edu")

# Column validation
all_needed <- c(demographic_vars, sacs_items, wemwbs_items, dass_items, cfq_items, atq_b_items)
missing_cols <- setdiff(all_needed, colnames(raw_data))
if (length(missing_cols)) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

selected_data <- raw_data %>%
  dplyr::select(dplyr::all_of(all_needed))

# Coerce items to numeric
item_vars <- c(sacs_items, wemwbs_items, dass_items, cfq_items, atq_b_items)
suppressWarnings({
  data_numeric <- selected_data %>%
    mutate(across(all_of(item_vars), ~ as.numeric(.)))
})
cat("Data subsetted and coerced to numeric for item variables.\n")

# ---------------- 1.3) Descriptives and ranges ---------------

item_descriptives <- psych::describe(dplyr::select(data_numeric, -all_of(demographic_vars)))
save_output(as.data.frame(item_descriptives), "1_item_descriptives")

# ------------------ 1.4) Outlier detection -------------------

numeric_items_only <- dplyr::select(data_numeric, -all_of(demographic_vars))
cc <- stats::na.omit(numeric_items_only)
if (nrow(cc) > 2) {
  md <- stats::mahalanobis(cc, colMeans(cc), stats::cov(cc))
  mahalanobis_results <- tibble::tibble(
    row_id_complete_cases = as.integer(attr(cc, "na.action")) %||% seq_len(nrow(cc)),  # robust index
    mahalanobis_distance = as.numeric(md),
    p_value = stats::pchisq(md, df = ncol(cc), lower.tail = FALSE),
    is_outlier = ifelse(stats::pchisq(md, df = ncol(cc), lower.tail = FALSE) < 0.001, "Yes", "No")
  )
  save_output(as.data.frame(mahalanobis_results), "2_multivariate_outliers_mahalanobis")
  cat(sum(mahalanobis_results$is_outlier == "Yes"), "multivariate outliers at p<.001.\n")
} else {
  cat("Insufficient complete cases for Mahalanobis outlier detection.\n")
}

# ------------------ 1.5) Missing data tests ------------------

miss_plot <- naniar::vis_miss(numeric_items_only)
ggplot2::ggsave(filename = file.path(output_dir, "3_missing_data_plot.png"), plot = miss_plot, width = 12, height = 8)

# Little's MCAR (BaylorEdPsych)
# Why: naniar::mcar_test dispatches BaylorEdPsych::LittleMCAR; call directly for clearer output.
mcar_data <- as.data.frame(numeric_items_only)
mcar_res <- try(BaylorEdPsych::LittleMCAR(mcar_data), silent = TRUE)
cat("Little's MCAR Test:\n")
print(mcar_res)

# ------------------ 2) Assumption testing --------------------

cc_mvn <- stats::na.omit(numeric_items_only)
if (nrow(cc_mvn) > 2) {
  mvn_results <- MVN::mvn(data = cc_mvn, mvnTest = "mardia")
  cat("Mardia multivariate normality:\n")
  print(mvn_results$multivariateNormality)
} else {
  cat("Insufficient complete cases for MVN test.\n")
}

# ------------------ 3) Frequentist CFA -----------------------

cat("\n--- Section 3.1: CFA of established measures ---\n")

cfq_model <- '
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7
'
cfq_fit <- lavaan::cfa(cfq_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(cfq_fit, fit.measures = TRUE, standardized = TRUE)

atq_b_model <- paste0('ATQB_latent =~ ', paste(atq_b_items, collapse = " + "))
atq_b_fit <- lavaan::cfa(atq_b_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(atq_b_fit, fit.measures = TRUE, standardized = TRUE)

wemwbs_model <- paste0('WEMWBS_latent =~ ', paste(wemwbs_items, collapse = " + "))
wemwbs_fit <- lavaan::cfa(wemwbs_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(wemwbs_fit, fit.measures = TRUE, standardized = TRUE)

dass_model_3factor <- '
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
'
dass_fit_3factor <- lavaan::cfa(dass_model_3factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(dass_fit_3factor, fit.measures = TRUE, standardized = TRUE)

dass_model_bifactor <- '
  GeneralDistress =~ DASS.21_1 + DASS.21_2 + DASS.21_3 + DASS.21_4 + DASS.21_5 + DASS.21_6 + DASS.21_7 + DASS.21_8 + DASS.21_9 + DASS.21_10 + DASS.21_11 + DASS.21_12 + DASS.21_13 + DASS.21_14 + DASS.21_15 + DASS.21_16 + DASS.21_17 + DASS.21_18 + DASS.21_19 + DASS.21_20 + DASS.21_21
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
'
dass_fit_bifactor <- lavaan::cfa(dass_model_bifactor, data = data_numeric, estimator = "MLR", missing = "fiml", orthogonal = TRUE)
summary(dass_fit_bifactor, fit.measures = TRUE, standardized = TRUE)

cat("\n--- Section 3.2: SACS model comparison ---\n")

sacs_model_1factor <- '
  SACS_latent =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'
sacs_fit_1factor <- lavaan::cfa(sacs_model_1factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(sacs_fit_1factor, fit.measures = TRUE, standardized = TRUE)

sacs_model_2factor <- '
  Centering    =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'
sacs_fit_2factor <- lavaan::cfa(sacs_model_2factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(sacs_fit_2factor, fit.measures = TRUE, standardized = TRUE)

fit_indices <- c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")
fit_summary_table <- rbind(
  CFQ_1F       = lavaan::fitMeasures(cfq_fit, fit_indices),
  ATQB_1F      = lavaan::fitMeasures(atq_b_fit, fit_indices),
  WEMWBS_1F    = lavaan::fitMeasures(wemwbs_fit, fit_indices),
  DASS_3F      = lavaan::fitMeasures(dass_fit_3factor, fit_indices),
  DASS_Bifactor= lavaan::fitMeasures(dass_fit_bifactor, fit_indices),
  SACS_1F      = lavaan::fitMeasures(sacs_fit_1factor, fit_indices),
  SACS_2F      = lavaan::fitMeasures(sacs_fit_2factor, fit_indices)
)
save_output(as.data.frame(fit_summary_table), "4_cfa_fit_summary")
cat("CFA fit summary saved.\n")

cat("\n--- Section 3.4: Full measurement model ---\n")

full_measurement_model_syntax <- '
  Centering    =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10

  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18

  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7

  ATQB_latent =~ ATQ_1b + ATQ_2b + ATQ_3b + ATQ_4b + ATQ_5b + ATQ_6b + ATQ_7b + ATQ_8b + ATQ_9b + ATQ_10b + ATQ_11b + ATQ_12b + ATQ_13b + ATQ_14b + ATQ_15b

  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14
'
full_measurement_fit <- lavaan::cfa(full_measurement_model_syntax, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(full_measurement_fit, fit.measures = TRUE)

latent_correlations <- lavaan::lavInspect(full_measurement_fit, "cor.lv")
save_output(as.data.frame(latent_correlations), "5_latent_variable_correlations")
cat("Latent variable correlations saved. Inspect values > 0.85.\n")

reliability_results <- semTools::reliability(full_measurement_fit)
save_output(as.data.frame(reliability_results), "6_scale_reliabilities")
cat("Scale reliabilities saved.\n")

# ------------------ 4) Bayesian SEM via blavaan --------------

cat("\n--- Section 4: Bayesian SEM ---\n")

bsem_model_syntax <- '
  Centering    =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10

  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18

  CFQ_latent  =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7
  ATQB_latent =~ ATQ_1b + ATQ_2b + ATQ_3b + ATQ_4b + ATQ_5b + ATQ_6b + ATQ_7b + ATQ_8b + ATQ_9b + ATQ_10b + ATQ_11b + ATQ_12b + ATQ_13b + ATQ_14b + ATQ_15b
  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14

  CFQ_latent  ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  ATQB_latent ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending

  WEMWBS_latent ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Depression    ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Anxiety       ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Stress        ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent

  WEMWBS_latent ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Depression    ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Anxiety       ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Stress        ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
'

bsem_fit <- blavaan::blavaan(
  bsem_model_syntax,
  data = data_numeric,
  auto.var = TRUE,
  auto.fix.first = TRUE,
  auto.cov.lv.x = TRUE,
  n.chains = 3,
  burnin = 2000,
  sample = 2000,
  dp = blavaan::dpriors(lambda = "dnorm(0, 0.5)"),
  seed = 1234
)

bsem_summary <- summary(bsem_fit, fit.measures = TRUE, ci = .95, standardized = TRUE)
cat("BSEM Model Summary:\n"); print(bsem_summary)

sink(file.path(output_dir, "7_bsem_full_summary.txt")); print(bsem_summary); sink()

# Convergence checks
psrf_vals <- try(blavaan::blavInspect(bsem_fit, "psrf"), silent = TRUE)
neff_vals <- try(blavaan::blavInspect(bsem_fit, "neff"), silent = TRUE)
if (!inherits(psrf_vals, "try-error")) {
  max_rhat <- max(psrf_vals, na.rm = TRUE)
  min_neff <- min(neff_vals, na.rm = TRUE)
  cat("\n--- BSEM Convergence Check ---\n")
  cat("Maximum R-hat (PSRF):", round(max_rhat, 4), "\n")
  cat("Minimum Effective Sample Size (n_eff):", min_neff, "\n")
  if (is.finite(max_rhat) && max_rhat > 1.1) {
    cat("WARNING: Potential non-convergence. Consider more iterations or re-specification.\n")
  } else {
    cat("R-hat ≤ 1.1 suggests convergence.\n")
  }
}

# Trace plots for first 10 parameters if available
tp <- try(plot(bsem_fit, pars = 1:10, plot.type = "trace"), silent = TRUE)
if (!inherits(tp, "try-error")) {
  ggplot2::ggsave(file.path(output_dir, "8_bsem_trace_plots.png"), tp, width = 12, height = 10)
  cat("Trace plots saved.\n")
} else {
  cat("Trace plotting not available in this session.\n")
}

cat("\nAnalysis complete. Outputs:", output_dir, "\n")
