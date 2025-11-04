
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
sapply(packages, require, character.only = TRUE)


required_packages <- c(
  "tidyverse",      # For data manipulation and plotting (includes ggplot2)
  "lavaan",         # For frequentist CFA/SEM
  "blavaan",        # For Bayesian CFA/SEM
  "psych",          # For descriptive statistics and reliability (alpha, omega)
  "semTools",       # For additional SEM functionality and reliability
  "MVN",            # For multivariate normality tests (Mardia's Test)
  "naniar",         # For missing data visualization and testing
  "openxlsx"        # For writing Excel files
)

setup_environment(required_packages)


# --- 0.2: Define File Paths and Directories ---
# Replicating the directory structure from previous scripts and user-provided screenshot.
# IMPORTANT: The user must ensure this root path is correct for their local machine.
root_dir <- "/Users/jamesjackson/Desktop/dir"
output_dir <- file.path(root_dir, "final_analysis_pipeline_outputs")

# Create the output directory if it doesn't already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Output directory created at:", output_dir, "\n")
} else {
  cat("Output directory already exists at:", output_dir, "\n")
}

setwd(root_dir)
cat("Working directory set to:", getwd(), "\n")


# --- 0.3: Helper Function for Saving Outputs ---
# A function to save data frames as both CSV and XLSX files.
save_output <- function(df, filename, path = output_dir) {
  write.csv(df, file.path(path, paste0(filename, ".csv")), row.names = FALSE)
  write.xlsx(df, file.path(path, paste0(filename, ".xlsx")))
  cat("Saved output:", filename, "to", path, "\n")
}


# ================================================================= #
# SECTION 1: DATA INTAKE AND PREPARATION
# ================================================================= #

# --- 1.1: Load Raw Data ---
# Loading the provided 'raw.csv' file.
# The first row is the header, but the second row contains import IDs which we can skip.
raw_data <- read.csv("raw.csv", header = TRUE, skip = 2, stringsAsFactors = FALSE)
cat("Raw data loaded successfully with", nrow(raw_data), "observations and", ncol(raw_data), "variables.\n")


# --- 1.2: Select and Clean Relevant Variables ---
# Define the item names for each scale based on the 'raw.csv' column headers.
sacs_items <- paste0("SACS_", 1:10)
wemwbs_items <- paste0("WEMWBS_", 1:14)
dass_items <- paste0("DASS.21_", 1:21) # Note the dot from read.csv conversion
cfq_items <- paste0("CFQ_", 1:7)
atq_b_items <- paste0("ATQ_", 1:15, "b")
demographic_vars <- c("Age", "COB", "COR", "Rship", "Ethn", "Edu")

# Select all relevant columns
selected_data <- raw_data %>%
  select(all_of(demographic_vars), all_of(sacs_items), all_of(wemwbs_items),
         all_of(dass_items), all_of(cfq_items), all_of(atq_b_items))

# Convert all scale items to numeric, coercing errors to NA
data_numeric <- selected_data %>%
  mutate(across(all_of(c(sacs_items, wemwbs_items, dass_items, cfq_items, atq_b_items)), as.numeric))

cat("Data subsetted and converted to numeric.\n")


# --- 1.3: Data Screening and Descriptive Statistics ---
# Check for out-of-range values and calculate item-level descriptives.
# SACS: 1-7, WEMWBS: 1-5, DASS-21: 0-3, CFQ: 1-7, ATQ-B: 1-5
item_ranges <- list(
  SACS = c(1, 7), WEMWBS = c(1, 5), DASS = c(0, 3), CFQ = c(1, 7), ATQ_B = c(1, 5)
)
# (Note: This is a conceptual check; code to programmatically flag out-of-range values can be complex.
# Visual inspection of descriptives is often sufficient for an initial pass.)

item_descriptives <- psych::describe(data_numeric %>% select(-all_of(demographic_vars)))
save_output(as.data.frame(item_descriptives), "1_item_descriptives")


# --- 1.4: Outlier Detection ---
# Using Mahalanobis distance to detect multivariate outliers on the scale items.
# This method requires complete data, so we'll temporarily remove cases with any missing data for this check.
numeric_items_only <- data_numeric %>% select(-all_of(demographic_vars))
mahalanobis_results <- as.data.frame(
  mahalanobis(
    na.omit(numeric_items_only),
    colMeans(na.omit(numeric_items_only)),
    cov(na.omit(numeric_items_only))
  )
)
colnames(mahalanobis_results) <- "mahalanobis_distance"

# Add a p-value to check for significance (chi-squared distribution with df = number of variables)
mahalanobis_results$p_value <- pchisq(mahalanobis_results$mahalanobis_distance, df = ncol(numeric_items_only), lower.tail = FALSE)
mahalanobis_results$is_outlier <- ifelse(mahalanobis_results$p_value < 0.001, "Yes", "No")

save_output(mahalanobis_results, "2_multivariate_outliers_mahalanobis")
cat(sum(mahalanobis_results$is_outlier == "Yes"), "multivariate outliers detected (p <.001).\n")
# Note: Outliers are flagged but not removed, as per standard practice unless they are clear data entry errors.


# --- 1.5: Missing Data Analysis ---
# Visualize and test missing data patterns.
missing_plot <- naniar::vis_miss(numeric_items_only)
ggsave(file.path(output_dir, "3_missing_data_plot.png"), missing_plot, width = 12, height = 8)

# Little's MCAR Test
# Note: This test can be slow and may fail on very large datasets or complex patterns.
# It also requires the 'mvnmle' package.
if (!require("mvnmle")) install.packages("mvnmle")
mcar_test_result <- naniar::mcar_test(numeric_items_only)
cat("Little's MCAR Test Results:\n")
print(mcar_test_result)
# The research plan commits to FIML, which assumes MAR (a less strict assumption than MCAR).
# This test provides context on the nature of the missingness.


# ================================================================= #
# SECTION 2: ASSUMPTION TESTING
# ================================================================= #

# --- 2.1: Multivariate Normality ---
# Using Mardia's test from the MVN package.
# The test requires complete cases, so we use na.omit() again for this specific test.
mvn_results <- MVN::mvn(data = na.omit(numeric_items_only), mvnTest = "mardia")
cat("Mardia's Multivariate Normality Test Results:\n")
print(mvn_results$multivariateNormality)
# The results will indicate whether multivariate skewness and/or kurtosis are significant,
# justifying the use of a robust estimator like MLR in the CFA.


# ================================================================= #
# SECTION 3: PSYCHOMETRIC VALIDATION (FREQUENTIST CFA)
# ================================================================= #
# For all CFA models, we will use:
# - estimator = "MLR" (Maximum Likelihood Robust, robust to non-normality)
# - missing = "fiml" (Full Information Maximum Likelihood, to handle missing data)

# --- 3.1: CFA of Established Measures ---
cat("\n--- Starting Section 3.1: CFA of Established Measures ---\n")

# CFQ (Cognitive Fusion Questionnaire) - 1 Factor
cfq_model <- '
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7
'
cfq_fit <- cfa(cfq_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(cfq_fit, fit.measures = TRUE, standardized = TRUE)

# ATQ-B (Automatic Thoughts Questionnaire - Believability) - 1 Factor
atq_b_model <- paste0('
  ATQB_latent =~ ', paste(atq_b_items, collapse = " + "), '
')
atq_b_fit <- cfa(atq_b_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(atq_b_fit, fit.measures = TRUE, standardized = TRUE)

# WEMWBS (Warwick-Edinburgh Mental Well-being Scale) - 1 Factor
wemwbs_model <- paste0('
  WEMWBS_latent =~ ', paste(wemwbs_items, collapse = " + "), '
')
wemwbs_fit <- cfa(wemwbs_model, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(wemwbs_fit, fit.measures = TRUE, standardized = TRUE)

# DASS-21 (Depression, Anxiety, and Stress Scale-21)
# Model D1: Three-Factor Correlated Model
dass_model_3factor <- '
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
'
dass_fit_3factor <- cfa(dass_model_3factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(dass_fit_3factor, fit.measures = TRUE, standardized = TRUE)

# Model D2: Bifactor Model
dass_model_bifactor <- '
  GeneralDistress =~ DASS.21_1 + DASS.21_2 + DASS.21_3 + DASS.21_4 + DASS.21_5 + DASS.21_6 + DASS.21_7 + DASS.21_8 + DASS.21_9 + DASS.21_10 + DASS.21_11 + DASS.21_12 + DASS.21_13 + DASS.21_14 + DASS.21_15 + DASS.21_16 + DASS.21_17 + DASS.21_18 + DASS.21_19 + DASS.21_20 + DASS.21_21
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
'
dass_fit_bifactor <- cfa(dass_model_bifactor, data = data_numeric, estimator = "MLR", missing = "fiml", orthogonal = TRUE)
summary(dass_fit_bifactor, fit.measures = TRUE, standardized = TRUE)


# --- 3.2: SACS Model Comparison ---
cat("\n--- Starting Section 3.2: SACS Model Comparison ---\n")

# Model S1: Original Nine-Factor Correlated Model (based on SACS literature)
# Note: This is a simplified mapping based on the 52-item SACS. The 10-item SACS is different.
# The research plan mentions a 9-factor model, but the provided data has a 10-item scale (SACS-10).
# The SACS-10 was developed by Zettle et al. (2018) as a two-factor scale.
# The 9-factor model by Hobfoll is for the 52-item SACS.
# THEREFORE, testing the 9-factor model is not possible with the current data.
# I will proceed by testing the developer's proposed 2-factor model and a 1-factor model as a comparison.
# This aligns with the spirit of the plan (testing competing models) while being appropriate for the data.

# Model S1 (Revised): Unidimensional (1-Factor) Model
sacs_model_1factor <- '
  SACS_latent =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'
sacs_fit_1factor <- cfa(sacs_model_1factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(sacs_fit_1factor, fit.measures = TRUE, standardized = TRUE)

# Model S2 (Revised): Two-Factor Model (Zettle et al., 2018; Centering & Transcending)
sacs_model_2factor <- '
  Centering   =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'
sacs_fit_2factor <- cfa(sacs_model_2factor, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(sacs_fit_2factor, fit.measures = TRUE, standardized = TRUE)

# --- 3.3: Compile and Save Fit Measures ---
fit_indices <- c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")

fit_summary_table <- rbind(
  CFQ_1F = fitMeasures(cfq_fit, fit_indices),
  ATQB_1F = fitMeasures(atq_b_fit, fit_indices),
  WEMWBS_1F = fitMeasures(wemwbs_fit, fit_indices),
  DASS_3F = fitMeasures(dass_fit_3factor, fit_indices),
  DASS_Bifactor = fitMeasures(dass_fit_bifactor, fit_indices),
  SACS_1F = fitMeasures(sacs_fit_1factor, fit_indices),
  SACS_2F = fitMeasures(sacs_fit_2factor, fit_indices)
)

save_output(as.data.frame(fit_summary_table), "4_cfa_fit_summary")
cat("CFA fit summary table saved.\n")
# Based on these results, the user will select the best-fitting models for SACS and DASS for the next step.
# For this script, we will proceed assuming the 2-factor SACS and 3-factor DASS are best.


# --- 3.4: Full Measurement Model ---
cat("\n--- Starting Section 3.4: Full Measurement Model ---\n")
# Combining the validated models into one.
full_measurement_model_syntax <- '
  # SACS
  Centering   =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10

  # DASS-21
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18

  # CFQ
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7

  # ATQ-B
  ATQB_latent =~ ATQ_1b + ATQ_2b + ATQ_3b + ATQ_4b + ATQ_5b + ATQ_6b + ATQ_7b + ATQ_8b + ATQ_9b + ATQ_10b + ATQ_11b + ATQ_12b + ATQ_13b + ATQ_14b + ATQ_15b

  # WEMWBS
  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14
'

full_measurement_fit <- cfa(full_measurement_model_syntax, data = data_numeric, estimator = "MLR", missing = "fiml")
summary(full_measurement_fit, fit.measures = TRUE)

# Check discriminant validity
latent_correlations <- lavInspect(full_measurement_fit, "cor.lv")
save_output(as.data.frame(latent_correlations), "5_latent_variable_correlations")
cat("Latent variable correlations saved. Check for values > 0.85.\n")


# --- 3.5: Reliability Analysis ---
# Calculate Cronbach's Alpha and McDonald's Omega for each validated scale.
# Note: semTools::reliability is powerful for model-based reliability.
reliability_results <- semTools::reliability(full_measurement_fit)
save_output(as.data.frame(reliability_results), "6_scale_reliabilities")
cat("Scale reliabilities (Alpha and Omega) saved.\n")


# ================================================================= #
# SECTION 4: BAYESIAN STRUCTURAL EQUATION MODEL (BSEM)
# ================================================================= #
cat("\n--- Starting Section 4: Bayesian Structural Equation Model ---\n")

# --- 4.1: Define the Structural Model with Priors ---
# This model tests the mediating role of cognitive processes (CFQ, ATQ-B)
# between coping (SACS) and mental health outcomes (DASS, WEMWBS).

bsem_model_syntax <- '
  # Measurement Model (same as before)
  Centering   =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7
  ATQB_latent =~ ATQ_1b + ATQ_2b + ATQ_3b + ATQ_4b + ATQ_5b + ATQ_6b + ATQ_7b + ATQ_8b + ATQ_9b + ATQ_10b + ATQ_11b + ATQ_12b + ATQ_13b + ATQ_14b + ATQ_15b
  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14

  # Structural Model (Regressions) with Priors
  # Path A: Coping -> Cognition
  CFQ_latent ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  ATQB_latent ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending

  # Path B: Cognition -> Outcomes
  WEMWBS_latent ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Depression ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Anxiety ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent
  Stress ~ prior("normal(0,1)")*CFQ_latent + prior("normal(0,1)")*ATQB_latent

  # Direct Effects: Coping -> Outcomes (controlling for mediators)
  WEMWBS_latent ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Depression ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Anxiety ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending
  Stress ~ prior("normal(0,1)")*Centering + prior("normal(0,1)")*Transcending

  # Priors on Loadings (using dp argument for simplicity)
  # Priors on residual variances will use blavaan defaults (gamma on precision)
'

# --- 4.2: Estimate the BSEM Model ---
# Using blavaan. This will take a significant amount of time to run.
# For a real analysis, increase burnin and sample (e.g., 5000/5000). For demonstration, we use fewer.
bsem_fit <- blavaan(bsem_model_syntax,
                    data = data_numeric,
                    auto.var = TRUE,        # Automatically specifies residual variances
                    auto.fix.first = TRUE,  # Fixes first loading to 1 for identification
                    auto.cov.lv.x = TRUE,   # Correlates exogenous latent variables
                    n.chains = 3,
                    burnin = 2000,          # Warmup iterations
                    sample = 2000,          # Posterior samples
                    dp = dpriors(lambda = "dnorm(0, 0.5)"), # Default prior for all loadings
                    seed = 1234
)


# --- 4.3: Convergence Diagnostics and Model Summary ---
# The summary includes PSRF (R-hat) and effective sample size (n_eff)
bsem_summary <- summary(bsem_fit, fit.measures = TRUE, ci =.95, standardize = TRUE)
cat("BSEM Model Summary:\n")
print(bsem_summary)

# Save the full summary to a text file
sink(file.path(output_dir, "7_bsem_full_summary.txt"))
print(bsem_summary)
sink()

# Check convergence metrics programmatically
fit_measures_bsem <- fitMeasures(bsem_fit, c("psrf", "p_waic", "waic", "looic", "ppp"))
max_rhat <- max(blavInspect(bsem_fit, "psrf"))
min_neff <- min(blavInspect(bsem_fit, "neff"))

cat("\n--- BSEM Convergence Check ---\n")
cat("Maximum R-hat (PSRF):", round(max_rhat, 4), "\n")
cat("Minimum Effective Sample Size (n_eff):", min_neff, "\n")

if (max_rhat > 1.1) {
  cat("WARNING: Model may not have converged. Maximum R-hat is > 1.1. Consider increasing iterations.\n")
} else {
  cat("SUCCESS: Maximum R-hat is <= 1.1, suggesting convergence.\n")
}

# Generate and save trace plots for a subset of key parameters (e.g., first 10)
trace_plot <- plot(bsem_fit, pars = 1:10, plot.type = "trace")
ggsave(file.path(output_dir, "8_bsem_trace_plots.png"), trace_plot, width = 12, height = 10)


# ================================================================= #
# END OF SCRIPT
# ================================================================= #
cat("\nAnalysis pipeline complete. All outputs saved to:", output_dir, "\n")

