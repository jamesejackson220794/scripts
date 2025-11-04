# ================================================================= #
# FINAL R SCRIPT FOR MASTER'S THESIS ANALYSIS (v7)
#
# This script tests the final, most defensible model based on
# previous psychometric checks.
#
# MODEL:
# 1. SACS (9-item) is a SINGLE unidimensional factor.
# 2. ATQ-B is REMOVED to solve multicollinearity.
# 3. DASS is the superior BIFACTOR model.
# 4. PATH: SACS (Uni) -> CFQ -> (DASS_General_Distress + WEMWBS)
# 5. BSEM is run in PARALLEL (3 cores).
# ================================================================= #

# --- 0.1: Environment Setup ---
setup_environment <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) {
    cat("Installing missing packages:", new_packages, "\n")
    install.packages(new_packages, dependencies = TRUE)
  }
  cat("Loading packages...\n")
  sapply(packages, require, character.only = TRUE)
}

# Added 'future' for parallel processing
required_packages <- c(
  "tidyverse", "dplyr", "lavaan", "blavaan", "psych", 
  "semTools", "MVN", "naniar", "openxlsx", "mvnmle", "future"
)
setup_environment(required_packages)


# --- 0.2: Define File Paths and Directories ---
root_dir <- "/Users/jamesjackson/Desktop/dir"
output_dir <- file.path(root_dir, "final_analysis_pipeline_outputs")
if (!dir.exists(output_dir)) dir.create(output_dir) else cat("Output directory already exists at:", output_dir, "\n")
setwd(root_dir)
cat("Working directory set to:", getwd(), "\n")


# --- 0.3: Helper Function for Saving Outputs ---
save_output <- function(df, filename, path = output_dir) {
  tryCatch({
    df_to_save <- as.data.frame(df)
    csv_path <- file.path(path, paste0(filename, ".csv"))
    xlsx_path <- file.path(path, paste0(filename, ".xlsx"))
    write.csv(df_to_save, csv_path, row.names = FALSE)
    openxlsx::write.xlsx(df_to_save, xlsx_path)
    cat("Saved output:", filename, "to", path, "\n")
  }, error = function(e) {
    cat("ERROR saving", filename, ":", e$message, "\n")
  })
}


# ================================================================= #
# SECTION 1: DATA INTAKE
# ================================================================= #

# --- 1.1: Load Raw Data ---
cat("Loading raw.csv...\n")
raw_data <- read.csv("raw.csv", header = TRUE, stringsAsFactors = FALSE)
raw_data <- raw_data[-1, ] # Remove Qualtrics metadata row
cat("Raw data loaded and metadata row removed.\n")

# --- 1.2: Select and Clean Relevant Variables for THIS model ---
# ** NOTE: SACS_9 and ATQ-B items are now excluded entirely **
sacs_items_9item <- paste0("SACS_", c(1:8, 10)) # SACS_9 is excluded
wemwbs_items <- paste0("WEMWBS_", 1:14)
dass_items <- paste0("DASS.21_", 1:21)
cfq_items <- paste0("CFQ_", 1:7)
demographic_vars <- c("Age", "COB", "COR", "Rship", "Ethn", "Edu")

# Select only the variables needed for the final model
selected_data_v7 <- raw_data %>%
  dplyr::select(all_of(demographic_vars), all_of(sacs_items_9item), all_of(wemwbs_items),
                all_of(dass_items), all_of(cfq_items))

data_numeric_v7 <- selected_data_v7 %>%
  dplyr::mutate(dplyr::across(all_of(c(sacs_items_9item, wemwbs_items, dass_items, cfq_items)), as.numeric))
cat("Data subsetted and converted to numeric for final model.\n")


# ================================================================= #
# SECTION 2: FINAL MEASUREMENT MODEL VALIDATION
# ================================================================= #
# This section confirms the validity of the *entire* measurement model
# *before* we add the structural (mediation) paths.

cat("\n--- Starting Section 2: Final Full Measurement Model (v7) ---\n")

# Define all factors for the final model
final_measurement_model_syntax <- '
  # SACS (9-item, Unidimensional)
  SACS_Uni_9item =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_10

  # DASS-21 (Bifactor)
  GeneralDistress =~ DASS.21_1 + DASS.21_2 + DASS.21_3 + DASS.21_4 + DASS.21_5 + DASS.21_6 + DASS.21_7 + DASS.21_8 + DASS.21_9 + DASS.21_10 + DASS.21_11 + DASS.21_12 + DASS.21_13 + DASS.21_14 + DASS.21_15 + DASS.21_16 + DASS.21_17 + DASS.21_18 + DASS.21_19 + DASS.21_20 + DASS.21_21
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18

  # CFQ (1-Factor)
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7

  # WEMWBS (1-Factor)
  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14
'

# Run the full measurement model CFA
final_measurement_fit <- cfa(final_measurement_model_syntax, data = data_numeric_v7, 
                             estimator = "MLR", missing = "fiml", orthogonal = TRUE)

cat("--- Full Measurement Model (v7) Summary ---\n")
print(summary(final_measurement_fit, fit.measures = TRUE))

# Save the critical fit and correlation files
fit_indices <- c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")
final_fit_summary <- fitMeasures(final_measurement_fit, fit_indices)
save_output(as.data.frame(final_fit_summary), "5_cfa_fit_summary_v7")

final_latent_correlations <- lavInspect(final_measurement_fit, "cor.lv")
save_output(as.data.frame(final_latent_correlations), "6_latent_correlations_v7")
cat("Final measurement model correlations saved. Check this file for multicollinearity (e.g., > 0.85).\n")

final_reliability <- semTools::compRelSEM(final_measurement_fit)
save_output(as.data.frame(final_reliability), "7_scale_reliabilities_v7")
cat("Final measurement model reliabilities saved.\n")


# ================================================================= #
# SECTION 3: FINAL BAYESIAN STRUCTURAL EQUATION MODEL (BSEM)
# ================================================================= #
cat("\n--- Starting Section 3: Final Bayesian Structural Equation Model (v7) ---\n")

# --- 3.1: Define the Final Structural Model ---
# This syntax combines the measurement model from above with the
# structural (mediation) paths.
final_bsem_syntax <- '
  # Measurement Model (Same as Section 2)
  SACS_Uni_9item =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_10
  GeneralDistress =~ DASS.21_1 + DASS.21_2 + DASS.21_3 + DASS.21_4 + DASS.21_5 + DASS.21_6 + DASS.21_7 + DASS.21_8 + DASS.21_9 + DASS.21_10 + DASS.21_11 + DASS.21_12 + DASS.21_13 + DASS.21_14 + DASS.21_15 + DASS.21_16 + DASS.21_17 + DASS.21_18 + DASS.21_19 + DASS.21_20 + DASS.21_21
  Depression =~ DASS.21_3 + DASS.21_5 + DASS.21_10 + DASS.21_13 + DASS.21_16 + DASS.21_17 + DASS.21_21
  Anxiety    =~ DASS.21_2 + DASS.21_4 + DASS.21_7 + DASS.21_9 + DASS.21_15 + DASS.21_19 + DASS.21_20
  Stress     =~ DASS.21_1 + DASS.21_6 + DASS.21_8 + DASS.21_11 + DASS.21_12 + DASS.21_14 + DASS.21_18
  CFQ_latent =~ CFQ_1 + CFQ_2 + CFQ_3 + CFQ_4 + CFQ_5 + CFQ_6 + CFQ_7
  WEMWBS_latent =~ WEMWBS_1 + WEMWBS_2 + WEMWBS_3 + WEMWBS_4 + WEMWBS_5 + WEMWBS_6 + WEMWBS_7 + WEMWBS_8 + WEMWBS_9 + WEMWBS_10 + WEMWBS_11 + WEMWBS_12 + WEMWBS_13 + WEMWBS_14

  # Structural (Mediation) Model with Priors
  
  # Path A: SACS -> Mediator (CFQ)
  CFQ_latent ~ a*SACS_Uni_9item
  
  # Path B: Mediator (CFQ) -> Outcomes (Distress & Wellbeing)
  GeneralDistress ~ b1*CFQ_latent
  WEMWBS_latent ~ b2*CFQ_latent
  
  # Path C' (Direct): SACS -> Outcomes
GeneralDistress ~ c1*SACS_Uni_9item
WEMWBS_latent ~ c2*SACS_Uni_9item

# Priors on structural paths
a ~ prior("normal(0,0.5)")
b1 ~ prior("normal(0,0.5)")
b2 ~ prior("normal(0,0.5)")
c1 ~ prior("normal(0,0.5)")
c2 ~ prior("normal(0,0.5)")

# Define Indirect and Total Effects for mediation
indirect_distress := a * b1
indirect_wellbeing := a * b2
total_distress := c1 + (a * b1)
total_wellbeing := c2 + (a * b2)
'

# --- 3.2: Estimate the BSEM Model (with Parallel Processing) ---
cat("Starting BSEM estimation... THIS WILL TAKE A LONG TIME.\n")
cat("Setting up parallel processing on 3 cores (for 3 chains)...\n")

# Enable parallel computation on 3 cores
future::plan("multisession", workers = 3)
options(mc.cores = 3)

final_bsem_fit <- blavaan(final_bsem_syntax,
                    data = data_numeric_v7,
                    auto.var = TRUE,
                    auto.fix.first = TRUE,
                    auto.cov.lv.x = TRUE,
                    orthogonal = TRUE, # For the bifactor DASS model
                    n.chains = 3,
                    burnin = 2000,
                    sample = 2000,
                    seed = 1234
)
cat("BSEM estimation complete.\n")

# Reset parallel plan to sequential
future::plan("sequential")

# --- 3.3: Convergence Diagnostics and Model Summary ---
final_bsem_summary <- summary(final_bsem_fit, fit.measures = TRUE, ci =.95, standardize = TRUE)
cat("--- Final BSEM Model Summary (v7) (check console) ---\n")
print(final_bsem_summary)

# Save the critical summary file
sink(file.path(output_dir, "8_bsem_full_summary_v7.txt"))
print(final_bsem_summary)
sink()

# Check convergence and fit
fit_measures_bsem_v7 <- fitMeasures(final_bsem_fit, c("psrf", "p_waic", "waic", "looic", "ppp"))
max_rhat_v7 <- max(blavInspect(final_bsem_fit, "psrf"))
min_neff_v7 <- min(blavInspect(final_bsem_fit, "neff"))

cat("\n--- BSEM Convergence Check (v7) ---\n")
cat("Maximum R-hat (PSRF):", round(max_rhat_v7, 4), "\n")
cat("Minimum Effective Sample Size (n_eff):", min_neff_v7, "\n")
cat("Posterior Predictive P-value (PPP):", fit_measures_bsem_v7["ppp"], "\n")

if (max_rhat_v7 > 1.1) {
  cat("WARNING: Model may not have converged. Maximum R-hat is > 1.1. Check trace plots.\n")
} else {
  cat("SUCCESS: Maximum R-hat is <= 1.1, suggesting convergence.\n")
}
if (fit_measures_bsem_v7["ppp"] < 0.05 || fit_measures_bsem_v7["ppp"] > 0.95) {
  cat("WARNING: Model fit is poor (PPP = ", round(fit_measures_bsem_v7["ppp"], 3), "). Interpret results with caution.\n")
} else {
  cat("SUCCESS: Model fit appears adequate (PPP = ", round(fit_measures_bsem_v7["ppp"], 3), ").\n")
}

# Generate and save trace plots
cat("Generating trace plots (v7)...\n")
trace_plot_g_v7 <- plot(final_bsem_fit, pars = 1:10, plot.type = "trace")
ggsave(file.path(output_dir, "9_bsem_trace_plots_v7.png"), trace_plot_g_v7, width = 12, height = 10)

# ================================================================= #
# END OF SCRIPT
# ================================================================= #
cat("\nAnalysis pipeline complete (v7). All outputs saved to:", output_dir, "\n")
