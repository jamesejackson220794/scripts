# ================================================================= #
# FINAL R SCRIPT FOR MASTER'S THESIS ANALYSIS (v9)
#
# PIVOT:
# Previous models (v1-v8) failed to converge, showing that a
# complex latent variable (SEM) model is not supported by the data.
# The scales (SACS, WEMWBS, etc.) show poor individual fit.
#
# NEW, ROBUST PLAN:
# 1. Acknowledge the poor CFA fit as a finding.
# 2. Create composite (mean) scores for each scale.
# 3. Use regression-based mediation (Baron & Kenny method)
#    with the `mediation` package to test the hypothesis.
#    This is a robust, valid, and defensible approach.
# 4. We will use the 9-item SACS (unidimensional) as the predictor.
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

# **NEW (v9)**: Added 'mediation' package. Removed 'blavaan', 'lavaan', 'semTools', 'future'.
required_packages <- c(
  "tidyverse", "dplyr", "psych", 
  "naniar", "openxlsx", "mvnmle", "mediation"
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
# SECTION 1: DATA INTAKE AND SCORING
# ================================================================= #

# --- 1.1: Load Raw Data ---
cat("Loading raw.csv...\n")
raw_data <- read.csv("raw.csv", header = TRUE, stringsAsFactors = FALSE)
raw_data <- raw_data[-1, ] # Remove Qualtrics metadata row
cat("Raw data loaded and metadata row removed.\n")

# --- 1.2: Select and Clean Relevant Variables ---
sacs_items_9item <- paste0("SACS_", c(1:8, 10)) # SACS_9 is excluded
wemwbs_items <- paste0("WEMWBS_", 1:14)
dass_d_items <- paste0("DASS.21_", c(3, 5, 10, 13, 16, 17, 21))
dass_a_items <- paste0("DASS.21_", c(2, 4, 7, 9, 15, 19, 20))
dass_s_items <- paste0("DASS.21_", c(1, 6, 8, 11, 12, 14, 18))
cfq_items <- paste0("CFQ_", 1:7)

# Select all items for scoring
all_items_to_keep <- c(sacs_items_9item, wemwbs_items, dass_d_items, dass_a_items, dass_s_items, cfq_items)
data_numeric_v9 <- raw_data %>%
  dplyr::select(all_of(all_items_to_keep)) %>%
  dplyr::mutate(dplyr::across(everything(), as.numeric))
cat("Data selected and converted to numeric.\n")

# --- 1.3: Create Composite Scale Scores ---
# We use rowMeans(..., na.rm = TRUE). This is a valid way
# to handle missing data at the score level (pro-rating).
scored_data <- data_numeric_v9 %>%
  dplyr::mutate(
    SACS_Score = rowMeans(dplyr::select(., all_of(sacs_items_9item)), na.rm = TRUE),
    CFQ_Score = rowMeans(dplyr::select(., all_of(cfq_items)), na.rm = TRUE),
    WEMWBS_Score = rowMeans(dplyr::select(., all_of(wemwbs_items)), na.rm = TRUE),
    DASS_D_Score = rowMeans(dplyr::select(., all_of(dass_d_items)), na.rm = TRUE),
    DASS_A_Score = rowMeans(dplyr::select(., all_of(dass_a_items)), na.rm = TRUE),
    DASS_S_Score = rowMeans(dplyr::select(., all_of(dass_s_items)), na.rm = TRUE)
  ) %>%
  # Create a total distress score for simplicity in one model
  dplyr::mutate(
    DASS_Total = (DASS_D_Score + DASS_A_Score + DASS_S_Score) / 3
  )

# Select only the final scale scores for analysis
analysis_data <- scored_data %>%
  dplyr::select(SACS_Score, CFQ_Score, WEMWBS_Score, DASS_Total, DASS_D_Score, DASS_A_Score, DASS_S_Score)

cat("Scale scores calculated.\n")

# Save descriptive statistics for the final scores
descriptives_v9 <- psych::describe(analysis_data)
save_output(as.data.frame(descriptives_v9), "1_scale_score_descriptives_v9")

# Save a correlation matrix of the final scores
cor_matrix_v9 <- cor(analysis_data, use = "pairwise.complete.obs")
save_output(as.data.frame(cor_matrix_v9), "2_scale_score_correlations_v9")


# ================================================================= #
# SECTION 2: REGRESSION-BASED MEDIATION
# ================================================================= #
# This is the final, valid analysis.
# We will test two main mediation models:
# 1. SACS -> CFQ -> Total Distress (DASS_Total)
# 2. SACS -> CFQ -> Wellbeing (WEMWBS_Score)

set.seed(1234) # for reproducible bootstrapping

# --- 2.1: Mediation Model 1 (Distress) ---
cat("\n--- Running Mediation Model 1 (SACS -> CFQ -> Total Distress) ---\n")

# Model 1a: The 'a' path (SACS predicting CFQ)
model_a <- lm(CFQ_Score ~ SACS_Score, data = analysis_data)

# Model 1b: The 'b' and 'c' paths
# (SACS and CFQ predicting Total Distress)
model_b_c1 <- lm(DASS_Total ~ SACS_Score + CFQ_Score, data = analysis_data)

# Run the mediation analysis
# (This runs 1000 bootstrapped simulations to get robust CIs)
mediation_fit_1 <- mediate(model_a, model_b_c1,
                           treat = "SACS_Score",
                           mediator = "CFQ_Score",
                           boot = TRUE, sims = 1000)

cat("--- Mediation 1 (Distress) Results ---\n")
print(summary(mediation_fit_1))
sink(file.path(output_dir, "3_mediation_summary_DASS_TOTAL_v9.txt"))
print(summary(mediation_fit_1))
sink()


# --- 2.2: Mediation Model 2 (Wellbeing) ---
cat("\n--- Running Mediation Model 2 (SACS -> CFQ -> Wellbeing) ---\n")

# Model 2a: The 'a' path (same as before)
# model_a <- lm(CFQ_Score ~ SACS_Score, data = analysis_data)

# Model 2b: The 'b' and 'c' paths
# (SACS and CFQ predicting Wellbeing)
model_b_c2 <- lm(WEMWBS_Score ~ SACS_Score + CFQ_Score, data = analysis_data)

# Run the mediation analysis
mediation_fit_2 <- mediate(model_a, model_b_c2,
                           treat = "SACS_Score",
                           mediator = "CFQ_Score",
                           boot = TRUE, sims = 1000)

cat("--- Mediation 2 (Wellbeing) Results ---\n")
print(summary(mediation_fit_2))
sink(file.path(output_dir, "4_mediation_summary_WEMWBS_v9.txt"))
print(summary(mediation_fit_2))
sink()

# ================================================================= #
# END OF SCRIPT
# ================================================================= #
cat("\nAnalysis pipeline complete (v9). All outputs saved to:", output_dir, "\n")
