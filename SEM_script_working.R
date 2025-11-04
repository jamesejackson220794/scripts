# ================================
# SEM_bayes_improved.R - Robust Bayesian SEM Analysis
# ================================

# ---- 0) Clean environment and package setup ----
# Clear any conflicting objects
rm(list = ls())

# Function to safely load/install packages
load_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load required packages
packages <- c("blavaan", "lavaan", "psych", "openxlsx", "tidyverse")
invisible(lapply(packages, load_package))

# Check blavaan version and backend availability
blav_version <- packageVersion("blavaan")
cat("blavaan version:", as.character(blav_version), "\n")

# Test which backends are available
test_backends <- function() {
  backends <- list()
  
  # Test Stan
  backends$stan <- tryCatch({
    suppressMessages(blavaan::blavInspect(NULL, "stan"))
    TRUE
  }, error = function(e) FALSE)
  
  # Test JAGS
  backends$jags <- tryCatch({
    suppressMessages(blavaan::blavInspect(NULL, "jags"))
    TRUE
  }, error = function(e) FALSE)
  
  return(backends)
}

available_backends <- test_backends()
cat("\nAvailable backends:\n")
print(available_backends)

# ---- 1) Setup paths and functions ----
root <- "/Users/jamesjackson/Desktop/dir"
phase_final <- file.path(root, "phase_final_analysis")
dir.create(phase_final, recursive = TRUE, showWarnings = FALSE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(phase_final, paste0("sem_bayes_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Output functions
write_results <- function(df, filename) {
  csv_path <- file.path(out_dir, paste0(filename, ".csv"))
  xlsx_path <- file.path(out_dir, paste0(filename, ".xlsx"))
  
  write.csv(df, csv_path, row.names = FALSE)
  
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(df, xlsx_path, overwrite = TRUE)
  }
  
  cat("Saved:", filename, "\n")
}

# ---- 2) Load and prepare data ----
# Find data file
data_paths <- c(
  file.path(root, "phase03_scoring", "scored_primary.csv"),
  file.path(root, "data", "raw copy.csv"),
  file.path(root, "data", "raw copy.xlsx")
)

df <- NULL
for (path in data_paths) {
  if (file.exists(path)) {
    if (grepl("\\.csv$", path)) {
      df <- read.csv(path, stringsAsFactors = FALSE)
    } else if (grepl("\\.xlsx$", path)) {
      df <- openxlsx::read.xlsx(path, sheet = 1)
    }
    cat("Loaded data from:", path, "\n")
    break
  }
}

if (is.null(df)) {
  stop("No data file found. Please check paths.")
}

# Normalize column names
names(df) <- gsub("^DASS-21_", "DASS.21_", names(df))

# ---- 3) Define items and compute scales ----
sacs_items <- paste0("SACS_", 1:10)
cfq_items <- paste0("CFQ_", 1:7)
atq_f_items <- paste0("ATQ_", 1:15, "f")
atq_b_items <- paste0("ATQ_", 1:15, "b")
wemwbs_items <- paste0("WEMWBS_", 1:14)
dass_items <- paste0("DASS.21_", 1:21)

# Helper function for scale computation
compute_scale <- function(data, items, scale_name) {
  if (all(items %in% names(data))) {
    data[[scale_name]] <- rowSums(data[items], na.rm = FALSE)  # Use FALSE to preserve NA
  }
  return(data)
}

# Compute scales if missing
df <- compute_scale(df, sacs_items, "SACS_Total")
df <- compute_scale(df, c("SACS_1", "SACS_2", "SACS_5", "SACS_6"), "SACS_Centering")
df <- compute_scale(df, c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_10"), "SACS_Transcending")
df <- compute_scale(df, cfq_items, "CFQ_Total")
df <- compute_scale(df, atq_b_items, "ATQ_Believability")
df <- compute_scale(df, atq_f_items, "ATQ_Frequency")
df <- compute_scale(df, wemwbs_items, "WEMWBS_Total")

# Compute DASS subscales
if (all(dass_items %in% names(df))) {
  dass_dep <- c(3, 5, 10, 13, 16, 17, 21)
  dass_anx <- c(2, 4, 7, 9, 15, 19, 20)
  dass_str <- c(1, 6, 8, 11, 12, 14, 18)
  
  df$DASS_Depression_x2 <- 2 * rowSums(df[paste0("DASS.21_", dass_dep)], na.rm = FALSE)
  df$DASS_Anxiety_x2 <- 2 * rowSums(df[paste0("DASS.21_", dass_anx)], na.rm = FALSE)
  df$DASS_Stress_x2 <- 2 * rowSums(df[paste0("DASS.21_", dass_str)], na.rm = FALSE)
  df$DASS_Total_x2 <- 2 * rowSums(df[dass_items], na.rm = FALSE)
}

# ---- 4) Prepare analysis dataset ----
# Define items for analysis (excluding SACS_9)
centering_items <- c("SACS_1", "SACS_2", "SACS_5", "SACS_6")
transcending_items <- c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_10")
all_sacs_items <- c(centering_items, transcending_items)

# Convert SACS items to ordered factors
for (item in all_sacs_items) {
  if (item %in% names(df)) {
    df[[item]] <- as.ordered(df[[item]])
  }
}

# Required variables for SEM
required_vars <- c("CFQ_Total", "ATQ_Believability", "WEMWBS_Total", "DASS_Total_x2")
missing_vars <- setdiff(required_vars, names(df))

if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

# Select complete cases for analysis
analysis_vars <- c(all_sacs_items, required_vars)
dat <- df[, analysis_vars]

# Check for missing data
n_complete <- sum(complete.cases(dat))
n_total <- nrow(dat)
cat("\nComplete cases:", n_complete, "out of", n_total, 
    "(", round(100 * n_complete/n_total, 1), "%)\n")

# Handle missing data
if (n_complete < n_total) {
  cat("Note: Using listwise deletion for missing data\n")
  dat <- dat[complete.cases(dat), ]
}

# ---- 5) Descriptive statistics ----
desc_items <- psych::describe(df[, all_sacs_items])
desc_scales <- psych::describe(df[, required_vars])

write_results(as.data.frame(desc_items), "01_descriptives_items")
write_results(as.data.frame(desc_scales), "02_descriptives_scales")

# ---- 6) Define SEM model ----
model_spec <- paste0("
  # Measurement model
  SACSc =~ ", paste(centering_items, collapse = " + "), "
  SACSt =~ ", paste(transcending_items, collapse = " + "), "
  
  # Structural paths to mediators
  CFQ_Total ~ a_c*SACSc + a_t*SACSt
  ATQ_Believability ~ g_c*SACSc + g_t*SACSt
  
  # Structural paths to outcomes
  DASS_Total_x2 ~ b1*CFQ_Total + b2*ATQ_Believability + c_c*SACSc + c_t*SACSt
  WEMWBS_Total ~ d1*CFQ_Total + d2*ATQ_Believability + e_c*SACSc + e_t*SACSt
  
  # Covariances
  CFQ_Total ~~ ATQ_Believability
  SACSc ~~ SACSt
  
  # Indirect effects via CFQ
  ind_cfq_dass_c := a_c*b1
  ind_cfq_dass_t := a_t*b1
  ind_cfq_well_c := a_c*d1
  ind_cfq_well_t := a_t*d1
  
  # Indirect effects via ATQ
  ind_atq_dass_c := g_c*b2
  ind_atq_dass_t := g_t*b2
  ind_atq_well_c := g_c*d2
  ind_atq_well_t := g_t*d2
  
  # Total indirect effects
  ind_total_dass_c := ind_cfq_dass_c + ind_atq_dass_c
  ind_total_dass_t := ind_cfq_dass_t + ind_atq_dass_t
  ind_total_well_c := ind_cfq_well_c + ind_atq_well_c
  ind_total_well_t := ind_cfq_well_t + ind_atq_well_t
")

# ---- 7) Run Bayesian SEM ----
cat("\n========== Running Bayesian SEM ==========\n")

# Function to run bsem with proper error handling
run_bayesian_sem <- function(model, data, ordered_vars, backend = "stan") {
  
  # Set MCMC parameters
  mcmc_params <- list(
    burnin = 4000,    # Increased from 3000
    sample = 10000,   # Increased from 6000
    n.chains = 3,     # Reduced from 4 for stability
    adapt = 1000      # Add adaptation phase
  )
  
  cat("\nAttempting Bayesian SEM with backend:", backend, "\n")
  
  fit <- tryCatch({
    
    # Choose backend-specific settings
    if (backend == "stan") {
      blavaan::bsem(
        model = model,
        data = data,
        ordered = ordered_vars,
        target = "stan",
        mcmcfile = TRUE,
        n.chains = mcmc_params$n.chains,
        burnin = mcmc_params$burnin,
        sample = mcmc_params$sample,
        seed = 12345,
        bcontrol = list(cores = 2)  # Limit cores for stability
      )
    } else if (backend == "jags") {
      blavaan::bsem(
        model = model,
        data = data,
        ordered = ordered_vars,
        target = "jags",
        n.chains = mcmc_params$n.chains,
        burnin = mcmc_params$burnin,
        sample = mcmc_params$sample,
        seed = 12345
      )
    } else {
      # Default (auto-select)
      blavaan::bsem(
        model = model,
        data = data,
        ordered = ordered_vars,
        n.chains = mcmc_params$n.chains,
        burnin = mcmc_params$burnin,
        sample = mcmc_params$sample,
        seed = 12345
      )
    }
    
  }, error = function(e) {
    cat("Error with", backend, ":", e$message, "\n")
    return(NULL)
  })
  
  return(fit)
}

# Try different backends
fit_bayes <- NULL
backends_to_try <- c("stan", "jags", "auto")

for (backend in backends_to_try) {
  if (backend == "stan" && !available_backends$stan) next
  if (backend == "jags" && !available_backends$jags) next
  
  fit_bayes <- run_bayesian_sem(model_spec, dat, all_sacs_items, backend)
  
  if (!is.null(fit_bayes)) {
    cat("\nSuccessfully fitted model with backend:", backend, "\n")
    break
  }
}

# ---- 8) Fallback to ML estimation if Bayesian fails ----
if (is.null(fit_bayes)) {
  cat("\n========== Bayesian estimation failed, using WLSMV ==========\n")
  
  fit_ml <- lavaan::sem(
    model = model_spec,
    data = dat,
    ordered = all_sacs_items,
    estimator = "WLSMV",
    se = "robust",
    missing = "listwise"
  )
  
  # Save ML results
  summary_ml <- summary(fit_ml, fit.measures = TRUE, standardized = TRUE)
  capture.output(summary_ml, file = file.path(out_dir, "03_ML_results.txt"))
  
  # Extract parameter estimates
  params_ml <- parameterEstimates(fit_ml, standardized = TRUE, ci = TRUE)
  write_results(params_ml, "04_ML_parameters")
  
  cat("\nML estimation completed successfully\n")
  
} else {
  # ---- 9) Extract and save Bayesian results ----
  cat("\n========== Extracting Bayesian results ==========\n")
  
  # Save summary
  summary_bayes <- summary(fit_bayes, standardized = TRUE)
  capture.output(summary_bayes, file = file.path(out_dir, "03_Bayes_results.txt"))
  
  # Check convergence
  psrf <- blavaan::blavInspect(fit_bayes, "psrf")
  if (!is.null(psrf)) {
    write_results(as.data.frame(psrf), "05_convergence_psrf")
    max_psrf <- max(psrf[is.finite(psrf)])
    cat("\nMax PSRF:", round(max_psrf, 3), 
        ifelse(max_psrf < 1.1, " (Good convergence)", " (Poor convergence)"), "\n")
  }
  
  # Posterior predictive check
  ppp <- tryCatch(blavaan::blavInspect(fit_bayes, "ppp"), error = function(e) NULL)
  if (!is.null(ppp)) {
    cat("Posterior Predictive P-value:", round(ppp, 3), "\n")
  }
  
  # Extract parameter estimates
  params_bayes <- parameterEstimates(fit_bayes, standardized = TRUE)
  write_results(params_bayes, "04_Bayes_parameters")
  
  # Extract MCMC chains for indirect effects
  mcmc_draws <- tryCatch({
    blavaan::blavInspect(fit_bayes, "mcmc")
  }, error = function(e) NULL)
  
  if (!is.null(mcmc_draws)) {
    # Combine chains
    if (is.list(mcmc_draws)) {
      combined_draws <- do.call(rbind, lapply(mcmc_draws, as.matrix))
    } else {
      combined_draws <- as.matrix(mcmc_draws)
    }
    
    # Extract indirect effects
    indirect_cols <- grep("^ind_", colnames(combined_draws), value = TRUE)
    
    if (length(indirect_cols) > 0) {
      indirect_summary <- data.frame(
        parameter = indirect_cols,
        mean = apply(combined_draws[, indirect_cols, drop = FALSE], 2, mean),
        median = apply(combined_draws[, indirect_cols, drop = FALSE], 2, median),
        sd = apply(combined_draws[, indirect_cols, drop = FALSE], 2, sd),
        lower_95 = apply(combined_draws[, indirect_cols, drop = FALSE], 2, quantile, 0.025),
        upper_95 = apply(combined_draws[, indirect_cols, drop = FALSE], 2, quantile, 0.975),
        row.names = NULL
      )
      
      write_results(indirect_summary, "06_indirect_effects")
      
      cat("\nIndirect effects summary:\n")
      print(round(indirect_summary[, -1], 3))
    }
  }
  
  cat("\nBayesian analysis completed successfully\n")
}

# ---- 10) Model comparison (if applicable) ----
if (!is.null(fit_bayes)) {
  cat("\n========== Model comparisons ==========\n")
  
  # Define alternative models
  model_direct <- paste0("
    SACSc =~ ", paste(centering_items, collapse = " + "), "
    SACSt =~ ", paste(transcending_items, collapse = " + "), "
    
    DASS_Total_x2 ~ SACSc + SACSt
    WEMWBS_Total ~ SACSc + SACSt
    
    SACSc ~~ SACSt
  ")
  
  # Fit alternative model
  fit_direct <- tryCatch({
    blavaan::bsem(
      model = model_direct,
      data = dat,
      ordered = all_sacs_items,
      target = ifelse(available_backends$stan, "stan", "jags"),
      n.chains = 2,
      burnin = 2000,
      sample = 4000,
      seed = 12345
    )
  }, error = function(e) NULL)
  
  if (!is.null(fit_direct)) {
    comparison <- tryCatch({
      blavaan::blavCompare(fit_bayes, fit_direct)
    }, error = function(e) NULL)
    
    if (!is.null(comparison)) {
      capture.output(comparison, file = file.path(out_dir, "07_model_comparison.txt"))
      cat("Model comparison saved\n")
    }
  }
}

cat("\n========== Analysis complete ==========\n")
cat("Results saved in:", out_dir, "\n")

# Print session info for reproducibility
sessionInfo_output <- capture.output(sessionInfo())
writeLines(sessionInfo_output, file.path(out_dir, "99_session_info.txt"))