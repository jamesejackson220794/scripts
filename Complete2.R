# ---- 0. Setup environment and paths ----

# Required packages for the analysis
required_packages <- c("lavaan", "psych", "mice", "boot", "openxlsx")
missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])
if(length(missing_packages) > 0){
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Define root directory and output directory for this analysis
root <- "/Users/jamesjackson/Desktop/dir"
out_dir <- file.path(root, "phase_new_analysis")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Utility: dual writer to CSV and XLSX (if openxlsx is available)
xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)
if(!xlsx_ok){
  warning("Package 'openxlsx' not installed; XLSX outputs will be skipped.")
}
write_dual <- function(df, stem, out_dir = NULL){
  if (is.null(out_dir)) out_dir <- get("out_dir", envir = .GlobalEnv)
  
  csv_path  <- file.path(out_dir, paste0(stem, ".csv"))
  xlsx_path <- file.path(out_dir, paste0(stem, ".xlsx"))
  
  # Write CSV
  tryCatch(
    {
      utils::write.csv(df, file = csv_path, row.names = FALSE)
    },
    error = function(e){
      write.csv(df, file = csv_path, row.names = FALSE)
    }
  )
  
  # Write XLSX if available
  if (xlsx_ok){
    openxlsx::write.xlsx(df, file = xlsx_path, overwrite = TRUE)
  }
  
  return(invisible(NULL))
}

# ---- 1. Load dataset ----

# Prefer the scored (cleaned) dataset from Phase 03, otherwise fallback to raw
scored_file <- file.path(root, "phase03_scoring", "scored_primary.csv")
if(file.exists(scored_file)){
  df <- read.csv(scored_file, stringsAsFactors = FALSE, check.names = TRUE)
} else {
  # Fallback: use raw.csv (assuming lexicon and previous scripts aren't available in this context)
  raw_file <- file.path(root, "raw.csv")
  if(!file.exists(raw_file)){
    stop("No input data found. Ensure 'scored_primary.csv' or 'raw.csv' is available in the project directory.")
  }
  df <- read.csv(raw_file, stringsAsFactors = FALSE, check.names = TRUE)
  # If raw is loaded, you may need to perform data cleaning and scoring here 
  # (not shown to avoid duplication; assumed to have been done in earlier scripts).
}

# Convert any demographic codes to appropriate types
if("Gender_label" %in% names(df)){
  df$Gender_label <- as.factor(df$Gender_label)  # treat gender as factor
}
# (Assume Age_clean is already numeric; if not, convert to numeric)

# ---- 2. Confirmatory Factor Analysis for SACS ----

# Define SACS item vectors (10 items total, 4 centering, 6 transcending)
sacs_items        <- sprintf("SACS_%d", 1:10)
sacs_center_items <- c("SACS_1","SACS_2","SACS_5","SACS_6")      # Centering items
sacs_trans_items  <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")  # Transcending items

# Ensure all SACS item columns are present
stopifnot(all(sacs_items %in% names(df)))

# Specify CFA models
model_2factor <- "
  Centering =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
  Centering ~~ Transcending  # allow factors to correlate
"
model_1factor <- "
  SAC_single =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
"

# Fit CFA models with WLSMV estimator for ordinal data (7-point Likert)
cfa_fit_2f <- lavaan::cfa(model_2factor, data = df, estimator = "WLSMV", ordered = sacs_items)
cfa_fit_1f <- lavaan::cfa(model_1factor, data = df, estimator = "WLSMV", ordered = sacs_items)

# Extract fit measures for both models
fit_measures <- c("chisq","df","pvalue","cfi","tli","rmsea","srmr")
fit_2f_vals <- lavaan::fitMeasures(cfa_fit_2f, fit_measures)
fit_1f_vals <- lavaan::fitMeasures(cfa_fit_1f, fit_measures)

# Prepare a comparative fit table
cfa_comparison <- data.frame(
  Model = c("Two-factor (Centering+Transcending)", "One-factor (SACS total)"),
  ChiSq = c(fit_2f_vals["chisq"], fit_1f_vals["chisq"]),
  df    = c(fit_2f_vals["df"],    fit_1f_vals["df"]),
  pvalue= c(fit_2f_vals["pvalue"],fit_1f_vals["pvalue"]),
  CFI   = c(fit_2f_vals["cfi"],   fit_1f_vals["cfi"]),
  TLI   = c(fit_2f_vals["tli"],   fit_1f_vals["tli"]),
  RMSEA = c(fit_2f_vals["rmsea"], fit_1f_vals["rmsea"]),
  SRMR  = c(fit_2f_vals["srmr"],  fit_1f_vals["srmr"]),
  stringsAsFactors = FALSE
)
# Add factor correlation for 2-factor model (standardized)
std_fit2f <- lavaan::inspect(cfa_fit_2f, "std")  # standardized solution
if(!is.null(std_fit2f$psi)){ 
  factor_corr <- std_fit2f$psi[1,2]  # correlation between Centering and Transcending
  cfa_comparison$Factor_Corr <- c(sprintf("%.2f", factor_corr), NA)
}

# Save fit comparison
write_dual(cfa_comparison, "CFA_model_fit_comparison")

# Extract standardized loadings for 2-factor model
std_loadings <- std_fit2f$lambda  # matrix of standardized loadings (items x factors)
loading_df <- data.frame(
  Item = rownames(std_loadings),
  Loading_Centering = std_loadings[,"Centering"],
  Loading_Transcending = std_loadings[,"Transcending"],
  row.names = NULL
)
write_dual(loading_df, "CFA_standardized_loadings")

# (If needed, one could also conduct a scaled chi-square difference test:
# lavTestLRT(cfa_fit_1f, cfa_fit_2f)  to statistically compare models.
# Here, we focus on fit indices differences (e.g., CFI change) as per guidelines.)

# ---- 3. Reliability (Cronbach's Alpha) for each scale ----

# Define item sets for each scale
cfq_items    <- sprintf("CFQ_%d", 1:7)
atq_f_items  <- sprintf("ATQ_%df", 1:15)
atq_b_items  <- sprintf("ATQ_%db", 1:15)
wemwbs_items <- sprintf("WEMWBS_%d", 1:14)
dass_items   <- sprintf("DASS.21_%d", 1:21)  # all DASS-21 item columns

# Ensure these item columns exist
stopifnot(all(c(cfq_items, atq_f_items, atq_b_items, wemwbs_items, dass_items) %in% names(df)))

# Function to compute alpha for a set of items (using only complete-case on those items)
alpha_for_items <- function(item_vec){
  sub_df <- df[ , item_vec]
  # Use only rows with no missing in these columns (complete cases for that scale)
  sub_df <- sub_df[complete.cases(sub_df), , drop=FALSE]
  if(nrow(sub_df) < 2){
    return(NA)  # not enough data
  }
  # Compute Cronbach's alpha using 'psych'
  a <- psych::alpha(sub_df, check.keys = FALSE)
  return(a$total$raw_alpha)
}

# Compute alphas
alpha_sacs_total      <- alpha_for_items(sacs_items)
alpha_sacs_center     <- alpha_for_items(sacs_center_items)
alpha_sacs_trans      <- alpha_for_items(sacs_trans_items)
alpha_cfq             <- alpha_for_items(cfq_items)
alpha_atq_frequency   <- alpha_for_items(atq_f_items)
alpha_atq_belief      <- alpha_for_items(atq_b_items)
alpha_wemwbs          <- alpha_for_items(wemwbs_items)
alpha_dass_total      <- alpha_for_items(dass_items)  # all 21 DASS items for total

# Prepare table of results
reliability_table <- data.frame(
  Scale               = c("SACS_Total", "SACS_Centering", "SACS_Transcending", 
                          "CFQ_Total", "ATQ_Frequency", "ATQ_Believability", 
                          "WEMWBS_Total", "DASS_21_Total"),
  Items               = c(length(sacs_items), length(sacs_center_items), length(sacs_trans_items),
                          length(cfq_items), length(atq_f_items), length(atq_b_items),
                          length(wemwbs_items), length(dass_items)),
  Cronbach_alpha      = c(alpha_sacs_total, alpha_sacs_center, alpha_sacs_trans, 
                          alpha_cfq, alpha_atq_frequency, alpha_atq_belief, 
                          alpha_wemwbs, alpha_dass_total)
)
# Round alpha values for neatness
reliability_table$Cronbach_alpha <- sprintf("%.3f", reliability_table$Cronbach_alpha)

# Save reliability results
write_dual(reliability_table, "Scale_Reliability_Alphas")

# ---- 4. Descriptive statistics (N, mean, SD, skew, kurtosis) ----

# List of main score variables (as per analysis plan)
score_vars <- c("SACS_Total","SACS_Centering","SACS_Transcending",
                "CFQ_Total",
                "ATQ_Frequency","ATQ_Believability","ATQ_Combined",
                "WEMWBS_Total",
                "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2")
score_vars <- intersect(score_vars, names(df))  # ensure present in data

# Function to get descriptive stats for one vector
desc_one <- function(x){
  x <- x[is.finite(x)]  # drop NA/NaN
  n <- length(x)
  res <- c(N = n, Missing = sum(!is.finite(df[[v]])), Mean = NA, SD = NA, Skewness = NA, Kurtosis_excess = NA)
  if(n > 0){
    m <- mean(x)
    s <- sd(x)
    res["Mean"] <- m
    res["SD"] <- s
    if(n > 2 && s > 0){
      res["Skewness"] <- mean(((x - m)/s)^3)
    }
    if(n > 3 && s > 0){
      res["Kurtosis_excess"] <- mean(((x - m)/s)^4) - 3
    }
  }
  return(res)
}

# Compute descriptives for each variable
desc_list <- lapply(score_vars, function(v) {
  stats <- desc_one(df[[v]])
  data.frame(Variable = v, t(stats), row.names = NULL)
})
descriptives <- do.call(rbind, desc_list)
# Round numeric columns for readability
numeric_cols <- c("Mean","SD","Skewness","Kurtosis_excess")
descriptives[numeric_cols] <- lapply(descriptives[numeric_cols], function(col) sprintf("%.3f", as.numeric(col)))
write_dual(descriptives, "Descriptive_Stats")

# ---- Normality tests (Shapiro-Wilk) ----
norm_tests <- data.frame(Variable = character(), W = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
for(v in score_vars){
  if(is.numeric(df[[v]]) || is.integer(df[[v]])){
    x <- df[[v]][is.finite(df[[v]])]
    if(length(x) >= 3){
      sw <- try(shapiro.test(x), silent = TRUE)
      if(inherits(sw, "htest")){
        norm_tests <- rbind(norm_tests, data.frame(Variable=v, W=unname(sw$statistic), p_value=sw$p.value))
      }
    }
  }
}
norm_tests$p_value <- sprintf("%.4f", norm_tests$p_value)
write_dual(norm_tests, "Normality_Tests_Shapiro")

# ---- Multicollinearity check (VIF for facets regression) ----
if(!requireNamespace("car", quietly = TRUE)){
  install.packages("car", dependencies = TRUE)
}
library(car)
# Build a regression model with SACS facets predicting distress (DASS_Total_x2) including covariates
predictors <- c("SACS_Centering","SACS_Transcending")
if("Age_clean" %in% names(df))      predictors <- c(predictors, "Age_clean")
if("Gender_label" %in% names(df))   predictors <- c(predictors, "Gender_label")
reg_formula <- as.formula(paste("DASS_Total_x2 ~", paste(predictors, collapse = " + ")))
vif_model <- try(lm(reg_formula, data = df), silent = TRUE)
vif_results <- data.frame(Variable = character(), VIF = numeric())
if(!inherits(vif_model, "try-error")){
  vifs <- car::vif(vif_model)
  # car::vif returns named vector (or matrix if factors), ensure it's in data frame form
  vif_df <- if(is.matrix(vifs)) {
    # If factor present, car::vif returns a matrix with GVIF etc. We take GVIF^(1/(2*Df)) for generalized VIF.
    gvif <- vifs[,"GVIF"]
    df_factor <- vifs[,"Df"]
    adj_vif <- gvif^(1/(2*df_factor))  # adjusted GVIF
    data.frame(Variable = rownames(vifs), VIF = adj_vif)
  } else {
    data.frame(Variable = names(vifs), VIF = as.numeric(vifs))
  }
  vif_results <- vif_df
}
write_dual(vif_results, "VIF_Facets_Model")


# ---- 5. Correlation matrix and significance tests ----

# Subset numeric data for correlation computations
numeric_data <- df[ , score_vars]
# Compute Pearson correlation matrix (pairwise complete observations)
cor_matrix <- suppressWarnings(cor(numeric_data, use = "pairwise.complete.obs", method = "pearson"))
cor_matrix_df <- data.frame(Variable = rownames(cor_matrix), cor_matrix, check.names = FALSE)
write_dual(cor_matrix_df, "Correlation_Matrix_Pearson")

# Compute Pearson correlations with significance (p-values) in long format
pairs <- t(combn(score_vars, 2))  # all unique variable pairs
cor_long <- data.frame(var1=character(), var2=character(), N=integer(), r=numeric(), p_value=numeric(), stringsAsFactors=FALSE)
for(i in 1:nrow(pairs)){
  v1 <- pairs[i,1]; v2 <- pairs[i,2]
  x <- df[[v1]]; y <- df[[v2]]
  ok <- is.finite(x) & is.finite(y)
  n_pair <- sum(ok)
  if(n_pair < 3) next  # skip if insufficient data
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method="pearson"))
  cor_long <- rbind(cor_long, data.frame(var1=v1, var2=v2, N=n_pair, r=unname(ct$estimate), p_value=ct$p.value, stringsAsFactors=FALSE))
}
# Format p-values for readability
cor_long$p_value <- sprintf("%.4f", cor_long$p_value)
write_dual(cor_long, "Correlation_LongFormat")

# ----  Bootstrapped 95% CIs for correlations ----
set.seed(123)
R <- 5000  # number of bootstrap resamples
boot_cor_results <- data.frame(var1=character(), var2=character(), boot_r=numeric(), CI_low=numeric(), CI_high=numeric(), stringsAsFactors=FALSE)
# Define a function to compute correlation on a resampled dataset for boot
boot_cor_fn <- function(data, idx, var1, var2){
  d <- data[idx, , drop=FALSE]
  x <- d[[var1]]; y <- d[[var2]]
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok) < 3) return(NA_real_)
  return(cor(x[ok], y[ok], method="pearson"))
}
for(i in 1:nrow(pairs)){
  v1 <- pairs[i,1]; v2 <- pairs[i,2]
  # Perform bootstrap on rows for this pair
  bt <- boot::boot(data = df, statistic = function(d, idx) boot_cor_fn(d, idx, v1, v2), R = R)
  t_vals <- bt$t[is.finite(bt$t)]
  if(length(t_vals) > 0){
    boot_mean <- mean(t_vals)
    ci <- quantile(t_vals, c(0.025, 0.975), na.rm = TRUE)
    boot_cor_results <- rbind(boot_cor_results, data.frame(
      var1 = v1, var2 = v2, 
      boot_r = boot_mean, 
      CI_low = ci[[1]], CI_high = ci[[2]], 
      stringsAsFactors = FALSE
    ))
  }
}
# Round results for output
boot_cor_results$boot_r  <- sprintf("%.3f", boot_cor_results$boot_r)
boot_cor_results$CI_low  <- sprintf("%.3f", boot_cor_results$CI_low)
boot_cor_results$CI_high <- sprintf("%.3f", boot_cor_results$CI_high)
write_dual(boot_cor_results, "Correlation_Bootstrap_CIs")


# ---- 6. Regression analyses for H2 (facets predicting outcomes) ----

# Helper function: get standardized beta (using lm.beta from package "lm.beta" if available, else manual)
standardize_coefs <- function(fit){
  if(requireNamespace("lm.beta", quietly = TRUE)){
    beta_fit <- lm.beta::lm.beta(fit)
    coefs <- summary(beta_fit)$coefficients
  } else {
    coefs <- summary(fit)$coefficients
    # Manually add standardized estimate if lm.beta is not installed:
    # beta = b * (sd(X)/sd(Y)) for each predictor; for factor, skip standardization.
    numeric_preds <- sapply(model.matrix(fit)[, -1, drop=FALSE], is.numeric)
    y_sd <- sd(model.response(model.frame(fit)), na.rm=TRUE)
    for(pred in names(numeric_preds)[numeric_preds]){
      x_sd <- sd(fit$model[[pred]], na.rm=TRUE)
      if(x_sd > 0){
        coefs[pred,"Std..Beta"] <- coefs[pred,"Estimate"] * (x_sd / y_sd)
      }
    }
  }
  return(as.data.frame(coefs))
}

# Define outcomes of interest
outcomes_H2 <- c("DASS_Total_x2", "WEMWBS_Total", "CFQ_Total", "ATQ_Believability")
# Ensure these outcomes exist in data
outcomes_H2 <- outcomes_H2[outcomes_H2 %in% names(df)]

# Run regressions for each outcome
for(outcome in outcomes_H2){
  formula_H2 <- as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending",
                                 if("Age_clean" %in% names(df)) "+ Age_clean" else "",
                                 if("Gender_label" %in% names(df)) "+ Gender_label" else ""))
  fit <- lm(formula_H2, data = df)
  # Coefficient table with estimates, std. error, t, p
  coef_tab <- coef(summary(fit))
  coef_df <- data.frame(Variable = rownames(coef_tab),
                        Estimate = coef_tab[,"Estimate"],
                        Std_Error= coef_tab[,"Std. Error"],
                        t_value = coef_tab[,"t value"],
                        p_value = coef_tab[,"Pr(>|t|)"],
                        stringsAsFactors = FALSE)
  # Add standardized beta
  coef_std_df <- standardize_coefs(fit)
  if("Std..Beta" %in% names(coef_std_df)){
    coef_df$Std_Beta <- coef_std_df$"Std..Beta"
  }
  # Round and format
  num_cols <- c("Estimate","Std_Error","t_value","p_value","Std_Beta")
  for(col in num_cols){
    if(col %in% names(coef_df)){
      coef_df[[col]] <- sprintf("%.4f", as.numeric(coef_df[[col]]))
    }
  }
  # Save coefficients
  out_stem <- paste0("Regression_", outcome, "_facets_coeffs")
  write_dual(coef_df, out_stem)
  
  # Hierarchical R-squared analysis:
  # R2 of model with both facets, R2 with only Transcending, R2 with only Centering
  base_data <- df[complete.cases(df[, all.vars(formula_H2)]), ]  # complete cases on all predictors & outcome
  if(nrow(base_data) > 0){
    # Model with both
    full_model <- lm(formula_H2, data = base_data)
    R2_full <- summary(full_model)$r.squared
    # Model with Transcending only (plus covars)
    form_trans <- update(formula_H2, paste(". ~ SACS_Transcending",
                                           if("Age_clean" %in% names(df)) "+ Age_clean" else "",
                                           if("Gender_label" %in% names(df)) "+ Gender_label" else ""))
    model_trans <- lm(form_trans, data = base_data)
    R2_trans <- summary(model_trans)$r.squared
    # Model with Centering only (plus covars)
    form_cent <- update(formula_H2, paste(". ~ SACS_Centering",
                                          if("Age_clean" %in% names(df)) "+ Age_clean" else "",
                                          if("Gender_label" %in% names(df)) "+ Gender_label" else ""))
    model_cent <- lm(form_cent, data = base_data)
    R2_cent <- summary(model_cent)$r.squared
    # Compile R2 results
    r2_table <- data.frame(
      Model = c("Transcending_only", "Centering_only", "Both_facets"),
      R_squared = c(R2_trans, R2_cent, R2_full),
      stringsAsFactors = FALSE
    )
    r2_table$R_squared <- sprintf("%.4f", r2_table$R_squared)
    write_dual(r2_table, paste0("Regression_", outcome, "_R2Comparison"))
  }
}


# ---- 7. Mediation analysis via bootstrapping ----

# Helper to run mediation bootstrap for a given X, M, Y (with covariates)
bootstrap_mediation <- function(X, M, Y, covars = NULL, n_boot = 5000){
  # Prepare data frame with complete cases on X, M, Y, and covariates
  vars <- c(X, M, Y, covars)
  med_data <- df[complete.cases(df[, vars]), vars, drop=FALSE]
  n <- nrow(med_data)
  if(n < 5){
    warning("Not enough data for mediation of ", X, " via ", M, " to ", Y)
    return(NULL)
  }
  # Fit models: M ~ X + covars, Y ~ X + M + covars, and Y ~ X + covars (for total effect)
  form_a <- as.formula(paste(M, "~", X, if(length(covars)) paste("+", paste(covars, collapse="+"))))
  form_b <- as.formula(paste(Y, "~", X, "+", M, if(length(covars)) paste("+", paste(covars, collapse="+"))))
  form_c <- as.formula(paste(Y, "~", X, if(length(covars)) paste("+", paste(covars, collapse="+"))))
  fit_a <- lm(form_a, data = med_data)
  fit_b <- lm(form_b, data = med_data)
  fit_c <- lm(form_c, data = med_data)
  
  # Extract coefficients
  a_coef <- coef(fit_a)[X]       # effect of X on M
  b_coef <- coef(fit_b)[M]       # effect of M on Y (X and covars controlled)
  c_prime <- coef(fit_b)[X]      # direct effect of X on Y controlling M
  c_total <- coef(fit_c)[X]      # total effect of X on Y (no M)
  
  # Bootstrap the indirect effect a*b
  set.seed(123)
  ab_vals <- numeric(n_boot)
  idx <- 1:n  # indices of data
  for(i in 1:n_boot){
    samp <- sample(idx, size = n, replace = TRUE)
    d <- med_data[samp, , drop=FALSE]
    # compute a and b on bootstrap sample
    fit_a_s <- lm(form_a, data = d)
    fit_b_s <- lm(form_b, data = d)
    a_s <- tryCatch(coef(fit_a_s)[X], error=function(e) NA)
    b_s <- tryCatch(coef(fit_b_s)[M], error=function(e) NA)
    ab_vals[i] <- a_s * b_s
  }
  ab_vals <- ab_vals[is.finite(ab_vals)]
  # Bootstrap CI for indirect
  ci <- quantile(ab_vals, c(0.025, 0.975), na.rm = TRUE)
  indirect_mean <- mean(ab_vals, na.rm = TRUE)
  prop_med <- (a_coef * b_coef) / c_total  # proportion mediated (using point estimates)
  
  # Return results in a list
  return(list(a = a_coef, b = b_coef, total_effect = c_total, direct_effect = c_prime,
              indirect_effect = a_coef * b_coef,  # point estimate
              indirect_boot_mean = indirect_mean,
              CI_low = ci[1], CI_high = ci[2],
              prop_mediated = prop_med,
              N = n))
}

# Define mediator-outcome combinations to test
mediators <- c("CFQ_Total", "ATQ_Believability")
outcomes_med <- c("DASS_Total_x2", "WEMWBS_Total")
covariates <- c()
if("Age_clean" %in% names(df)) covariates <- c(covariates, "Age_clean")
if("Gender_label" %in% names(df)) covariates <- c(covariates, "Gender_label")

med_results <- data.frame(Predictor=character(), Mediator=character(), Outcome=character(),
                          Indirect_Effect=numeric(), CI_lower=numeric(), CI_upper=numeric(),
                          Total_Effect=numeric(), Direct_Effect=numeric(), Prop_Mediated=numeric(),
                          stringsAsFactors=FALSE)
for(M in mediators){
  for(Y in outcomes_med){
    res <- bootstrap_mediation(X="SACS_Total", M=M, Y=Y, covars=covariates, n_boot=5000)
    if(!is.null(res)){
      med_results <- rbind(med_results, data.frame(
        Predictor = "SACS_Total",
        Mediator  = M,
        Outcome   = Y,
        Indirect_Effect = res$indirect_effect,
        CI_lower  = res$CI_low,
        CI_upper  = res$CI_high,
        Total_Effect = res$total_effect,
        Direct_Effect= res$direct_effect,
        Prop_Mediated = res$prop_mediated,
        stringsAsFactors=FALSE
      ))
    }
  }
}
# Format numeric results
num_cols <- c("Indirect_Effect","CI_lower","CI_upper","Total_Effect","Direct_Effect","Prop_Mediated")
med_results[num_cols] <- lapply(med_results[num_cols], function(x) sprintf("%.4f", as.numeric(x)))
write_dual(med_results, "Mediation_Indirect_Effects")

# ---- 8. Multiple Imputation for sensitivity ----

# Define variables to impute (exclude redundant composites like SACS_Total, ATQ_Combined, DASS_Total_x2 to avoid linear dependencies)
vars_to_impute <- c("SACS_Centering","SACS_Transcending",
                    "CFQ_Total",
                    "ATQ_Frequency","ATQ_Believability",
                    "WEMWBS_Total",
                    "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2",
                    "Age_clean","Gender_label")
# Ensure all are present
vars_to_impute <- vars_to_impute[vars_to_impute %in% names(df)]

# Initialize mice imputation (5 multiple imputations)
imp <- mice::mice(df[ , vars_to_impute], m = 5, seed = 100, printFlag = FALSE)

# Create one completed dataset for correlation analysis (using first imputation as representative)
complete_data1 <- mice::complete(imp, action = 1)

# Ensure score_vars are present in imputed dataset
score_vars_imp <- intersect(score_vars, names(complete_data1))

# Recompute correlation matrix with imputed data (full N)
cor_matrix_imp <- cor(complete_data1[ , score_vars_imp], use = "pairwise.complete.obs", method = "pearson")
cor_matrix_imp_df <- data.frame(Variable = rownames(cor_matrix_imp), cor_matrix_imp, check.names = FALSE)
write_dual(cor_matrix_imp_df, "Correlation_Matrix_MI")

# Re-run facet regressions with pooled MI data
for(outcome in outcomes_H2){
  formula_H2 <- as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending",
                                 if("Age_clean" %in% names(df)) "+ Age_clean" else "",
                                 if("Gender_label" %in% names(df)) "+ Gender_label" else ""))
  # Fit model on each imputed set and pool
  fit_imp <- try(mice::with(imp, lm(formula_H2)), silent = TRUE)
  if(inherits(fit_imp, "try-error")){
    next  # skip if model cannot be fitted (e.g., no variability in outcome in some imp)
  }
  pooled <- mice::pool(fit_imp)
  pool_summary <- summary(pooled)  # pooled coefficients
  # Prepare output dataframe
  coef_df <- data.frame(Variable = pool_summary$term,
                        Estimate = pool_summary$estimate,
                        Std_Error = pool_summary$std.error,
                        df = pool_summary$df,
                        p_value = pool_summary$p.value,
                        stringsAsFactors = FALSE)
  # Format numeric columns
  coef_df$Estimate <- sprintf("%.4f", coef_df$Estimate)
  coef_df$Std_Error<- sprintf("%.4f", coef_df$Std_Error)
  coef_df$p_value  <- sprintf("%.4f", coef_df$p_value)
  # Save
  write_dual(coef_df, paste0("Regression_", outcome, "_facets_MI_coeffs"))
}


