# file: scripts/sac_analysis_fixed.R

# ---------------- Boot ----------------
options(stringsAsFactors = FALSE, scipen = 999, digits = 3, warn = 1)

log_step <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste(..., collapse = "")))
  flush.console()
}

install_and_load <- function(pkgs) {
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))
}

# Required packages
install_and_load(c(
  "tidyverse","ggplot2","lavaan","semTools","psych","ppcor","boot","car",
  "lmtest","sandwich","mgcv","mice","MVN","interactions","effectsize","domir",
  "BaylorEdPsych"
))

set.seed(123)

# --------------- Data load ---------------
# Why: Make sourcing predictable whether 'data' exists or not.
DATA_PATH <- Sys.getenv("SAC_DATA_CSV", unset = "your_data_file.csv")
if (!exists("data", inherits = TRUE)) {
  if (file.exists(DATA_PATH)) {
    log_step("Reading data from: ", DATA_PATH)
    data <- read.csv(DATA_PATH, check.names = TRUE)
  } else {
    stop("No data found. Supply a data.frame named `data` or set SAC_DATA_CSV to a CSV path.")
  }
}
stopifnot(is.data.frame(data))

# --------------- Name normalization ---------------
# Why: Support both DASS_# and DASS.21_#; ATQ believability '_Believability' vs 'b'.
cn <- colnames(data)

# DASS prefix resolve
dass_prefix <- if (any(startsWith(cn, "DASS_"))) "DASS_" else if (any(startsWith(cn, "DASS.21_"))) "DASS.21_" else NULL
if (is.null(dass_prefix)) warning("No DASS columns found. Downstream steps may skip.")

# ATQ believability columns: prefer *_Believability, else ATQ_#b
has_atq_bel <- any(grepl("^ATQ_[0-9]+_Believability$", cn))
has_atq_b   <- any(grepl("^ATQ_[0-9]+b$", cn))

# Age variable
if (!"Age_clean" %in% cn) {
  if ("Age" %in% cn) {
    data$Age_clean <- suppressWarnings(as.numeric(data$Age))
  } else {
    data$Age_clean <- NA_real_
  }
}

# --------------- PART 1: Data quality ---------------
log_step("=== PART 1: DATA QUALITY ===")

# Missingness summary
missing_summary <- data.frame(
  Variable   = names(data),
  Missing_N  = colSums(is.na(data)),
  Missing_Pct= round(colSums(is.na(data))/nrow(data)*100, 2)
) |> dplyr::filter(Missing_Pct > 0)
print(missing_summary)

# --------------- Scale scoring ---------------
# SACS (drop 9)
centering_items    <- c("SACS_1","SACS_2","SACS_5","SACS_6")
transcending_items <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_10")

# Validate presence softly
needed <- c(centering_items, transcending_items, paste0("CFQ_",1:7), paste0("WEMWBS_",1:14))
missing_needed <- setdiff(needed, colnames(data))
if (length(missing_needed)) log_step("Note: missing expected items: ", paste(missing_needed, collapse = ", "))

data$SACS_Centering    <- rowMeans(dplyr::select(data, dplyr::any_of(centering_items)), na.rm = TRUE)
data$SACS_Transcending <- rowMeans(dplyr::select(data, dplyr::any_of(transcending_items)), na.rm = TRUE)
data$SACS_Total        <- rowMeans(dplyr::select(data, dplyr::any_of(c(centering_items, transcending_items))), na.rm = TRUE)

# CFQ
data$CFQ_Total <- rowMeans(dplyr::select(data, dplyr::any_of(paste0("CFQ_",1:7))), na.rm = TRUE)

# ATQ believability and frequency
if (has_atq_bel) {
  data$ATQ_Believability <- rowMeans(dplyr::select(data, dplyr::any_of(paste0("ATQ_",1:15,"_Believability"))), na.rm = TRUE)
  data$ATQ_Frequency     <- rowMeans(dplyr::select(data, dplyr::any_of(paste0("ATQ_",1:15,"_Frequency"))), na.rm = TRUE)
} else if (has_atq_b) {
  data$ATQ_Believability <- rowMeans(dplyr::select(data, dplyr::any_of(paste0("ATQ_",1:15,"b"))), na.rm = TRUE)
} else {
  data$ATQ_Believability <- NA_real_
}

# WEMWBS
data$WEMWBS_Total <- rowSums(dplyr::select(data, dplyr::any_of(paste0("WEMWBS_",1:14))), na.rm = TRUE)

# DASS x2
d_dep <- c(3,5,10,13,16,17,21); d_anx <- c(2,4,7,9,15,19,20); d_str <- c(1,6,8,11,12,14,18)
mk <- function(idx) paste0(dass_prefix, idx)
if (!is.null(dass_prefix)) {
  data$DASS_Depression_x2 <- rowSums(dplyr::select(data, dplyr::any_of(mk(d_dep))), na.rm = FALSE) * 2
  data$DASS_Anxiety_x2    <- rowSums(dplyr::select(data, dplyr::any_of(mk(d_anx))), na.rm = FALSE) * 2
  data$DASS_Stress_x2     <- rowSums(dplyr::select(data, dplyr::any_of(mk(d_str))), na.rm = FALSE) * 2
  data$DASS_Total_x2      <- with(data, DASS_Depression_x2 + DASS_Anxiety_x2 + DASS_Stress_x2)
}

usable_n <- sum(stats::complete.cases(dplyr::select(data, dplyr::any_of(
  c("SACS_Centering","SACS_Transcending","CFQ_Total","ATQ_Believability","WEMWBS_Total","DASS_Total_x2")
)))))
log_step("Usable N for primary analyses: ", usable_n)

# Descriptives
key_vars <- c("SACS_Centering","SACS_Transcending","SACS_Total","CFQ_Total","ATQ_Believability",
              "WEMWBS_Total","DASS_Total_x2","DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2")
key_vars <- intersect(key_vars, names(data))
desc_stats <- psych::describe(data[, key_vars, drop = FALSE])
print(round(desc_stats[, c("n","mean","sd","min","max","skew","kurtosis")], 2))

# Robust outliers via MAD and winsorize
winsorize <- function(x, probs = c(0.01, 0.99)) {
  qu <- stats::quantile(x, probs = probs, na.rm = TRUE)
  x[x < qu[1]] <- qu[1]; x[x > qu[2]] <- qu[2]; x
}
for (v in key_vars) data[[paste0(v,"_wins")]] <- winsorize(data[[v]])

# Little’s MCAR with BaylorEdPsych
mcar_df <- dplyr::select(data, dplyr::any_of(c(paste0("SACS_",1:10), paste0("CFQ_",1:7))))
if (ncol(mcar_df) > 0) {
  log_step("Little's MCAR test (BaylorEdPsych::LittleMCAR)")
  print(try(BaylorEdPsych::LittleMCAR(mcar_df), silent = TRUE))
}

# MVN via Mardia
complete_key <- stats::na.omit(data[, key_vars, drop = FALSE])
if (nrow(complete_key) > 4) {
  mvn_result <- MVN::mvn(complete_key, mvnTest = "mardia")
  print(mvn_result$multivariateNormality)
}

# --------------- PART 2: CFA of SACS ---------------
log_step("=== PART 2: CFA ===")

model_1factor <- '
  SAC =~ SACS_1 + SACS_2 + SACS_5 + SACS_6 + SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
'

model_2factor <- '
  Centering =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  Transcending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
'

fit_1factor <- lavaan::cfa(model_1factor, data = data,
                           estimator = "WLSMV", missing = "pairwise",
                           std.lv = TRUE, ordered = TRUE)

fit_2factor <- lavaan::cfa(model_2factor, data = data,
                           estimator = "WLSMV", missing = "pairwise",
                           std.lv = TRUE, ordered = TRUE)

fit_comparison <- data.frame(
  Model = c("One-Factor","Two-Factor"),
  ChiSq = c(lavaan::fitMeasures(fit_1factor)["chisq"], lavaan::fitMeasures(fit_2factor)["chisq"]),
  df    = c(lavaan::fitMeasures(fit_1factor)["df"],    lavaan::fitMeasures(fit_2factor)["df"]),
  CFI   = c(lavaan::fitMeasures(fit_1factor)["cfi"],   lavaan::fitMeasures(fit_2factor)["cfi"]),
  TLI   = c(lavaan::fitMeasures(fit_1factor)["tli"],   lavaan::fitMeasures(fit_2factor)["tli"]),
  RMSEA = c(lavaan::fitMeasures(fit_1factor)["rmsea"], lavaan::fitMeasures(fit_2factor)["rmsea"]),
  SRMR  = c(lavaan::fitMeasures(fit_1factor)["srmr"],  lavaan::fitMeasures(fit_2factor)["srmr"]),
  BIC   = c(lavaan::fitMeasures(fit_1factor)["bic"],   lavaan::fitMeasures(fit_2factor)["bic"])
)
print(round(fit_comparison, 3))

delta_cfi   <- fit_comparison$CFI[2]   - fit_comparison$CFI[1]
delta_rmsea <- fit_comparison$RMSEA[2] - fit_comparison$RMSEA[1]
delta_bic   <- fit_comparison$BIC[1]   - fit_comparison$BIC[2]  # + favors 2F

cat("\nΔCFI:", round(delta_cfi,3),
    " ΔRMSEA:", round(delta_rmsea,3),
    " ΔBIC:", round(delta_bic,1), "\n")  # :contentReference[oaicite:5]{index=5}

# HTMT
htmt_result <- semTools::htmt(model_2factor, data = data)
cat("HTMT(Centering,Transcending):", round(htmt_result[1,1], 3), "\n")  # :contentReference[oaicite:6]{index=6}

# Bifactor
model_bifactor <- '
  g      =~ SACS_1 + SACS_2 + SACS_5 + SACS_6 + SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
  s_cent =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  s_trans=~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
  g ~~ 0*s_cent
  g ~~ 0*s_trans
  s_cent ~~ 0*s_trans
'
fit_bifactor <- lavaan::cfa(model_bifactor, data = data, estimator = "WLSMV",
                            missing = "pairwise", std.lv = TRUE, orthogonal = TRUE, ordered = TRUE)

# Reliability: alpha and omega_h for Transcending
alpha_centering    <- psych::alpha(dplyr::select(data, dplyr::any_of(centering_items)))
alpha_transcending <- psych::alpha(dplyr::select(data, dplyr::any_of(transcending_items)))
cat("\nAlpha Centering:", round(alpha_centering$total$raw_alpha,3),
    " Alpha Transcending:", round(alpha_transcending$total$raw_alpha,3), "\n")  # :contentReference[oaicite:7]{index=7}
omega_trans <- try(psych::omegah(dplyr::select(data, dplyr::any_of(transcending_items)), nfactors = 1), silent = TRUE)
if (!inherits(omega_trans,"try-error")) cat("Omega_h Transcending:", round(omega_trans$omega_h,3), "\n")  # :contentReference[oaicite:8]{index=8}

# --------------- PART 3: Correlations & suppression ---------------
log_step("=== PART 3: CORRELATIONS & SUPPRESSION ===")
outcome_vars   <- intersect(c("WEMWBS_Total","DASS_Total_x2","DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2"), names(data))
predictor_vars <- intersect(c("SACS_Centering","SACS_Transcending","CFQ_Total","ATQ_Believability"), names(data))

if (length(c(outcome_vars,predictor_vars)) > 1) {
  cor_matrix <- stats::cor(dplyr::select(data, dplyr::any_of(c(outcome_vars,predictor_vars))),
                           use = "pairwise.complete.obs")
  print(round(cor_matrix[outcome_vars, intersect(c("SACS_Centering","SACS_Transcending"), colnames(cor_matrix)), drop=FALSE], 3))
}

detect_suppression <- function(outcome, data) {
  sub <- data[, c(outcome,"SACS_Centering","SACS_Transcending"), drop = FALSE]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (nrow(sub) < 5) return(NULL)
  
  r_cent  <- stats::cor(sub$SACS_Centering,  sub[[outcome]])
  r_trans <- stats::cor(sub$SACS_Transcending,sub[[outcome]])
  r_ct    <- stats::cor(sub$SACS_Centering,  sub$SACS_Transcending)
  
  m <- stats::lm(stats::as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending")), data = sub)
  # std betas
  sdy <- stats::sd(sub[[outcome]]); sdc <- stats::sd(sub$SACS_Centering); sdt <- stats::sd(sub$SACS_Transcending)
  beta <- stats::coef(m)[c("SACS_Centering","SACS_Transcending")]
  std_beta <- c(Centering = beta[1] * sdc/sdy, Transcending = beta[2] * sdt/sdy)
  
  # semi-partial (ppcor)
  sp <- ppcor::spcor(sub)$estimate
  sr_cent  <- sp["SACS_Centering", outcome]
  sr_trans <- sp["SACS_Transcending", outcome]
  
  tibble::tibble(
    Outcome = outcome,
    r_Centering = r_cent, r_Transcending = r_trans, r_Cent_Trans = r_ct,
    Beta_Centering = std_beta["Centering"], Beta_Transcending = std_beta["Transcending"],
    sr_Centering = sr_cent, sr_Transcending = sr_trans,
    Sign_Reversal_Trans = sign(r_trans) != sign(std_beta["Transcending"]),
    Suppression_Trans   = abs(sr_trans) > abs(r_trans) * 1.1
  )
}
suppression_results <- purrr::compact(lapply(outcome_vars, detect_suppression, data = data)) |> dplyr::bind_rows()
if (nrow(suppression_results)) print(round(suppression_results, 3))

# --------------- PART 4: Non-linearity (GAM) ---------------
log_step("=== PART 4: NON-LINEARITY (GAM) ===")
if (all(c("DASS_Total_x2","SACS_Transcending","SACS_Centering","CFQ_Total","Age_clean") %in% names(data))) {
  gam_dass <- mgcv::gam(DASS_Total_x2 ~ s(SACS_Transcending, k = 7) + SACS_Centering + CFQ_Total + Age_clean,
                        data = data, method = "REML", na.action = na.exclude)
  print(summary(gam_dass))
  edf_trans <- summary(gam_dass)$s.table["s(SACS_Transcending)", "edf"]
  pval_trans<- summary(gam_dass)$s.table["s(SACS_Transcending)", "p-value"]
  cat("EDF:", round(edf_trans,2), " p:", round(pval_trans,4), "\n")  # :contentReference[oaicite:9]{index=9}
  
  if (!is.na(edf_trans) && edf_trans > 1.5) {
    data$SACS_Trans_c  <- scale(data$SACS_Transcending, scale = FALSE)
    data$SACS_Trans_c2 <- data$SACS_Trans_c^2
    m_lin  <- lm(DASS_Total_x2 ~ SACS_Trans_c + SACS_Centering + CFQ_Total + Age_clean, data = data)
    m_quad <- lm(DASS_Total_x2 ~ SACS_Trans_c + SACS_Trans_c2 + SACS_Centering + CFQ_Total + Age_clean, data = data)
    print(anova(m_lin, m_quad))  # :contentReference[oaicite:10]{index=10}
  }
}

# --------------- PART 5: Incremental validity ---------------
log_step("=== PART 5: INCREMENTAL VALIDITY ===")
test_incremental_validity <- function(outcome, data) {
  keep <- stats::complete.cases(dplyr::select(data, dplyr::any_of(c(outcome,"SACS_Centering","SACS_Transcending","CFQ_Total","Age_clean"))))
  d <- data[keep, , drop = FALSE]; if (nrow(d) < 10) return(NULL)
  
  m1 <- lm(stats::as.formula(paste(outcome, "~ SACS_Centering + CFQ_Total + Age_clean")), data = d)
  m2 <- lm(stats::as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending + CFQ_Total + Age_clean")), data = d)
  an <- anova(m1, m2)
  r2_base <- summary(m1)$r.squared; r2_full <- summary(m2)$r.squared
  delta_r2 <- r2_full - r2_base; f2 <- delta_r2 / max(1 - r2_full, .Machine$double.eps)
  
  boot_delta_r2 <- function(data, idx) {
    dd <- data[idx, , drop = FALSE]
    a <- lm(stats::as.formula(paste(outcome, "~ SACS_Centering + CFQ_Total + Age_clean")), data = dd)
    b <- lm(stats::as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending + CFQ_Total + Age_clean")), data = dd)
    summary(b)$r.squared - summary(a)$r.squared
  }
  bt <- boot::boot(d, boot_delta_r2, R = 2000)
  ci <- try(boot::boot.ci(bt, type = "bca", conf = 0.95)$bca[4:5], silent = TRUE)
  if (inherits(ci,"try-error")) ci <- c(NA_real_, NA_real_)
  
  list(delta_r2 = delta_r2, f2 = f2, F = an$F[2], p = an$`Pr(>F)`[2],
       ci_lo = ci[1], ci_hi = ci[2],
       vif = car::vif(m2))
}
iv <- purrr::map(outcome_vars, ~ test_incremental_validity(.x, data)) |> stats::setNames(outcome_vars)
print(purrr::map(iv, ~ if (is.null(.x)) NULL else unlist(.x)))  # :contentReference[oaicite:11]{index=11}

# Dominance analysis (guarded)
if (requireNamespace("domir", quietly = TRUE)) {
  for (outcome in outcome_vars) {
    keep <- stats::complete.cases(dplyr::select(data, dplyr::any_of(c(outcome,"SACS_Centering","SACS_Transcending","CFQ_Total"))))
    d <- data[keep, , drop = FALSE]; if (nrow(d) < 10) next
    dom <- domir::domin(stats::as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending + CFQ_Total")),
                        reg = "lm", data = d)
    cat("\nDominance for", outcome, "\n"); print(round(dom$General_Dominance, 3))
  }
}

# --------------- PART 6: Moderation ---------------
log_step("=== PART 6: MODERATION ===")
if (all(c("DASS_Total_x2","SACS_Centering","SACS_Transcending","CFQ_Total") %in% names(data))) {
  data$Cent_c  <- scale(data$SACS_Centering,    scale = FALSE)
  data$Trans_c <- scale(data$SACS_Transcending, scale = FALSE)
  
  mod_model <- lm(DASS_Total_x2 ~ Cent_c + Trans_c + Cent_c:Trans_c + CFQ_Total + Age_clean, data = data)
  print(summary(mod_model))
  cat("Interaction Cent×Trans:", coef(summary(mod_model))["Cent_c:Trans_c","Pr(>|t|)"], "\n")  # :contentReference[oaicite:12]{index=12}
  if (coef(summary(mod_model))["Cent_c:Trans_c","Pr(>|t|)"] < 0.05) {
    interactions::sim_slopes(mod_model, pred = Trans_c, modx = Cent_c, jnplot = TRUE)
    suppressWarnings(print(interactions::johnson_neyman(mod_model, pred = Trans_c, modx = Cent_c)))
  }
  
  data$CFQ_c <- scale(data$CFQ_Total, scale = FALSE)
  mod_cfq <- lm(DASS_Total_x2 ~ Trans_c + CFQ_c + Trans_c:CFQ_c + Cent_c + Age_clean, data = data)
  print(summary(mod_cfq))
}

# --------------- PART 7: Mediation ---------------
log_step("=== PART 7: MEDIATION ===")
if (all(c("SACS_Centering","SACS_Transcending","CFQ_Total","DASS_Total_x2") %in% names(data))) {
  mediation_cfq <- '
    CFQ_Total ~ a1*SACS_Centering + a2*SACS_Transcending + Age_clean
    DASS_Total_x2 ~ b*CFQ_Total + cp1*SACS_Centering + cp2*SACS_Transcending + Age_clean
    SACS_Centering ~~ SACS_Transcending
    indirect_cent := a1*b
    indirect_trans := a2*b
    total_cent := cp1 + (a1*b)
    total_trans := cp2 + (a2*b)
    diff_indirect := indirect_cent - indirect_trans
  '
  fit_med_cfq <- lavaan::sem(mediation_cfq, data = data, missing = "fiml",
                             estimator = "MLM", se = "bootstrap", bootstrap = 2000)
  cat("CFQ mediation fit: CFI=", round(lavaan::fitMeasures(fit_med_cfq)["cfi"],3),
      " RMSEA=", round(lavaan::fitMeasures(fit_med_cfq)["rmsea"],3), "\n")  # :contentReference[oaicite:13]{index=13}
}

if (all(c("SACS_Centering","SACS_Transcending","ATQ_Believability","DASS_Total_x2") %in% names(data))) {
  mediation_atq <- '
    ATQ_Believability ~ a1*SACS_Centering + a2*SACS_Transcending + Age_clean
    DASS_Total_x2 ~ b*ATQ_Believability + cp1*SACS_Centering + cp2*SACS_Transcending + Age_clean
    SACS_Centering ~~ SACS_Transcending
    indirect_cent := a1*b
    indirect_trans := a2*b
  '
  fit_med_atq <- lavaan::sem(mediation_atq, data = data, missing = "fiml",
                             estimator = "MLM", se = "bootstrap", bootstrap = 2000)
  cat("ATQ mediation fitted.\n")
}

if (all(c("SACS_Centering","SACS_Transcending","CFQ_Total","ATQ_Believability","DASS_Total_x2") %in% names(data))) {
  multi_med <- '
    CFQ_Total ~ a1_cfq*SACS_Centering + a2_cfq*SACS_Transcending + Age_clean
    ATQ_Believability ~ a1_atq*SACS_Centering + a2_atq*SACS_Transcending + Age_clean
    DASS_Total_x2 ~ b_cfq*CFQ_Total + b_atq*ATQ_Believability + cp1*SACS_Centering + cp2*SACS_Transcending + Age_clean
    SACS_Centering ~~ SACS_Transcending
    CFQ_Total ~~ ATQ_Believability
    ind_cent_cfq := a1_cfq * b_cfq
    ind_cent_atq := a1_atq * b_atq
    ind_trans_cfq := a2_cfq * b_cfq
    ind_trans_atq := a2_atq * b_atq
    total_ind_cent := ind_cent_cfq + ind_cent_atq
    total_ind_trans := ind_trans_cfq + ind_trans_atq
    cent_cfq_vs_atq := ind_cent_cfq - ind_cent_atq
    trans_cfq_vs_atq := ind_trans_cfq - ind_trans_atq
  '
  fit_multi_med <- lavaan::sem(multi_med, data = data, missing = "fiml",
                               estimator = "MLM", se = "bootstrap", bootstrap = 2000)
  param_multi <- lavaan::parameterEstimates(fit_multi_med, boot.ci.type = "bca.simple")
  print(param_multi[param_multi$label %in% c("ind_cent_cfq","ind_cent_atq","ind_trans_cfq","ind_trans_atq",
                                             "cent_cfq_vs_atq","trans_cfq_vs_atq"),
                    c("label","est","ci.lower","ci.upper","pvalue")])  # :contentReference[oaicite:14]{index=14}
}

# --------------- PART 9: Summary ---------------
log_step("=== SUMMARY ===")
cat("Usable N:", usable_n, "\n")
if (exists("delta_cfi")) {
  cat("SACS 2F vs 1F ΔCFI:", round(delta_cfi,3),
      " ΔRMSEA:", round(delta_rmsea,3),
      " ΔBIC:", round(delta_bic,1), "\n")  # :contentReference[oaicite:15]{index=15}
}
if (exists("htmt_result")) cat("HTMT:", round(htmt_result[1,1],3), "\n")
if (!inherits(omega_trans,"try-error")) cat("Omega_h Transcending:", round(omega_trans$omega_h,3), "\n")

if (length(iv)) {
  for (nm in names(iv)) if (!is.null(iv[[nm]]))
    cat("ΔR²", nm, ":", round(iv[[nm]]$delta_r2,3), "CI[", round(iv[[nm]]$ci_lo,3), ",", round(iv[[nm]]$ci_hi,3), "]\n")
}

log_step("Analysis complete.")
