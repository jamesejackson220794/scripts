## tVSc6.R — SACS structure, ΔR², Wald, MIMIC, WEMWBS MIMIC, HTMT (ordinal-safe)
suppressPackageStartupMessages({
  library(lavaan)       # 0.6-19+
  library(semTools)     # for lavCor (polychorics)
  library(sandwich)     # robust vcov for OLS Wald
  library(lmtest)       # coeftest
  library(car)          # linearHypothesis
  library(dplyr)
  library(stringr)
})

## -----------------------------
## 0) EXPECTED DATA & VARIABLES
## -----------------------------
## Assumes 'dat' exists in the parent environment after source()
## Required columns:
##   SACS_1 ... SACS_10 (ordered Likert), Age_clean, Gender_label (factor with Woman/Non-Binary/Man),
##   WEMWBS_Total, DASS_Total_x2, DASS_Depression_x2, DASS_Anxiety_x2, DASS_Stress_x2

stopifnot(exists("dat"))
req <- c(paste0("SACS_",1:10),"Age_clean","Gender_label",
         "WEMWBS_Total","DASS_Total_x2","DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2")
miss <- setdiff(req, names(dat))
if(length(miss)){
  stop(sprintf("Missing required variables in 'dat': %s", paste(miss, collapse=", ")))
}

## ensure factor/coding
dat <- dat %>%
  mutate(Gender_label = factor(Gender_label))

## item sets
c_items <- c("SACS_1","SACS_2","SACS_5","SACS_6")
t_items <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_10")
all_items <- c(c_items,t_items)
all10_items <- c("SACS_1","SACS_2","SACS_3","SACS_4","SACS_5","SACS_6","SACS_7","SACS_8","SACS_9","SACS_10")

ordered_items10 <- all10_items
ordered_items_trim <- all_items

## -----------------------------
## 1) CFA: 1-factor vs 2-factor (10 items)
## -----------------------------
model_1f_10 <- '
  SACSg =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'
model_2f_10 <- '
  SACSc =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  SACSt =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
  SACSc ~~ SACSt
'

fit_1f_10 <- cfa(model_1f_10, data=dat, ordered=ordered_items10, estimator="WLSMV", std.lv=TRUE)
fit_2f_10 <- cfa(model_2f_10, data=dat, ordered=ordered_items10, estimator="WLSMV", std.lv=TRUE)

cat("\n--- 1-factor (10 items) ---\n")
print(fitMeasures(fit_1f_10, c("cfi","tli","rmsea","srmr","bic")))
cat("\n--- FULL 2-factor (10 items) ---\n")
print(fitMeasures(fit_2f_10, c("cfi","tli","rmsea","srmr","bic")))
r_full <- lavInspect(fit_2f_10, "cor.lv")[1,2]
cat("Latent r(SACSc,SACSt) full = ", round(r_full,3), "\n")

## delta fit indices
fm1 <- fitMeasures(fit_1f_10, c("cfi","tli","rmsea","bic"))
fm2 <- fitMeasures(fit_2f_10, c("cfi","tli","rmsea","bic"))
delta <- c(dCFI.cfi = fm2["cfi"]-fm1["cfi"],
           dTLI.tli = fm2["tli"]-fm1["tli"],
           dRMSEA.rmsea = fm2["rmsea"]-fm1["rmsea"],
           dBIC.bic = fm2["bic"]-fm1["bic"])
cat("\nDelta (2f - 1f):\n"); print(delta)

## --- Heywood check helper (no 'heywood' in lavInspect) ---
heywood_check <- function(fit){
  PE <- parameterEstimates(fit, standardized=TRUE)
  # Negative variances (observed & latent)
  ov_var <- subset(PE, op=="~~" & lhs %in% all10_items & rhs==lhs)
  lv_var <- subset(PE, op=="~~" & lhs %in% c("SACSc","SACSt","SACSg") & rhs==lhs)
  any_ov_neg_raw <- any(ov_var$est < 0, na.rm=TRUE)
  any_lv_neg_raw <- any(lv_var$est < 0, na.rm=TRUE)
  ov_var_std <- subset(PE, op=="~~" & lhs %in% all10_items & rhs==lhs)
  lv_var_std <- subset(PE, op=="~~" & lhs %in% c("SACSc","SACSt","SACSg") & rhs==lhs)
  any_ov_neg_std <- any(ov_var_std$std.all < 0, na.rm=TRUE)
  any_lv_neg_std <- any(lv_var_std$std.all < 0, na.rm=TRUE)
  # Loadings >1 (offenders)
  load <- subset(PE, op=="=~")
  ov_off <- subset(load, std.all > 1 | std.all < -1, select=c("lhs","est","std.all"))
  lv_off <- subset(lv_var_std, std.all > 1 | std.all < -1, select=c("lhs","est","std.all"))
  list(any_ov_neg_raw=any_ov_neg_raw,
       any_ov_neg_std=any_ov_neg_std,
       any_lv_neg_raw=any_lv_neg_raw,
       any_lv_neg_std=any_lv_neg_std,
       ov_offenders=ov_off,
       lv_offenders=lv_off)
}
cat("\nHeywood check (FULL 2f):\n"); print(heywood_check(fit_2f_10))

## --------------------------------
## 2) Trim: drop item 9, re-fit 2f
## --------------------------------
model_2f_trim <- '
  SACSc =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  SACSt =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
  SACSc ~~ SACSt
'
fit_2f_trim <- cfa(model_2f_trim, data=dat, ordered=ordered_items_trim, estimator="WLSMV", std.lv=TRUE)

cat("\n--- TRIM 2-factor (drop SACS_9) ---\n")
print(fitMeasures(fit_2f_trim, c("cfi","tli","rmsea","srmr","bic")))
r_trim <- lavInspect(fit_2f_trim, "cor.lv")[1,2]
cat("Latent r(SACSc,SACSt) trim = ", round(r_trim,3), "\n")
cat("\nHeywood check (TRIM 2f):\n"); print(heywood_check(fit_2f_trim))

## --------------------------------
## 3) HTMT from polychorics (2 blocks)
## --------------------------------
htmt_from_cor <- function(S, setA, setB){
  # numerator: mean |cor| across A-B pairs
  num <- mean(abs(S[setA, setB]), na.rm=TRUE)
  # denom: geometric mean of mean |cor| within A and within B (off-diagonal)
  within_mean <- function(ix){
    M <- abs(S[ix, ix])
    if(length(ix) < 2) return(NA_real_)
    m <- M[upper.tri(M)]
    mean(m, na.rm=TRUE)
  }
  den <- sqrt(within_mean(setA) * within_mean(setB))
  num/den
}
# polychoric correlation among trimmed items
S_poly <- lavCor(dat[,ordered_items_trim], ordered=ordered_items_trim, output="cor")
htmt_CT <- htmt_from_cor(S_poly, setA=c_items, setB=t_items)
cat("\nHTMT (C vs T, trim; polychoric): ", round(htmt_CT,3), "\n")

## --------------------------------
## 4) Manifest ΔR² + robust Wald
## --------------------------------
mk_formula <- function(y, full=TRUE){
  if(full){
    as.formula(sprintf("%s ~ Age_clean + Gender_label + SACS_Centering_trim + SACS_Transcending_trim", y))
  } else {
    as.formula(sprintf("%s ~ Age_clean + Gender_label + SACS_Centering_trim", y))
  }
}
## Create trimmed composites (means) for OLS
dat$SACS_Centering_trim   <- rowMeans(dat[,c_items], na.rm=TRUE)
dat$SACS_Transcending_trim<- rowMeans(dat[,t_items], na.rm=TRUE)

run_ols_block <- function(y){
  form_base <- mk_formula(y, full=FALSE)
  form_full <- mk_formula(y, full=TRUE)
  m1 <- lm(form_base, data=dat)
  m2 <- lm(form_full, data=dat)
  cat("\n--- OLS (", y, "): base ---\n", sep="")
  print(summary(m1))
  cat("\n--- OLS (", y, "): full ---\n", sep="")
  print(summary(m2))
  cat("\n\nΔR² F-test:\n"); print(anova(m1,m2))
  # Robust Wald test: βC = βT
  Vrob <- vcovHC(m2, type="HC3")
  cat("\nWald test (robust) of β_SACSc = β_SACSt:\n\n")
  print(car::linearHypothesis(m2, "SACS_Centering_trim - SACS_Transcending_trim = 0", vcov.=Vrob))
}

run_ols_block("DASS_Depression_x2")
run_ols_block("DASS_Anxiety_x2")
run_ols_block("DASS_Stress_x2")

## --------------------------------
## 5) MIMIC models (trim): DASS & WEMWBS -> latent C/T
## --------------------------------
mimic_template <- function(exog){
  sprintf('
    SACSc =~ %s
    SACSt =~ %s
    SACSc ~ %s + Age_clean + Gender_label
    SACSt ~ %s + Age_clean + Gender_label
    SACSc ~~ SACSt
  ',
          paste(c_items, collapse=" + "),
          paste(t_items, collapse=" + "),
          exog, exog)
}

fit_mimic <- function(exog, label){
  mod <- mimic_template(exog)
  fit <- sem(mod, data=dat,
             ordered=ordered_items_trim, estimator="WLSMV",
             std.lv=TRUE, se="robust.sem")
  fm <- fitMeasures(fit, c("cfi","tli","rmsea","srmr","bic"))
  cat("\n--- MIMIC (", exog, " -> latent SACSc/SACSt; WLSMV, TRIM) ---\n", sep="")
  print(summary(fit, fit.measures=TRUE, standardized=TRUE))
  invisible(fit)
}

fit_mimic("WEMWBS_Total", "WEMWBS")
fit_mimic("DASS_Total_x2", "DASS_Total")
fit_mimic("DASS_Depression_x2", "DASS_Depression")
fit_mimic("DASS_Anxiety_x2", "DASS_Anxiety")

## Note: DASS_Stress_x2 can be added similarly if desired:
# fit_mimic("DASS_Stress_x2", "DASS_Stress")

## --------------------------------
## 6) Diagnostics: VCOV eigenvalues (informative, not fatal)
## --------------------------------
safe_smallest_eig <- function(fit){
  V <- try(lavInspect(fit, "vcov"), silent=TRUE)
  if(inherits(V,"try-error")) return(NA_real_)
  min(eigen(V, symmetric=TRUE, only.values=TRUE)$values)
}
eigs <- c(
  one_f_10 = safe_smallest_eig(fit_1f_10),
  two_f_10 = safe_smallest_eig(fit_2f_10),
  two_f_trim = safe_smallest_eig(fit_2f_trim)
)
cat("\nSmallest VCOV eigenvalues (informative):\n"); print(eigs)
if(any(!is.na(eigs) & eigs < 0)){
  cat("Note: near-zero/negative eigenvalues are common with WLSMV+thresholds at small N; use robust SEs and inspect parameters.\n")
}
