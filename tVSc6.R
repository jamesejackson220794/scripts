## =================== Packages & options ===================
suppressPackageStartupMessages({
  library(lavaan)
  library(semTools)   # for htmt() & reliability()
  library(car)        # for linearHypothesis (robust Wald)
  library(sandwich)   # robust vcov for OLS (HC3)
})

options(lavaan.warn = "short")

## =================== Data expectations ====================
stopifnot(exists("dat"))

## Ensure gender dummies exist for MIMIC (ref = "Man")
if (!all(c("Gender_label__Non.Binary","Gender_label__Woman") %in% names(dat))) {
  if (!"Gender_label" %in% names(dat)) stop("dat$Gender_label not found.")
  dat$Gender_label <- factor(dat$Gender_label)
  if ("Man" %in% levels(dat$Gender_label)) dat$Gender_label <- relevel(dat$Gender_label, ref="Man")
  Xg <- model.matrix(~ Gender_label, data = dat)
  if ("Gender_labelNon-Binary" %in% colnames(Xg))
    dat$Gender_label__Non.Binary <- as.numeric(Xg[,"Gender_labelNon-Binary"])
  if ("Gender_labelWoman" %in% colnames(Xg))
    dat$Gender_label__Woman <- as.numeric(Xg[,"Gender_labelWoman"])
  # fallbacks for alternate contrast names
  if (!"Gender_label__Non.Binary" %in% names(dat) && "Gender_labelNon.Binary" %in% colnames(Xg))
    dat$Gender_label__Non.Binary <- as.numeric(Xg[,"Gender_labelNon.Binary"])
  if (!"Gender_label__Woman" %in% names(dat) && "Gender_label.Woman" %in% colnames(Xg))
    dat$Gender_label__Woman <- as.numeric(Xg[,"Gender_label.Woman"])
}

## =================== Item sets ============================
c_items_full   <- c("SACS_1","SACS_2","SACS_5","SACS_6")
t_items_full   <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")
ordered_full   <- c(c_items_full, t_items_full)

c_items_trim   <- c("SACS_1","SACS_2","SACS_5","SACS_6")
t_items_trim   <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_10") # drop SACS_9
ordered_trim   <- c(c_items_trim, t_items_trim)

## Coerce items to ordered
coerce_to_ordered <- function(dd, items) {
  for (v in items) {
    if (!is.ordered(dd[[v]])) {
      if (is.numeric(dd[[v]]) || is.integer(dd[[v]]) || is.factor(dd[[v]])) {
        dd[[v]] <- ordered(dd[[v]])
      } else stop(sprintf("Item %s is not numeric/factor; cannot coerce to ordered.", v))
    }
  }
  dd
}
dat <- coerce_to_ordered(dat, unique(c(ordered_full, ordered_trim)))

## =================== CFA models ===========================
model_1f <- "
  SACSg =~ SACS_1 + SACS_2 + SACS_5 + SACS_6 + SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
"
model_2f_full <- "
  SACSc =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  SACSt =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10
  SACSc ~~ SACSt
"
model_2f_trim <- "
  SACSc =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
  SACSt =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
  SACSc ~~ SACSt
"

## Fit 1-factor
fit_1f <- cfa(model_1f, data=dat, ordered=ordered_full,
              estimator="WLSMV", parameterization="theta", std.lv=TRUE)
fm_1f <- fitMeasures(fit_1f, c("cfi","tli","rmsea","srmr","bic"))
cat("\n--- 1-factor (10 items) ---\n"); print(fm_1f)

## Fit FULL 2-factor (10 items)
fit_2f_full <- cfa(model_2f_full, data=dat, ordered=ordered_full,
                   estimator="WLSMV", parameterization="theta", std.lv=TRUE)
fm_full <- fitMeasures(fit_2f_full, c("cfi","tli","rmsea","srmr","bic"))
cat("\n--- FULL 2-factor (10 items) ---\n"); print(fm_full)
fc_full <- lavaan::inspect(fit_2f_full, "std")$psi["SACSc","SACSt"]
cat("Latent r(SACSc,SACSt) full = ", round(fc_full,3), "\n")

## Δ indices (2f - 1f)
delta <- c(
  dCFI   = fm_full["cfi"]   - fm_1f["cfi"],
  dTLI   = fm_full["tli"]   - fm_1f["tli"],
  dRMSEA = fm_full["rmsea"] - fm_1f["rmsea"],
  dBIC   = fm_1f["bic"]     - fm_full["bic"]
)
cat("\nDelta (2f - 1f):\n"); print(delta)

## Heywood checker
check_heywood <- function(fit) {
  pe <- parameterEstimates(fit, standardized=TRUE)
  ov <- subset(pe, op=="~~" & lhs==rhs & grepl("^SACS_", lhs))
  lv <- subset(pe, op=="~~" & lhs==rhs & !grepl("^SACS_", lhs))
  list(
    any_ov_neg_raw = any(ov$est < 0, na.rm=TRUE),
    any_ov_neg_std = any(ov$std.all < 0, na.rm=TRUE),
    any_lv_neg_raw = any(lv$est < 0, na.rm=TRUE),
    any_lv_neg_std = any(lv$std.all < 0, na.rm=TRUE),
    ov_offenders = ov[which(ov$est<0 | ov$std.all<0), c("lhs","est","std.all")],
    lv_offenders = lv[which(lv$est<0 | lv$std.all<0), c("lhs","est","std.all")]
  )
}
cat("\nHeywood check (FULL 2f):\n"); print(check_heywood(fit_2f_full))

## TRIM 2-factor (drop SACS_9)
fit_2f_trim <- cfa(model_2f_trim, data=dat, ordered=ordered_trim,
                   estimator="WLSMV", parameterization="theta", std.lv=TRUE)
fm_trim <- fitMeasures(fit_2f_trim, c("cfi","tli","rmsea","srmr","bic"))
cat("\n--- TRIM 2-factor (drop SACS_9) ---\n"); print(fm_trim)
fc_trim <- lavaan::inspect(fit_2f_trim, "std")$psi["SACSc","SACSt"]
cat("Latent r(SACSc,SACSt) trim = ", round(fc_trim,3), "\n")
cat("\nHeywood check (TRIM 2f):\n"); print(check_heywood(fit_2f_trim))

## HTMT (explicit semTools namespace + fallback)
get_htmt <- function(fit) {
  out <- try(semTools::htmt(fit), silent=TRUE)
  if (inherits(out,"try-error")) {
    cat("\n[Note] semTools::htmt failed; computing simple latent correlation proxy instead.\n")
    return(lavaan::inspect(fit, "std")$psi)
  } else return(out)
}
cat("\nHTMT (trim):\n"); print(get_htmt(fit_2f_trim))

## Reliability (try; not all estimators output every index)
rel_trim <- try(semTools::reliability(fit_2f_trim), silent=TRUE)
if (!inherits(rel_trim,"try-error")) { cat("\nReliability (trim):\n"); print(rel_trim) }

## =================== Factor scores (trim) & OLS ======================
lv_trim <- lavPredict(fit_2f_trim, method="EBM", type="lv")
dat$SACS_Centering_trim    <- as.numeric(scale(lv_trim[,"SACSc"]))
dat$SACS_Transcending_trim <- as.numeric(scale(lv_trim[,"SACSt"]))

ols_with_tests <- function(form_base, form_full, data, label) {
  m1 <- lm(as.formula(form_base), data=data)
  m2 <- lm(as.formula(form_full), data=data)
  cat(sprintf("\n--- OLS (%s): base ---\n", label)); print(summary(m1))
  cat(sprintf("\n--- OLS (%s): full ---\n", label)); print(summary(m2))
  cat("\nΔR² F-test:\n"); print(anova(m1, m2))
  V <- vcovHC(m2, type="HC3")
  ht <- car::linearHypothesis(m2, "SACS_Centering_trim = SACS_Transcending_trim", vcov.=V)
  cat("\nWald test (robust) of β_SACSc = β_SACSt:\n"); print(ht)
}

# WEMWBS
ols_with_tests(
  "WEMWBS_Total ~ Age_clean + Gender_label + SACS_Centering_trim",
  "WEMWBS_Total ~ Age_clean + Gender_label + SACS_Centering_trim + SACS_Transcending_trim",
  dat, "WEMWBS"
)

# DASS total
ols_with_tests(
  "DASS_Total_x2 ~ Age_clean + Gender_label + SACS_Centering_trim",
  "DASS_Total_x2 ~ Age_clean + Gender_label + SACS_Centering_trim + SACS_Transcending_trim",
  dat, "DASS Total"
)

# DASS subscales
for (lab in c("Depression","Anxiety","Stress")) {
  y <- switch(lab,
              "Depression"="DASS_Depression_x2",
              "Anxiety"   ="DASS_Anxiety_x2",
              "Stress"    ="DASS_Stress_x2")
  ols_with_tests(
    sprintf("%s ~ Age_clean + Gender_label + SACS_Centering_trim", y),
    sprintf("%s ~ Age_clean + Gender_label + SACS_Centering_trim + SACS_Transcending_trim", y),
    dat, paste0("DASS—",lab)
  )
}

## =================== MIMIC models (trim; WLSMV) ======================
run_mimic <- function(pred) {
  model <- sprintf('
    SACSc =~ %s
    SACSt =~ %s
    SACSc ~ %s + Age_clean + Gender_label__Non.Binary + Gender_label__Woman
    SACSt ~ %s + Age_clean + Gender_label__Non.Binary + Gender_label__Woman
    SACSc ~~ SACSt
  ',
                   paste(c_items_trim, collapse=" + "),
                   paste(t_items_trim, collapse=" + "),
                   pred, pred)
  
  fit <- sem(model, data=dat, ordered=ordered_trim,
             estimator="WLSMV", parameterization="theta", std.lv=TRUE)
  cat(sprintf("\n--- MIMIC (%s -> latent SACSc/SACSt; WLSMV, TRIM) ---\n", pred))
  print(summary(fit, fit.measures=TRUE, standardized=TRUE))
  invisible(fit)
}

fit_mimic_wem <- run_mimic("WEMWBS_Total")
fit_mimic_dt  <- run_mimic("DASS_Total_x2")
fit_mimic_dd  <- run_mimic("DASS_Depression_x2")
fit_mimic_da  <- run_mimic("DASS_Anxiety_x2")
fit_mimic_ds  <- run_mimic("DASS_Stress_x2")

## =================== Decision prompts ====================
cat("\nDecision guide:\n",
    "- Prefer 2-factor SACS if ΔCFI ≥ .01, ΔRMSEA ≤ -0.015, ΔBIC > 10, and HTMT < .85–.90.\n",
    "- Independent impacts supported if (manifest OLS): each β ≠ 0 AND ΔR² ≳ .01, with robust Wald rejecting β_SACSc = β_SACSt.\n",
    "- For MIMIC, prioritise standardized paths under WLSMV (theta) and inspect residuals.\n")
