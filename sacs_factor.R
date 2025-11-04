## =============================
## SACS one- vs two-factor check
## =============================

## 0) Packages (align with prior conventions)
pkgs <- c("lavaan","psych","openxlsx")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(lavaan)
library(psych)
xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)

## 1) Paths and item names (as in previous scripts)
root         <- "/Users/jamesjackson/Desktop/dir"
out_dir      <- file.path(root, "phase_final_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sacs_items <- sprintf("SACS_%d", 1:10)
# Centring items per your scoring: 1,2,5,6
sacs_centering_items    <- c("SACS_1","SACS_2","SACS_5","SACS_6")
# Transcendence items per your scoring: 3,4,7,8,9,10
sacs_transcendence_items<- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")

# Helper to write CSV/XLSX like your scripts
write_dual <- function(df, stem, directory = out_dir){
  csv_path <- file.path(directory, paste0(stem, ".csv"))
  utils::write.csv(df, csv_path, row.names = FALSE)
  if (xlsx_ok){
    xlsx_path <- file.path(directory, paste0(stem, ".xlsx"))
    openxlsx::write.xlsx(df, file = xlsx_path, overwrite = TRUE)
  }
  invisible(csv_path)
}

## 2) Data source
## Use df if it already exists in your session (e.g., from Complete.R).
## Otherwise, fall back to your standard load (scored_primary.csv or raw copy.xlsx).
if (!exists("df")) {
  scored_file <- file.path(root, "phase03_scoring", "scored_primary.csv")
  if (file.exists(scored_file)) {
    df <- utils::read.csv(scored_file, check.names = TRUE, stringsAsFactors = FALSE)
  } else {
    raw_xlsx <- file.path(root, "data", "raw copy.xlsx")
    if (!file.exists(raw_xlsx)) stop("Data not found: provide df, or ensure scored_primary.csv/raw copy.xlsx exist.")
    if (!xlsx_ok) install.packages("openxlsx")
    library(openxlsx)
    df <- as.data.frame(read.xlsx(raw_xlsx, sheet = 1))
    names(df) <- make.names(names(df))  # just in case
  }
}

## 3) Keep complete SACS cases; coerce to ordered factors for WLSMV CFA
stopifnot(all(sacs_items %in% names(df)))
df_sacs <- df[, sacs_items]
# Coerce numeric if needed
df_sacs[] <- lapply(df_sacs, function(x) suppressWarnings(as.numeric(x)))
# Retain rows with all 10 SACS present
cc <- stats::complete.cases(df_sacs)
df_sacs <- df_sacs[cc, , drop = FALSE]

# Treat Likert 1–7 as ordered categories
df_ord <- as.data.frame(lapply(df_sacs, function(x){
  if (all(is.na(x))) return(ordered(x))
  ordered(x, levels = sort(na.omit(unique(x))))
}))

cat("SACS CFA sample size (complete SACS): ", nrow(df_ord), "\n")

## 4) Specify models
# One-factor SACS
mod_1f <- '
SACS =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'

# Two-factor SACS (correlated factors): Centring vs Transcendence
mod_2f <- sprintf('
Centring =~ %s
Transcendence =~ %s
Centring ~~ Transcendence
',
                  paste(sacs_centering_items, collapse=" + "),
                  paste(sacs_transcendence_items, collapse=" + ")
)

## 5) Fit CFA with WLSMV (robust for ordinal Likert)
fit_1f <- cfa(model = mod_1f, data = df_ord, estimator = "WLSMV", ordered = sacs_items, std.lv = TRUE)
fit_2f <- cfa(model = mod_2f, data = df_ord, estimator = "WLSMV", ordered = sacs_items, std.lv = TRUE)

## 6) Extract and save fit indices & loadings
get_fit <- function(fit){
  fm <- fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
  as.data.frame(as.list(round(fm, 3)))
}

fit_tab <- rbind(
  cbind(Model = "SACS_1Factor",        get_fit(fit_1f)),
  cbind(Model = "SACS_2Factor_C+T",    get_fit(fit_2f))
)

write_dual(fit_tab, "SACS_CFA_fit")

# Standardised loadings tables
ld_1f <- parameterEstimates(fit_1f, standardized = TRUE)
ld_2f <- parameterEstimates(fit_2f, standardized = TRUE)
ld_1f_std <- subset(ld_1f, op == "=~", select = c("lhs","rhs","est","se","z","pvalue","std.all"))
ld_2f_std <- subset(ld_2f, op == "=~", select = c("lhs","rhs","est","se","z","pvalue","std.all"))

write_dual(ld_1f_std, "SACS_CFA_loadings_1factor")
write_dual(ld_2f_std, "SACS_CFA_loadings_2factor")

## 7) Model comparison
# With WLSMV, use DIFFTEST (scaled difference test) via anova() on lavaan fits
cmp <- tryCatch(anova(fit_1f, fit_2f), error = function(e) NULL)

# Also capture factor correlation from the 2-factor model
fcorr <- subset(parameterEstimates(fit_2f, standardized = TRUE),
                lhs == "Centring" & rhs == "Transcendence" & op == "~~")

## 8) Reliability by factor (omega, using polychorics)
# Polychoric matrix
poly <- psych::polychoric(df_sacs)$rho
# Omega for all 10 as single factor
om_1f <- psych::omega(poly, nfactors = 1, plot = FALSE)
# Omega by subscale (Centring and Transcendence)
poly_C  <- psych::polychoric(df_sacs[, sacs_centering_items])$rho
poly_T  <- psych::polychoric(df_sacs[, sacs_transcendence_items])$rho
om_C <- psych::omega(poly_C, nfactors = 1, plot = FALSE)
om_T <- psych::omega(poly_T, nfactors = 1, plot = FALSE)

rel_tab <- data.frame(
  Scale = c("SACS_1Factor_All10","SACS_Centring_4","SACS_Transcendence_6"),
  Omega_Total = round(c(om_1f$omega.tot, om_C$omega.tot, om_T$omega.tot), 3),
  Omega_H     = round(c(om_1f$omega.h,   om_C$omega.h,   om_T$omega.h), 3)
)
write_dual(rel_tab, "SACS_omega_reliability")

## 9) Optional exploratory check (EFA) to mirror “recent paper” claim
## Parallel analysis + oblimin 1- vs 2-factor MINRES, on polychorics (saved for transparency)
pa <- fa.parallel(poly, n.obs = nrow(df_ord), fa = "fa", fm = "minres", ylabel = FALSE, plot = FALSE)
efa1 <- fa(poly, nfactors = 1, n.obs = nrow(df_ord), fm = "minres", rotate = "oblimin")
efa2 <- fa(poly, nfactors = 2, n.obs = nrow(df_ord), fm = "minres", rotate = "oblimin")

efa1_tab <- data.frame(Item = rownames(efa1$loadings[]), Loading = round(efa1$loadings[,1], 3))
efa2_tab <- as.data.frame(round(efa2$loadings[,], 3)); efa2_tab$Item <- rownames(efa2_tab)
efa2_tab <- efa2_tab[, c("Item", colnames(efa2$loadings))]

write_dual(efa1_tab, "SACS_EFA_1factor_loadings")
write_dual(efa2_tab, "SACS_EFA_2factor_loadings")

## 10) Write a compact text summary
summ_lines <- c(
  "SACS CFA one- vs two-factor comparison",
  paste0("N (complete SACS): ", nrow(df_ord)),
  "\n-- Fit indices --",
  capture.output(print(fitMeasures(fit_1f, c("cfi","tli","rmsea","srmr","chisq","df","pvalue")), digits = 3)),
  capture.output(print(fitMeasures(fit_2f, c("cfi","tli","rmsea","srmr","chisq","df","pvalue")), digits = 3)),
  "\n-- Scaled DIFFTEST (WLSMV) --",
  if (is.null(cmp)) "DIFFTEST not available" else capture.output(print(cmp)),
  "\n-- Factor correlation (2-factor) --",
  if (nrow(fcorr)==0) "NA" else paste0("r(Centring, Transcendence) = ", round(fcorr$est, 3),
                                       " [std = ", round(fcorr$std.all, 3), "]"),
  "\n-- Reliability (omega) --",
  capture.output(print(rel_tab, row.names = FALSE))
)
writeLines(summ_lines, file.path(out_dir, "SACS_factor_models_summary.txt"))

cat("Done. See files in:", out_dir, "\n",
    "- SACS_CFA_fit.(csv/xlsx)\n",
    "- SACS_CFA_loadings_1factor.(csv/xlsx)\n",
    "- SACS_CFA_loadings_2factor.(csv/xlsx)\n",
    "- SACS_omega_reliability.(csv/xlsx)\n",
    "- SACS_EFA_1factor_loadings.(csv/xlsx)\n",
    "- SACS_EFA_2factor_loadings.(csv/xlsx)\n",
    "- SACS_factor_models_summary.txt\n", sep = "")
## ============================
## SACS CFA: diagnose → refit
## ============================

suppressPackageStartupMessages({
  library(lavaan); library(psych); library(dplyr)
})

## --- Assumes df with SACS_1..SACS_10 is in memory (from your earlier run) ---
sacs_items <- sprintf("SACS_%d", 1:10)
stopifnot(all(sacs_items %in% names(df)))
df_sacs <- df[, sacs_items]

## Coerce numeric; keep complete SACS cases
df_sacs[] <- lapply(df_sacs, function(x) suppressWarnings(as.numeric(x)))
df_sacs <- df_sacs[complete.cases(df_sacs), , drop = FALSE]

message("N complete SACS: ", nrow(df_sacs))

## ---------- 2.1 Category diagnostics ----------
cat_counts <- lapply(df_sacs, function(x) table(factor(x, levels = sort(unique(na.omit(x))))))
cat_overview <- do.call(rbind, lapply(names(cat_counts), function(v){
  cc <- cat_counts[[v]]
  data.frame(Item=v, Cat=as.numeric(names(cc)), n=as.integer(cc),
             stringsAsFactors = FALSE)
}))
print(cat_overview %>% group_by(Item) %>% summarise(min_cat_n=min(n), used_cats=n()), n = 50)

## Flag items with extremely sparse categories (min cell < 3)
sparse_items <- cat_overview %>% group_by(Item) %>% summarise(min_cat_n=min(n)) %>%
  filter(min_cat_n < 3) %>% pull(Item)

if (length(sparse_items)){
  message("Sparse categories detected in: ", paste(sparse_items, collapse=", "))
}

## ---------- 2.2 Optional collapsing of extreme sparse categories ----------
## Collapses the outermost categories until all cells >= min_n (default 5).
## Conservative; only collapses if needed.
collapse_sparse <- function(x, min_n = 5){
  x <- as.numeric(x)
  lev <- sort(unique(x))
  if (length(lev) < 3) return(x)
  repeat {
    tt <- table(factor(x, levels = lev))
    if (all(tt >= min_n) || length(lev) <= 3) break
    ## collapse the smallest extreme side
    left  <- lev[1]; right <- lev[length(lev)]
    if (tt[1] <= tt[length(tt)]) {
      x[x == left] <- lev[2]
    } else {
      x[x == right] <- lev[length(lev)-1]
    }
    lev <- sort(unique(x))
  }
  x
}

df_sacs_re <- df_sacs
if (length(sparse_items)){
  for (v in sparse_items) df_sacs_re[[v]] <- collapse_sparse(df_sacs_re[[v]], min_n = 5)
  message("Collapsed extreme categories for sparse items (min cell target = 5).")
}

## ---------- 2.3 Build ordered data frames ----------
ordify <- function(d){
  as.data.frame(lapply(d, function(x){
    x <- as.numeric(x)
    ordered(x, levels = sort(unique(na.omit(x))))
  }))
}
df_ord_raw <- ordify(df_sacs)
df_ord_re  <- ordify(df_sacs_re)

## ---------- 2.4 Model specs ----------
sacs_centering_items     <- c("SACS_1","SACS_2","SACS_5","SACS_6")
sacs_transcendence_items <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")

mod_1f <- '
SACS =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10
'

mod_2f <- sprintf('
Centring =~ %s
Transcendence =~ %s
Centring ~~ Transcendence
',
                  paste(sacs_centering_items, collapse=" + "),
                  paste(sacs_transcendence_items, collapse=" + ")
)

## ---------- 2.5 Fit with WLSMV + theta parameterization ----------
fit_safe <- function(model, dat){
  cfa(model = model, data = dat, estimator = "WLSMV",
      ordered = names(dat), parameterization = "theta", std.lv = TRUE)
}

fit1_raw <- fit_safe(mod_1f, df_ord_raw)
fit2_raw <- fit_safe(mod_2f, df_ord_raw)

## If warnings about vcov persist, try the collapsed data
need_refit <- any(grepl("not identified|not positive definite", fitMeasures(fit2_raw, "cfi"),
                        ignore.case = TRUE)) || !lavInspect(fit2_raw, "converged")

if (!lavInspect(fit1_raw, "converged") || !lavInspect(fit2_raw, "converged")) need_refit <- TRUE

if (need_refit && length(sparse_items)){
  message("Refitting on category-collapsed data due to instability.")
  fit1 <- fit_safe(mod_1f, df_ord_re)
  fit2 <- fit_safe(mod_2f, df_ord_re)
  dat_used <- "collapsed"
} else {
  fit1 <- fit1_raw; fit2 <- fit2_raw
  dat_used <- "raw"
}

message("Data used for final comparison: ", dat_used)

## ---------- 2.6 Quick stability checks ----------
heywood <- function(fit){
  pe <- parameterEstimates(fit, standardized = TRUE)
  any(pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0, na.rm = TRUE)
}
cat("Heywood (negative residual variances) 1-factor: ", heywood(fit1), "\n")
cat("Heywood (negative residual variances) 2-factor: ", heywood(fit2), "\n")

## Factor correlation
fc <- subset(parameterEstimates(fit2, standardized = TRUE),
             lhs == "Centring" & rhs == "Transcendence" & op == "~~")
if (nrow(fc)) cat("r(Centring, Transcendence): ", round(fc$std.all, 3), "\n")

## ---------- 2.7 Fit indices + DIFFTEST ----------
get_fit <- function(fit)
  round(fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr")), 3)

fi1 <- get_fit(fit1); fi2 <- get_fit(fit2)
print(rbind(`SACS_1Factor`=fi1, `SACS_2Factor_C+T`=fi2))

## DIFFTEST (proper for WLSMV)
cmp <- try(anova(fit1, fit2), silent = TRUE)
if (!inherits(cmp, "try-error")) print(cmp) else message("DIFFTEST unavailable.")

## ---------- 2.8 Loadings (standardised) ----------
ld1 <- subset(parameterEstimates(fit1, standardized = TRUE), op == "=~",
              select = c("lhs","rhs","est","se","z","pvalue","std.all"))
ld2 <- subset(parameterEstimates(fit2, standardized = TRUE), op == "=~",
              select = c("lhs","rhs","est","se","z","pvalue","std.all"))
cat("\nStd. loadings (1-factor):\n"); print(ld1, row.names = FALSE)
cat("\nStd. loadings (2-factor):\n"); print(ld2, row.names = FALSE)

## ---------- 2.9 Reliability (omega_total only; skip omega_h for 1-factor) ----------
poly_raw <- try(psych::polychoric(df_sacs)$rho, silent = TRUE)
if (!inherits(poly_raw, "try-error")) {
  om_all <- psych::omega(poly_raw, nfactors = 1, plot = FALSE)
  cat("\nOmega_total (all 10 items, 1-factor): ", round(om_all$omega.tot, 3), "\n")
  poly_C <- try(psych::polychoric(df_sacs[, sacs_centering_items])$rho, silent = TRUE)
  poly_T <- try(psych::polychoric(df_sacs[, sacs_transcendence_items])$rho, silent = TRUE)
  if (!inherits(poly_C, "try-error")) {
    om_C <- psych::omega(poly_C, nfactors=1, plot=FALSE)
    cat("Omega_total (Centring-4): ", round(om_C$omega.tot, 3), "\n")
  }
  if (!inherits(poly_T, "try-error")) {
    om_T <- psych::omega(poly_T, nfactors=1, plot=FALSE)
    cat("Omega_total (Transcendence-6): ", round(om_T$omega.tot, 3), "\n")
  }
} else {
  message("Polychoric failed; consider continuous MLR sensitivity.")
}

## ---------- 2.10 Sensitivity (continuous, MLR) ----------
fit1_mlr <- cfa(mod_1f, data = df_sacs, estimator = "MLR", std.lv = TRUE)
fit2_mlr <- cfa(mod_2f, data = df_sacs, estimator = "MLR", std.lv = TRUE)
cat("\n[MLR sensitivity] CFI/TLI/RMSEA/SRMR (1f vs 2f):\n")
print(rbind(`1f`=round(fitMeasures(fit1_mlr, c("cfi","tli","rmsea","srmr")),3),
            `2f`=round(fitMeasures(fit2_mlr, c("cfi","tli","rmsea","srmr")),3)))
