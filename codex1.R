# ================================================================
# Master SAC Analysis Pipeline
# ---------------------------------------------------------------
# Consolidates the strongest elements across historical scripts to
# deliver a single, reproducible workflow for the SAC project.
#
# Major capabilities:
#   * Data ingestion with flexible scoring fallbacks
#   * Missing-data diagnostics (Little's MCAR, MVN, Shapiro)
#   * Reliability, descriptives, and assumption checks
#   * SACS EFA (polychoric) and CFA comparisons
#   * Correlations, bootstraps, regression, and mediation tests
#   * Multiple imputation sensitivity analyses
#   * Bayesian SEM (blavaan) with WLSMV fallback
#   * Structured exports to timestamped output folders (CSV + XLSX)
# ================================================================

## ---- 0. Setup ------------------------------------------------------------
# Default project root (override by defining `root` before sourcing)
if (!exists("root")) {
  root <- "/Users/jamesjackson/Desktop/dir"
}

analysis_root <- file.path(root, "phase_final_analysis")
if (!dir.exists(analysis_root)) dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(analysis_root, paste0("master_analysis_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper to ensure packages are installed and loaded
ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing missing package: ", pkg)
      try(install.packages(pkg, dependencies = TRUE), silent = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

core_packages <- c(
  "tidyverse", "openxlsx", "psych", "lavaan", "blavaan", "mice",
  "boot", "car", "lm.beta", "GPArotation", "naniar", "MissMech",
  "MVN", "energy", "moments"
)
ensure_packages(core_packages)

set.seed(2025)

xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)

# Unified writer (CSV always, XLSX when available)
write_dual <- function(df, stem, directory = out_dir) {
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  csv_path  <- file.path(directory, paste0(stem, ".csv"))
  utils::write.csv(df, csv_path, row.names = FALSE)
  if (xlsx_ok) {
    openxlsx::write.xlsx(df, file.path(directory, paste0(stem, ".xlsx")), overwrite = TRUE)
  }
  invisible(csv_path)
}

## ---- 1. Data ingestion & scoring ----------------------------------------
message("\n[1] Data ingestion")

# Candidate paths in order of preference
scored_paths <- c(
  file.path(root, "phase03_scoring", "scored_primary.csv"),
  file.path(root, "scored_primary.csv")
)
raw_paths <- c(
  file.path(root, "raw.csv"),
  file.path(root, "data", "raw copy.csv"),
  file.path(root, "data", "raw copy.xlsx")
)

load_pre_scored <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) {
      message("Loaded scored data from: ", p)
      if (grepl("\\.xlsx$", p, ignore.case = TRUE)) {
        return(as.data.frame(openxlsx::read.xlsx(p, sheet = 1)))
      } else {
        return(utils::read.csv(p, check.names = TRUE, stringsAsFactors = FALSE))
      }
    }
  }
  NULL
}

# Scoring helper -----------------------------------------------------------
score_scale <- function(data, items, min_prop = 0.8, multiplier = 1) {
  miss_items <- setdiff(items, names(data))
  if (length(miss_items)) return(rep(NA_real_, nrow(data)))
  sub <- data[, items, drop = FALSE]
  sub[] <- lapply(sub, function(x) suppressWarnings(as.numeric(x)))
  available <- rowSums(!is.na(sub))
  scores <- rowSums(sub, na.rm = TRUE) * multiplier
  threshold <- ceiling(min_prop * length(items))
  scores[available < threshold] <- NA_real_
  scores
}

score_all_scales <- function(df) {
  sacs_total_items     <- sprintf("SACS_%d", 1:10)
  sacs_centering_items <- c("SACS_1", "SACS_2", "SACS_5", "SACS_6")
  sacs_trans_items     <- c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_9", "SACS_10")
  cfq_items            <- sprintf("CFQ_%d", 1:7)
  atq_f_items          <- sprintf("ATQ_%df", 1:15)
  atq_b_items          <- sprintf("ATQ_%db", 1:15)
  wemwbs_items         <- sprintf("WEMWBS_%d", 1:14)
  dass_items           <- sprintf("DASS.21_%d", 1:21)
  dass_dep             <- c(3, 5, 10, 13, 16, 17, 21)
  dass_anx             <- c(2, 4, 7, 9, 15, 19, 20)
  dass_str             <- c(1, 6, 8, 11, 12, 14, 18)
  
  df$SACS_Total        <- score_scale(df, sacs_total_items)
  df$SACS_Centering    <- score_scale(df, sacs_centering_items)
  df$SACS_Transcending <- score_scale(df, sacs_trans_items)
  df$CFQ_Total         <- score_scale(df, cfq_items)
  df$ATQ_Frequency     <- score_scale(df, atq_f_items)
  df$ATQ_Believability <- score_scale(df, atq_b_items)
  df$ATQ_Combined      <- score_scale(df, c(atq_f_items, atq_b_items))
  df$WEMWBS_Total      <- score_scale(df, wemwbs_items)
  df$DASS_Depression_x2<- score_scale(df, paste0("DASS.21_", dass_dep), multiplier = 2)
  df$DASS_Anxiety_x2   <- score_scale(df, paste0("DASS.21_", dass_anx), multiplier = 2)
  df$DASS_Stress_x2    <- score_scale(df, paste0("DASS.21_", dass_str), multiplier = 2)
  df$DASS_Total_x2     <- score_scale(df, dass_items, multiplier = 2)
  
  score_vars <- c(
    "SACS_Total","SACS_Centering","SACS_Transcending",
    "CFQ_Total","ATQ_Frequency","ATQ_Believability","ATQ_Combined",
    "WEMWBS_Total",
    "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2"
  )
  for (sv in intersect(score_vars, names(df))) {
    zname <- paste0(sv, "_z")
    vals <- as.numeric(df[[sv]])
    df[[zname]] <- scale(vals)
  }
  df
}

# Load data ---------------------------------------------------------------
df <- NULL
if (exists("df") && is.data.frame(df)) {
  message("Using data.frame `df` already in environment.")
} else {
  df <- load_pre_scored(scored_paths)
}

if (is.null(df)) {
  message("No scored dataset found. Attempting to load raw data and score.")
  for (p in raw_paths) {
    if (!file.exists(p)) next
    message("Loaded raw data from: ", p)
    if (grepl("\\.xlsx$", p, ignore.case = TRUE)) {
      df <- as.data.frame(openxlsx::read.xlsx(p, sheet = 1))
    } else {
      df <- utils::read.csv(p, check.names = TRUE, stringsAsFactors = FALSE)
    }
    break
  }
  if (is.null(df)) stop("No data source found. Please ensure scored_primary.csv or raw.csv is available.")
  names(df) <- gsub("^DASS-21_", "DASS.21_", names(df))
  df <- score_all_scales(df)
} else {
  names(df) <- gsub("^DASS-21_", "DASS.21_", names(df))
  required_scores <- c("SACS_Total","SACS_Centering","SACS_Transcending","CFQ_Total",
                       "ATQ_Believability","WEMWBS_Total","DASS_Total_x2")
  missing_scores <- setdiff(required_scores, names(df))
  if (length(missing_scores)) {
    message("Scored dataset missing composites; computing: ", paste(missing_scores, collapse = ", "))
    df <- score_all_scales(df)
  }
}

analysis_n <- nrow(df)
message("Rows in analysis dataset: ", analysis_n)

write_dual(head(df, 10), "preview_first10")

## ---- 2. Missingness & normality diagnostics ------------------------------
message("\n[2] Missingness diagnostics")
miss_dir <- file.path(out_dir, "missingness")

default_scores <- intersect(c(
  "SACS_Total","SACS_Centering","SACS_Transcending",
  "CFQ_Total","ATQ_Believability","ATQ_Frequency","ATQ_Combined",
  "WEMWBS_Total",
  "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2",
  "Age_clean","Gender_label"
), names(df))

miss_summary <- data.frame(
  variable = default_scores,
  n_missing = sapply(default_scores, function(v) sum(is.na(df[[v]]))),
  pct_missing = round(100 * sapply(default_scores, function(v) mean(is.na(df[[v]]))), 2)
)
write_dual(miss_summary, "M0_missing_summary", miss_dir)

md_pattern <- tryCatch(mice::md.pattern(df[default_scores], plot = FALSE), error = function(e) NULL)
if (is.null(md_pattern)) {
  write_dual(data.frame(note = "md.pattern failed (no missingness or singular pattern)."),
             "M1_missing_patterns", miss_dir)
} else {
  write_dual(as.data.frame(md_pattern, check.names = FALSE), "M1_missing_patterns", miss_dir)
}

# Little's MCAR tests
mcar_df <- df[default_scores]
if ("Gender_label" %in% names(mcar_df)) mcar_df$Gender_label <- as.numeric(as.factor(mcar_df$Gender_label))
keep_cols <- vapply(mcar_df, function(x) sum(!is.na(x)) >= 3 && length(unique(stats::na.omit(x))) >= 2, logical(1))
mcar_df <- mcar_df[, keep_cols, drop = FALSE]

mcar_rows <- list()
if (ncol(mcar_df) > 1) {
  res_naniar <- try(naniar::mcar_test(mcar_df), silent = TRUE)
  if (!inherits(res_naniar, "try-error")) {
    stat <- suppressWarnings(as.numeric(res_naniar$statistic))[1]
    df_stat <- if (!is.null(res_naniar$parameter)) as.numeric(res_naniar$parameter)[1] else NA_real_
    p_val <- suppressWarnings(as.numeric(res_naniar$p.value))[1]
    mcar_rows[[length(mcar_rows) + 1]] <- data.frame(method = "naniar::mcar_test",
                                                     statistic = stat, df = df_stat, p_value = p_val)
  }
  res_missmech <- try(MissMech::TestMCARNormality(as.matrix(mcar_df)), silent = TRUE)
  if (!inherits(res_missmech, "try-error")) {
    mcar_rows[[length(mcar_rows) + 1]] <- data.frame(method = "MissMech::TestMCARNormality",
                                                     statistic = res_missmech$chi.square,
                                                     df = res_missmech$df,
                                                     p_value = res_missmech$p.value)
  }
}
if (length(mcar_rows)) {
  write_dual(do.call(rbind, mcar_rows), "M2_LittleMCAR_results", miss_dir)
} else {
  write_dual(data.frame(note = "MCAR tests unavailable (insufficient variables or package error)."),
             "M2_LittleMCAR_results", miss_dir)
}

# Multivariate normality (Henze-Zirkler, Mardia, Energy)
num_df <- mcar_df[, vapply(mcar_df, is.numeric, logical(1)), drop = FALSE]
num_cc <- stats::na.omit(num_df)
normality_rows <- list()
if (nrow(num_cc) >= 10 && ncol(num_cc) >= 2) {
  hz <- try(MVN::mvn(num_cc, mvnTest = "hz", covariance = TRUE, desc = FALSE), silent = TRUE)
  if (!inherits(hz, "try-error") && is.list(hz$HZ)) {
    normality_rows[[length(normality_rows) + 1]] <- data.frame(test = "Henze-Zirkler",
                                                               statistic = hz$HZ$HZ,
                                                               p_value = hz$HZ$p.value,
                                                               n_complete = nrow(num_cc),
                                                               p_vars = ncol(num_cc))
  }
  md <- try(MVN::mvn(num_cc, mvnTest = "mardia", covariance = TRUE, desc = FALSE), silent = TRUE)
  if (!inherits(md, "try-error") && is.data.frame(md$multivariateNormality)) {
    mtab <- md$multivariateNormality
    for (i in seq_len(nrow(mtab))) {
      normality_rows[[length(normality_rows) + 1]] <- data.frame(test = as.character(mtab[i, 1]),
                                                                 statistic = as.numeric(mtab[i, 2]),
                                                                 p_value = as.numeric(mtab[i, 4]),
                                                                 n_complete = nrow(num_cc),
                                                                 p_vars = ncol(num_cc))
    }
  }
  et <- try(energy::mvnorm.etest(as.matrix(num_cc), R = 1000), silent = TRUE)
  if (!inherits(et, "try-error")) {
    normality_rows[[length(normality_rows) + 1]] <- data.frame(test = "Energy E-test",
                                                               statistic = et$statistic,
                                                               p_value = et$p.value,
                                                               n_complete = nrow(num_cc),
                                                               p_vars = ncol(num_cc))
  }
}
if (length(normality_rows)) {
  write_dual(do.call(rbind, normality_rows), "M3_multivariate_normality", miss_dir)
} else {
  write_dual(data.frame(note = "MVN tests skipped (need >=10 complete numeric cases)."),
             "M3_multivariate_normality", miss_dir)
}

# Univariate Shapiro + moments
shapiro_rows <- list()
for (v in colnames(num_df)) {
  xv <- stats::na.omit(num_df[[v]])
  if (length(xv) >= 3 && length(xv) <= 5000) {
    sw <- try(stats::shapiro.test(xv), silent = TRUE)
    if (!inherits(sw, "try-error")) {
      skew <- try(moments::skewness(xv), silent = TRUE)
      kurt <- try(moments::kurtosis(xv) - 3, silent = TRUE)
      shapiro_rows[[length(shapiro_rows) + 1]] <- data.frame(
        variable = v,
        n = length(xv),
        shapiro_W = unname(sw$statistic),
        shapiro_p = sw$p.value,
        skewness = if (!inherits(skew, "try-error")) as.numeric(skew) else NA_real_,
        excess_kurtosis = if (!inherits(kurt, "try-error")) as.numeric(kurt) else NA_real_
      )
    }
  }
}
if (length(shapiro_rows)) {
  write_dual(do.call(rbind, shapiro_rows), "M4_univariate_normality", miss_dir)
}

## ---- 3. Reliability & descriptives --------------------------------------
message("\n[3] Reliability and descriptives")
reliability_dir <- file.path(out_dir, "reliability_descriptives")

alpha_items <- list(
  SACS_Total        = sprintf("SACS_%d", 1:10),
  SACS_Centering    = c("SACS_1", "SACS_2", "SACS_5", "SACS_6"),
  SACS_Transcending = c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_9", "SACS_10"),
  CFQ_Total         = sprintf("CFQ_%d", 1:7),
  ATQ_Frequency     = sprintf("ATQ_%df", 1:15),
  ATQ_Believability = sprintf("ATQ_%db", 1:15),
  WEMWBS_Total      = sprintf("WEMWBS_%d", 1:14),
  DASS_Total        = sprintf("DASS.21_%d", 1:21)
)

alpha_table <- lapply(names(alpha_items), function(scale) {
  items <- alpha_items[[scale]]
  if (!all(items %in% names(df))) return(NULL)
  sub <- df[, items]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (nrow(sub) < 3) return(NULL)
  a <- psych::alpha(sub, check.keys = FALSE)
  data.frame(scale = scale,
             n_items = length(items),
             alpha_raw = a$total$raw_alpha,
             alpha_std = a$total$std.alpha,
             stringsAsFactors = FALSE)
})
alpha_table <- do.call(rbind, alpha_table)
if (!is.null(alpha_table)) write_dual(alpha_table, "R1_cronbach_alpha", reliability_dir)

score_vars <- intersect(c(
  "SACS_Total","SACS_Centering","SACS_Transcending",
  "CFQ_Total","ATQ_Frequency","ATQ_Believability","ATQ_Combined",
  "WEMWBS_Total",
  "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2"
), names(df))

desc_stats <- lapply(score_vars, function(v) {
  x <- df[[v]]
  finite_x <- x[is.finite(x)]
  n <- length(finite_x)
  data.frame(
    variable = v,
    N = n,
    Missing = sum(!is.finite(x)),
    Mean = if (n) mean(finite_x) else NA_real_,
    SD = if (n > 1) stats::sd(finite_x) else NA_real_,
    Skewness = if (n > 2) mean(((finite_x - mean(finite_x))/stats::sd(finite_x))^3) else NA_real_,
    Kurtosis_excess = if (n > 3) mean(((finite_x - mean(finite_x))/stats::sd(finite_x))^4) - 3 else NA_real_
  )
})
write_dual(do.call(rbind, desc_stats), "R2_descriptive_stats", reliability_dir)

# Shapiro-Wilk on composites
shapiro_scores <- lapply(score_vars, function(v) {
  x <- stats::na.omit(df[[v]])
  if (length(x) < 3 || length(x) > 5000) return(NULL)
  sw <- try(stats::shapiro.test(x), silent = TRUE)
  if (inherits(sw, "try-error")) return(NULL)
  data.frame(variable = v, W = unname(sw$statistic), p_value = sw$p.value)
})
shapiro_scores <- do.call(rbind, shapiro_scores)
if (!is.null(shapiro_scores)) write_dual(shapiro_scores, "R3_shapiro_scores", reliability_dir)

## ---- 4. SACS EFA (polychoric) -------------------------------------------
message("\n[4] SACS EFA")
efa_dir <- file.path(out_dir, "sacs_efa")

sacs_items <- sprintf("SACS_%d", 1:10)
if (all(sacs_items %in% names(df))) {
  sacs_dat <- df[, sacs_items]
  sacs_dat[] <- lapply(sacs_dat, function(x) suppressWarnings(as.numeric(x)))
  sacs_dat <- sacs_dat[complete.cases(sacs_dat), , drop = FALSE]
  if (nrow(sacs_dat) >= 30) {
    tab_min <- sapply(sacs_dat, function(x) min(table(factor(x, levels = sort(unique(x))))) )
    sparse <- names(tab_min)[tab_min < 3]
    if (length(sparse)) {
      collapse_sparse <- function(x, min_n = 5) {
        x <- as.numeric(x); lev <- sort(unique(x))
        if (length(lev) <= 3) return(x)
        repeat {
          tt <- table(factor(x, levels = lev))
          if (all(tt >= min_n) || length(lev) <= 3) break
          if (tt[1] <= tt[length(tt)]) {
            x[x == lev[1]] <- lev[2]
          } else {
            x[x == lev[length(lev)]] <- lev[length(lev) - 1]
          }
          lev <- sort(unique(x))
        }
        x
      }
      for (v in sparse) sacs_dat[[v]] <- collapse_sparse(sacs_dat[[v]])
    }
    poly <- psych::polychoric(sacs_dat)
    R <- poly$rho
    KMO_res <- psych::KMO(R)
    bart <- psych::cortest.bartlett(R, n = nrow(sacs_dat))
    write_dual(as.data.frame(R), "E0_polychoric", efa_dir)
    write_dual(data.frame(Measure = "KMO_overall", Value = KMO_res$MSA), "E1_KMO_overall", efa_dir)
    write_dual(data.frame(Item = names(KMO_res$MSAi), MSA = KMO_res$MSAi), "E2_KMO_items", efa_dir)
    write_dual(data.frame(Chi2 = bart$chisq, df = bart$df, p = bart$p.value), "E3_Bartlett", efa_dir)
    
    fa_par <- psych::fa.parallel(R, n.obs = nrow(sacs_dat), fa = "fa", fm = "minres", plot = FALSE)
    vss <- psych::VSS(R, n = 5, rotate = "none", fm = "minres", plot = FALSE)
    write_dual(data.frame(Parallel_nf = fa_par$nfact,
                          MAP_choice = which.min(vss$map)),
               "E4_factor_recommendations", efa_dir)
    
    efa1 <- psych::fa(R, nfactors = 1, n.obs = nrow(sacs_dat), fm = "minres", rotate = "oblimin")
    efa2 <- psych::fa(R, nfactors = 2, n.obs = nrow(sacs_dat), fm = "minres", rotate = "oblimin")
    
    load1 <- data.frame(Item = rownames(efa1$loadings), g = round(efa1$loadings[, 1], 3))
    load2 <- data.frame(Item = rownames(efa2$loadings), F1 = round(efa2$loadings[, 1], 3), F2 = round(efa2$loadings[, 2], 3))
    fit_tab <- data.frame(
      Model = c("EFA_1f", "EFA_2f"),
      RMSEA = c(efa1$RMSEA[1], efa2$RMSEA[1]),
      TLI = c(efa1$TLI, efa2$TLI),
      BIC = c(efa1$BIC, efa2$BIC)
    )
    write_dual(load1, "E5_loadings_1factor", efa_dir)
    write_dual(load2, "E6_loadings_2factor", efa_dir)
    write_dual(fit_tab, "E7_fit", efa_dir)
    write_dual(as.data.frame(round(efa2$Phi, 3)), "E8_phi", efa_dir)
  } else {
    write_dual(data.frame(note = "Insufficient complete SACS cases for EFA."), "E_warning", efa_dir)
  }
} else {
  write_dual(data.frame(note = "SACS items missing; EFA skipped."), "E_warning", efa_dir)
}

## ---- 5. CFA comparisons --------------------------------------------------
message("\n[5] SACS CFA")
cfa_dir <- file.path(out_dir, "cfa")

if (all(sacs_items %in% names(df))) {
  ordered_items <- sacs_items
  cfa_data <- df
  cfa_fit_2f <- lavaan::cfa(
    "Centering =~ SACS_1 + SACS_2 + SACS_5 + SACS_6\nTranscending =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10\nCentering ~~ Transcending",
    data = cfa_data,
    estimator = "WLSMV",
    ordered = ordered_items
  )
  cfa_fit_1f <- lavaan::cfa(
    "SAC_single =~ SACS_1 + SACS_2 + SACS_3 + SACS_4 + SACS_5 + SACS_6 + SACS_7 + SACS_8 + SACS_9 + SACS_10",
    data = cfa_data,
    estimator = "WLSMV",
    ordered = ordered_items
  )
  fit_measures <- c("chisq","df","pvalue","cfi","tli","rmsea","srmr")
  fit_2 <- lavaan::fitMeasures(cfa_fit_2f, fit_measures)
  fit_1 <- lavaan::fitMeasures(cfa_fit_1f, fit_measures)
  fit_tab <- data.frame(
    Model = c("Two-factor", "One-factor"),
    ChiSq = c(fit_2["chisq"], fit_1["chisq"]),
    df = c(fit_2["df"], fit_1["df"]),
    pvalue = c(fit_2["pvalue"], fit_1["pvalue"]),
    CFI = c(fit_2["cfi"], fit_1["cfi"]),
    TLI = c(fit_2["tli"], fit_1["tli"]),
    RMSEA = c(fit_2["rmsea"], fit_1["rmsea"]),
    SRMR = c(fit_2["srmr"], fit_1["srmr"])
  )
  write_dual(fit_tab, "C1_fit_comparison", cfa_dir)
  
  std_2f <- lavaan::inspect(cfa_fit_2f, "std")
  load_df <- data.frame(
    Item = rownames(std_2f$lambda),
    Centering = std_2f$lambda[, "Centering"],
    Transcending = std_2f$lambda[, "Transcending"]
  )
  if (!is.null(std_2f$psi)) {
    factor_corr <- std_2f$psi["Centering", "Transcending"]
    write_dual(data.frame(Factor_Correlation = factor_corr), "C3_factor_correlation", cfa_dir)
  }
  write_dual(load_df, "C2_standardized_loadings", cfa_dir)
} else {
  write_dual(data.frame(note = "SACS items missing; CFA skipped."), "C_warning", cfa_dir)
}

## ---- 6. Correlations, VIF, regressions, bootstraps ----------------------
message("\n[6] Correlations and regressions")
cor_dir <- file.path(out_dir, "correlations_regressions")

if (length(score_vars) >= 2) {
  numeric_data <- df[, score_vars]
  cor_mat <- suppressWarnings(stats::cor(numeric_data, use = "pairwise.complete.obs"))
  write_dual(data.frame(Variable = rownames(cor_mat), cor_mat, check.names = FALSE),
             "CR1_correlation_matrix", cor_dir)
  
  combs <- t(combn(score_vars, 2))
  cor_long <- lapply(seq_len(nrow(combs)), function(i) {
    v1 <- combs[i, 1]; v2 <- combs[i, 2]
    x <- df[[v1]]; y <- df[[v2]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NULL)
    ct <- stats::cor.test(x[ok], y[ok])
    data.frame(var1 = v1, var2 = v2, N = sum(ok), r = ct$estimate, p_value = ct$p.value)
  })
  cor_long <- do.call(rbind, cor_long)
  if (!is.null(cor_long)) write_dual(cor_long, "CR2_correlation_long", cor_dir)
  
  # Bootstrap correlations
  boot_fn <- function(data, indices, v1, v2) {
    d <- data[indices, , drop = FALSE]
    x <- d[[v1]]; y <- d[[v2]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NA_real_)
    stats::cor(x[ok], y[ok])
  }
  boot_rows <- list()
  for (i in seq_len(nrow(combs))) {
    v1 <- combs[i, 1]; v2 <- combs[i, 2]
    bt <- boot::boot(df, statistic = function(d, idx) boot_fn(d, idx, v1, v2), R = 5000)
    vals <- bt$t[is.finite(bt$t)]
    if (!length(vals)) next
    ci <- stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE)
    boot_rows[[length(boot_rows) + 1]] <- data.frame(var1 = v1, var2 = v2,
                                                     mean_r = mean(vals),
                                                     ci_lower = ci[[1]], ci_upper = ci[[2]])
  }
  if (length(boot_rows)) write_dual(do.call(rbind, boot_rows), "CR3_correlation_bootstrap", cor_dir)
}

# VIF check for SACS facets predicting DASS_Total_x2 (with covariates)
if (all(c("SACS_Centering","SACS_Transcending","DASS_Total_x2") %in% names(df))) {
  preds <- c("SACS_Centering", "SACS_Transcending")
  if ("Age_clean" %in% names(df)) preds <- c(preds, "Age_clean")
  if ("Gender_label" %in% names(df)) preds <- c(preds, "Gender_label")
  formula_vif <- as.formula(paste("DASS_Total_x2 ~", paste(preds, collapse = " + ")))
  vif_model <- try(stats::lm(formula_vif, data = df), silent = TRUE)
  if (!inherits(vif_model, "try-error")) {
    vif_vals <- car::vif(vif_model)
    if (is.matrix(vif_vals)) {
      adj <- vif_vals[, "GVIF"]^(1/(2 * vif_vals[, "Df"]))
      vif_df <- data.frame(Variable = rownames(vif_vals), VIF = adj)
    } else {
      vif_df <- data.frame(Variable = names(vif_vals), VIF = vif_vals)
    }
    write_dual(vif_df, "CR4_vif_facets", cor_dir)
  }
}

# Regression models -------------------------------------------------------
outcomes <- intersect(c("DASS_Total_x2", "WEMWBS_Total", "CFQ_Total", "ATQ_Believability"), names(df))
for (outcome in outcomes) {
  formula_reg <- as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending",
                                  if ("Age_clean" %in% names(df)) " + Age_clean" else "",
                                  if ("Gender_label" %in% names(df)) " + Gender_label" else ""))
  fit <- stats::lm(formula_reg, data = df)
  coef_tab <- summary(fit)$coefficients
  coef_df <- data.frame(term = rownames(coef_tab), coef_tab, row.names = NULL)
  if (requireNamespace("lm.beta", quietly = TRUE)) {
    beta_fit <- try(lm.beta::lm.beta(fit), silent = TRUE)
    if (!inherits(beta_fit, "try-error")) {
      std_coefs <- summary(beta_fit)$coefficients
      coef_df$Std_Beta <- std_coefs[rownames(coef_tab), "Standardized"]
    }
  }
  write_dual(coef_df, paste0("CR5_regression_", outcome), cor_dir)
  
  cc <- df[complete.cases(df[, all.vars(formula_reg)]), ]
  if (nrow(cc) > 10) {
    full_model <- stats::lm(formula_reg, data = cc)
    form_trans <- update(formula_reg, paste(". ~ SACS_Transcending",
                                            if ("Age_clean" %in% names(df)) " + Age_clean" else "",
                                            if ("Gender_label" %in% names(df)) " + Gender_label" else ""))
    form_center <- update(formula_reg, paste(". ~ SACS_Centering",
                                             if ("Age_clean" %in% names(df)) " + Age_clean" else "",
                                             if ("Gender_label" %in% names(df)) " + Gender_label" else ""))
    model_trans <- stats::lm(form_trans, data = cc)
    model_center <- stats::lm(form_center, data = cc)
    r2_table <- data.frame(
      Model = c("Transcending_only", "Centering_only", "Both_facets"),
      R_squared = c(summary(model_trans)$r.squared,
                    summary(model_center)$r.squared,
                    summary(full_model)$r.squared)
    )
    write_dual(r2_table, paste0("CR6_r2_", outcome), cor_dir)
  }
}

## ---- 7. Mediation analyses ----------------------------------------------
message("\n[7] Mediation bootstraps")
med_dir <- file.path(out_dir, "mediation")

doit_mediation <- function(X, M, Y, covars = NULL, n_boot = 5000) {
  vars <- c(X, M, Y, covars)
  med_data <- df[complete.cases(df[, vars]), vars, drop = FALSE]
  if (nrow(med_data) < 20) return(NULL)
  form_a <- as.formula(paste(M, "~", paste(c(X, covars), collapse = " + ")))
  form_b <- as.formula(paste(Y, "~", paste(c(X, M, covars), collapse = " + ")))
  form_c <- as.formula(paste(Y, "~", paste(c(X, covars), collapse = " + ")))
  fit_a <- stats::lm(form_a, data = med_data)
  fit_b <- stats::lm(form_b, data = med_data)
  fit_c <- stats::lm(form_c, data = med_data)
  a <- coef(fit_a)[X]
  b <- coef(fit_b)[M]
  c_prime <- coef(fit_b)[X]
  c_total <- coef(fit_c)[X]
  idx <- seq_len(nrow(med_data))
  ab_vals <- replicate(n_boot, {
    samp <- sample(idx, replace = TRUE)
    d <- med_data[samp, , drop = FALSE]
    fa <- stats::lm(form_a, data = d)
    fb <- stats::lm(form_b, data = d)
    a_s <- try(coef(fa)[X], silent = TRUE)
    b_s <- try(coef(fb)[M], silent = TRUE)
    if (inherits(a_s, "try-error") || inherits(b_s, "try-error")) return(NA_real_)
    as.numeric(a_s) * as.numeric(b_s)
  })
  ab_vals <- ab_vals[is.finite(ab_vals)]
  if (!length(ab_vals)) return(NULL)
  ci <- stats::quantile(ab_vals, c(0.025, 0.975))
  data.frame(
    Predictor = X,
    Mediator = M,
    Outcome = Y,
    Indirect = a * b,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Total = c_total,
    Direct = c_prime,
    Prop_Mediated = if (!is.na(c_total) && c_total != 0) (a * b)/c_total else NA_real_,
    N = nrow(med_data)
  )
}

covariates <- intersect(c("Age_clean", "Gender_label"), names(df))
med_rows <- list()
for (M in intersect(c("CFQ_Total", "ATQ_Believability"), names(df))) {
  for (Y in intersect(c("DASS_Total_x2", "WEMWBS_Total"), names(df))) {
    res <- doit_mediation("SACS_Total", M, Y, covars = covariates)
    if (!is.null(res)) med_rows[[length(med_rows) + 1]] <- res
  }
}
if (length(med_rows)) write_dual(do.call(rbind, med_rows), "MED1_results", med_dir)

## ---- 8. Multiple imputation ---------------------------------------------
message("\n[8] Multiple imputation sensitivity")
mi_dir <- file.path(out_dir, "multiple_imputation")

mi_vars <- intersect(c(
  "SACS_Centering","SACS_Transcending",
  "CFQ_Total","ATQ_Believability","ATQ_Frequency","ATQ_Combined",
  "WEMWBS_Total",
  "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2",
  "Age_clean","Gender_label"
), names(df))

if (length(mi_vars) >= 3) {
  imp <- try(mice::mice(df[, mi_vars], m = 5, seed = 100, printFlag = FALSE), silent = TRUE)
  if (!inherits(imp, "try-error")) {
    complete1 <- mice::complete(imp, 1)
    score_vars_imp <- intersect(score_vars, names(complete1))
    if (length(score_vars_imp) >= 2) {
      cor_imp <- stats::cor(complete1[, score_vars_imp], use = "pairwise.complete.obs")
      write_dual(data.frame(Variable = rownames(cor_imp), cor_imp, check.names = FALSE),
                 "MI1_correlation_matrix", mi_dir)
    }
    for (outcome in outcomes) {
      form <- as.formula(paste(outcome, "~ SACS_Centering + SACS_Transcending",
                               if ("Age_clean" %in% names(df)) " + Age_clean" else "",
                               if ("Gender_label" %in% names(df)) " + Gender_label" else ""))
      fit_imp <- try(mice::with(imp, stats::lm(form)), silent = TRUE)
      if (!inherits(fit_imp, "try-error")) {
        pooled <- mice::pool(fit_imp)
        write_dual(summary(pooled), paste0("MI2_regression_", outcome), mi_dir)
      }
    }
  } else {
    write_dual(data.frame(note = "mice failed; check data for perfect collinearity."), "MI_error", mi_dir)
  }
} else {
  write_dual(data.frame(note = "Not enough variables for MI."), "MI_warning", mi_dir)
}

## ---- 9. Bayesian SEM with fallback --------------------------------------
message("\n[9] Bayesian SEM")
sem_dir <- file.path(out_dir, "sem")
dir.create(sem_dir, showWarnings = FALSE)

sem_items_center <- c("SACS_1", "SACS_2", "SACS_5", "SACS_6")
sem_items_trans  <- c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_9", "SACS_10")
sem_required <- c(sem_items_center, sem_items_trans,
                  "CFQ_Total", "ATQ_Believability", "WEMWBS_Total", "DASS_Total_x2")
if (all(sem_required %in% names(df))) {
  sem_data <- df[, sem_required]
  for (item in c(sem_items_center, sem_items_trans)) sem_data[[item]] <- as.ordered(sem_data[[item]])
  sem_complete <- sem_data[complete.cases(sem_data), ]
  if (nrow(sem_complete) >= 80) {
    sem_model <- "
      SAC_Center =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
      SAC_Trans  =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_9 + SACS_10

      CFQ_Total          ~ a_c*SAC_Center + a_t*SAC_Trans
      ATQ_Believability  ~ g_c*SAC_Center + g_t*SAC_Trans

      DASS_Total_x2 ~ b1*CFQ_Total + b2*ATQ_Believability + c_c*SAC_Center + c_t*SAC_Trans
      WEMWBS_Total  ~ d1*CFQ_Total + d2*ATQ_Believability + e_c*SAC_Center + e_t*SAC_Trans

      CFQ_Total ~~ ATQ_Believability
      SAC_Center ~~ SAC_Trans

      ind_cfq_dass_c := a_c*b1
      ind_cfq_dass_t := a_t*b1
      ind_atq_dass_c := g_c*b2
      ind_atq_dass_t := g_t*b2
      ind_cfq_well_c := a_c*d1
      ind_cfq_well_t := a_t*d1
      ind_atq_well_c := g_c*d2
      ind_atq_well_t := g_t*d2
      ind_total_dass_c := ind_cfq_dass_c + ind_atq_dass_c
      ind_total_dass_t := ind_cfq_dass_t + ind_atq_dass_t
      ind_total_well_c := ind_cfq_well_c + ind_atq_well_c
      ind_total_well_t := ind_cfq_well_t + ind_atq_well_t
    "
    
    run_bayes <- function(target) {
      try(blavaan::bsem(
        model = sem_model,
        data = sem_complete,
        ordered = c(sem_items_center, sem_items_trans),
        target = target,
        n.chains = 3,
        burnin = 4000,
        sample = 8000,
        seed = 1234,
        bcontrol = list(cores = 2)
      ), silent = TRUE)
    }
    
    fit_bayes <- run_bayes("stan")
    if (inherits(fit_bayes, "try-error") || is.null(fit_bayes)) fit_bayes <- run_bayes("jags")
    
    if (!inherits(fit_bayes, "try-error") && !is.null(fit_bayes)) {
      capture.output(summary(fit_bayes, standardized = TRUE),
                     file = file.path(sem_dir, "SEM1_bayes_summary.txt"))
      params <- lavaan::parameterEstimates(fit_bayes, standardized = TRUE)
      write_dual(params, "SEM2_bayes_parameters", sem_dir)
      psrf <- try(blavaan::blavInspect(fit_bayes, "psrf"), silent = TRUE)
      if (!inherits(psrf, "try-error") && !is.null(psrf)) {
        write_dual(as.data.frame(psrf), "SEM3_psrf", sem_dir)
      }
      ppp <- try(blavaan::blavInspect(fit_bayes, "ppp"), silent = TRUE)
      if (!inherits(ppp, "try-error")) {
        write_dual(data.frame(ppp = ppp), "SEM4_ppp", sem_dir)
      }
      draws <- try(blavaan::blavInspect(fit_bayes, "mcmc"), silent = TRUE)
      if (!inherits(draws, "try-error") && !is.null(draws)) {
        chain <- if (is.list(draws)) do.call(rbind, lapply(draws, as.matrix)) else as.matrix(draws)
        ind_cols <- grep("^ind_", colnames(chain), value = TRUE)
        if (length(ind_cols)) {
          ind_summary <- data.frame(
            parameter = ind_cols,
            mean = apply(chain[, ind_cols, drop = FALSE], 2, mean),
            median = apply(chain[, ind_cols, drop = FALSE], 2, median),
            sd = apply(chain[, ind_cols, drop = FALSE], 2, sd),
            lower_95 = apply(chain[, ind_cols, drop = FALSE], 2, quantile, 0.025),
            upper_95 = apply(chain[, ind_cols, drop = FALSE], 2, quantile, 0.975)
          )
          write_dual(ind_summary, "SEM5_indirect_effects", sem_dir)
        }
      }
    } else {
      fit_ml <- lavaan::sem(
        sem_model,
        data = sem_complete,
        ordered = c(sem_items_center, sem_items_trans),
        estimator = "WLSMV",
        se = "robust"
      )
      capture.output(summary(fit_ml, fit.measures = TRUE, standardized = TRUE),
                     file = file.path(sem_dir, "SEM_ML_summary.txt"))
      write_dual(lavaan::parameterEstimates(fit_ml, standardized = TRUE),
                 "SEM_ML_parameters", sem_dir)
    }
  } else {
    write_dual(data.frame(note = "Insufficient complete cases for SEM."), "SEM_warning", sem_dir)
  }
} else {
  write_dual(data.frame(note = "Required SEM variables missing."), "SEM_warning", sem_dir)
}

## ---- 10. Session information -------------------------------------------
message("\n[10] Session info saved")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("\nMaster SAC analysis complete. Outputs in: ", out_dir)