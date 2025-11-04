## =========================================================
## Stage M: Missing-data diagnostics (MCAR / MAR / MNAR)
## DSUR workspace compatible; robust to separation/collinearity
## + Multivariate normality tests (Henze–Zirkler, Mardia), univariate Shapiro
## + NEW: Item-level MCAR (Little) per scale with md.pattern outputs
## =========================================================

# ---- 0) Paths, packages, helpers ----
root <- if (!exists("root")) "/Users/jamesjackson/Desktop/dir" else root
phase_final <- file.path(root, "phase_final_analysis")
if (!dir.exists(phase_final)) dir.create(phase_final, recursive = TRUE, showWarnings = FALSE)

# write_dual: reuse if present, else define
if (!exists("write_dual")) {
  xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)
  write_dual <- function(df, stem, out_dir = phase_final) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    csv_path <- file.path(out_dir, paste0(stem, ".csv"))
    utils::write.csv(df, csv_path, row.names = FALSE)
    if (xlsx_ok) openxlsx::write.xlsx(df, file.path(out_dir, paste0(stem, ".xlsx")), overwrite = TRUE)
    invisible(csv_path)
  }
} else {
  xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)
}

need <- function(pkgs) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) try(utils::install.packages(p), silent = TRUE)
}
need(c("mice","naniar","MissMech","brglm2","glmnet","MVN","energy","moments"))

safelib <- function(pkg) { if (requireNamespace(pkg, quietly = TRUE)) suppressPackageStartupMessages(library(pkg, character.only = TRUE)) }
safelib("mice")  # md.pattern, imputation utilities

grab1 <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else unname(as.numeric(x)[1])

# ---- 1) Load data using your precedence rules ----
sacs_items   <- sprintf("SACS_%d", 1:10)
cfq_items    <- sprintf("CFQ_%d", 1:7)
atq_f_items  <- sprintf("ATQ_%df", 1:15)
atq_b_items  <- sprintf("ATQ_%db", 1:15)
wemwbs_items <- sprintf("WEMWBS_%d", 1:14)
dass_items   <- sprintf("DASS.21_%d", 1:21)

scored_file <- file.path(root, "phase03_scoring", "scored_primary.csv")
if (file.exists(scored_file)) {
  df <- utils::read.csv(scored_file, check.names = TRUE, stringsAsFactors = FALSE)
} else {
  raw_xlsx <- file.path(root, "data", "raw copy.xlsx")
  if (!file.exists(raw_xlsx)) stop("No data file found for missingness diagnostics.")
  if (!xlsx_ok) stop("openxlsx is required to read the raw Excel file.")
  df <- as.data.frame(openxlsx::read.xlsx(raw_xlsx, sheet = 1))
  names(df) <- sub("^DASS-21_", "DASS.21_", names(df))
}

# analysis variables: totals + demographics if present
score_vars <- intersect(c(
  "SACS_Total","CFQ_Total","ATQ_Believability","ATQ_Frequency","ATQ_Combined",
  "WEMWBS_Total","DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2",
  "Age_clean","Gender_label"
), names(df))

miss_dir <- file.path(phase_final, "missingness")
if (!dir.exists(miss_dir)) dir.create(miss_dir, showWarnings = FALSE)

# ---- 2) Descriptive missingness summaries (totals) ----
miss_summ <- data.frame(
  variable = score_vars,
  n_missing = sapply(score_vars, function(v) sum(is.na(df[[v]]))),
  pct_missing = round(100 * sapply(score_vars, function(v) mean(is.na(df[[v]]))), 2),
  stringsAsFactors = FALSE
)
write_dual(miss_summ, "M0_missing_summary", miss_dir)

md_pat <- tryCatch(mice::md.pattern(df[score_vars], plot = FALSE), error = function(e) NULL)
if (is.null(md_pat)) {
  write_dual(data.frame(note = "md.pattern failed; likely no missingness or singular patterns."),
             "M1_missing_patterns", miss_dir)
} else {
  write_dual(data.frame(md_pat, check.names = FALSE), "M1_missing_patterns", miss_dir)
}

# ---- 3) MCAR test (Little-style) on totals w/ robust extraction ----
X_mcar <- df[score_vars]
if ("Gender_label" %in% names(X_mcar)) X_mcar$Gender_label <- as.numeric(as.factor(X_mcar$Gender_label))
keep <- vapply(X_mcar, function(x) sum(!is.na(x)) >= 3 && length(unique(na.omit(x))) >= 2, logical(1))
X_mcar <- X_mcar[, keep, drop = FALSE]

mcar_written <- FALSE
if (requireNamespace("naniar", quietly = TRUE)) {
  res <- try(naniar::mcar_test(X_mcar), silent = TRUE)
  if (!inherits(res, "try-error")) {
    if (is.data.frame(res)) {
      stat <- grab1(res$statistic); df1 <- grab1(if (!is.null(res$df)) res$df else res$parameter); p <- grab1(res$p.value)
    } else {
      stat <- grab1(res$statistic); df1 <- grab1(res$parameter); p <- grab1(res$p.value)
    }
    out <- data.frame(method = "naniar::mcar_test (Little)", stat = stat, df = df1, p_value = p)
    write_dual(out, "M2_LittleMCAR_result", miss_dir)
    mcar_written <- TRUE
  }
}
if (!mcar_written && requireNamespace("MissMech", quietly = TRUE)) {
  res <- try(MissMech::TestMCARNormality(as.matrix(X_mcar)), silent = TRUE)
  if (!inherits(res, "try-error")) {
    out <- data.frame(method = "MissMech::TestMCARNormality",
                      stat = grab1(res$chi.square), df = grab1(res$df), p_value = grab1(res$p.value))
    write_dual(out, "M2_LittleMCAR_result", miss_dir)
    mcar_written <- TRUE
  }
}
if (!mcar_written) {
  write_dual(data.frame(note = "No MCAR test ran. Install/update `naniar` or `MissMech`."),
             "M2_LittleMCAR_ERROR", miss_dir)
}

# ---- 3b) Multivariate and univariate normality diagnostics (totals) ----
mvn_written <- FALSE
if (ncol(X_mcar) > 0) {
  X_num <- X_mcar[, vapply(X_mcar, is.numeric, logical(1)), drop = FALSE]
  X_cc  <- stats::na.omit(X_num)
  if (nrow(X_cc) >= 10 && ncol(X_cc) >= 2) {
    if (requireNamespace("MVN", quietly = TRUE)) {
      hz <- try(MVN::mvn(X_cc, mvnTest = "hz", covariance = TRUE, desc = FALSE), silent = TRUE)
      md <- try(MVN::mvn(X_cc, mvnTest = "mardia", covariance = TRUE, desc = FALSE), silent = TRUE)
      rows <- list()
      if (!inherits(hz, "try-error") && is.list(hz$HZ)) {
        rows[["HZ"]] <- data.frame(test = "Henze-Zirkler",
                                   stat = grab1(hz$HZ$HZ),
                                   df = NA_real_,
                                   p_value = grab1(hz$HZ$p.value),
                                   n_complete = nrow(X_cc),
                                   p_vars = ncol(X_cc))
      }
      if (!inherits(md, "try-error") && is.list(md$multivariateNormality)) {
        mtab <- md$multivariateNormality
        for (i in seq_len(nrow(mtab))) {
          rows[[paste0("Mardia_", i)]] <- data.frame(
            test = as.character(mtab[i, 1]),
            stat = as.numeric(mtab[i, 2]),
            df = NA_real_,
            p_value = as.numeric(mtab[i, 4]),
            n_complete = nrow(X_cc),
            p_vars = ncol(X_cc)
          )
        }
      }
      if (length(rows)) {
        mvn_df <- do.call(rbind, rows)
        write_dual(mvn_df, "M2a_MVN_multivariate_normality", miss_dir)
        mvn_written <- TRUE
      }
    }
    if (!mvn_written && requireNamespace("energy", quietly = TRUE)) {
      et <- try(energy::mvnorm.etest(as.matrix(X_cc), R = 1000), silent = TRUE)
      if (!inherits(et, "try-error")) {
        mvn_df <- data.frame(test = "Energy E-test",
                             stat = grab1(et$statistic),
                             df = NA_real_,
                             p_value = grab1(et$p.value),
                             n_complete = nrow(X_cc),
                             p_vars = ncol(X_cc))
        write_dual(mvn_df, "M2a_MVN_multivariate_normality", miss_dir)
        mvn_written <- TRUE
      }
    }
  } else {
    write_dual(data.frame(note = "MVN tests skipped: need >=10 complete cases and >=2 numeric vars.",
                          n_complete = nrow(X_cc), p_vars = ncol(X_cc)),
               "M2a_MVN_multivariate_normality", miss_dir)
  }
  
  # Univariate Shapiro–Wilk per variable
  uni_rows <- list()
  for (v in colnames(X_num)) {
    xv <- stats::na.omit(X_num[[v]])
    if (length(xv) >= 3 && length(xv) <= 5000) {
      sw <- try(stats::shapiro.test(xv), silent = TRUE)
      if (!inherits(sw, "try-error")) {
        sk <- try(moments::skewness(xv), silent = TRUE); if (inherits(sk, "try-error")) sk <- NA_real_
        kt <- try(moments::kurtosis(xv) - 3, silent = TRUE); if (inherits(kt, "try-error")) kt <- NA_real_
        uni_rows[[v]] <- data.frame(variable = v,
                                    n = length(xv),
                                    shapiro_W = grab1(sw$statistic),
                                    shapiro_p = grab1(sw$p.value),
                                    skewness = as.numeric(sk),
                                    excess_kurtosis = as.numeric(kt))
      }
    }
  }
  if (length(uni_rows)) {
    write_dual(do.call(rbind, uni_rows), "M2b_univariate_normality", miss_dir)
  } else {
    write_dual(data.frame(note = "Univariate normality skipped: insufficient sample sizes."),
               "M2b_univariate_normality", miss_dir)
  }
} else {
  write_dual(data.frame(note = "No variables available for MVN diagnostics."), "M2a_MVN_multivariate_normality", miss_dir)
  write_dual(data.frame(note = "No variables available for univariate normality."), "M2b_univariate_normality", miss_dir)
}

# ---- 3c) NEW: Item-level MCAR (per scale) + md.pattern ----
# Helper to run Little's test on an item set and write results
run_mcar_items <- function(data, items, scale_name, out_dir = miss_dir) {
  items <- items[items %in% names(data)]
  if (length(items) < 2) {
    write_dual(data.frame(note = paste0("Not enough items present for ", scale_name)), 
               paste0("M1i_", scale_name, "_mdpattern"), out_dir)
    write_dual(data.frame(note = paste0("MCAR skipped for ", scale_name)), 
               paste0("M2i_", scale_name, "_LittleMCAR_ERROR"), out_dir)
    return(NULL)
  }
  # Write md.pattern for the item block
  mdp <- tryCatch(mice::md.pattern(data[items], plot = FALSE), error = function(e) NULL)
  if (is.null(mdp)) {
    write_dual(data.frame(note = paste0("md.pattern failed for ", scale_name)), 
               paste0("M1i_", scale_name, "_mdpattern"), out_dir)
  } else {
    write_dual(data.frame(mdp, check.names = FALSE), 
               paste0("M1i_", scale_name, "_mdpattern"), out_dir)
  }
  
  # Little's test using naniar, fallback to MissMech
  X <- data[items]
  # coerce factors/characters to numeric where needed (ordinal Likert assumed already numeric in scored data)
  for (nm in names(X)) if (is.factor(X[[nm]])) X[[nm]] <- as.numeric(X[[nm]])
  mcar_ok <- FALSE
  if (requireNamespace("naniar", quietly = TRUE)) {
    res <- try(naniar::mcar_test(X), silent = TRUE)
    if (!inherits(res, "try-error")) {
      if (is.data.frame(res)) {
        stat <- grab1(res$statistic); df1 <- grab1(if (!is.null(res$df)) res$df else res$parameter); p <- grab1(res$p.value)
      } else {
        stat <- grab1(res$statistic); df1 <- grab1(res$parameter); p <- grab1(res$p.value)
      }
      write_dual(data.frame(scale = scale_name, method = "naniar::mcar_test (Little)", stat = stat, df = df1, p_value = p),
                 paste0("M2i_", scale_name, "_LittleMCAR_result"), out_dir)
      mcar_ok <- TRUE
    }
  }
  if (!mcar_ok && requireNamespace("MissMech", quietly = TRUE)) {
    res <- try(MissMech::TestMCARNormality(as.matrix(X)), silent = TRUE)
    if (!inherits(res, "try-error")) {
      write_dual(data.frame(scale = scale_name, method = "MissMech::TestMCARNormality",
                            stat = grab1(res$chi.square), df = grab1(res$df), p_value = grab1(res$p.value)),
                 paste0("M2i_", scale_name, "_LittleMCAR_result"), out_dir)
      mcar_ok <- TRUE
    }
  }
  if (!mcar_ok) {
    write_dual(data.frame(scale = scale_name, note = "No item-level MCAR test ran (package/shape issue)."),
               paste0("M2i_", scale_name, "_LittleMCAR_ERROR"), out_dir)
  }
  invisible(NULL)
}

# Run per-scale item MCAR & patterns (only if at least some items present)
scale_sets <- list(
  SACS = sacs_items,
  CFQ = cfq_items,
  ATQ_Frequency = atq_f_items,
  ATQ_Believability = atq_b_items,
  WEMWBS = wemwbs_items,
  DASS21 = dass_items
)
for (nm in names(scale_sets)) run_mcar_items(df, scale_sets[[nm]], nm, miss_dir)

# Build a compact index of per-scale item-level Little’s results
idx_rows <- list()
for (nm in names(scale_sets)) {
  res_path <- file.path(miss_dir, paste0("M2i_", nm, "_LittleMCAR_result.csv"))
  err_path <- file.path(miss_dir, paste0("M2i_", nm, "_LittleMCAR_ERROR.csv"))
  if (file.exists(res_path)) {
    tb <- try(utils::read.csv(res_path, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(tb, "try-error") && nrow(tb) >= 1) {
      idx_rows[[nm]] <- data.frame(scale = nm, stat = tb$stat[1], df = tb$df[1], p_value = tb$p_value[1])
    }
  } else if (file.exists(err_path)) {
    idx_rows[[nm]] <- data.frame(scale = nm, stat = NA_real_, df = NA_real_, p_value = NA_real_, note = "ERROR")
  } else {
    idx_rows[[nm]] <- data.frame(scale = nm, stat = NA_real_, df = NA_real_, p_value = NA_real_, note = "NOT_RUN")
  }
}
if (length(idx_rows)) {
  write_dual(do.call(rbind, idx_rows), "M2i_item_level_Little_index", miss_dir)
}

# ---- 4) MAR diagnostics with robust fitting pipeline (totals) ----
drop_constant <- function(X) {
  keep <- vapply(X, function(x) length(unique(na.omit(x))) > 1, logical(1))
  X[, keep, drop = FALSE]
}
drop_duplicate <- function(X) {
  if (ncol(X) < 2) return(X)
  X[, !duplicated(lapply(seq_len(ncol(X)), function(i) paste(X[, i], collapse = "\r"))), drop = FALSE]
}
drop_high_corr <- function(X, thr = 0.98) {
  num <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
  if (ncol(num) < 2) return(X)
  R <- suppressWarnings(stats::cor(num, use = "pairwise.complete.obs"))
  R[is.na(R)] <- 0
  remove <- character(0)
  for (j in seq_len(ncol(R))) {
    if (colnames(R)[j] %in% remove) next
    high <- which(abs(R[, j]) > thr)
    high <- setdiff(high, j)
    if (length(high)) remove <- union(remove, colnames(R)[high])
  }
  X[, setdiff(colnames(X), remove), drop = FALSE]
}

vars_for_mar <- setdiff(score_vars, "Gender_label")
mar_results <- list()
have_brglm2 <- requireNamespace("brglm2", quietly = TRUE)
have_glmnet <- requireNamespace("glmnet", quietly = TRUE)

for (V in vars_for_mar) {
  Rv <- as.integer(is.na(df[[V]]))
  if (length(unique(na.omit(Rv))) < 2) next
  
  preds <- setdiff(score_vars, V)
  preds <- preds[sapply(preds, function(p) mean(is.na(df[[p]])) < 0.5)]
  if (!length(preds)) next
  
  P0 <- df[, preds, drop = FALSE]
  if ("Gender_label" %in% colnames(P0)) P0$Gender_label <- as.factor(P0$Gender_label)
  
  imp <- suppressWarnings(mice::mice(P0, m = 1, maxit = 5, printFlag = FALSE))
  X <- mice::complete(imp, 1)
  
  X <- drop_constant(X)
  X <- drop_duplicate(X)
  X <- drop_high_corr(X, thr = 0.98)
  if (!ncol(X)) next
  
  mm <- model.matrix(~ . , data = X)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  
  fit_ok <- FALSE
  fit <- try(glm(Rv ~ ., data = X, family = binomial(), control = glm.control(maxit = 200)), silent = TRUE)
  if (!inherits(fit, "try-error") && isTRUE(fit$converged) && all(is.finite(coef(fit)))) fit_ok <- TRUE
  
  if (fit_ok) {
    coefs <- summary(fit)$coef
    out <- data.frame(
      variable_missing = V,
      predictor = rownames(coefs),
      estimate = coefs[, "Estimate"],
      se = coefs[, "Std. Error"],
      z = coefs[, "z value"],
      p = coefs[, "Pr(>|z|)"],
      row.names = NULL
    )
    mar_results[[V]] <- out[out$predictor != "(Intercept)", ]
    next
  }
  
  if (have_brglm2) {
    fit2 <- try(brglm2::brglm(Rv ~ ., data = X, family = stats::binomial("meanBR")), silent = TRUE)
    if (!inherits(fit2, "try-error")) {
      s2 <- summary(fit2)
      coefs <- s2$coef
      out <- data.frame(
        variable_missing = V,
        predictor = rownames(coefs),
        estimate = coefs[, "Estimate"],
        se = coefs[, "Std. Error"],
        z = coefs[, "z value"],
        p = coefs[, "Pr(>|z|)"],
        row.names = NULL
      )
      mar_results[[V]] <- out[out$predictor != "(Intercept)", ]
      next
    }
  }
  
  if (have_glmnet) {
    y <- Rv
    cv <- glmnet::cv.glmnet(mm, y, family = "binomial", alpha = 0, nfolds = 5)
    b  <- as.matrix(glmnet::coef.glmnet(cv$glmnet.fit, s = cv$lambda.1se))
    est <- data.frame(
      variable_missing = V,
      predictor = rownames(b)[-1],
      estimate = as.numeric(b[-1, 1]),
      se = NA_real_, z = NA_real_, p = NA_real_,
      row.names = NULL
    )
    mar_results[[V]] <- est
  }
}

if (length(mar_results)) {
  mar_df <- do.call(rbind, mar_results)
  if (!all(is.na(mar_df$p))) {
    mar_df$q_fdr <- p.adjust(mar_df$p, method = "fdr")
    mar_df$evidence_FOR_not_MCAR <- !is.na(mar_df$q_fdr) & mar_df$q_fdr < 0.10
  } else {
    mar_df$q_fdr <- NA_real_
    mar_df$evidence_FOR_not_MCAR <- NA
    mar_df$note <- "Ridge fallback used; no p-values. Interpret magnitudes only."
  }
  write_dual(mar_df, "M3_MAR_logit_missingness_on_observed", miss_dir)
} else {
  write_dual(data.frame(note = "Insufficient variability to fit MAR models."),
             "M3_MAR_logit_missingness_on_observed", miss_dir)
}

# ---- 5) MNAR plausibility: pattern-mixture sensitivity (totals) ----
pmix_vars <- score_vars[sapply(score_vars, function(v) mean(is.na(df[[v]])) > 0 & mean(is.na(df[[v]])) < 0.9)]
if (length(pmix_vars) > 0) {
  meth <- mice::make.method(df[score_vars]); meth[meth != ""] <- "pmm"
  pred <- mice::make.predictorMatrix(df[score_vars]); diag(pred) <- 0
  imp_all <- suppressWarnings(mice(df[score_vars], m = 5, maxit = 10, method = meth,
                                   predictorMatrix = pred, printFlag = FALSE))
  mids <- mice::complete(imp_all, "all")
  
  sens_rows <- list()
  for (V in pmix_vars) {
    obs_vals <- df[[V]][!is.na(df[[V]])]
    miss_idx <- which(is.na(df[[V]]))
    if (length(obs_vals) < 5 || length(miss_idx) < 5) next
    imps <- unlist(lapply(mids, function(d) d[[V]][miss_idx]))
    sens_rows[[V]] <- data.frame(
      variable = V,
      n_obs = sum(!is.na(df[[V]])),
      n_miss = sum(is.na(df[[V]])),
      mean_obs = mean(obs_vals),
      mean_imp_for_miss = mean(imps),
      delta_imp_minus_obs = mean(imps) - mean(obs_vals),
      sd_obs = stats::sd(obs_vals),
      sd_imp_for_miss = stats::sd(imps)
    )
  }
  if (length(sens_rows)) {
    sens_df <- do.call(rbind, sens_rows)
    sens_df$mnar_flag_delta <- with(sens_df, abs(delta_imp_minus_obs) > 0.3 * sd_obs)
    write_dual(sens_df, "M4_pattern_mixture_sensitivity", miss_dir)
  } else {
    write_dual(data.frame(note = "Not enough missingness to compute pattern-mixture deltas."),
               "M4_pattern_mixture_sensitivity", miss_dir)
  }
} else {
  write_dual(data.frame(note = "No variables with both observed and missing values among analysis vars."),
             "M4_pattern_mixture_sensitivity", miss_dir)
}

# ---- 6) Output index ----
conclusion <- data.frame(
  key_files = c(
    "M0_missing_summary.csv",
    "M1_missing_patterns.csv",
    if (file.exists(file.path(miss_dir,"M2_LittleMCAR_result.csv")))
      "M2_LittleMCAR_result.csv" else "M2_LittleMCAR_ERROR.csv",
    "M2a_MVN_multivariate_normality.csv",
    "M2b_univariate_normality.csv",
    # New item-level products:
    "M1i_*_mdpattern.csv",
    "M2i_item_level_Little_index.csv",
    "M3_MAR_logit_missingness_on_observed.csv",
    "M4_pattern_mixture_sensitivity.csv"
  ),
  interpretation_hint = c(
    "Scan pct_missing. Focus >5%.",
    "Pattern structure: monotone vs block.",
    "MCAR supported if p>.05; else reject MCAR (totals).",
    "MVN tests inform MCAR’s normality assumption.",
    "Per-variable normality context.",
    "Item-level md.patterns per scale.",
    "Per-scale Little’s p-values (items).",
    "Significant predictors of R_V support MAR (use q<.10).",
    "Large |delta| vs SD suggests MNAR risk."
  )
)
write_dual(conclusion, "M5_missingness_conclusion_index", miss_dir)
cat("Missingness diagnostics complete. See:", miss_dir, "\n")
