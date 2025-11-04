# ==========================================================
# med_SACS_as_mediator.R
# Runs mediation models where SACS (Centering/Transcending/Total) is the mediator
# Exports one CSV with all results:
#   /Users/jamesjackson/Desktop/dir/mediation_results.csv
# Paths use working.csv at the same directory.
# ==========================================================

# 0) Packages ---------------------------------------------------------------
required_packages <- c("readr","dplyr","tidyr","tibble","purrr","stringr","lavaan")
missing_packages <- setdiff(required_packages, installed.packages()[, "Package"])
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# 1) Paths -----------------------------------------------------------------
in_csv  <- "/Users/jamesjackson/Desktop/dir/working.csv"
out_csv <- "/Users/jamesjackson/Desktop/dir/mediation_results.csv"

if (!file.exists(in_csv)) stop("Input not found: ", in_csv)

# 2) Load data -------------------------------------------------------------
d_all <- readr::read_csv(in_csv, show_col_types = FALSE)

# 3) Helper to run (parallel) mediation and extract results ----------------
run_med <- function(data, x, m_list, y, label, cvars = character(), n_boot = 5000) {
  needed <- unique(c(x, m_list, y, cvars))
  if (!all(needed %in% names(data))) {
    msg <- paste("Skipping:", label, "| missing vars:", paste(setdiff(needed, names(data)), collapse=", "))
    warning(msg)
    return(tibble(model_label = label, x = x, mediator = paste(m_list, collapse = "+"),
                  y = y, N = NA_integer_, a = NA_real_, a_p = NA_real_, b = NA_real_,
                  b_p = NA_real_, cprime = NA_real_, cprime_p = NA_real_, ab = NA_real_,
                  ab_ci_lower = NA_real_, ab_ci_upper = NA_real_, ab_sig = NA,
                  total = NA_real_, total_ci_lower = NA_real_, total_ci_upper = NA_real_,
                  effect = "skipped"))
  }
  
  # Keep only vars of interest (retain NAs for FIML)
  d <- data[, needed, drop = FALSE]
  N <- nrow(d)
  
  # rename to x, y, m1..mk for lavaan convenience
  names_map <- setNames(names(d), names(d))
  names_map[names_map == x] <- "x"
  names_map[names_map == y] <- "y"
  if (length(m_list)) {
    for (i in seq_along(m_list)) names_map[names_map == m_list[i]] <- paste0("m", i)
  }
  colnames(d) <- unname(names_map)
  
  # Build lavaan model
  paths_m <- if (length(m_list)) {
    paste0("m", seq_along(m_list), " ~ a", seq_along(m_list), " * x",
           if (length(cvars)) paste0(" + ", paste(cvars, collapse=" + ")) else "")
  } else character(0)
  
  path_y  <- paste0(
    "y ~ cprime * x",
    if (length(m_list)) paste0(" + ", paste0("b", seq_along(m_list), " * m", seq_along(m_list), collapse = " + ")) else "",
    if (length(cvars)) paste0(" + ", paste(cvars, collapse=" + ")) else ""
  )
  
  effects <- if (length(m_list)) paste0("ab", seq_along(m_list), " := a", seq_along(m_list), " * b", seq_along(m_list)) else character(0)
  total_def <- if (length(m_list)) paste0("total := cprime + ", paste0("ab", seq_along(m_list), collapse = " + ")) else "total := cprime"
  diff_def  <- if (length(m_list) == 2) "diff_ab := ab1 - ab2" else ""
  
  model_syntax <- paste(c(paths_m, path_y, effects, total_def, diff_def), collapse = "\n")
  
  # Robust paths + bootstrap CIs for indirects
  fit_mlr  <- lavaan::sem(model_syntax, data = d, estimator = "MLR", missing = "fiml")
  fit_boot <- lavaan::sem(model_syntax, data = d, se = "bootstrap", bootstrap = n_boot, missing = "fiml")
  
  pe_mlr  <- parameterEstimates(fit_mlr, standardized = TRUE)
  pe_boot <- parameterEstimates(fit_boot, standardized = TRUE, ci = TRUE, boot.ci.type = "bca.simple")
  
  # Grab path coefficients (a's, b's, c')
  get_lab <- function(pe, lbl) pe %>% dplyr::filter(.data$label == lbl)
  get_val <- function(pe, lbl, field) { r <- get_lab(pe, lbl); if (nrow(r)) r[[field]][1] else NA_real_ }
  
  cprime_est <- get_val(pe_mlr, "cprime", "est")
  cprime_p   <- get_val(pe_mlr, "cprime", "pvalue")
  
  total_est  <- get_val(pe_boot, "total", "est")
  total_lo   <- get_val(pe_boot, "total", "ci.lower")
  total_hi   <- get_val(pe_boot, "total", "ci.upper")
  
  # Build rows per mediator
  rows <- list()
  if (length(m_list)) {
    for (k in seq_along(m_list)) {
      a_est <- get_val(pe_mlr, paste0("a", k), "est")
      a_p   <- get_val(pe_mlr, paste0("a", k), "pvalue")
      b_est <- get_val(pe_mlr, paste0("b", k), "est")
      b_p   <- get_val(pe_mlr, paste0("b", k), "pvalue")
      
      ab_row <- get_lab(pe_boot, paste0("ab", k))
      ab_est <- if (nrow(ab_row)) ab_row$est[1] else NA_real_
      ab_lo  <- if (nrow(ab_row)) ab_row$`ci.lower`[1] else NA_real_
      ab_hi  <- if (nrow(ab_row)) ab_row$`ci.upper`[1] else NA_real_
      ab_sig <- if (is.na(ab_lo) || is.na(ab_hi)) NA else (ab_lo * ab_hi > 0)
      
      rows[[k]] <- tibble(
        model_label = label,
        x = x,
        mediator = m_list[k],
        y = y,
        N = N,
        a = a_est, a_p = a_p,
        b = b_est, b_p = b_p,
        cprime = cprime_est, cprime_p = cprime_p,
        ab = ab_est, ab_ci_lower = ab_lo, ab_ci_upper = ab_hi, ab_sig = ab_sig,
        total = total_est, total_ci_lower = total_lo, total_ci_upper = total_hi,
        effect = "indirect"
      )
    }
  } else {
    # No mediator case (rare here), still return total/cprime
    rows[[1]] <- tibble(
      model_label = label, x = x, mediator = "", y = y, N = N,
      a = NA_real_, a_p = NA_real_, b = NA_real_, b_p = NA_real_,
      cprime = cprime_est, cprime_p = cprime_p,
      ab = NA_real_, ab_ci_lower = NA_real_, ab_ci_upper = NA_real_, ab_sig = NA,
      total = total_est, total_ci_lower = total_lo, total_ci_upper = total_hi,
      effect = "total_only"
    )
  }
  
  out <- dplyr::bind_rows(rows)
  
  # add diff_ab row for parallel (2 mediators)
  if (length(m_list) == 2) {
    diff_row <- get_lab(pe_boot, "diff_ab")
    out <- dplyr::bind_rows(
      out,
      tibble(
        model_label = label, x = x, mediator = paste(m_list, collapse = " - "),
        y = y, N = N,
        a = NA_real_, a_p = NA_real_, b = NA_real_, b_p = NA_real_,
        cprime = cprime_est, cprime_p = cprime_p,
        ab = if (nrow(diff_row)) diff_row$est[1] else NA_real_,
        ab_ci_lower = if (nrow(diff_row)) diff_row$`ci.lower`[1] else NA_real_,
        ab_ci_upper = if (nrow(diff_row)) diff_row$`ci.upper`[1] else NA_real_,
        ab_sig = if (nrow(diff_row)) (diff_row$`ci.lower`[1] * diff_row$`ci.upper`[1] > 0) else NA,
        total = total_est, total_ci_lower = total_lo, total_ci_upper = total_hi,
        effect = "diff_indirect"
      )
    )
  }
  
  out
}

# 4) Define models to run ---------------------------------------------------
models <- list(
  # Primary Centering mediator
  list(x = "atq_f_total", m = c("sacs_centering"), y = "dass_dep",   label = "ATQ_F -> Centering -> Dep"),
  list(x = "t_cfq",       m = c("sacs_centering"), y = "dass_dep",   label = "CFQ -> Centering -> Dep"),
  list(x = "atq_f_total", m = c("sacs_centering"), y = "t_wemwbs",   label = "ATQ_F -> Centering -> WEMWBS"),
  
  # Parallel facets Centering + Transcending
  list(x = "atq_f_total", m = c("sacs_centering","sacs_transcending"), y = "dass_dep",
       label = "ATQ_F -> (Centering, Transcending) -> Dep"),
  list(x = "atq_f_total", m = c("sacs_centering","sacs_transcending"), y = "t_wemwbs",
       label = "ATQ_F -> (Centering, Transcending) -> WEMWBS"),
  
  # Optional: Total SACS as mediator (include if present in your data)
  list(x = "atq_f_total", m = c("t_sacs"), y = "dass_dep",   label = "ATQ_F -> SACS_total -> Dep"),
  list(x = "t_cfq",       m = c("t_sacs"), y = "dass_dep",   label = "CFQ -> SACS_total -> Dep"),
  list(x = "atq_f_total", m = c("t_sacs"), y = "t_wemwbs",   label = "ATQ_F -> SACS_total -> WEMWBS")
)

# 5) Run models and export -------------------------------------------------
results <- purrr::map_dfr(models, function(mod) {
  run_med(d_all, x = mod$x, m_list = mod$m, y = mod$y, label = mod$label, cvars = character(), n_boot = 5000)
})

# order columns nicely
results <- results %>%
  dplyr::mutate(
    ab_sig = dplyr::case_when(isTRUE(ab_sig) ~ "yes", isFALSE(ab_sig) ~ "no", TRUE ~ NA_character_)
  ) %>%
  dplyr::select(model_label, x, mediator, y, N,
                a, a_p, b, b_p, cprime, cprime_p,
                ab, ab_ci_lower, ab_ci_upper, ab_sig,
                total, total_ci_lower, total_ci_upper,
                effect)

readr::write_csv(results, out_csv)
cat("✓ Saved mediation results to:\n  ", out_csv, "\n", sep = "")
