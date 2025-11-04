# ============================================================
# Phase 02 — Clean & QC (DSUR-aligned; base R + utils only)
# Outputs -> /Users/jamesjackson/Desktop/dir/phase02_clean_qc
# Keeps ALL cases with complete SACS; others retained but flagged
# ============================================================

## ---- 0. Paths -----------------------------------------------------------
root <- "/Users/jamesjackson/Desktop/dir"
p01   <- file.path(root, "phase01_intake")
p02   <- file.path(root, "phase02_clean_qc")
if (!dir.exists(p02)) dir.create(p02, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load data -------------------------------------------------------
# Prefer the raw snapshot created in Phase 01, else fall back to raw.csv
raw_csv <- file.path(p01, "raw_snapshot.csv")
if (!file.exists(raw_csv)) raw_csv <- file.path(root, "raw.csv")
df <- utils::read.csv(raw_csv, check.names = TRUE, stringsAsFactors = FALSE)

# Small helper to write simple CSVs/TXT
wcsv <- function(x, path) utils::write.csv(x, path, row.names = FALSE)
wtxt <- function(x, path) writeLines(x, con = path)

## ---- 2. Identify key columns (robust to dots vs hyphens) ----------------
has_all <- function(nms, pool) all(nms %in% pool)

cols <- names(df)

# SACS items
sacs_items <- sprintf("SACS_%d", 1:10)
stopifnot(has_all(sacs_items, cols))

# CFQ items
cfq_items <- sprintf("CFQ_%d", 1:7)
stopifnot(has_all(cfq_items, cols))

# ATQ items (frequency f and believability b)
atq_f <- sprintf("ATQ_%df", 1:15)
atq_b <- sprintf("ATQ_%db", 1:15)
stopifnot(has_all(atq_f, cols), has_all(atq_b, cols))

# DASS items: names appear as DASS.21_1 .. DASS.21_21 in this dataset
# (handle either dotted or hyphenated just in case)
dass_dot <- paste0("DASS.21_", 1:21)
dass_hyp <- paste0("DASS-21_", 1:21)
if (has_all(dass_dot, cols)) {
  dass_items <- dass_dot
} else if (has_all(dass_hyp, cols)) {
  dass_items <- dass_hyp
} else {
  stop("Could not locate DASS-21 item columns.")
}

# WEMWBS items
wemwbs_items <- sprintf("WEMWBS_%d", 1:14)
stopifnot(has_all(wemwbs_items, cols))

# Demographics (Age; Gender is Q1 with optional Q1_4_TEXT for 'Other')
age_var <- "Age"
gender_code <- "Q1"
gender_other_text <- "Q1_4_TEXT"
stopifnot(age_var %in% cols, gender_code %in% cols)

# Timing columns (Qualtrics often uses dots)
dur_total   <- "Duration..in.seconds."
time_sacs   <- "SACStime_Page.Submit"
time_wemwbs <- "WEMWBStime_Page.Submit"
time_dass   <- "DASS.time_Page.Submit"
time_cfq    <- "CFQtime_Page.Submit"
time_atq    <- "ATQtime_Page.Submit"

timing_cols <- c(dur_total, time_sacs, time_wemwbs, time_dass, time_cfq, time_atq)
timing_cols <- timing_cols[timing_cols %in% cols]  # keep only those that exist

## ---- 3. Standardise missing & coerce types ------------------------------
# Treat "" as NA for text-esque fields
for (v in names(df)) {
  if (is.character(df[[v]])) df[[v]][df[[v]] == ""] <- NA
}

# Coerce item responses to numeric
to_numeric <- c(sacs_items, cfq_items, atq_f, atq_b, dass_items, wemwbs_items)
for (v in to_numeric) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))

# Age clean
df$Age_clean <- suppressWarnings(as.numeric(df[[age_var]]))
df$age_flag_out_of_range <- ifelse(is.na(df$Age_clean), FALSE,
                                   df$Age_clean < 16 | df$Age_clean > 100)

# Gender recode from Q1 (codes 1..5 per README; present data show 1..3)
gender_lab <- c(
  "1" = "Woman",
  "2" = "Man",
  "3" = "Non-Binary",
  "4" = "Other",
  "5" = "Prefer not to say"
)
df$Gender_code  <- suppressWarnings(as.integer(df[[gender_code]]))
df$Gender_label <- factor(ifelse(is.na(df$Gender_code), NA,
                                 gender_lab[as.character(df$Gender_code)]),
                          levels = gender_lab)

# Preserve free-text for 'Other' where available
if (gender_other_text %in% cols) {
  df$Gender_other_text <- df[[gender_other_text]]
}

## ---- 4. Completeness flags by scale -------------------------------------
is_complete_rowset <- function(d, vars) {
  # TRUE only if all required vars are non-missing in that row
  rowSums(is.na(d[vars])) == 0
}

df$sacs_complete   <- is_complete_rowset(df, sacs_items)
df$cfq_complete    <- is_complete_rowset(df, cfq_items)
df$atq_f_complete  <- is_complete_rowset(df, atq_f)
df$atq_b_complete  <- is_complete_rowset(df, atq_b)
df$wemwbs_complete <- is_complete_rowset(df, wemwbs_items)
df$dass_complete   <- is_complete_rowset(df, dass_items)

# Retention rule: keep ALL rows with complete SACS
df$retain_primary_pool <- df$sacs_complete

# ≤10% missing flags (contingency for possible MI later; not used now)
lt10 <- function(vars) rowMeans(is.na(df[vars])) <= 0.10
df$sacs_le10miss   <- lt10(sacs_items)
df$cfq_le10miss    <- lt10(cfq_items)
df$atq_f_le10miss  <- lt10(atq_f)
df$atq_b_le10miss  <- lt10(atq_b)
df$wemwbs_le10miss <- lt10(wemwbs_items)
df$dass_le10miss   <- lt10(dass_items)

## ---- 5. Timing QC (flag only; do not exclude) ---------------------------
# DSUR-friendly approach: log-transform skewed times; IQR fences for flags.
flag_timing <- function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  lg   <- log(pmax(xnum, 1)) # avoid log(0)
  q1   <- stats::quantile(lg, 0.25, na.rm = TRUE)
  q3   <- stats::quantile(lg, 0.75, na.rm = TRUE)
  iqr  <- q3 - q1
  low  <- q1 - 1.5 * iqr
  high <- q3 + 1.5 * iqr
  list(
    low_threshold_log  = as.numeric(low),
    high_threshold_log = as.numeric(high),
    fast  = lg < low,
    slow  = lg > high
  )
}

timing_summary <- list()
for (v in timing_cols) {
  out <- flag_timing(df[[v]])
  df[[paste0(v, "_fast_flag")]] <- isTRUE(out$fast)  # vectorized below
  df[[paste0(v, "_slow_flag")]] <- isTRUE(out$slow)
  # (Replace with vectorized assignment:)
  df[[paste0(v, "_fast_flag")]] <- out$fast
  df[[paste0(v, "_slow_flag")]] <- out$slow
  timing_summary[[v]] <- c(low_log = out$low_threshold_log,
                           high_log = out$high_threshold_log)
}

# Convenience overall flags
fast_cols <- grep("_fast_flag$", names(df), value = TRUE)
slow_cols <- grep("_slow_flag$", names(df), value = TRUE)
df$any_fast_time_flag <- rowSums(df[fast_cols], na.rm = TRUE) > 0
df$any_slow_time_flag <- rowSums(df[slow_cols], na.rm = TRUE) > 0

## ---- 6. Straight-lining / low-variance flags ----------------------------
wip_sd_flag <- function(vars, tol = 0.25) {
  # within-person SD across listed items; flag if very low variance
  sd_vals <- apply(df[vars], 1, function(r) stats::sd(r, na.rm = TRUE))
  is.na(sd_vals) <- FALSE
  sd_vals <= tol
}
df$sacs_lowvar_flag   <- wip_sd_flag(sacs_items, tol = 0.25)
df$cfq_lowvar_flag    <- wip_sd_flag(cfq_items,  tol = 0.25)
df$atq_f_lowvar_flag  <- wip_sd_flag(atq_f,      tol = 0.25)
df$atq_b_lowvar_flag  <- wip_sd_flag(atq_b,      tol = 0.25)
df$wemwbs_lowvar_flag <- wip_sd_flag(wemwbs_items, tol = 0.25)
df$dass_lowvar_flag   <- wip_sd_flag(dass_items, tol = 0.25)

## ---- 7. Out-of-range checks (basic) -------------------------------------
# Most items are Likert 1..7 (SACS/CFQ/ATQ/WEMWBS) and 0..3 for DASS.
oor_flag <- function(vars, minv, maxv) {
  apply(df[vars], 2, function(x) ifelse(is.na(x), FALSE, x < minv | x > maxv))
}
oor <- list(
  SACS   = oor_flag(sacs_items,   1, 7),
  CFQ    = oor_flag(cfq_items,    1, 7),
  ATQ_F  = oor_flag(atq_f,        1, 7),
  ATQ_B  = oor_flag(atq_b,        1, 7),
  WEMWBS = oor_flag(wemwbs_items, 1, 5),  # WEMWBS is 1..5
  DASS   = oor_flag(dass_items,   0, 3)   # DASS-21 items 0..3
)

# Count any out-of-range by row to help QC
df$any_oor_items <- Reduce(`|`, lapply(oor, function(m) rowSums(m) > 0))

## ---- 8. Export main cleaned base + QC flags -----------------------------
# Keep core cleaned fields + flags; retain all original columns as well
clean_cols <- unique(c(
  names(df), "Age_clean", "age_flag_out_of_range",
  "Gender_code", "Gender_label", "Gender_other_text",
  "sacs_complete","cfq_complete","atq_f_complete","atq_b_complete",
  "wemwbs_complete","dass_complete",
  "sacs_le10miss","cfq_le10miss","atq_f_le10miss","atq_b_le10miss",
  "wemwbs_le10miss","dass_le10miss",
  "retain_primary_pool","any_fast_time_flag","any_slow_time_flag",
  "sacs_lowvar_flag","cfq_lowvar_flag","atq_f_lowvar_flag","atq_b_lowvar_flag",
  "wemwbs_lowvar_flag","dass_lowvar_flag","any_oor_items"
))

clean_base <- df[clean_cols]
wcsv(clean_base, file.path(p02, "clean_base.csv"))

# A slim QC-only table (one row per participant)
qc_keep <- c("Age_clean","Gender_code","Gender_label","retain_primary_pool",
             "sacs_complete","cfq_complete","atq_f_complete","atq_b_complete",
             "wemwbs_complete","dass_complete",
             "any_fast_time_flag","any_slow_time_flag",
             "sacs_lowvar_flag","cfq_lowvar_flag","atq_f_lowvar_flag",
             "atq_b_lowvar_flag","wemwbs_lowvar_flag","dass_lowvar_flag",
             "any_oor_items")
qc_keep <- qc_keep[qc_keep %in% names(df)]
qc_flags <- df[, qc_keep, drop = FALSE]
wcsv(qc_flags, file.path(p02, "qc_flags.csv"))

## ---- 9. Decisions log (human-readable) ----------------------------------
decisions <- c(
  "Phase 02 — Cleaning & QC decisions",
  paste0("Data source: ", basename(raw_csv)),
  "",
  "Retention rule:",
  "- Retain all cases with complete SACS (retain_primary_pool = TRUE). Others retained for sensitivity.",
  "",
  "Timing flags:",
  paste0("- Columns considered: ", paste(timing_cols, collapse = ", ")),
  "- Log-transform times; IQR fences define 'fast' and 'slow' flags; no exclusions applied.",
  "",
  "Straight-lining flags:",
  "- Within-person SD threshold (<= 0.25) for each scale; flagged only (no deletions).",
  "",
  "Out-of-range checks:",
  "- SACS/CFQ/ATQ/WEMWBS assumed 1..7 (WEMWBS is 1..5), DASS items 0..3; flagged only.",
  "",
  "Completeness rules:",
  "- Scores will be computed only when required items are all present.",
  "- ≤10% missing flags saved for potential MI later (not used in primary analyses)."
)

# Add per-timing variable thresholds (log-scale) if available
if (length(timing_summary)) {
  decisions <- c(decisions, "", "Timing thresholds (log-scale):")
  for (nm in names(timing_summary)) {
    th <- timing_summary[[nm]]
    decisions <- c(decisions,
                   sprintf("  * %s: low=%.3f, high=%.3f", nm, th["low_log"], th["high_log"]))
  }
}
wtxt(decisions, file.path(p02, "decisions.txt"))

## ---- 10. Console summary -------------------------------------------------
cat("\nPhase 02 complete.\n",
    "- Clean base: ", file.path(p02, 'clean_base.csv'), "\n",
    "- QC flags  : ", file.path(p02, 'qc_flags.csv'), "\n",
    "- Decisions : ", file.path(p02, 'decisions.txt'), "\n", sep = "")
