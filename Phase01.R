## ============================================================
## Phase 01 — Intake and Structure Check
##
## This script reads the raw dataset into R, verifies the
## presence of expected columns according to the lexicon, and
## produces basic structure and missingness reports. It saves
## intermediate artifacts to the designated phase folder.
##
## Conventions follow Andy Field's DSUR text as closely as
## practical: base R functions are used wherever possible,
## avoiding additional packages unless explicitly mentioned in
## DSUR. All outputs are written to
## '/Users/jamesjackson/Desktop/dir/phase01_intake'.
## ============================================================

## ---- 0. Define paths ----------------------------------------------------
root <- "/Users/jamesjackson/Desktop/dir"
phase01 <- file.path(root, "phase01_intake")
if (!dir.exists(phase01)) dir.create(phase01, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Read lexicon and define expected columns ------------------------
lex_path <- file.path(root, "lexicon.csv")
if (!file.exists(lex_path)) {
  stop("lexicon.csv not found in the project root. Please ensure it exists.")
}

lexicon <- utils::read.csv(lex_path, stringsAsFactors = FALSE)

## Extract the unique raw variable names from the lexicon
expected_cols <- unique(lexicon$raw_variable)

## ---- 2. Read raw data ---------------------------------------------------
raw_path <- file.path(root, "raw.csv")
if (!file.exists(raw_path)) {
  stop("raw.csv not found in the project root. Please ensure it exists.")
}

# read raw CSV using base functions
df_raw <- utils::read.csv(raw_path, stringsAsFactors = FALSE)

## ---- 3. Verify presence of expected columns ----------------------------
missing_columns <- setdiff(expected_cols, names(df_raw))

verification_report <- data.frame(
  expected_variable = expected_cols,
  present = expected_cols %in% names(df_raw),
  stringsAsFactors = FALSE
)

utils::write.csv(
  verification_report,
  file = file.path(phase01, "verification_report.csv"),
  row.names = FALSE
)

if (length(missing_columns) > 0) {
  warning(
    paste0(
      "The following expected columns are missing in raw.csv: ",
      paste(missing_columns, collapse = ", "),
      "\nCheck lexicon.csv and raw.csv for consistency."
    )
  )
}

## ---- 4. Structure snapshot ---------------------------------------------
structure_lines <- capture.output(str(df_raw))
writeLines(structure_lines, file.path(phase01, "structure_snapshot.txt"))

# Save the first 10 rows to a CSV for quick inspection
utils::write.csv(
  head(df_raw, 10),
  file = file.path(phase01, "first10rows.csv"),
  row.names = FALSE
)

## ---- 5. Missingness report ---------------------------------------------
missing_counts <- sapply(df_raw, function(x) sum(is.na(x)))
missing_report <- data.frame(
  variable = names(missing_counts),
  missing = as.integer(missing_counts),
  stringsAsFactors = FALSE
)

utils::write.csv(
  missing_report,
  file = file.path(phase01, "missingness_report.csv"),
  row.names = FALSE
)

## ---- 6. Rename map (identity) ------------------------------------------
rename_map <- data.frame(
  original_name = names(df_raw),
  new_name = names(df_raw),
  stringsAsFactors = FALSE
)
utils::write.csv(
  rename_map,
  file = file.path(phase01, "rename_map.csv"),
  row.names = FALSE
)

## ---- 7. Save snapshot of raw data --------------------------------------
utils::write.csv(
  df_raw,
  file = file.path(phase01, "raw_snapshot.csv"),
  row.names = FALSE
)

saveRDS(df_raw, file = file.path(phase01, "raw_snapshot.rds"))

## ---- 8. Console message -------------------------------------------------
cat(
  "\nPhase 01 intake and structure check complete.\n",
  "- Verification report: ", file.path(phase01, "verification_report.csv"), "\n",
  "- Structure snapshot: ", file.path(phase01, "structure_snapshot.txt"), "\n",
  "- Missingness report: ", file.path(phase01, "missingness_report.csv"), "\n",
  "- Rename map: ", file.path(phase01, "rename_map.csv"), "\n",
  "- Raw data snapshot saved in CSV and RDS forms.\n",
  sep = ""
)