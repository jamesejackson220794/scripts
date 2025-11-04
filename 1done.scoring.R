# Score Psychometric Scales - 1done Dataset
# ================================================================================
# Applies validated scoring algorithms to create scale totals and subscales
# Based on scoring documentation and maintains naming conventions
# ================================================================================

# 1. SETUP -----------------------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(here)

# 2. LOAD 1DONE DATASET ---------------------------------------------------------
cat("Loading 1done dataset for scoring...\n")

# Check if 1done dataset exists, if not load the main extracted dataset
if (file.exists("df_1done.csv")) {
  df_1done <- read_csv("df_1done.csv", col_types = cols(.default = "c"))
  cat("Loaded df_1done.csv with", nrow(df_1done), "participants\n")
} else if (file.exists("df_extracted.csv")) {
  df_1done <- read_csv("df_extracted.csv", col_types = cols(.default = "c"))
  cat("Loaded df_extracted.csv with", nrow(df_1done), "participants\n")
  cat("Note: Consider running the 1done filter script first\n")
} else {
  stop("Error: No input dataset found. Please ensure df_1done.csv or df_extracted.csv exists.")
}

# 3. DATA PREPARATION ------------------------------------------------------------
cat("Preparing data for scoring...\n")

# Convert scale items to numeric (preserve non-scale columns as character)
scale_item_patterns <- c(
  "^WEMWBS_[0-9]+$", "^DASS21_[0-9]+$", "^CSES_[0-9]+$", 
  "^SACS_[0-9]+$", "^ABS_[0-9]+$", "^SWLS_[0-9]+$",
  "^CFQ_[0-9]+$", "^BRS_[0-9]+$", "^ATQ_freq_[0-9]+$", "^ATQ_belief_[0-9]+$"
)

# Identify scale item columns
scale_columns <- names(df_1done)[grepl(paste(scale_item_patterns, collapse = "|"), names(df_1done))]

# Convert scale items to numeric
df_scored <- df_1done %>%
  mutate(across(all_of(scale_columns), ~ suppressWarnings(as.numeric(.))))

# 4. SCORING FUNCTIONS -----------------------------------------------------------

# Helper function: Score scale with missing data handling
score_scale <- function(data, item_cols, reverse_items = NULL, reverse_max = NULL) {
  scale_data <- data %>% select(all_of(item_cols))
  
  # Apply reverse scoring if specified
  if (!is.null(reverse_items) && !is.null(reverse_max)) {
    reverse_cols <- paste0(names(scale_data)[1] %>% str_extract("^[A-Z]+"), "_", reverse_items)
    reverse_cols <- intersect(reverse_cols, names(scale_data))
    scale_data <- scale_data %>%
      mutate(across(all_of(reverse_cols), ~ reverse_max + 1 - .))
  }
  
  # Calculate sum with missing data handling
  scale_data %>%
    pmap_dbl(~ {
      values <- c(...)
      if (any(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
    })
}

# 5. SCALE SCORING ---------------------------------------------------------------
cat("Computing scale scores...\n")

# 5.1 WEMWBS (14 items, range 14-70)
cat("  • WEMWBS...\n")
wemwbs_items <- paste0("WEMWBS_", 1:14)
df_scored$t.wemwbs <- score_scale(df_scored, wemwbs_items)

# 5.2 DASS-21 (21 items total, 3 subscales of 7 items each)
cat("  • DASS-21...\n")
dass_all_items <- paste0("DASS21_", 1:21)
dass_dep_items <- paste0("DASS21_", c(3,5,10,13,16,17,21))    # Depression
dass_anx_items <- paste0("DASS21_", c(2,4,7,9,15,19,20))     # Anxiety  
dass_str_items <- paste0("DASS21_", c(1,6,8,11,12,14,18))    # Stress

df_scored$t.dass <- score_scale(df_scored, dass_all_items)
df_scored$t.dep <- score_scale(df_scored, dass_dep_items) 
df_scored$t.anx <- score_scale(df_scored, dass_anx_items)
df_scored$t.stress <- score_scale(df_scored, dass_str_items)

# 5.3 CSES (12 items, reverse items 2,4,6,8,10,12 on 1-5 scale)
cat("  • CSES...\n")
cses_items <- paste0("CSES_", 1:12)
cses_reverse <- c(2,4,6,8,10,12)
df_scored$t.cses <- score_scale(df_scored, cses_items, cses_reverse, 5)

# 5.4 SACS (10 items, range 10-70)
cat("  • SACS...\n")
sacs_items <- paste0("SACS_", 1:10)
df_scored$t.sacs <- score_scale(df_scored, sacs_items)

# 5.5 SWLS (5 items, range 5-35)
cat("  • SWLS...\n")
swls_items <- paste0("SWLS_", 1:5)
df_scored$t.swls <- score_scale(df_scored, swls_items)

# 5.6 CFQ-7 (7 items, range 7-49)
cat("  • CFQ-7...\n")
cfq_items <- paste0("CFQ_", 1:7)
df_scored$t.cfq <- score_scale(df_scored, cfq_items)

# 5.7 BRS (6 items, range 6-30)
cat("  • BRS...\n")
brs_items <- paste0("BRS_", 1:6)
df_scored$t.brs <- score_scale(df_scored, brs_items)

# 5.8 ABS-SF (24 items: 12 irrational + 12 rational)
cat("  • ABS-SF...\n")
# Irrational items: 1,2,4,7,9,11,13,14,16,18,20,22
# Rational items: 3,5,6,8,10,12,15,17,19,21,23,24
abs_irr_items <- paste0("ABS_", c(1,2,4,7,9,11,13,14,16,18,20,22))
abs_rat_items <- paste0("ABS_", c(3,5,6,8,10,12,15,17,19,21,23,24))

df_scored$t.abs.irr <- score_scale(df_scored, abs_irr_items)
df_scored$t.abs.rat <- score_scale(df_scored, abs_rat_items)

# 5.9 ATQ-15 (15 items x 2 dimensions: frequency + believability)
cat("  • ATQ-15...\n")
atq_freq_items <- paste0("ATQ_freq_", 1:15)
atq_belief_items <- paste0("ATQ_belief_", 1:15)

df_scored$f.atq <- score_scale(df_scored, atq_freq_items)
df_scored$b.atq <- score_scale(df_scored, atq_belief_items)

# ATQ total (frequency + believability)
df_scored$t.atq <- ifelse(
  is.na(df_scored$f.atq) | is.na(df_scored$b.atq),
  NA_real_,
  df_scored$f.atq + df_scored$b.atq
)

# 6. QUALITY CHECKS --------------------------------------------------------------
cat("\nPerforming quality checks...\n")

# Check score distributions and ranges
score_vars <- c("t.wemwbs", "t.dass", "t.dep", "t.anx", "t.stress", 
                "t.cses", "t.sacs", "t.swls", "t.cfq", "t.brs",
                "t.abs.irr", "t.abs.rat", "f.atq", "b.atq", "t.atq")

score_summary <- df_scored %>%
  select(all_of(score_vars)) %>%
  summarise(across(everything(), list(
    n = ~ sum(!is.na(.)),
    min = ~ min(., na.rm = TRUE),
    max = ~ max(., na.rm = TRUE), 
    mean = ~ round(mean(., na.rm = TRUE), 2),
    sd = ~ round(sd(., na.rm = TRUE), 2)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
  separate(stat, into = c("scale", "statistic"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = statistic, values_from = value)

cat("Score summary statistics:\n")
print(score_summary, n = Inf)

# Check for any impossible values
cat("\nChecking for out-of-range values:\n")

expected_ranges <- list(
  t.wemwbs = c(14, 70),
  t.dass = c(0, 63),
  t.dep = c(0, 21), 
  t.anx = c(0, 21),
  t.stress = c(0, 21),
  t.cses = c(12, 60),
  t.sacs = c(10, 70),
  t.swls = c(5, 35),
  t.cfq = c(7, 49),
  t.brs = c(6, 30),
  t.abs.irr = c(0, 48),
  t.abs.rat = c(0, 48),
  f.atq = c(15, 75),
  b.atq = c(15, 75),
  t.atq = c(30, 150)
)

for (scale in names(expected_ranges)) {
  if (scale %in% names(df_scored)) {
    values <- df_scored[[scale]][!is.na(df_scored[[scale]])]
    range_expected <- expected_ranges[[scale]]
    out_of_range <- sum(values < range_expected[1] | values > range_expected[2])
    if (out_of_range > 0) {
      cat("WARNING:", scale, "has", out_of_range, "values outside expected range", 
          paste(range_expected, collapse = "-"), "\n")
    }
  }
}

# 7. COMPLETION ANALYSIS ---------------------------------------------------------
cat("\nAnalyzing scale completion rates:\n")

completion_rates <- df_scored %>%
  select(all_of(score_vars)) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "scale", values_to = "completed") %>%
  mutate(
    total = nrow(df_scored),
    completion_rate = round((completed / total) * 100, 1)
  ) %>%
  arrange(desc(completion_rate))

print(completion_rates)

# 8. CREATE NEW SEQUENTIAL PARTICIPANT IDs -------------------------------------
cat("\nCreating new sequential participant IDs...\n")

# Create new sequential IDs for the scored dataset
df_scored <- df_scored %>%
  arrange(participant_id_num) %>%  # Sort by existing ID for consistency
  mutate(
    new_participant_id = sprintf("%02d", row_number()),  # 01, 02, 03, etc.
    .before = 1  # Place as first column
  )

cat("Assigned new participant IDs:", min(df_scored$new_participant_id), "to", max(df_scored$new_participant_id), "\n")

# Create mapping table for reference
id_mapping <- df_scored %>%
  select(new_participant_id, participant_id_num, external_data_ref) %>%
  arrange(new_participant_id)

cat("ID mapping created for", nrow(id_mapping), "participants\n")

# 9. SAVE SCORED DATASET --------------------------------------------------------
cat("\nSaving scored dataset...\n")

# Create output directory
processed_dir <- here("data", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Save scored dataset
saveRDS(df_scored, file.path(processed_dir, "df_1done_scored.rds"))
write_csv(df_scored, file.path(processed_dir, "df_1done_scored.csv"))
write_csv(df_scored, "df_1done_scored.csv")  # Also in project root

# Create scores-only dataset for easy analysis
df_scores_only <- df_scored %>%
  select(
    # New sequential ID
    new_participant_id,
    # Original identifiers (for mapping back if needed)
    participant_id_num, external_data_ref,
    # Demographics  
    Age, CountryBorn, CountryRes, Relationship, Ethnicity, Education,
    # Journal completion
    total_journals,
    # Scale scores
    all_of(score_vars)
  )

write_csv(df_scores_only, file.path(processed_dir, "scale_scores_summary.csv"))
write_csv(df_scores_only, "scale_scores_summary.csv")

# Save ID mapping table
write_csv(id_mapping, file.path(processed_dir, "participant_id_mapping.csv"))
write_csv(id_mapping, "participant_id_mapping.csv")

# 10. CREATE SCORING REPORT ------------------------------------------------------
cat("\nGenerating scoring report...\n")

scoring_report <- list(
  dataset_info = list(
    input_file = ifelse(file.exists("df_1done.csv"), "df_1done.csv", "df_extracted.csv"),
    participants = nrow(df_scored),
    scales_scored = length(score_vars),
    id_range = paste(min(df_scored$new_participant_id), "to", max(df_scored$new_participant_id))
  ),
  scale_descriptions = list(
    "t.wemwbs" = "Warwick-Edinburgh Mental Well-being Scale (14 items, range 14-70)",
    "t.dass" = "DASS-21 Total (21 items, range 0-63)",
    "t.dep" = "DASS-21 Depression subscale (7 items, range 0-21)",
    "t.anx" = "DASS-21 Anxiety subscale (7 items, range 0-21)", 
    "t.stress" = "DASS-21 Stress subscale (7 items, range 0-21)",
    "t.cses" = "Core Self-Evaluations Scale (12 items, 6 reversed, range 12-60)",
    "t.sacs" = "Self-as-Context Scale (10 items, range 10-70)",
    "t.swls" = "Satisfaction with Life Scale (5 items, range 5-35)",
    "t.cfq" = "Cognitive Fusion Questionnaire (7 items, range 7-49)",
    "t.brs" = "Brief Resilience Scale (6 items, range 6-30)",
    "t.abs.irr" = "ABS-SF Irrational Beliefs (12 items, range 0-48)",
    "t.abs.rat" = "ABS-SF Rational Beliefs (12 items, range 0-48)",
    "f.atq" = "ATQ-15 Frequency (15 items, range 15-75)",
    "b.atq" = "ATQ-15 Believability (15 items, range 15-75)", 
    "t.atq" = "ATQ-15 Total (frequency + believability, range 30-150)"
  ),
  completion_rates = completion_rates
)

# Save report as JSON for reference
jsonlite::write_json(scoring_report, file.path(processed_dir, "scoring_report.json"), 
                     pretty = TRUE, auto_unbox = TRUE)

# 11. FINAL SUMMARY -------------------------------------------------------------
cat("\n✓ Scoring completed successfully!\n")
cat("\nFiles created:\n")
cat("• df_1done_scored.csv (full dataset with scores)\n")
cat("• scale_scores_summary.csv (demographics + scores only)\n")
cat("• participant_id_mapping.csv (ID mapping table)\n")
cat("• scoring_report.json (detailed scoring information)\n")

cat("\nDataset summary:\n")
cat("• Participants with scores:", nrow(df_scored), "\n")
cat("• New participant IDs:", min(df_scored$new_participant_id), "to", max(df_scored$new_participant_id), "\n")
cat("• Scales scored:", length(score_vars), "\n")
cat("• Variables in full dataset:", ncol(df_scored), "\n")
cat("• Average completion rate:", round(mean(completion_rates$completion_rate), 1), "%\n")

cat("\nScale score variables added:\n")
for (var in score_vars) {
  cat("• ", var, "\n")
}

cat("\nNew participant ID system:\n")
cat("• Sequential IDs from 01 to", sprintf("%02d", nrow(df_scored)), "\n")
cat("• Original IDs preserved for reference\n")
cat("• Mapping table created for de-identification\n")

cat("\nReady for statistical analysis! 📊\n")