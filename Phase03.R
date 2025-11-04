# ===============================================================
# Phase 03 — Scoring (DSUR-aligned)
# 
# Reads phase02_clean_qc/clean_base.csv, computes scale scores,
# and writes scored data + a scale-key mapping.
# Outputs: scored_primary.csv, scale_key.csv
# ===============================================================

## ---- 0. Paths -----------------------------------------------------------
root  <- "/Users/jamesjackson/Desktop/dir"
p02   <- file.path(root, "phase02_clean_qc")
p03   <- file.path(root, "phase03_scoring")
if (!dir.exists(p03)) dir.create(p03, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load cleaned data -----------------------------------------------
clean_path <- file.path(p02, "clean_base.csv")
df <- utils::read.csv(clean_path, check.names = TRUE, stringsAsFactors = FALSE)

## ---- 2. Helper to compute scores safely ---------------------------------
compute_sum <- function(d, vars, complete_flag) {
  # Returns sum across vars if complete_flag is TRUE, otherwise NA
  vals <- rowSums(d[vars], na.rm = FALSE)
  if (!missing(complete_flag)) vals[!complete_flag] <- NA
  return(vals)
}

## ---- 3. Define item sets -------------------------------------------------
sacs_total_items     <- sprintf("SACS_%d", 1:10)
sacs_centering_items <- c("SACS_1", "SACS_2", "SACS_5", "SACS_6")
sacs_transcend_items <- c("SACS_3", "SACS_4", "SACS_7", "SACS_8", "SACS_9", "SACS_10")

cfq_items <- sprintf("CFQ_%d", 1:7)

atq_f_items <- sprintf("ATQ_%df", 1:15)
atq_b_items <- sprintf("ATQ_%db", 1:15)

dass_items <- paste0("DASS.21_", 1:21)
# DASS subscales: Depression (3,5,10,13,16,17,21), Anxiety (2,4,7,9,15,19,20), Stress (1,6,8,11,12,14,18)
dass_depr <- c(3, 5, 10, 13, 16, 17, 21)
dass_anxi <- c(2, 4, 7, 9, 15, 19, 20)
dass_stre <- c(1, 6, 8, 11, 12, 14, 18)
dass_depr_items <- paste0("DASS.21_", dass_depr)
dass_anxi_items <- paste0("DASS.21_", dass_anxi)
dass_stre_items <- paste0("DASS.21_", dass_stre)

wemwbs_items <- sprintf("WEMWBS_%d", 1:14)

## ---- 4. Compute SACS scores ---------------------------------------------
df$SACS_Total       <- compute_sum(df, sacs_total_items, df$sacs_complete)
df$SACS_Centering   <- compute_sum(df, sacs_centering_items, df$sacs_complete)
df$SACS_Transcending<- compute_sum(df, sacs_transcend_items, df$sacs_complete)

## ---- 5. Compute CFQ score -----------------------------------------------
df$CFQ_Total <- compute_sum(df, cfq_items, df$cfq_complete)

## ---- 6. Compute ATQ scores ----------------------------------------------
df$ATQ_Frequency   <- compute_sum(df, atq_f_items, df$atq_f_complete)
df$ATQ_Believability <- compute_sum(df, atq_b_items, df$atq_b_complete)
# Combined = sum of all 30 items (frequency + believability)
# Only valid when both halves are complete
df$ATQ_Combined    <- compute_sum(df, c(atq_f_items, atq_b_items),
                                  df$atq_f_complete & df$atq_b_complete)

## ---- 7. Compute WEMWBS score --------------------------------------------
df$WEMWBS_Total <- compute_sum(df, wemwbs_items, df$wemwbs_complete)

## ---- 8. Compute DASS scores (×2) ----------------------------------------
df$DASS_Depression_x2 <- 2 * compute_sum(df, dass_depr_items, df$dass_complete)
df$DASS_Anxiety_x2    <- 2 * compute_sum(df, dass_anxi_items, df$dass_complete)
df$DASS_Stress_x2     <- 2 * compute_sum(df, dass_stre_items, df$dass_complete)
df$DASS_Total_x2      <- 2 * compute_sum(df, dass_items, df$dass_complete)

## ---- 9. (Optional) z-scores for comparability ---------------------------
score_vars <- c("SACS_Total","SACS_Centering","SACS_Transcending",
                "CFQ_Total","ATQ_Frequency","ATQ_Believability","ATQ_Combined",
                "WEMWBS_Total",
                "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2")
for (sv in score_vars) {
  sv_z <- paste0(sv, "_z")
  vals <- df[[sv]]
  df[[sv_z]] <- (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
}

## ---- 10. Write scored dataset ------------------------------------------
scored_file <- file.path(p03, "scored_primary.csv")
utils::write.csv(df, scored_file, row.names = FALSE)

## ---- 11. Build scale key ------------------------------------------------
scale_key <- data.frame(
  variable_name = c("SACS_Total","SACS_Centering","SACS_Transcending",
                    "CFQ_Total","ATQ_Frequency","ATQ_Believability","ATQ_Combined",
                    "WEMWBS_Total",
                    "DASS_Depression_x2","DASS_Anxiety_x2","DASS_Stress_x2","DASS_Total_x2"),
  scale    = c("SACS","SACS","SACS",
               "CFQ-7","ATQ-15","ATQ-15","ATQ-15",
               "WEMWBS-14",
               "DASS-21","DASS-21","DASS-21","DASS-21"),
  subscale = c("Total","Centering","Transcending",
               "Total","Frequency","Believability","Combined",
               "Total",
               "Depression","Anxiety","Stress","Total"),
  items = c(paste(sacs_total_items, collapse = ", "),
            paste(sacs_centering_items, collapse = ", "),
            paste(sacs_transcend_items, collapse = ", "),
            paste(cfq_items, collapse = ", "),
            paste(atq_f_items, collapse = ", "),
            paste(atq_b_items, collapse = ", "),
            paste(c(atq_f_items, atq_b_items), collapse = ", "),
            paste(wemwbs_items, collapse = ", "),
            paste(dass_depr_items, collapse = ", "),
            paste(dass_anxi_items, collapse = ", "),
            paste(dass_stre_items, collapse = ", "),
            paste(dass_items, collapse = ", ")),
  notes = c(
    "Sum of 10 SACS items; valid only when sacs_complete",
    "Sum of SACS_1,2,5,6; valid only when sacs_complete",
    "Sum of SACS_3,4,7,8,9,10; valid only when sacs_complete",
    "Sum of 7 CFQ items; valid only when cfq_complete",
    "Sum of 15 ATQ frequency items; valid only when atq_f_complete",
    "Sum of 15 ATQ believability items; valid only when atq_b_complete",
    "Sum of all 30 ATQ items; valid only when both atq_f_complete & atq_b_complete",
    "Sum of 14 WEMWBS items; valid only when wemwbs_complete",
    "Sum of 7 Depression items ×2; valid only when dass_complete",
    "Sum of 7 Anxiety items ×2; valid only when dass_complete",
    "Sum of 7 Stress items ×2; valid only when dass_complete",
    "Sum of 21 DASS items ×2; valid only when dass_complete"
  ),
  stringsAsFactors = FALSE
)
scale_key_file <- file.path(p03, "scale_key.csv")
utils::write.csv(scale_key, scale_key_file, row.names = FALSE)

## ---- 12. Console summary -----------------------------------------------
cat("\nPhase 03 complete.\n",
    "- Scored data saved as: ", scored_file, "\n",
    "- Scale key saved as: ", scale_key_file, "\n", sep = "")
