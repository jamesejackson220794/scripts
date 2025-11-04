# ===============================================================
# Phase 04 — Exploratory Data Analysis (DSUR-aligned)
# Inputs : /Users/jamesjackson/Desktop/dir/phase03_scoring/scored_primary.csv
# Outputs: /Users/jamesjackson/Desktop/dir/phase04_eda/
#   • descriptives.csv           : mean, SD, skew, kurtosis, n, missing
#   • qc_overview.csv            : counts of complete cases and timing flags
#   • correlations.csv / cor_long.csv : r matrix and tidy pairwise r & p
#   • histogram/density plots    : hist_<var>.png
#   • boxplots by gender (if Gender_label exists)
#   • age trends (if Age_clean exists)
#   • cor_tiles.png              : heatmap of pairwise r
#   • console summary printed at end
# ===============================================================

## ---- 0. Paths -----------------------------------------------------------
root <- "/Users/jamesjackson/Desktop/dir"
p03  <- file.path(root, "phase03_scoring")
p04  <- file.path(root, "phase04_eda")
if (!dir.exists(p04)) dir.create(p04, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load data -------------------------------------------------------
scored_path <- file.path(p03, "scored_primary.csv")
df <- utils::read.csv(scored_path, check.names = TRUE, stringsAsFactors = FALSE)

## Convenience writer
wcsv <- function(x, path) utils::write.csv(x, path, row.names = FALSE)

## ---- 2. Identify variables ---------------------------------------------
# Primary scale scores
score_vars <- c(
  "SACS_Total", "SACS_Centering", "SACS_Transcending",
  "CFQ_Total",
  "ATQ_Frequency", "ATQ_Believability", "ATQ_Combined",
  "WEMWBS_Total",
  "DASS_Depression_x2", "DASS_Anxiety_x2", "DASS_Stress_x2", "DASS_Total_x2"
)
score_vars <- intersect(score_vars, names(df))

# Demographics if available
age_var    <- if ("Age_clean"    %in% names(df)) "Age_clean" else NA
gender_var <- if ("Gender_label" %in% names(df)) "Gender_label" else NA

## ---- 3. Descriptive statistics -----------------------------------------
desc_one <- function(x) {
  # Remove non-finite values for descriptive stats
  ok <- is.finite(x)
  xx <- x[ok]
  n <- length(xx)
  out <- list(n = n, missing = sum(!ok))
  if (n == 0) {
    out$mean <- NA_real_; out$sd <- NA_real_; out$skew <- NA_real_; out$kurtosis_excess <- NA_real_
  } else {
    m <- mean(xx); s <- stats::sd(xx)
    out$mean <- m; out$sd <- s
    # Skewness: mean((x - mean)^3) / sd^3
    out$skew <- if (n > 2 && s > 0) mean(((xx - m) / s)^3) else NA_real_
    # Excess kurtosis: mean((x - mean)^4) / sd^4 - 3
    out$kurtosis_excess <- if (n > 3 && s > 0) mean(((xx - m) / s)^4) - 3 else NA_real_
  }
  return(out)
}

desclist <- lapply(score_vars, function(v) {
  cbind(variable = v, as.data.frame(desc_one(df[[v]]), stringsAsFactors = FALSE))
})
descriptives <- do.call(rbind, desclist)
wcsv(descriptives, file.path(p04, "descriptives.csv"))

## ---- 4. Minimal ggplot theme -------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

theme_dsur <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

## ---- 5. Plot distributions ---------------------------------------------
plot_hist <- function(df, v, outdir) {
  g <- ggplot(df, aes(x = .data[[v]])) +
    geom_histogram(bins = 30, color = "black", fill = "#CCCCCC") +
    geom_density(aes(y = ..density..), color = "darkblue") +
    labs(title = paste0(v, " distribution"),
         x = v, y = "Count / Density") +
    theme_dsur()
  ggsave(filename = file.path(outdir, paste0("hist_", v, ".png")),
         plot = g, width = 6, height = 4, dpi = 150)
}
invisible(lapply(score_vars, function(v) plot_hist(df, v, p04)))

## ---- 6. Boxplots by gender (if present) --------------------------------
if (!is.na(gender_var)) {
  df[[gender_var]] <- factor(df[[gender_var]])
  for (v in score_vars) {
    g <- ggplot(df, aes(x = .data[[gender_var]], y = .data[[v]])) +
      geom_boxplot(outlier.alpha = 0.5) +
      labs(title = paste0(v, " by ", gender_var),
           x = gender_var, y = v) +
      theme_dsur()
    ggsave(filename = file.path(p04, paste0("box_", v, "_by_", gender_var, ".png")),
           plot = g, width = 6, height = 4, dpi = 150)
  }
}

## ---- 7. Age trends (if present) ----------------------------------------
if (!is.na(age_var)) {
  for (v in score_vars) {
    g <- ggplot(df, aes(x = .data[[age_var]], y = .data[[v]])) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(title = paste0(v, " vs ", age_var),
           x = age_var, y = v) +
      theme_dsur()
    ggsave(filename = file.path(p04, paste0("age_trend_", v, ".png")),
           plot = g, width = 6, height = 4, dpi = 150)
  }
}

## ---- 8. Correlations ----------------------------------------------------
# Pearson correlation matrix (pairwise complete)
numdat <- df[score_vars]
cor_mat <- suppressWarnings(stats::cor(numdat, use = "pairwise.complete.obs", method = "pearson"))
# Save wide matrix
wcsv(data.frame(var = rownames(cor_mat), cor_mat, check.names = FALSE),
     file.path(p04, "correlations.csv"))

# Tidy long form with r and p
pairs <- t(combn(score_vars, 2))
get_rp <- function(a, b) {
  x <- df[[a]]; y <- df[[b]]
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3) return(data.frame(var1 = a, var2 = b, n = n, r = NA, p = NA))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
  data.frame(var1 = a, var2 = b, n = n, r = unname(ct$estimate), p = ct$p.value)
}
cor_long <- do.call(rbind, apply(pairs, 1, function(ix) get_rp(ix[1], ix[2])))
wcsv(cor_long, file.path(p04, "cor_long.csv"))

## ---- 9. Correlation heatmap --------------------------------------------
if (length(score_vars) >= 3) {
  # Create long table for plotting (upper triangle only)
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("Var1", "Var2", "r")
  # Keep only Var1 < Var2 to avoid duplicates
  cor_df <- subset(cor_df, match(Var1, rownames(cor_mat)) < match(Var2, colnames(cor_mat)))
  g <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
    labs(title = "Pearson correlations (pairwise-complete)", fill = "r") +
    theme_dsur() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(filename = file.path(p04, "cor_tiles.png"),
         plot = g, width = 7.5, height = 6, dpi = 150)
}

## ---- 10. QC overview summary -------------------------------------------
qc_summary <- data.frame(
  n_rows            = nrow(df),
  n_complete_SACS   = sum(df$sacs_complete,   na.rm = TRUE),
  n_complete_CFQ    = sum(df$cfq_complete,    na.rm = TRUE),
  n_complete_ATQF   = sum(df$atq_f_complete,  na.rm = TRUE),
  n_complete_ATQB   = sum(df$atq_b_complete,  na.rm = TRUE),
  n_complete_WEMWBS = sum(df$wemwbs_complete, na.rm = TRUE),
  n_complete_DASS   = sum(df$dass_complete,   na.rm = TRUE),
  any_fast_time     = sum(df$any_fast_time_flag, na.rm = TRUE),
  any_slow_time     = sum(df$any_slow_time_flag, na.rm = TRUE)
)
wcsv(qc_summary, file.path(p04, "qc_overview.csv"))

## ---- 11. Console summary -----------------------------------------------
cat("\nPhase 04 complete.\n",
    "- Descriptives saved to: ", file.path(p04, 'descriptives.csv'), "\n",
    "- Correlation matrix (wide) to: ", file.path(p04, 'correlations.csv'), "\n",
    "- Correlation long form to: ", file.path(p04, 'cor_long.csv'), "\n",
    "- Histogram/density plots saved as hist_<var>.png\n",
    if (!is.na(gender_var)) paste0("- Boxplots by gender saved as box_<var>_by_", gender_var, ".png\n") else "",
    if (!is.na(age_var)) paste0("- Age trends saved as age_trend_<var>.png\n") else "",
    "- Correlation heatmap saved as cor_tiles.png\n",
    "- QC overview summary saved as qc_overview.csv\n",
    sep = "")
