# --- Setup ---
suppressPackageStartupMessages({
  if (!requireNamespace("GGally", quietly = TRUE)) install.packages("GGally")
  library(psych)       # for polychoric and describe
  library(MVN)         # for mvn()
  library(corrplot)    # for correlation plots
  library(ggplot2)
  library(GGally)
})

# --- Select SACS items ---
sacs_items <- paste0("SACS_", 1:10)
if (!exists("dat")) stop("Data object `dat` not found in global environment.")
if (!all(sacs_items %in% names(dat))) stop("Missing SACS items in `dat`.")

dat_sub <- dat[, sacs_items]

# --- Coerce to numeric for continuous approximation ---
dat_cont <- as.data.frame(lapply(dat_sub, function(x) as.numeric(as.character(x))))

# --- 1. Univariate descriptive stats ---
desc <- psych::describe(dat_cont)
print(desc[, c("mean", "sd", "skew", "kurtosis")])
cat("\n[INFO] Rule-of-thumb: |skew| > 2 or |kurtosis| > 7 = problematic.\n")

# --- 2. Multivariate normality via MVN::mvn() ---
mvn_result <- MVN::mvn(data = dat_cont)
print(mvn_result$multivariateNormality)

# Optional QQ Plot for multivariate normality
# (Current MVN versions do NOT accept multivariatePlot as argument to mvn())
plot(mvn_result)  # Produces Mardia QQ and histogram plots if interactive

# --- 3. Pearson vs Polychoric correlations ---
R_pearson  <- cor(dat_cont, use = "pairwise.complete.obs")
R_polychor <- psych::polychoric(dat_sub, correct = 0)$rho

par(mfrow = c(1,2))
corrplot(R_pearson,  title = "Pearson Correlation",  mar = c(0,0,2,0))
corrplot(R_polychor, title = "Polychoric Correlation", mar = c(0,0,2,0))
par(mfrow = c(1,1))

# --- 4. Bivariate contingency sparsity check ---
item_pairs <- combn(sacs_items, 2, simplify = FALSE)
sparse_flags <- lapply(item_pairs, function(pair) {
  tbl <- table(dat_sub[[pair[1]]], dat_sub[[pair[2]]])
  zero_cells <- sum(tbl == 0)
  list(pair = pair, zeros = zero_cells, total = length(tbl))
})
flagged <- Filter(function(x) x$zeros / x$total > 0.2, sparse_flags)

if (length(flagged) > 0) {
  cat("\n[WARNING] Pairs with >20% empty cells (may bias polychoric estimates):\n")
  print(flagged)
} else {
  cat("\n[INFO] No major sparsity detected in bivariate ordinal tables.\n")
}

# --- 5. Optional: Pairwise ordinal plots (jittered) ---
try(
  GGally::ggpairs(dat_sub, 
                  mapping = ggplot2::aes(colour = NULL), 
                  lower = list(continuous = "points")),
  silent = TRUE
)
