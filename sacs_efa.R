## =========================================================
## SACS EFA (polychoric) — 1 vs 2 factors, robust exports
## =========================================================

suppressPackageStartupMessages({
  library(psych)         # fa(), fa.parallel, polychoric, KMO, VSS
  library(GPArotation)   # oblimin rotation
})

## ---- 0) Paths, helpers, reproducibility ----
root    <- "/Users/jamesjackson/Desktop/dir"
out_dir <- file.path(root, "phase_final_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)
write_dual <- function(df, stem, directory = out_dir){
  csv_path <- file.path(directory, paste0(stem, ".csv"))
  utils::write.csv(df, csv_path, row.names = FALSE)
  if (xlsx_ok){
    openxlsx::write.xlsx(df, file = file.path(directory, paste0(stem, ".xlsx")), overwrite = TRUE)
  }
  invisible(csv_path)
}
set.seed(2025)

## ---- 1) Get data (use df from environment or your standard files) ----
sacs_items <- paste0("SACS_", 1:10)

if (!exists("df")) {
  scored_file <- file.path(root, "phase03_scoring", "scored_primary.csv")
  if (file.exists(scored_file)) {
    df <- utils::read.csv(scored_file, check.names = TRUE, stringsAsFactors = FALSE)
  } else {
    raw_xlsx <- file.path(root, "data", "raw copy.xlsx")
    if (!file.exists(raw_xlsx)) stop("Data not found: load df, or ensure scored_primary.csv/raw copy.xlsx exist.")
    if (!xlsx_ok) install.packages("openxlsx")
    library(openxlsx)
    df <- as.data.frame(read.xlsx(raw_xlsx, sheet = 1))
    names(df) <- make.names(names(df))
  }
}
stopifnot(all(sacs_items %in% names(df)))

dat <- df[, sacs_items, drop = FALSE]
dat[] <- lapply(dat, function(x) suppressWarnings(as.numeric(x)))
dat  <- dat[stats::complete.cases(dat), , drop = FALSE]
cat("EFA N (complete SACS):", nrow(dat), "\n")

## ---- 2) Category diagnostics and conservative collapsing (if needed) ----
tab_min <- sapply(dat, function(x) min(table(factor(x, levels = sort(unique(x))))))
sparse  <- names(tab_min)[tab_min < 3]

collapse_sparse <- function(x, min_n = 5){
  x <- as.numeric(x); lev <- sort(unique(x))
  if (length(lev) <= 3) return(x)
  repeat {
    tt <- table(factor(x, levels = lev))
    if (all(tt >= min_n) || length(lev) <= 3) break
    if (tt[1] <= tt[length(tt)]) x[x == lev[1]] <- lev[2] else x[x == lev[length(lev)]] <- lev[length(lev)-1]
    lev <- sort(unique(x))
  }
  x
}

dat2 <- dat
if (length(sparse)){
  for (v in sparse) dat2[[v]] <- collapse_sparse(dat2[[v]], min_n = 5)
  cat("Collapsed extreme categories for:", paste(sparse, collapse = ", "), "\n")
}

## ---- 3) Polychoric correlations + KMO & Bartlett (with safe extraction) ----
pc <- psych::polychoric(dat2)
R  <- pc$rho

KMO_res <- psych::KMO(R)
bart    <- psych::cortest.bartlett(R, n = nrow(dat2))

# Safe KMO extraction (overall + per-item)
kmo_overall <- if (!is.null(KMO_res$MSA)) KMO_res$MSA else NA_real_
kmo_items   <- if (!is.null(KMO_res$MSAi)) KMO_res$MSAi else rep(NA_real_, ncol(R))
kmo_overall_df <- data.frame(Measure = "KMO_overall_MSA", Value = round(unname(kmo_overall), 3))
kmo_items_df   <- data.frame(Item = colnames(R), MSA = round(as.numeric(kmo_items), 3))

# Bartlett tidy
bart_df <- data.frame(Chi2 = round(bart$chisq, 3), df = bart$df, p = signif(bart$p.value, 3))

# Exports
write_dual(as.data.frame(R),                      "SACS_EFA_polychoric")
write_dual(kmo_overall_df,                        "SACS_EFA_KMO_overall")
write_dual(kmo_items_df,                          "SACS_EFA_KMO_items")
write_dual(bart_df,                               "SACS_EFA_Bartlett")

cat("KMO overall:", round(kmo_overall, 3), " | Bartlett p:", signif(bart$p.value, 3), "\n")

## ---- 4) Factor number: Parallel analysis + MAP ----
fa_par <- psych::fa.parallel(R, n.obs = nrow(dat2), fa = "fa", fm = "minres", plot = FALSE)
nf_pa  <- fa_par$nfact

# MAP from VSS (pick index of minimum MAP value)
vss_obj <- psych::VSS(R, n = 5, rotate = "none", fm = "minres", plot = FALSE)
nf_map_choice <- which.min(vss_obj$map)

n_recommend <- max(1, min(3, nf_pa)) # bound to 1..3 for SACS

write_dual(data.frame(Parallel_n = nf_pa,
                      MAP_choice = nf_map_choice,
                      Rec_n = n_recommend),
           "SACS_EFA_nfactor_recommendation")

cat("Parallel suggests:", nf_pa, " | MAP suggests:", nf_map_choice,
    " | Proceeding to 1- and 2-factor EFAs.\n")

## ---- 5) EFA solutions (MINRES, oblimin), 1 vs 2 factors ----
efa1 <- psych::fa(R, nfactors = 1, n.obs = nrow(dat2), fm = "minres", rotate = "oblimin")
efa2 <- psych::fa(R, nfactors = 2, n.obs = nrow(dat2), fm = "minres", rotate = "oblimin")

# Tidy loadings
L1 <- as.data.frame(round(efa1$loadings[, 1, drop = FALSE], 3))
L1$Item <- rownames(L1); names(L1)[1] <- "g"
L1 <- L1[, c("Item","g")]

L2 <- as.data.frame(round(efa2$loadings[, 1:2, drop = FALSE], 3))
L2$Item <- rownames(L2); names(L2)[1:2] <- c("F1","F2")
L2 <- L2[, c("Item","F1","F2")]

# Factor correlation (Phi) for 2-factor
Phi <- round(efa2$Phi, 3)

# Fit/summary table
fit_tab <- data.frame(
  Model       = c("EFA_1f","EFA_2f"),
  RMSEA       = round(c(efa1$RMSEA[1],  efa2$RMSEA[1]), 3),
  RMSEA_L     = round(c(efa1$RMSEA[2],  efa2$RMSEA[2]), 3),
  RMSEA_U     = round(c(efa1$RMSEA[3],  efa2$RMSEA[3]), 3),
  TLI         = round(c(efa1$TLI,       efa2$TLI), 3),
  BIC         = round(c(efa1$BIC,       efa2$BIC), 1),
  SS_Loadings = round(c(sum(efa1$loadings[]^2),
                        sum(efa2$loadings[]^2)), 3)
)

# Exports
write_dual(fit_tab, "SACS_EFA_fit")
write_dual(L1,      "SACS_EFA_loadings_1factor")
write_dual(L2,      "SACS_EFA_loadings_2factor")
write_dual(as.data.frame(Phi), "SACS_EFA_phi_2factor")

## ---- 6) Console decision aids (concise) ----
cat("\n--- DECISION AIDS ---\n")
cat("1) Prefer the solution with lower BIC & RMSEA, acceptable TLI.\n")
cat("2) For a true 2-factor, items {1,2,5,6} should cluster together; {3,4,7,8,9,10} on the other.\n")
cat("3) Flag cross-loaders: any |secondary| ≥ .30.\n")
cat("4) If Phi ≥ .95, treat factors as empirically one.\n")
cat("---------------------\n")

print(fit_tab, row.names = FALSE)
cat("\n2-factor Phi (factor correlation):\n"); print(Phi)
cat("\nTop of 2-factor loadings (sorted by absolute max):\n")
print(L2[order(apply(abs(L2[,c('F1','F2')]), 1, max), decreasing = TRUE), ], row.names = FALSE)

cat("\nFiles written to:\n", out_dir, "\n")

# After your EFA objects exist:
efa2$uniquenesses  # check which item is ~0 or negative
efa2$communality   # check for >1 (rare) or ~1 (problematic)

# Re-run EFA with alternative extractions (better behaved with polychorics)
efa2_uls <- psych::fa(R, nfactors=2, n.obs=nrow(dat2), fm="uls",  rotate="oblimin")
efa2_pa  <- psych::fa(R, nfactors=2, n.obs=nrow(dat2), fm="pa",   rotate="oblimin")

# Compare quickly
data.frame(
  Model = c("MINRES","ULS","PA"),
  RMSEA = c(efa2$RMSEA[1], efa2_uls$RMSEA[1], efa2_pa$RMSEA[1]),
  TLI   = c(efa2$TLI,      efa2_uls$TLI,      efa2_pa$TLI)
)

