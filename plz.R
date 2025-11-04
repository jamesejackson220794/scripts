# ================================
# SEM_bayes.R  (fixed)
# ================================

# ---- 0) Packages ----
req <- c("blavaan","lavaan","psych","boot","car","openxlsx")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
library(blavaan); library(lavaan); library(psych); library(boot); library(car)

# ---- 1) Paths, output, helpers ----
root <- "/Users/jamesjackson/Desktop/dir"
phase_final <- file.path(root, "phase_final_analysis")
if (!dir.exists(phase_final)) dir.create(phase_final, recursive = TRUE, showWarnings = FALSE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(phase_final, paste0("sem_bayes_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_ok <- requireNamespace("openxlsx", quietly = TRUE)
write_dual <- function(df, stem, out_dir_local = out_dir){
  csv_path  <- file.path(out_dir_local, paste0(stem, ".csv"))
  utils::write.csv(df, csv_path, row.names = FALSE)
  if (xlsx_ok) openxlsx::write.xlsx(df, file.path(out_dir_local, paste0(stem, ".xlsx")), overwrite = TRUE)
  invisible(csv_path)
}
log_path <- file.path(out_dir, "00_log.txt")
log_append <- function(...) cat(paste0(..., "\n"), file = log_path, append = TRUE)

msg <- function(...) { m <- paste0(...); message(m); log_append(m) }

# ---- 2) Load data ----
sacs_items   <- sprintf("SACS_%d", 1:10)
cfq_items    <- sprintf("CFQ_%d", 1:7)
atq_f_items  <- sprintf("ATQ_%df", 1:15)
atq_b_items  <- sprintf("ATQ_%db", 1:15)
wemwbs_items <- sprintf("WEMWBS_%d", 1:14)
dass_items   <- sprintf("DASS.21_%d", 1:21)

scored_file <- file.path(root, "phase03_scoring", "scored_primary.csv")
raw_csv     <- file.path(root, "data", "raw copy.csv")
raw_xlsx    <- file.path(root, "data", "raw copy.xlsx")

if (file.exists(scored_file)) {
  df <- read.csv(scored_file, stringsAsFactors = FALSE, check.names = TRUE)
  msg("Loaded scored file: ", scored_file)
} else if (file.exists(raw_csv)) {
  df <- read.csv(raw_csv, stringsAsFactors = FALSE, check.names = TRUE)
  msg("Loaded raw CSV: ", raw_csv)
} else if (file.exists(raw_xlsx)) {
  if (!xlsx_ok) { install.packages("openxlsx", dependencies = TRUE); library(openxlsx) }
  df <- as.data.frame(openxlsx::read.xlsx(raw_xlsx, sheet = 1))
  msg("Loaded raw XLSX: ", raw_xlsx)
} else {
  stop("No input data found. Expecting scored_primary.csv or data/raw copy.[csv|xlsx].")
}

# Normalise DASS names (hyphen -> dot)
names(df) <- sub("^DASS-21_", "DASS.21_", names(df))

# Coerce plausible numeric vectors
num_cols <- grep("^(SACS_|CFQ_|ATQ_|WEMWBS_|DASS\\.21_|Age$)", names(df), value = TRUE)
if (length(num_cols)) df[num_cols] <- lapply(df[num_cols], function(x) suppressWarnings(as.numeric(x)))

# ---- 3) Compute/alias scale scores if missing ----
sum_if  <- function(d, items) if (all(items %in% names(d))) rowSums(d[items], na.rm = TRUE) else NULL

if (!"SACS_Total" %in% names(df)) {
  sc <- sum_if(df, sacs_items); if (!is.null(sc)) df$SACS_Total <- sc
}
if (!"SACS_Centering" %in% names(df) && all(c("SACS_1","SACS_2","SACS_5","SACS_6") %in% names(df)))
  df$SACS_Centering <- rowSums(df[, c("SACS_1","SACS_2","SACS_5","SACS_6")], na.rm = TRUE)
if (!"SACS_Transcending" %in% names(df)) {
  need <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")
  need <- need[need %in% names(df)]
  if (length(need) >= 4) df$SACS_Transcending <- rowSums(df[, need], na.rm = TRUE)
}

if (!"CFQ_Total" %in% names(df)) { sc <- sum_if(df, cfq_items); if (!is.null(sc)) df$CFQ_Total <- sc }
if (!"ATQ_Believability" %in% names(df)) { scb <- sum_if(df, atq_b_items); if (!is.null(scb)) df$ATQ_Believability <- scb }
if (!"ATQ_Frequency" %in% names(df))    { scf <- sum_if(df, atq_f_items); if (!is.null(scf)) df$ATQ_Frequency <- scf }
if ("ATQ_Bel" %in% names(df) && !"ATQ_Believability" %in% names(df)) df$ATQ_Believability <- df$ATQ_Bel
if (!"WEMWBS_Total" %in% names(df)) { sc <- sum_if(df, wemwbs_items); if (!is.null(sc)) df$WEMWBS_Total <- sc }

dass_idx <- list(dep=c(3,5,10,13,16,17,21), anx=c(2,4,7,9,15,19,20), str=c(1,6,8,11,12,14,18))
if (all(dass_items %in% names(df))) {
  df$DASS_Depression_x2 <- 2 * rowSums(df[paste0("DASS.21_", dass_idx$dep)], na.rm = TRUE)
  df$DASS_Anxiety_x2    <- 2 * rowSums(df[paste0("DASS.21_", dass_idx$anx)], na.rm = TRUE)
  df$DASS_Stress_x2     <- 2 * rowSums(df[paste0("DASS.21_", dass_idx$str)], na.rm = TRUE)
  df$DASS_Total_x2      <- 2 * rowSums(df[dass_items], na.rm = TRUE)
}
if (!"DASS_Total_x2" %in% names(df) && "DASS_Total" %in% names(df)) df$DASS_Total_x2 <- df$DASS_Total

# ---- 4) Define dataset and checks ----
dat <- df
drop9 <- TRUE
c_items <- c("SACS_1","SACS_2","SACS_5","SACS_6")
t_items <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10"); if (drop9) t_items <- setdiff(t_items, "SACS_9")
all_items <- c(c_items, t_items)

missing_items <- setdiff(all_items, names(dat))
if (length(missing_items)) stop("Missing SACS item columns: ", paste(missing_items, collapse = ", "))

# Convert SACS items to ordered
dat[all_items] <- lapply(dat[all_items], function(x) { if (is.numeric(x)) x <- round(x); as.ordered(x) })

needed <- c("CFQ_Total","ATQ_Believability","WEMWBS_Total","DASS_Total_x2")
missing_needed <- setdiff(needed, names(dat))
if (length(missing_needed)) stop("Missing observed variables required for SEM: ", paste(missing_needed, collapse = ", "))

msg("Rows: ", nrow(dat), " | Cols: ", ncol(dat))

# ---- 5) Robust missingness + descriptives ----
miss_summarise <- function(d){
  nm <- names(d)
  n  <- nrow(d)
  data.frame(
    variable = nm,
    n_miss = vapply(nm, function(v) sum(is.na(d[[v]])), integer(1)),
    pct_miss = vapply(nm, function(v) 100*mean(is.na(d[[v]])), numeric(1)),
    stringsAsFactors = FALSE
  )
}
miss_df <- miss_summarise(dat[, c(all_items, needed)])
write_dual(miss_df, "01_missingness")

write_dual(as.data.frame(psych::describe(dat[, all_items])), "02_desc_SACS_items")
write_dual(as.data.frame(psych::describe(dat[, needed])), "03_desc_outcomes")

msg("Missingness: nonzero vars = ", sum(miss_df$n_miss > 0))

# ---- 6) Bayesian SEM model ----
model_med <- paste0("
  SACSc =~ ", paste(c_items, collapse = " + "), "
  SACSt =~ ", paste(t_items, collapse = " + "), "

  CFQ_Total         ~ a_c*SACSc + a_t*SACSt
  ATQ_Believability ~ g_c*SACSc + g_t*SACSt

  DASS_Total_x2 ~ b1*CFQ_Total + b2*ATQ_Believability + c_c*SACSc + c_t*SACSt
  WEMWBS_Total  ~ d1*CFQ_Total + d2*ATQ_Believability + e_c*SACSc + e_t*SACSt

  CFQ_Total ~~ ATQ_Believability
  SACSc ~~ SACSt

  ind_cfq_dass_c := a_c*b1
  ind_cfq_dass_t := a_t*b1
  ind_atq_dass_c := g_c*b2
  ind_atq_dass_t := g_t*b2
  ind_total_dass_c := ind_cfq_dass_c + ind_atq_dass_c
  ind_total_dass_t := ind_cfq_dass_t + ind_atq_dass_t

  ind_cfq_well_c := a_c*d1
  ind_cfq_well_t := a_t*d1
  ind_atq_well_c := g_c*d2
  ind_atq_well_t := g_t*d2
  ind_total_well_c := ind_cfq_well_c + ind_atq_well_c
  ind_total_well_t := ind_cfq_well_t + ind_atq_well_t
")

fit_txt  <- file.path(out_dir, "04_bsem_fit_summary.txt")
cmp_txt  <- file.path(out_dir, "05_model_comparisons.txt")

fit_bsem <- try(
  bsem(model_med,
       data = dat,
       ordered = all_items,
       convergence = "psrf",
       target = "stan",
       mcmc = list(burnin = 3000, sample = 6000, n.chains = 4),
       seed = 123),
  silent = TRUE
)

if (inherits(fit_bsem, "try-error")) {
  msg("blavaan/Stan failed; using lavaan WLSMV fallback.")
  fit_wls <- cfa(model_med, data = dat, ordered = all_items, estimator = "WLSMV")
  smry <- capture.output(summary(fit_wls, fit.measures = TRUE, standardized = TRUE))
  writeLines(smry, fit_txt)
} else {
  smry <- capture.output(summary(fit_bsem, fit.measures = TRUE, standardized = TRUE))
  writeLines(smry, fit_txt)
  psrf <- try(blavInspect(fit_bsem, "psrf"), silent = TRUE)
  if (!inherits(psrf, "try-error")) write_dual(as.data.frame(round(psrf, 3)), "06_psrf")
  ppp  <- try(blavInspect(fit_bsem, "ppp"), silent = TRUE)
  if (!inherits(ppp, "try-error")) write_dual(data.frame(PPP = ppp), "07_ppp")
  msg("PSRF and PPP extracted if available.")
}

# ---- 7) Comparators ----
model_direct <- paste0("
  SACSc =~ ", paste(c_items, collapse = " + "), "
  SACSt =~ ", paste(t_items, collapse = " + "), "
  DASS_Total_x2 ~ SACSc + SACSt
  WEMWBS_Total  ~ SACSc + SACSt
  SACSc ~~ SACSt
")
model_1f <- paste0("
  SACSall =~ ", paste(all_items, collapse = " + "), "
  CFQ_Total         ~ SACSall
  ATQ_Believability ~ SACSall
  DASS_Total_x2 ~ SACSall + CFQ_Total + ATQ_Believability
  WEMWBS_Total  ~ SACSall + CFQ_Total + ATQ_Believability
  CFQ_Total ~~ ATQ_Believability
")

cmp_out <- list()
if (exists("fit_bsem") && !inherits(fit_bsem, "try-error")) {
  fit_direct <- bsem(model_direct, data = dat, ordered = all_items, target = "stan", mcmc = list(sample = 2000), seed = 123)
  fit_1f     <- bsem(model_1f,   data = dat, ordered = all_items, target = "stan", mcmc = list(sample = 2000), seed = 123)
  cmp_out$two_vs_direct <- capture.output( blavCompare(fit_bsem, fit_direct) )
  cmp_out$two_vs_one    <- capture.output( blavCompare(fit_bsem, fit_1f) )
} else {
  fit_direct <- cfa(model_direct, data = dat, ordered = all_items, estimator = "WLSMV")
  fit_1f     <- cfa(model_1f,     data = dat, ordered = all_items, estimator = "WLSMV")
  fit_med    <- cfa(model_med,    data = dat, ordered = all_items, estimator = "WLSMV")
  cmp_out$two_vs_direct <- capture.output( lavTestLRT(fit_direct, fit_med) )
  cmp_out$two_vs_one    <- capture.output( lavTestLRT(fit_1f,     fit_med) )
}
writeLines(unlist(cmp_out), cmp_txt)

# ---- 8) Defined parameters (indirects) ----
if (exists("fit_bsem") && !inherits(fit_bsem, "try-error")) {
  pe <- parameterEstimates(fit_bsem, standardized = TRUE, ci = TRUE)
} else {
  pe <- parameterEstimates(cfa(model_med, data = dat, ordered = all_items, estimator = "WLSMV"),
                           standardized = TRUE, ci = TRUE)
}
ind_rows <- !is.na(pe$label) & grepl("^ind_", pe$label)
if (any(ind_rows)) write_dual(pe[ind_rows, ], "08_indirects")

msg("Analysis complete. Outputs in: ", out_dir)
