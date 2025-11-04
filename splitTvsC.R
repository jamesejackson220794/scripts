## =========================================================
## SACS MG-CFA (ordinal, WLSMV) with automatic category pooling
## - Two-group splits by DASS_Total_x2 and WEMWBS_Total (median)
## - Configural & metric fits; Δ-fit; latent r(SACSc,SACSt) per group
## - HTMT via polychoric (psych); optional semTools::htmt if available
## =========================================================

suppressPackageStartupMessages({
  library(lavaan)
  ok_psych    <- requireNamespace("psych",     quietly = TRUE)
  ok_semTools <- requireNamespace("semTools",  quietly = TRUE)
})

## ---------- 0) Inputs & checks ----------
need <- c(paste0("SACS_", 1:10), "DASS_Total_x2", "WEMWBS_Total")
if (!exists("dat")) stop("`dat` not found in Global Env.")
miss <- setdiff(need, names(dat))
if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse=", "))

## ---------- 1) Options ----------
drop9            <- TRUE        # set FALSE to keep SACS_9
estimator        <- "WLSMV"
std.lv           <- TRUE
min_per_group    <- 1           # per-group minimum cell count per category (1 is usually enough)
report_fits      <- c("cfi","tli","rmsea","srmr","bic")

## ---------- 2) Item sets ----------
c_items_full <- c("SACS_1","SACS_2","SACS_5","SACS_6")
t_items_full <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_9","SACS_10")
if (drop9) {
  t_items <- setdiff(t_items_full, "SACS_9")
  c_items <- c_items_full
} else {
  t_items <- t_items_full
  c_items <- c_items_full
}
all_items <- c(c_items, t_items)

## ---------- 3) Utilities ----------
as_ordered_like <- function(x) {
  if (is.ordered(x)) return(x)
  if (is.numeric(x)) return(ordered(x))
  ordered(as.character(x))
}

median_split <- function(x, low_name="Low", high_name="High") {
  med <- median(x, na.rm = TRUE)
  factor(ifelse(x <= med, low_name, high_name), levels = c(low_name, high_name))
}

## Pool adjacent **edge** categories to remove zero/sparse cells by group
pool_adjacent_per_group <- function(df, item, group_var, min_per_group = 1) {
  x <- df[[item]]
  if (!is.ordered(x)) x <- as_ordered_like(x)
  g <- df[[group_var]]
  
  lev <- levels(x); xi <- as.integer(x); K <- length(lev)
  mapping <- seq_len(K); pooled <- FALSE
  
  tab_counts <- function(xi, g, mapping) {
    xm     <- mapping[xi]
    uniq   <- sort(unique(xm))
    renum  <- match(xm, uniq)
    T      <- table(g, renum)
    list(T=T, uniq=uniq, renum=renum)
  }
  
  repeat {
    tc <- tab_counts(xi, g, mapping); T <- tc$T
    bad <- which(T < min_per_group, arr.ind = TRUE)
    if (length(bad) == 0) break
    
    Kcurr <- ncol(T)
    if (Kcurr <= 2) {
      # Collapse smallest edge and stop
      col_tot <- colSums(T)
      edge <- if (col_tot[1] <= col_tot[Kcurr]) 1L else Kcurr
      uniq <- sort(unique(mapping))
      if (edge == 1L) { mapping[mapping == uniq[1]] <- uniq[2] }
      else            { mapping[mapping == uniq[Kcurr]] <- uniq[Kcurr-1L] }
      pooled <- TRUE
      break
    }
    
    col_tot <- colSums(T)
    edge <- if (col_tot[1] <= col_tot[Kcurr]) 1L else Kcurr
    uniq <- sort(unique(mapping))
    if (edge == 1L) { mapping[mapping == uniq[1]] <- uniq[2] }
    else            { mapping[mapping == uniq[length(uniq)]] <- uniq[length(uniq)-1L] }
    pooled <- TRUE
  }
  
  xm <- mapping[xi]; uniq <- sort(unique(xm)); xf_int <- match(xm, uniq)
  # Build human-readable merged labels
  orig_lab <- setNames(lev, seq_along(lev))
  groups <- lapply(uniq, function(code) which(mapping == code))
  new_labels <- vapply(groups, function(idx) paste(orig_lab[idx], collapse="|"), character(1))
  xf <- ordered(xf_int, labels = new_labels)
  attr(xf, "pooled") <- pooled
  attr(xf, "map")    <- list(mapping = mapping, new_labels = new_labels)
  xf
}

apply_pooling_for_split <- function(df, items, group_var, min_per_group = 1) {
  out <- df
  log_map <- list()
  for (it in items) {
    out[[it]] <- as_ordered_like(out[[it]])
    T0 <- table(out[[group_var]], out[[it]])
    if (any(T0 < min_per_group)) {
      pooled <- pool_adjacent_per_group(out, it, group_var, min_per_group = min_per_group)
      out[[it]] <- pooled
      log_map[[it]] <- attr(pooled, "map")
    } else {
      log_map[[it]] <- list(mapping = seq_along(levels(out[[it]])),
                            new_labels = levels(out[[it]]))
    }
  }
  attr(out, "pool_log") <- log_map
  out
}

## HTMT via polychoric (handles unequal category counts per item)
htmt_poly_twofactor <- function(df, c_items, t_items) {
  if (!ok_psych) return(NA_real_)
  sub <- df[, c(c_items, t_items), drop = FALSE]
  PC  <- psych::polychoric(sub, correct = 0)  # prints "Converted non-numeric..." sometimes; harmless
  R   <- PC$rho
  colnames(R) <- rownames(R) <- c(c_items, t_items)
  
  mean_within <- function(items) {
    if (length(items) < 2) return(NA_real_)
    M <- abs(R[items, items, drop = FALSE])
    mean(M[upper.tri(M, diag = FALSE)], na.rm = TRUE)
  }
  mean_between <- function(A, B) mean(abs(R[A, B, drop = FALSE]), na.rm = TRUE)
  
  num  <- mean_between(c_items, t_items)
  denC <- mean_within(c_items)
  denT <- mean_within(t_items)
  as.numeric(num / sqrt(denC * denT))
}

## Two-group CFA (configural & metric)
fit_two_group <- function(df, group_var, model, ord_vars, estimator="WLSMV", std.lv=TRUE) {
  if (!is.factor(df[[group_var]])) df[[group_var]] <- factor(df[[group_var]])
  if (nlevels(df[[group_var]]) != 2) stop("Grouping var must have exactly 2 levels.")
  
  fit_config <- cfa(model, data = df, ordered = ord_vars,
                    estimator = estimator, std.lv = std.lv, group = group_var)
  
  fit_metric <- cfa(model, data = df, ordered = ord_vars,
                    estimator = estimator, std.lv = std.lv, group = group_var,
                    group.equal = c("loadings"))
  
  fm_conf <- fitMeasures(fit_config, report_fits)
  fm_metr <- fitMeasures(fit_metric, report_fits)
  delta   <- setNames(fm_metr - fm_conf, paste0("d.", names(fm_conf)))
  
  cor_lv_conf <- lavInspect(fit_config, "cor.lv")
  get_r <- function(clist) {
    r1 <- clist[[1]]["SACSc","SACSt"]; r2 <- clist[[2]]["SACSc","SACSt"]
    c(low = as.numeric(r1), high = as.numeric(r2))
  }
  r_latent <- get_r(cor_lv_conf)
  
  list(
    fit_config = fit_config,
    fit_metric = fit_metric,
    fit_conf_meas = fm_conf,
    fit_metr_meas = fm_metr,
    fit_delta = delta,
    r_latent = r_latent
  )
}

## ---------- 4) Coerce items globally to ordered ----------
for (v in all_items) dat[[v]] <- as_ordered_like(dat[[v]])

## ---------- 5) Build grouping factors ----------
dat$grp_DASS   <- median_split(dat$DASS_Total_x2,  low_name="Low_DASS",   high_name="High_DASS")
dat$grp_WEMBWS <- median_split(dat$WEMWBS_Total,   low_name="Low_WEMWBS", high_name="High_WEMWBS")

## ---------- 6) 2-factor model ----------
model_2f <- paste0(
  "SACSc =~ ", paste(c_items, collapse=" + "), "\n",
  "SACSt =~ ", paste(t_items, collapse=" + "), "\n",
  "SACSc ~~ SACSt\n"
)

## ---------- 7) Run by DASS split (with pre-pooling) ----------
cat("\n===== Two-group MG-CFA by DASS_Total_x2 (", ifelse(drop9,"9 items; SACS_9 dropped","10 items"), ") =====\n", sep="")
dat_DASS <- dat[, c(all_items, "grp_DASS"), drop = FALSE]
dat_DASS_p <- apply_pooling_for_split(dat_DASS, all_items, "grp_DASS", min_per_group = min_per_group)

pool_log <- attr(dat_DASS_p, "pool_log")
changed <- vapply(pool_log, function(m) {
  # changed if mapping collapsed any codes OR labels show merged categories
  (length(unique(m$mapping)) < length(m$mapping)) || any(grepl("\\|", m$new_labels))
}, logical(1))
if (any(changed)) {
  cat("\n[INFO] Category pooling performed for DASS split on items:\n")
  print(names(which(changed)))
}

res_DASS <- fit_two_group(df = dat_DASS_p, group_var = "grp_DASS",
                          model = model_2f, ord_vars = all_items,
                          estimator = estimator, std.lv = std.lv)

print(rbind(CONFIGURAL = res_DASS$fit_conf_meas, METRIC = res_DASS$fit_metr_meas), digits = 3)
cat("\nΔ (METRIC - CONFIGURAL):\n"); print(res_DASS$fit_delta, digits = 3)

cat("\nLatent r(SACSc, SACSt) by group (CONFIGURAL):\n")
cat("  Low_DASS : ", sprintf("%.3f", res_DASS$r_latent["low"]),  "\n", sep="")
cat("  High_DASS: ", sprintf("%.3f", res_DASS$r_latent["high"]), "\n", sep="")

# HTMT (polychoric)
htmt_DASS_low  <- htmt_poly_twofactor(subset(dat_DASS_p, grp_DASS=="Low_DASS"),  c_items, t_items)
htmt_DASS_high <- htmt_poly_twofactor(subset(dat_DASS_p, grp_DASS=="High_DASS"), c_items, t_items)
cat("\nHTMT (polychoric) by DASS group:\n")
cat("  Low_DASS : ", sprintf("%.3f", htmt_DASS_low),  "\n", sep="")
cat("  High_DASS: ", sprintf("%.3f", htmt_DASS_high), "\n", sep="")

if (ok_semTools) {
  cat("\n[semTools::htmt] (configural model):\n")
  # Do NOT pass unsupported args; be version-agnostic
  try(print(semTools::htmt(res_DASS$fit_config)), silent = TRUE)
}

## ---------- 8) Run by WEMWBS split (with pre-pooling) ----------
cat("\n===== Two-group MG-CFA by WEMWBS_Total (", ifelse(drop9,"9 items; SACS_9 dropped","10 items"), ") =====\n", sep="")
dat_W <- dat[, c(all_items, "grp_WEMBWS"), drop = FALSE]   # <-- correct name
dat_W_p <- apply_pooling_for_split(dat_W, all_items, "grp_WEMBWS", min_per_group = min_per_group)

pool_log_W <- attr(dat_W_p, "pool_log")
changed_W <- vapply(pool_log_W, function(m) {
  (length(unique(m$mapping)) < length(m$mapping)) || any(grepl("\\|", m$new_labels))
}, logical(1))
if (any(changed_W)) {
  cat("\n[INFO] Category pooling performed for WEMWBS split on items:\n")
  print(names(which(changed_W)))
}

res_W <- fit_two_group(df = dat_W_p, group_var = "grp_WEMBWS",
                       model = model_2f, ord_vars = all_items,
                       estimator = estimator, std.lv = std.lv)

print(rbind(CONFIGURAL = res_W$fit_conf_meas, METRIC = res_W$fit_metr_meas), digits = 3)
cat("\nΔ (METRIC - CONFIGURAL):\n"); print(res_W$fit_delta, digits = 3)

cat("\nLatent r(SACSc, SACSt) by group (CONFIGURAL):\n")
cat("  Low_WEMWBS : ", sprintf("%.3f", res_W$r_latent["low"]),  "\n", sep="")
cat("  High_WEMWBS: ", sprintf("%.3f", res_W$r_latent["high"]), "\n", sep="")

# HTMT (polychoric)
htmt_W_low  <- htmt_poly_twofactor(subset(dat_W_p, grp_WEMBWS=="Low_WEMWBS"),  c_items, t_items)
htmt_W_high <- htmt_poly_twofactor(subset(dat_W_p, grp_WEMBWS=="High_WEMWBS"), c_items, t_items)
cat("\nHTMT (polychoric) by WEMWBS group:\n")
cat("  Low_WEMWBS : ", sprintf("%.3f", htmt_W_low),  "\n", sep="")
cat("  High_WEMWBS: ", sprintf("%.3f", htmt_W_high), "\n", sep="")

if (ok_semTools) {
  cat("\n[semTools::htmt] (configural model):\n")
  try(print(semTools::htmt(res_W$fit_config)), silent = TRUE)
}

## ---------- 9) Interpretation prompts ----------
cat("\n\n====== INTERPRETATION GUIDE ======\n",
    "- Sparse categories per group were repaired by conservative edge pooling.\n",
    "- Invariance: prefer ΔCFI ≤ .01 and ΔRMSEA ≤ .015 (configural → metric).\n",
    "- Latent r(SACSc,SACSt): lower values = better discriminant validity.\n",
    "- HTMT < .85–.90 supports discriminant validity; treat as descriptive with small N.\n",
    "- With n≈88 (≈44/group), focus on pattern-level evidence; CIs will be wide.\n", sep="")
