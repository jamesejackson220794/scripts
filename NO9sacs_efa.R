library(lavaan)

# Use your collapsed ordinal data 'dat2' from the EFA step:
dat_ord <- as.data.frame(lapply(dat2, function(x) ordered(x, levels = sort(unique(x)))))

center_items <- c("SACS_1","SACS_2","SACS_5","SACS_6")
transc_core  <- c("SACS_3","SACS_4","SACS_7","SACS_8","SACS_10")

mod2_09 <- '
Centring =~ SACS_1 + SACS_2 + SACS_5 + SACS_6
Transcendence =~ SACS_3 + SACS_4 + SACS_7 + SACS_8 + SACS_10
Centring ~~ Transcendence
'

fit2_09 <- cfa(mod2_09, data=dat_ord, estimator="WLSMV",
               ordered=names(dat_ord), parameterization="theta")

fi <- function(f) round(fitMeasures(f, c("cfi","tli","rmsea","srmr","chisq","df","pvalue")),3)
print(fi(fit2_09))

# Heywood safety check
pe <- parameterEstimates(fit2_09, standardized = TRUE)
any(pe$op=="~~" & pe$lhs==pe$rhs & pe$est < 0)  # should be FALSE
