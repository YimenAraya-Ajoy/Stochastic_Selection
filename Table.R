require(brms)
library(brms)
library(lme4)

source("Analysis.R")
# Fixed effects
fixL <- round(fixef(modL), 3)
fixR <- round(fixef(modR), 3)

# Posterior draws
drawsL <- as_draws_df(modL)
drawsR <- as_draws_df(modR)

# Optima (clutch size)
optL_draws <- drawsL$b_c_z / (-2 * drawsL$b_c_z2)
optR_draws <- drawsR$b_c_z / (-2 * drawsR$b_c_z2)

round(quantile(optL_draws, c(.025, .5, .975), na.rm = TRUE), 3)
round(quantile(optR_draws, c(.025, .5, .975), na.rm = TRUE), 3)

# Random effect variance components - log-linear model
V1_L  <- drawsL$sd_year__Intercept^2
V2_L  <- drawsL$sd_year__c_z^2
Cov_L <- drawsL$sd_year__Intercept * 
         drawsL$sd_year__c_z * 
         drawsL$cor_year__Intercept__c_z

V1_L_q   <- round(quantile(V1_L,  c(.025, .5, .975), na.rm = TRUE), 3)
V2_L_q   <- round(quantile(V2_L,  c(.025, .5, .975), na.rm = TRUE), 3)
Cov_L_q  <- round(quantile(Cov_L, c(.025, .5, .975), na.rm = TRUE), 3)

# Random effect variance components - relative fitness model
# Note: modR has no random intercept, so only random slope variance
V2_R  <- drawsR$sd_year__c_z^2
V2_R_q <- round(quantile(V2_R, c(.025, .5, .975), na.rm = TRUE), 3)

# Unexplained variance
shapeL <- exp(drawsL$b_Intercept) + 
          exp(drawsL$b_Intercept)^2 / drawsL$shape
shapeL_q <- round(quantile(shapeL, c(.025, .5, .975), na.rm = TRUE), 3)

sigma2_R  <- drawsR$sigma^2
sigma2_R_q <- round(quantile(sigma2_R, c(.025, .5, .975), na.rm = TRUE), 3)

# Helper function for formatting
fmt <- function(q) paste0(q[2], " (", q[1], ", ", q[3], ")")

# Build table
row_names <- c(
  "\\textit{Fixed effects}",
  "$\\beta_0$",
  "$\\beta_1$",
  "$\\beta_2$",
  "$N$",
  "\\textit{Random effects}",
  "$V_{\\mu_0}$",
  "$V_{\\mu_1}$",
  "$C_{\\mu_{01}}$",
  "\\textit{Unexplained}"
)

dR <- matrix("", nrow = length(row_names), ncol = 2)
colnames(dR) <- c("Log-linear", "Relative fitness")
rownames(dR) <- row_names

# Fixed effects - note row indices match β0, β1, β2, N
# Check that fixef row order matches: Intercept, c_z, c_z2, N
dR[2, 1] <- fmt(fixL["Intercept",   c("Q2.5","Estimate","Q97.5")])
dR[3, 1] <- fmt(fixL["c_z",         c("Q2.5","Estimate","Q97.5")])
dR[4, 1] <- fmt(fixL["c_z2",        c("Q2.5","Estimate","Q97.5")])
dR[5, 1] <- fmt(fixL["N",           c("Q2.5","Estimate","Q97.5")])

dR[2, 2] <- fmt(fixR["Intercept",   c("Q2.5","Estimate","Q97.5")])
dR[3, 2] <- fmt(fixR["c_z",  c("Q2.5","Estimate","Q97.5")])
dR[4, 2] <- fmt(fixR["c_z2", c("Q2.5","Estimate","Q97.5")])
dR[5, 2] <- fmt(fixR["N",    c("Q2.5","Estimate","Q97.5")])

# Random effects
dR[7, 1] <- fmt(V1_L_q)
dR[8, 1] <- fmt(V2_L_q)
dR[9, 1] <- fmt(Cov_L_q)
dR[10, 1] <- fmt(shapeL_q)

dR[7, 2] <- "—"          # no random intercept in modR
dR[8, 2] <- fmt(V2_R_q)
dR[9, 2] <- "—"          # no Cµ01 in modR
dR[10, 2] <- fmt(sigma2_R_q)

# Print table
library(xtable)

tab <- print(
  xtable(
    data.frame(dR),
    caption = "Posterior medians and 95\\% credible intervals for fixed and 
               random effects from Bayesian mixed-effects models of great tit 
               recruitment.",
    align = c("l", "l", "l")
  ),
  floating = FALSE,
  include.rownames = TRUE,
  sanitize.text.function = identity,
  sanitize.colnames.function = function(x) {
    paste0("\\multicolumn{1}{c}{", x, "}")
  }
)

writeLines(tab, "/home/yi/Dropbox/Apps/Overleaf/FluctuatingSelection/model_results.tex")
