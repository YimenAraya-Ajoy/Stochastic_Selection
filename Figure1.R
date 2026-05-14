source("Results.R")

# ── Bridge: map analysis objects to figure objects ────────────────────────────
n_years     <- nlevels(d$year)
year_levels <- levels(d$year)

# Fixed effect point estimates
estsL <- fixL[, "Estimate"]
estsR <- fixR[, "Estimate"]

# Variance component medians
V1  <- V1_L_q
V2  <- V2_L_q
C12 <- Cov_L_q

# Optimum point estimates
optL <- median(optL_draws, na.rm = TRUE)
optR <- median(optR_draws, na.rm = TRUE)

# ── Full posterior draws matrix (base R) ──────────────────────────────────────
draws_mat <- as.matrix(drawsL)

# Fixed effects
b0_draws  <- draws_mat[, "b_Intercept"]
b1_draws  <- draws_mat[, "b_c_z"]
b2_draws  <- draws_mat[, "b_c_z2"]

# Variance components (brms stores SDs, square for variances)
sd0_draws  <- draws_mat[, "sd_year__Intercept"]
sd1_draws  <- draws_mat[, "sd_year__c_z"]
cor_draws  <- draws_mat[, "cor_year__Intercept__c_z"]

V1_draws   <- sd0_draws^2
V2_draws   <- sd1_draws^2
C12_draws  <- cor_draws * sd0_draws * sd1_draws

# Sanity check
stopifnot(length(b0_draws) == length(V1_draws))

# ── Year-specific random effects ──────────────────────────────────────────────
int_names   <- paste0("r_year[", year_levels, ",Intercept]")
slope_names <- paste0("r_year[", year_levels, ",c_z]")

# Posterior means of random effects
coefs_L <- colMeans(draws_mat[, c(int_names, slope_names)])

alpha_t <- estsL["Intercept"] + coefs_L[int_names]
beta_t  <- estsL["c_z"]       + coefs_L[slope_names]
gamma_t <- rep(estsL["c_z2"], n_years)

# ── Grids ─────────────────────────────────────────────────────────────────────
z      <- seq(min(d$c_z, na.rm = TRUE), max(d$c_z, na.rm = TRUE), by = 0.01)
z_mean <- z
n.z    <- length(z)
z_mat  <- matrix(z, nrow = n.z, ncol = 1)

# ── Point-estimate fitness surfaces ───────────────────────────────────────────
r_hat     <- estsL["Intercept"] + estsL["c_z"] * z_mean + estsL["c_z2"] * z_mean^2
sigma_E_z <- V1[2] + V2[2] * z_mean^2 + 2 * C12[2] * z_mean
r_arith   <- r_hat + 0.5 * sigma_E_z
r_stoch   <- r_hat - 0.5 * sigma_E_z
w_z       <- estsR["Intercept"] + estsR["c_z"] * z_mean + estsR["c_z2"] * z_mean^2

# ── Annual fitness functions matrix (n.z x n_years) ───────────────────────────
Alpha <- matrix(alpha_t, nrow = n.z, ncol = n_years, byrow = TRUE)
Beta  <- matrix(beta_t,  nrow = 1,   ncol = n_years, byrow = TRUE)
Gamma <- matrix(gamma_t, nrow = 1,   ncol = n_years, byrow = TRUE)
eta   <- Alpha + z_mat %*% Beta + (z_mat^2) %*% Gamma

# ── Posterior CI for stochastic growth rate ───────────────────────────────────
r_stoch_draws <- sapply(z_mean, function(zi) {
  r_hat_i   <- b0_draws + b1_draws * zi + b2_draws * zi^2
  sigma_e_i <- V1_draws + V2_draws * zi^2 + 2 * C12_draws * zi
  r_hat_i - 0.5 * sigma_e_i
})
# r_stoch_draws: n_draws x n.z
r_stoch_lo <- apply(r_stoch_draws, 2, quantile, 0.025)
r_stoch_hi <- apply(r_stoch_draws, 2, quantile, 0.975)

# ── Helper: vertical line to maximum ─────────────────────────────────────────
vline_to_max <- function(x, y, lty = 1, col = "black") {
  i <- which.max(y)
  segments(x0 = x[i], y0 = -4, x1 = x[i], y1 = y[i], lty = lty, col = col)
}

# ── Figure 1 ──────────────────────────────────────────────────────────────────
out_file <- "/home/yi/Dropbox/Apps/Overleaf/FluctuatingSelection/Figures/Figure1.pdf"
pdf(out_file, width = 9, height = 5)
par(mfrow = c(1, 2))

# Panel A: annual fitness functions
eta_mean <- exp(rowMeans(eta))
plot(eta_mean ~ z,
     type = "l", ylim = c(0, 3),
     ylab = "Number of recruits",
     xlab = "Clutch size (sd)")
for (i in seq_len(n_years)) lines(z, exp(eta[, i]), col = "gray")
lines(z, eta_mean, lwd = 2)
segments(x0 = optL, y0 = -1,
         x1 = optL, y1 = exp(r_hat[which.max(r_hat)]))
mtext("A)", side = 3, adj = 0)

# Panel B: adaptive topographies
plot(r_arith ~ z_mean,
     type = "l", lty = 2,
     ylab = "Log fitness",
     xlab = "Mean clutch size (sd)",
     ylim = c(-1, 0.5), xlim = c(-2, 4))
lines(z_mean, r_hat,    lty = 4)
lines(z_mean, log(w_z), lty = 3)

# CI band then solid line for stochastic growth rate
polygon(c(z_mean, rev(z_mean)),
        c(r_stoch_lo, rev(r_stoch_hi)),
        col = adjustcolor("black", alpha.f = 0.15),
        border = NA)
lines(z_mean, r_stoch, lty = 1, lwd = 2)

vline_to_max(z_mean, r_arith,  lty = 2)
vline_to_max(z_mean, r_hat,    lty = 4)
vline_to_max(z_mean, r_stoch,  lty = 1)
vline_to_max(z_mean, log(w_z), lty = 3)

legend("topright",
       legend = c(
         expression(log ~ E[t](W[t](z))),
         expression(hat(r)(z)),
         expression(tilde(r)(bar(z))),
         expression(w(z))
       ),
       lty = c(2, 4, 1, 3),
       lwd = c(1, 1, 2, 1),
       bty = "n")
mtext("B)", side = 3, adj = 0)

dev.off()

