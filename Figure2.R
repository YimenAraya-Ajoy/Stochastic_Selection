source("Results.R")

# ── Extract posterior draws for variance components ───────────────────────────
sd0_draws  <- drawsL$sd_year__Intercept
sd1_draws  <- drawsL$sd_year__c_z
cor_draws  <- drawsL$cor_year__Intercept__c_z

Vmu0_draws  <- sd0_draws^2
Vmu1_draws  <- sd1_draws^2
Cmu01_draws <- cor_draws * sd0_draws * sd1_draws

# Compute σ²_e(z̄) across the full posterior at each z̄ value
n_z <- length(z_mean)
sigma_E_draws <- matrix(NA, nrow = length(sd0_draws), ncol = n_z)
for (i in seq_along(sd0_draws)) {
  sigma_E_draws[i, ] <- Vmu0_draws[i] +
                        Vmu1_draws[i]  * z_mean^2 +
                        2 * Cmu01_draws[i] * z_mean
}

# Posterior summaries at each z̄
sigma_E_mean <- apply(sigma_E_draws, 2, mean)
sigma_E_lo95 <- apply(sigma_E_draws, 2, quantile, 0.025)
sigma_E_hi95 <- apply(sigma_E_draws, 2, quantile, 0.975)

# ── Figure 2 ─────────────────────────────────────────────────────────────────
out_file2 <- "/home/yi/Dropbox/Apps/Overleaf/FluctuatingSelection/Figures/Figure2.pdf"
pdf(out_file2, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

# ── Panel A: covariance penalty c(z, z̄) ──────────────────────────────────────
z_vals    <- c(-1, 0, 1)
z_bar_seq <- z_mean
cov_by_z  <- sapply(z_vals, function(zi) {
  V1[2] + V2[2] * zi * z_bar_seq + C12[2] * (zi + z_bar_seq)
})

plot(cov_by_z[, 2] ~ z_bar_seq,
     type = "l", lwd = 2,
     ylab = expression("Environ. covariance " * c(z, bar(z))),
     xlab = expression("Mean clutch size " * bar(z) * " (sd)"),
     ylim = range(cov_by_z))
lines(z_bar_seq, cov_by_z[, 1], lwd = 2, lty = 2)
lines(z_bar_seq, cov_by_z[, 3], lwd = 2, lty = 3)
mtext("A)", 3, adj = 0)

zbar0 <- which(round(z_bar_seq, 1) == 2)[1]
text(2, cov_by_z[zbar0, 1], labels = expression(z == -1), pos = 1, cex = 0.9)
text(2, cov_by_z[zbar0, 2], labels = expression(z ==  0), pos = 1, cex = 0.9)
text(2, cov_by_z[zbar0, 3], labels = expression(z ==  1), pos = 1, cex = 0.9)

# ── Panel B: σ²_e(z̄) posterior mean + 95% credible band ─────────────────────
plot(sigma_E_z ~ z_mean, type = "n",
     ylab = expression("Environ. stochasticity " * sigma[e]^2),
     xlab = expression("Mean clutch size " * bar(z) * " (sd)"),
     ylim = c(0, 3.5))
polygon(c(z_mean, rev(z_mean)),
        c(sigma_E_lo95, rev(sigma_E_hi95)),
        col = adjustcolor("gray70", alpha.f = 0.4),
        border = NA)
lines(z_mean, sigma_E_mean, lwd = 2)
mtext("B)", 3, adj = 0)

dev.off()
