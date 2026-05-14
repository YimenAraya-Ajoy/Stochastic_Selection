source("Results.R")

# ── Grid ──────────────────────────────────────────────────────────────────────
z_grid    <- seq(-1, 1, length.out = 50)
zbar_grid <- seq(-1, 1, length.out = 50)

Z    <- outer(z_grid, zbar_grid, function(z, zb) z)
Zbar <- outer(z_grid, zbar_grid, function(z, zb) zb)

# ── Stochastic relative fitness m̃(z, z̄) ──────────────────────────────────────
r_hat_mat <- estsL["Intercept"] + estsL["c_z"] * Z + estsL["c_z2"] * Z^2
penalty   <- V1[2] + C12[2] * (Z + Zbar) + V2[2] * (Z * Zbar)
Mtilde    <- r_hat_mat - penalty

# ── Ridge z*(z̄) and equilibrium z̄* ───────────────────────────────────────────
z_star  <- (estsL["c_z"] - (C12[2] + V2[2] * zbar_grid)) / 
           (-2 * estsL["c_z2"])
z_eq    <- -(estsL["c_z"] - C12[2]) / (2 * estsL["c_z2"] - V2[2])
zbar_eq <- z_eq

# ── Colour palette ────────────────────────────────────────────────────────────
cols <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf",
                           "#fdae61", "#d7191c"))(50)

# ── Plot ──────────────────────────────────────────────────────────────────────
out_file <- "/home/yi/Dropbox/Apps/Overleaf/FluctuatingSelection/Figures/Figure3.pdf"
pdf(out_file, width = 8, height = 6)
par(mar = c(4.5, 5.5, 2.5, 1))

filled.contour(
  z_grid, zbar_grid, Mtilde,
  color.palette = function(n) cols,
  nlevels = 50,
  xlab = expression("Individual clutch size " ~ (z)),
  ylab = expression("Population mean clutch size " ~ (bar(z))),
  key.axes = {
    axis(4)
    text(x = par("usr")[2] + 2,
         y = mean(par("usr")[3:4]),
         labels = expression("Relative fitness:"~
                                  tilde(m)(z, bar(z))),
         srt = 270,
         xpd = NA,
         cex = 1)                                    # <-- moved outside expression()
  },
  plot.axes = {
    axis(1); axis(2)
    lines(z_star, zbar_grid, col = "black", lwd = 1)
    abline(0, 1, lty = 2, col = "gray30")
    points(z_eq, zbar_eq, pch = 19, cex = 1.3)
    text(z_eq, zbar_eq,
         labels = expression(bar(z)^"*"),
         pos = 4, cex = 1.2)
    text(-0.5, -0.6,
         labels = expression(z == bar(z)),
         cex = 1.2)
  }
)


dev.off()
