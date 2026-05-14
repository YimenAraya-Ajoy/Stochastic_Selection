library(brms)

# Load data and create a Section-Year identifier
d1 <- read.csv("data.csv")
d1$SY <- paste(d1$year, d1$Section)

# Drop rows with missing hatch date (required for analysis)
d <- d1[complete.cases(d1$hatchdate), ]

# --- Clutch size (cs): z-score standardisation ---
m_cs  <- mean(d$cs, na.rm = TRUE)
sd_cs <- sd(d$cs, na.rm = TRUE)
d$c_z  <- (d$cs - m_cs) / sd_cs   # global z-score
d$c_z2 <- d$c_z^2                  # quadratic term

# --- Lay date (layday): within-year centring after global z-scoring ---
d$l_z    <- as.vector(scale(d$layday))            # global z-score
l_zc_m   <- tapply(d$l_z, d$year, mean)           # per-year mean of z-score
d$l_zc_m <- as.vector(l_zc_m[match(d$year, names(l_zc_m))])  # map back to rows
d$l_zc   <- d$l_z - d$l_zc_m     # within-year centred lay date
d$l_zc2  <- d$l_zc^2              # quadratic term

# --- Year as factor + numeric index for tapply indexing ---
d$year  <- as.factor(d$year)
d$year2 <- as.numeric(d$year)

# --- Composite fitness: survival + recruitment ---
d$r2s <- d$surv + d$recruit

# --- Population density: scaled number of individuals per year ---
N    <- table(d$year)
d$N  <- scale(as.vector(scale(N[match(d$year, names(N))])))

# --- Relative recruitment and fledgling success (individual / annual mean) ---
m_rec       <- tapply(d$recruit,     d$year2, mean, na.rm = TRUE)
d$m_rec     <- m_rec[d$year2]
d$rel_rec   <- d$recruit / d$m_rec

# --- Model 1: absolute recruitment (counts) ---
# Negative binomial to handle overdispersion in count data.
# Random slopes for clutch size and lay date allowed to co-vary across years.
modL <- brm(
  recruit ~ c_z * l_zc + c_z2 + l_zc2 + N + (c_z + l_zc | year),
  family = negbinomial(link = "log"),
  data   = d,
  chains = 3, cores = 3, iter = 4000
)

# --- Model 2: relative recruitment (individual / annual mean) ---
# Gaussian model for ratio-scaled fitness.
# No random intercept (0 +) so only slopes vary by year,
# avoiding confounding with the within-year centering already applied.
modR <- brm(
  rel_rec ~ c_z * l_zc + c_z2 + l_zc2 + N + (0 + c_z + l_zc | year),
  family = gaussian(),
  data   = d,
  chains = 3, cores = 3, iter = 4000,
  seed   = 1
)
