## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dev = "svglite",
  fig.ext = "svg"
)

## ----simulate-----------------------------------------------------------------
library(spacc)
set.seed(42)
n_sites <- 80; n_species <- 40
coords <- data.frame(x = runif(n_sites, 0, 100), y = runif(n_sites, 0, 100))
species <- matrix(0L, n_sites, n_species)   # presence/absence, spatially clustered
for (sp in seq_len(n_species)) {
  cx <- runif(1, 20, 80); cy <- runif(1, 20, 80)
  prob <- exp(-0.001 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rbinom(n_sites, 1, prob)
}
colnames(species) <- paste0("sp", seq_len(n_species))

## ----basic-sac----------------------------------------------------------------
sac <- spacc(species, coords, n_seeds = 30, method = "knn", progress = FALSE)
sac

## ----curves-matrix------------------------------------------------------------
dim(sac$curves)
sac$curves[1:3, 1:6]

## ----summary-sac--------------------------------------------------------------
summary(sac)

## ----plot-sac, fig.cap = "Spatial species accumulation curve with 95% confidence ribbon."----
plot(sac)

## ----df-sac-------------------------------------------------------------------
sac_df <- as.data.frame(sac)
head(sac_df)
tail(sac_df, 3)

## ----custom-plot, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Manual ggplot from the summary data frame."----
library(ggplot2)
ggplot(sac_df, aes(sites, mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#4CAF50") +
  geom_line(linewidth = 1, colour = "#2E7D32") +
  labs(x = "Sites", y = "Cumulative species") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----compare-methods----------------------------------------------------------
sac_kncn <- spacc(species, coords, n_seeds = 30, method = "kncn", progress = FALSE)
sac_rand <- spacc(species, coords, n_seeds = 30, method = "random", progress = FALSE)

## ----combine, fig.cap = "kNN, kNCN, and random accumulation compared."--------
combined <- c(knn = sac, kncn = sac_kncn, random = sac_rand)
plot(combined)

## ----custom-order-------------------------------------------------------------
# Accumulate west-to-east (e.g. an elevation rank or survey date)
sweep_order <- order(coords$x)
sac_order <- spacc(species, coords, order = sweep_order, progress = FALSE)
sac_order

## ----compare-test-------------------------------------------------------------
sac_a <- spacc(species[, 1:20], coords, n_seeds = 30, progress = FALSE)
sac_b <- spacc(species[, 21:40], coords, n_seeds = 30, progress = FALSE)
comp <- compare(sac_a, sac_b, method = "permutation", n_perm = 199)
comp

## ----compare-df---------------------------------------------------------------
as.data.frame(comp)

## ----extrapolate--------------------------------------------------------------
fit <- extrapolate(sac, model = "lomolino")
fit

## ----plot-fit, fig.cap = "Fitted Lomolino curve with asymptote estimate."-----
plot(fit)

## ----predict-fit--------------------------------------------------------------
predict(fit, n = c(40, 80, 160, 320))

## ----coef-confint-------------------------------------------------------------
coef(fit)
confint(fit, parm = "a")

## ----compare-models-----------------------------------------------------------
cm <- compareModels(sac, models = c("michaelis-menten", "lomolino", "asymptotic"))
cm

## ----distances----------------------------------------------------------------
d <- distances(coords, method = "euclidean")
d

sac1 <- spacc(species[, 1:20], d, n_seeds = 30, progress = FALSE)
sac2 <- spacc(species[, 21:40], d, n_seeds = 30, progress = FALSE)

## ----hill---------------------------------------------------------------------
set.seed(7)
abund <- matrix(rpois(n_sites * n_species, lambda = 2), n_sites, n_species)
colnames(abund) <- colnames(species)
hill <- spaccHill(abund, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
tail(as.data.frame(hill), 3)

## ----plot-hill, fig.cap = "Hill number accumulation at q = 0, 1, 2.", eval = requireNamespace("ggplot2", quietly = TRUE)----
plot(hill)

## ----beta---------------------------------------------------------------------
beta <- spaccBeta(species, coords, n_seeds = 20, index = "sorensen", progress = FALSE)
tail(as.data.frame(beta), 3)

## ----plot-beta, fig.cap = "Beta diversity partitioned into turnover and nestedness.", eval = requireNamespace("ggplot2", quietly = TRUE)----
plot(beta)

## ----coverage-----------------------------------------------------------------
cov <- spaccCoverage(abund, coords, n_seeds = 20, coverage = "chiu", progress = FALSE)
tail(as.data.frame(cov), 3)

## ----interp-coverage----------------------------------------------------------
ic <- interpolateCoverage(cov, target = c(0.90, 0.95))
colMeans(ic)

## ----estimators---------------------------------------------------------------
chao1(abund)
chao2(species)

## ----estimators-df------------------------------------------------------------
rbind(
  as.data.frame(chao1(abund)),
  as.data.frame(chao2(species))
)

## ----beta-decay---------------------------------------------------------------
bd <- betaDecay(species, coords, index = "sorensen", model = "exponential", progress = FALSE)
bd

## ----metrics------------------------------------------------------------------
m <- spaccMetrics(species, coords, metrics = c("slope_10", "half_richness", "auc"),
                  method = "knn", progress = FALSE)
m
head(as.data.frame(m))

## ----s3-core------------------------------------------------------------------
print(sac)
head(as.data.frame(sac), 2)

## ----s3-combine---------------------------------------------------------------
grouped <- c(first = sac1, second = sac2)
print(grouped)
sub <- sac[1:5]
print(sub)

## ----s3-autoplot, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "autoplot returns a ggplot object."----
ggplot2::autoplot(sac)

## ----s3-assf, eval = requireNamespace("sf", quietly = TRUE)-------------------
m_sf <- as_sf(m)
class(m_sf)
m_sf

