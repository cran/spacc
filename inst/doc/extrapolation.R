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
n_sites <- 100
n_species <- 50

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  lambda <- 4 * exp(-0.0008 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rpois(n_sites, lambda)
}
colnames(species) <- paste0("sp", seq_len(n_species))
pa <- (species > 0) * 1L

## ----sac----------------------------------------------------------------------
sac <- spacc(pa, coords, n_seeds = 30, progress = FALSE)

## ----models-------------------------------------------------------------------
models <- c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic")
fits <- lapply(models, function(m) extrapolate(sac, model = m))
names(fits) <- models

# Compare AIC
data.frame(
  model = models,
  asymptote = sapply(fits, function(f) round(f$asymptote, 1)),
  AIC = sapply(fits, function(f) round(f$aic, 1))
)

## ----plot-best, fig.cap = "Best-fitting asymptotic model."--------------------
best <- fits[[which.min(sapply(fits, function(f) f$aic))]]
plot(best) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )

## ----evt----------------------------------------------------------------------
fit_evt <- extrapolate(sac, model = "evt")
fit_evt

## ----compare-models-----------------------------------------------------------
cm <- compareModels(sac, progress = FALSE)
cm

## ----compare-table------------------------------------------------------------
as.data.frame(cm)[, c("model", "AIC", "delta_AIC", "AIC_weight", "converged")]

## ----coef-confint-------------------------------------------------------------
coef(best)
suppressMessages(confint(best))

## ----residuals----------------------------------------------------------------
obs <- as.data.frame(best)
resid <- obs$observed - obs$predicted
round(c(rmse = sqrt(mean(resid^2)), max_abs = max(abs(resid))), 3)

## ----seeds--------------------------------------------------------------------
sapply(c(5, 20, 50), function(k) {
  s <- spacc(pa, coords, n_seeds = k, progress = FALSE)
  round(extrapolate(s, model = "michaelis-menten")$asymptote, 1)
})

## ----cov-extrap---------------------------------------------------------------
cov <- spaccCoverage(species, coords, n_seeds = 20, progress = FALSE)

ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99, 1.0), q = 0)
ext

## ----plot-cov-extrap, fig.cap = "Coverage-based extrapolation of species richness."----
plot(ext) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )

## ----dar, eval = requireNamespace("sf", quietly = TRUE)-----------------------
dar_result <- dar(species, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
dar_result

## ----plot-dar, fig.cap = "Diversity-area relationship for q = 0, 1, 2.", eval = requireNamespace("sf", quietly = TRUE)----
plot(dar_result) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )

## ----predict------------------------------------------------------------------
predict(best, n = c(50, 100, 200, 500))

