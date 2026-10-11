## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 5,
  dev = "svglite",
  fig.ext = "svg",
  dev.args = list(pointsize = 14)
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("spacc")            # CRAN
# 
# install.packages("pak")              # development version
# pak::pak("gcol33/spacc")

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
fits <- lapply(models, function(m) extrapolate(sac, model = m, interval = "none"))
names(fits) <- models

# Compare AIC
data.frame(
  model = models,
  asymptote = sapply(fits, function(f) round(f$asymptote, 1)),
  AIC = sapply(fits, function(f) round(f$aic, 1))
)

## ----plot-best, fig.cap = "Best-fitting asymptotic model with bootstrap band."----
best_name <- names(fits)[which.min(sapply(fits, function(f) f$aic))]
best <- extrapolate(sac, model = best_name)
plot(best) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )

## ----evt----------------------------------------------------------------------
fit_evt <- extrapolate(sac, model = "evt")
fit_evt

## ----compare-models-----------------------------------------------------------
cm <- compareModels(sac)
cm

## ----compare-table------------------------------------------------------------
as.data.frame(cm)[, c("model", "AIC", "delta_AIC", "AIC_weight", "converged")]

## ----coef-confint-------------------------------------------------------------
coef(best)
confint(best)
confint(best, method = "profile")

## ----residuals----------------------------------------------------------------
best$gof$rmse
best

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

## ----finite-units-------------------------------------------------------------
set.seed(1)
grid <- sapply(runif(300, 0.001, 0.03),
               function(p) rbinom(400, 1, p) * (1L + rpois(400, 1)))
grid <- grid[, colSums(grid) > 0]
cells <- data.frame(x = rep(1:20, 20), y = rep(1:20, each = 20))
surveyed <- sample(400, 320)
survey <- grid[surveyed, ]
detected <- colSums(survey) > 0

c(true = ncol(grid), observed = sum(detected),
  chao2 = chao2(survey)$estimate,
  MoRE = chao2(survey, n_units = 400)$estimate)

cov_wr <- spaccCoverage(survey, cells[surveyed, ], n_seeds = 20,
                        coverage = "chiu", progress = FALSE, seed = 1)
cov_grid <- spaccCoverage(survey, cells[surveyed, ], n_seeds = 20,
                          coverage = "chiu", n_units = 400,
                          progress = FALSE, seed = 1)
c(true = sum(grid[, detected]) / sum(grid),
  with_replacement = mean(cov_wr$coverage[, 320]),
  n_units = mean(cov_grid$coverage[, 320]))

## ----dar----------------------------------------------------------------------
dar_result <- dar(species, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
dar_result

## ----plot-dar, fig.cap = "Diversity-area relationship for q = 0, 1, 2."-------
plot(dar_result) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )

## ----dar-coef-----------------------------------------------------------------
coef(dar_result)

## ----predict------------------------------------------------------------------
predict(best, n = c(50, 100, 200))

## ----predict-interval---------------------------------------------------------
predict(best, n = c(50, 100, 200), interval = "bootstrap")

## ----area---------------------------------------------------------------------
area_fit <- extrapolateArea(pa, coords, n_subareas = 8,
                            model = "michaelis-menten", R = 100,
                            seed = 1, progress = FALSE)
area_fit

## ----area-predict-------------------------------------------------------------
predict(area_fit, area = c(area_fit$observed_area, 2 * area_fit$observed_area),
        interval = "bootstrap")

