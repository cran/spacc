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

## ----transparent, include = FALSE---------------------------------------------
transparent <- if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA)
  )
} else NULL

## ----eval = FALSE-------------------------------------------------------------
# install.packages("spacc")            # CRAN
# 
# install.packages("pak")              # development version
# pak::pak("gcol33/spacc")

## ----data---------------------------------------------------------------------
library(spacc)
set.seed(123)

n_sites <- 60
coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

## ----data-species-------------------------------------------------------------
n_species <- 25
species <- matrix(0, nrow = n_sites, ncol = n_species)
centers_x <- runif(n_species, 0, 100)
centers_y <- runif(n_species, 0, 100)

for (sp in seq_len(n_species)) {
  dists <- sqrt((coords$x - centers_x[sp])^2 + (coords$y - centers_y[sp])^2)
  prob <- exp(-dists / 30)
  species[, sp] <- rpois(n_sites, lambda = prob * 3)
}
colnames(species) <- paste0("sp", seq_len(n_species))

## ----data-sac-----------------------------------------------------------------
sac <- spacc(species, coords, n_seeds = 10, progress = FALSE)

## ----decay, eval = requireNamespace("ggplot2", quietly = TRUE)----------------
decay <- betaDecay(species, coords, progress = FALSE)
plot(decay) + transparent

## ----decay-summary------------------------------------------------------------
summary(decay)

## ----decay-coef---------------------------------------------------------------
coef(decay)
decay$half_life

## ----dist-decay, eval = requireNamespace("ggplot2", quietly = TRUE)-----------
dd <- distanceDecay(species, coords, n_seeds = 30, progress = FALSE, seed = 1)
plot(dd) + transparent

## ----zeta, eval = requireNamespace("ggplot2", quietly = TRUE)-----------------
zd <- zetaDiversity(species, coords, orders = 1:10, n_samples = 50,
                    progress = FALSE, seed = 1)
plot(zd) + transparent

## ----zeta-ratio, eval = requireNamespace("ggplot2", quietly = TRUE)-----------
plot(zd, type = "ratio") + transparent

## ----zeta-df------------------------------------------------------------------
zeta_tab <- as.data.frame(zd)
head(zeta_tab)

## ----ses, eval = requireNamespace("ggplot2", quietly = TRUE)------------------
result <- ses(sac, species, coords,
              null_model = "curveball", n_perm = 19,
              progress = FALSE, seed = 1)
print(result)

## ----ses-plot, eval = requireNamespace("ggplot2", quietly = TRUE)-------------
plot(result, type = "curve") + transparent

## ----ses-hist, eval = requireNamespace("ggplot2", quietly = TRUE)-------------
plot(result, type = "histogram") + transparent

## ----ses-compare--------------------------------------------------------------
res_freq <- ses(sac, species, coords,
                null_model = "frequency", n_perm = 19,
                progress = FALSE, seed = 1)
c(curveball = mean(result$ses), frequency = mean(res_freq$ses))

