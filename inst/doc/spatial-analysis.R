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
} else {
  NULL
}

## ----eval = FALSE-------------------------------------------------------------
# install.packages("spacc")            # CRAN
# 
# install.packages("pak")              # development version
# pak::pak("gcol33/spacc")

## ----simulate-----------------------------------------------------------------
library(spacc)

set.seed(42)
n_sites <- 80
n_species <- 40

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

# Species with varying range sizes (some endemic, some widespread)
species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  # First 10 species are narrow-ranged (endemics)
  spread <- if (sp <= 10) 0.005 else 0.001
  prob <- exp(-spread * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rbinom(n_sites, 1, prob)
}
colnames(species) <- paste0("sp", seq_len(n_species))

## ----range-sizes--------------------------------------------------------------
range_size <- colSums(species)
tapply(range_size, c(rep("endemic", 10), rep("widespread", 30)), mean)

## ----endemism-----------------------------------------------------------------
end <- spaccEndemism(species, coords, n_seeds = 20, progress = FALSE, seed = 1)
end

## ----endemism-summary---------------------------------------------------------
es <- summary(end)
head(es[, c("sites", "mean_richness", "mean_endemism", "endemism_proportion")], 4)

## ----plot-endemism, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Total richness and endemic richness as area expands."----
plot(end) + transparent

## ----sfar---------------------------------------------------------------------
sac <- spacc(species, coords, n_seeds = 20, progress = FALSE, seed = 1)

set.seed(123)
patches <- kmeans(coords, centers = 5)$cluster

sfar_fit <- sfar(sac, patches)
sfar_fit

## ----sfar-coef----------------------------------------------------------------
round(sfar_fit$coef, 3)

## ----plot-sfar, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "SFAR fit: observed richness and the area-plus-fragmentation prediction."----
plot(sfar_fit) + transparent

## ----sesars-------------------------------------------------------------------
set.seed(7)
effort <- rpois(n_sites, 10) + 1

sesars_fit <- sesars(sac, effort)
sesars_fit

## ----sesars-coef--------------------------------------------------------------
round(sesars_fit$coef, 3)
sesars_fit$r_squared

## ----plot-sesars, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "SESARS fit with the area-plus-effort prediction overlaid on observed richness."----
plot(sesars_fit) + transparent

## ----metrics------------------------------------------------------------------
met <- spaccMetrics(species, coords,
                    metrics = c("slope_10", "half_richness", "auc"),
                    progress = FALSE)
summary(met)

## ----metrics-df---------------------------------------------------------------
head(met$metrics[, c("site_id", "x", "y", "slope_10", "half_richness", "auc")], 4)

## ----plot-metrics, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Spatial heatmap of the initial accumulation slope per site."----
plot(met, metric = "slope_10", type = "heatmap") + transparent

## ----metrics-sf, eval = requireNamespace("sf", quietly = TRUE)----------------
met_sf <- as_sf(met, crs = 32631)
class(met_sf)

## ----mem----------------------------------------------------------------------
mem <- spatialEigenvectors(coords)
mem

## ----mem-summary--------------------------------------------------------------
head(summary(mem), 5)

## ----plot-mem, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "First spatial eigenvectors mapped over the landscape."----
plot(mem, type = "map", n_vectors = 4) + transparent

## ----partition----------------------------------------------------------------
slope <- met$metrics$slope_10
part <- spatialPartition(slope, mem)
part

## ----plot-partition, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Variance in per-site slope split into spatial and non-spatial components."----
plot(part) + transparent

## ----wavefront----------------------------------------------------------------
wf <- spaccWavefront(species, coords, n_seeds = 20, n_steps = 40,
                progress = FALSE, seed = 1)
wf

## ----plot-wavefront, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Species captured against expanding radius from seed points."----
plot(wf) + transparent

## ----subsample----------------------------------------------------------------
keep <- subsample(coords, method = "grid", cell_size = 20, seed = 1)
length(keep)

sac_full <- spacc(species, coords, n_seeds = 20, progress = FALSE, seed = 1)
sac_thin <- spacc(species[keep, ], coords[keep, ], n_seeds = 20,
                  progress = FALSE, seed = 1)
combined <- c(full = sac_full, thinned = sac_thin)

## ----plot-subsample, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Accumulation curve before and after grid thinning."----
plot(combined) + transparent

