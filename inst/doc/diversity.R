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

set.seed(123)
n_sites <- 60
n_species <- 30

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

# Abundance matrix with spatial clustering
species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  lambda <- 5 * exp(-0.001 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rpois(n_sites, lambda)
}
colnames(species) <- paste0("sp", seq_len(n_species))

## ----sim-peek-----------------------------------------------------------------
dim(species)
species[1:4, 1:8]
sum(species > 0) / length(species)   # fraction of nonzero cells

## ----hill---------------------------------------------------------------------
hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
hill

## ----plot-hill, fig.cap = "Hill number accumulation for q = 0, 1, 2."---------
library(ggplot2)
plot(hill) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----hill-df------------------------------------------------------------------
hill_df <- as.data.frame(hill)
head(hill_df, 3)
# Final-site values per q
hill_df[hill_df$sites == max(hill_df$sites), c("q", "mean", "lower", "upper")]

## ----hill-summary-------------------------------------------------------------
hs <- summary(hill, ci_level = 0.90)
tail(hs, 3)

## ----abg----------------------------------------------------------------------
# Alpha diversity: per-site Hill numbers
alpha <- alphaDiversity(species, q = c(0, 1, 2))
colMeans(alpha)

# Gamma diversity: pooled regional diversity
gamma <- gammaDiversity(species, q = c(0, 1, 2))
gamma

# Full partition: gamma = alpha * beta
partition <- diversityPartition(species, q = c(0, 1, 2))
partition

## ----partition-summary--------------------------------------------------------
summary(partition)

## ----ci-band------------------------------------------------------------------
ci <- as.data.frame(hill)
ci0 <- ci[ci$q == 0, ]
# Width of the 95% band at a few site counts
data.frame(sites = c(5, 20, 40, 60),
           width = ci0$upper[c(5, 20, 40, 60)] - ci0$lower[c(5, 20, 40, 60)])

## ----evenness-----------------------------------------------------------------
ev <- evenness(species, q = seq(0.1, 3, by = 0.1))
ev

## ----plot-evenness, fig.cap = "Hill evenness profile across orders."----------
plot(ev) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----evenness-pielou----------------------------------------------------------
evenness(species, type = "pielou")

## ----profile------------------------------------------------------------------
prof <- diversityProfile(species, q = seq(0, 3, by = 0.2))
prof

## ----plot-profile, fig.cap = "Diversity profile: per-site mean (solid) and regional gamma (dashed)."----
plot(prof) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----beta---------------------------------------------------------------------
pa <- (species > 0) * 1L
beta <- spaccBeta(pa, coords, n_seeds = 20, progress = FALSE)
beta

## ----plot-beta, fig.cap = "Spatial beta diversity accumulation with turnover/nestedness."----
plot(beta) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----beta-df------------------------------------------------------------------
tail(as.data.frame(beta), 3)

## ----beta-pool----------------------------------------------------------------
band <- cut(coords$x, breaks = c(-Inf, 33, 66, Inf), labels = c("west", "centre", "east"))
beta_pool <- spaccBeta(pa, coords, n_seeds = 20, progress = FALSE, pool = band)
beta_pool

## ----plot-beta-pool, fig.cap = "Beta partition between sites and between pooled bands."----
plot(beta_pool, resolution = "both") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----hillbeta-----------------------------------------------------------------
hb <- spaccHillBeta(species, coords, q = c(0, 1, 2), n_seeds = 15, progress = FALSE)
hb

## ----plot-hillbeta, fig.cap = "Hill beta (effective number of communities) along accumulation."----
plot(hb, component = "beta") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----hillbeta-homog-----------------------------------------------------------
set.seed(42)
lambda <- rlnorm(80, -0.5, 1.3)
homog <- matrix(rpois(60 * 80, rep(lambda, each = 60)), nrow = 60)
homog_xy <- data.frame(x = runif(60, 0, 100), y = runif(60, 0, 100))

hb_obs <- spaccHillBeta(homog, homog_xy, n_seeds = 15, progress = FALSE, seed = 1)
hb_std <- spaccHillBeta(homog, homog_xy, n_seeds = 15, progress = FALSE, seed = 1,
                        standardize = "coverage")
at <- c(2, 10, 30, 60)
drift <- rbind(observed = colMeans(hb_obs$beta$q0)[at],
               standardised = colMeans(hb_std$beta$q0)[at])
colnames(drift) <- paste(at, "sites")
drift

## ----beta-func----------------------------------------------------------------
# Simulate two continuous traits
traits <- data.frame(
  body_size = rnorm(n_species),
  wing_length = rnorm(n_species)
)
rownames(traits) <- colnames(species)

beta_func <- spaccBeta(pa, coords, traits = traits, n_seeds = 20, progress = FALSE)

## ----plot-beta-func, fig.cap = "Functional beta diversity accumulation."------
plot(beta_func) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----beta-phylo, eval = requireNamespace("ape", quietly = TRUE)---------------
library(ape)
tree <- rcoal(n_species, tip.label = colnames(species))

beta_phylo <- spaccBeta(pa, coords, tree = tree, n_seeds = 20, progress = FALSE)

## ----plot-beta-phylo, fig.cap = "Phylogenetic beta diversity accumulation.", eval = requireNamespace("ape", quietly = TRUE)----
plot(beta_phylo) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----phylo, eval = requireNamespace("ape", quietly = TRUE)--------------------
phylo_acc <- spaccPhylo(pa, coords, tree,
                        metric = c("mpd", "mntd", "pd"),
                        n_seeds = 20, progress = FALSE)

## ----plot-phylo, fig.cap = "Phylogenetic diversity accumulation.", eval = requireNamespace("ape", quietly = TRUE)----
plot(phylo_acc) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----evoheritage, eval = requireNamespace("ape", quietly = TRUE)--------------
evo_acc <- spaccPhylo(pa, coords, tree,
                      metric = "evoheritage", rho = c(0, 0.5, 5),
                      n_seeds = 20, progress = FALSE)
evo_df <- summary(evo_acc)
evo_df[evo_df$sites == max(evo_df$sites), c("metric", "mean", "lower", "upper")]
sum(colSums(pa) > 0)   # observed richness

## ----func---------------------------------------------------------------------
func_acc <- spaccFunc(species, coords, traits,
                      metric = c("fdis"),
                      n_seeds = 20, progress = FALSE)

## ----plot-func, fig.cap = "Functional diversity accumulation."----------------
plot(func_acc) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----func-k-------------------------------------------------------------------
func_k <- spaccFunc(species, coords, traits, metric = "K",
                    n_seeds = 20, progress = FALSE)
k_df <- summary(func_k)
k_df[k_df$sites %in% c(1, 10, max(k_df$sites)), c("metric", "sites", "mean")]

## ----profile-func-------------------------------------------------------------
fp <- diversityProfile(species, traits = traits, q = seq(0, 3, by = 0.5))
fp

## ----redundancy, eval = requireNamespace("ape", quietly = TRUE)---------------
red <- spaccRedundancy(species, coords, tree = tree, n_seeds = 20, progress = FALSE)
red

## ----plot-redundancy, fig.cap = "Taxonomic over phylogenetic Hill number along accumulation.", eval = requireNamespace("ape", quietly = TRUE)----
plot(red) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----rao-phylo, eval = requireNamespace("ape", quietly = TRUE)----------------
phylo_rao <- spaccPhylo(species, coords, tree,
                        metric = "rao",
                        n_seeds = 20, progress = FALSE)

## ----plot-rao-phylo, fig.cap = "Phylogenetic Rao's Q accumulation.", eval = requireNamespace("ape", quietly = TRUE)----
plot(phylo_rao) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----rao-func-----------------------------------------------------------------
func_rao <- spaccFunc(species, coords, traits,
                      metric = "rao",
                      n_seeds = 20, progress = FALSE)
tail(as.data.frame(func_rao), 3)

## ----rao-cover----------------------------------------------------------------
cover <- species / max(species)   # fractional values in [0, 1]
func_rao_cover <- spaccFunc(cover, coords, traits,
                            metric = "rao", n_seeds = 10, progress = FALSE)
tail(as.data.frame(func_rao_cover), 3)

## ----coverage-----------------------------------------------------------------
cov <- spaccCoverage(species, coords, n_seeds = 20, progress = FALSE)
cov

## ----plot-coverage, fig.cap = "Coverage-based spatial rarefaction."-----------
plot(cov) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----coverage-chiu------------------------------------------------------------
cov_chiu <- spaccCoverage(species, coords, coverage = "chiu",
                          n_seeds = 20, progress = FALSE)
tail(as.data.frame(cov_chiu), 3)

## ----interpolate--------------------------------------------------------------
interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
summary(interp)

## ----hillcov------------------------------------------------------------------
hc <- spaccHillCoverage(species, coords, q = c(0, 1, 2),
                        n_seeds = 15, progress = FALSE)

## ----plot-hillcov, fig.cap = "Hill numbers plotted against sample coverage."----
plot(hc, xaxis = "coverage") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----custom-metric------------------------------------------------------------
# Shannon entropy of the cumulative community
shannon <- function(comm) {
  p <- comm[comm > 0] / sum(comm)
  -sum(p * log(p))
}

div <- spaccDiversity(species, coords, shannon,
                      method = "knn", n_seeds = 20, progress = FALSE)

## ----plot-custom, fig.cap = "Custom (Shannon entropy) accumulation curve."----
plot(div, ylab = "Cumulative Shannon entropy") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

