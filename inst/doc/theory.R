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
library(spacc)
set.seed(1)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("spacc")            # CRAN
# 
# install.packages("pak")              # development version
# pak::pak("gcol33/spacc")

## ----toy-points---------------------------------------------------------------
pts <- data.frame(
  x = c(1, 2, 1.5, 3, 3.5, 5, 5.5, 6, 2, 4,
        7, 7.5, 8, 1, 6.5, 4.5, 2.5, 8.5, 5, 3),
  y = c(1, 1.2, 2, 1, 2.5, 1, 2, 3, 3.5, 4,
        1.5, 3, 2, 4.5, 4, 3.5, 4.8, 4, 5, 5.5)
)

## ----toy-knn-order------------------------------------------------------------
knn_order <- function(coords, focus) {
  order(sqrt((coords$x - focus[1])^2 + (coords$y - focus[2])^2))
}
focus <- c(3, 2.5)
ord <- knn_order(pts, focus)

## ----toy-knn-plot, fig.width=6, fig.height=5, fig.alt="Twenty points numbered by increasing distance from one fixed focal point."----
plot(pts$x, pts$y, pch = 19, col = "grey70", cex = 1.4,
     xlab = "x", ylab = "y", main = "Fixed-focus kNN ordering")
points(focus[1], focus[2], pch = 4, col = "#C62828", cex = 2, lwd = 2)
text(pts$x, pts$y, labels = match(seq_len(nrow(pts)), ord), pos = 3, cex = 1)

## ----method-schematic, fig.width=6.5, fig.height=6, fig.alt="Four panels showing the same twenty points coloured by visiting order under fixed-focus kNN, a nearest-neighbour walk, random, and collector rules."----
pal <- function(o) grDevices::hcl.colors(length(o), "Greens", rev = TRUE)[order(o)]
nn_walk_order <- function(coords, seed) {
  n <- nrow(coords); visited <- logical(n); result <- integer(n)
  current <- seed; visited[current] <- TRUE; result[1] <- current
  for (step in 2:n) {
    distances <- sqrt((coords$x - coords$x[current])^2 +
                      (coords$y - coords$y[current])^2)
    distances[visited] <- Inf
    current <- which.min(distances)
    visited[current] <- TRUE
    result[step] <- current
  }
  result
}

orders <- list(
  kNN       = knn_order(pts, focus),
  nn_walk   = nn_walk_order(pts, 1),
  random    = sample(nrow(pts)),
  collector = seq_len(nrow(pts))
)
op <- par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
for (nm in names(orders)) {
  rank_k <- match(seq_len(nrow(pts)), orders[[nm]])
  plot(pts$x, pts$y, pch = 19, cex = 1.6, col = pal(rank_k),
       xlab = "", ylab = "", main = nm, axes = FALSE); box()
}
par(op)

## ----impl-call----------------------------------------------------------------
coords  <- data.frame(x = runif(120), y = runif(120))
species <- matrix(rbinom(120 * 40, 1, 0.25), nrow = 120)

sac <- spacc(species, coords, method = "knn", n_seeds = 100,
             progress = FALSE, seed = 42)
sac

## ----impl-order---------------------------------------------------------------
ord_mat <- t(replicate(20, sample(nrow(species))))
sac_user <- spacc(species, coords, order = ord_mat, progress = FALSE)
sac_user$n_seeds

## ----backend-equiv------------------------------------------------------------
e <- spacc(species, coords, method = "nn_walk", backend = "exact",
           n_seeds = 30, progress = FALSE, seed = 7)
k <- spacc(species, coords, method = "nn_walk", backend = "kdtree",
           n_seeds = 30, progress = FALSE, seed = 7)
c(exact_mean_end = mean(e$curves[, ncol(e$curves)]),
  kdtree_mean_end = mean(k$curves[, ncol(k$curves)]))

## ----seed-band, fig.width=6.5, fig.height=4.5, fig.alt="A spatial accumulation curve with a shaded confidence band. The mean curve rises and saturates; the band is widest in the early, low-effort region where the choice of starting site matters most."----
plot(sac)

## ----seed-count---------------------------------------------------------------
b50  <- spacc(species, coords, n_seeds = 50,  progress = FALSE, seed = 1)
b300 <- spacc(species, coords, n_seeds = 300, progress = FALSE, seed = 1)
mid <- round(ncol(b50$curves) / 4)
c(seeds_50  = diff(quantile(b50$curves[,  mid], c(.025, .975))),
  seeds_300 = diff(quantile(b300$curves[, mid], c(.025, .975))))

## ----compare-random, fig.width=6.5, fig.height=4.5, fig.alt="Two accumulation curves on the same axes: the random-order curve rises faster and the spatial kNN curve lags behind it, with the shaded gap between them representing spatial turnover."----
sp <- spacc(species, coords, method = "knn",    n_seeds = 100, progress = FALSE, seed = 3)
rd <- spacc(species, coords, method = "random", n_seeds = 100, progress = FALSE, seed = 3)
cmp <- compare(sp, rd)
plot(cmp)

## ----downstream-taste, fig.width=6.5, fig.height=4.5, fig.alt="Hill-number accumulation curves for orders q=0, q=1, and q=2 along the spatial ordering, with the q=0 richness curve highest and the q=2 inverse-Simpson curve lowest."----
hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 50, progress = FALSE, seed = 9)
plot(hill)

## ----session-info-------------------------------------------------------------
sessionInfo()

