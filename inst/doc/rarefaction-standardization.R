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

## ----theme-helper, include = FALSE--------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  transparent <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA)
  )
}

## ----eval = FALSE-------------------------------------------------------------
# install.packages("spacc")            # CRAN
# 
# install.packages("pak")              # development version
# pak::pak("gcol33/spacc")

## ----data---------------------------------------------------------------------
library(spacc)
set.seed(42)

n_sites <- 60
coords <- data.frame(x = runif(n_sites), y = runif(n_sites))

# Variable total abundances across sites (realistic uneven sampling)
lambdas <- rep(c(1, 3, 5), each = 20)
species <- matrix(0, nrow = n_sites, ncol = 20)
for (i in seq_len(n_sites)) {
  species[i, ] <- rpois(20, lambda = lambdas[i])
}
colnames(species) <- paste0("sp", 1:20)

## ----effort-bias--------------------------------------------------------------
tier <- factor(lambdas)
indiv <- rowSums(species)
rich <- rowSums(species > 0)
aggregate(cbind(individuals = indiv, richness = rich), list(tier = tier), mean)

## ----rarefy, eval = requireNamespace("ggplot2", quietly = TRUE)---------------
# Rarefy to the pooled abundance grid (q = 0, richness)
rare <- rarefy(species)
print(rare)

## ----rarefy-extract-----------------------------------------------------------
rare_df <- as.data.frame(rare)
head(rare_df, 4)
tail(rare_df, 2)

## ----rarefy-plot, eval = requireNamespace("ggplot2", quietly = TRUE)----------
plot(rare) + transparent

## ----rarefy-hill--------------------------------------------------------------
rare_q1 <- rarefy(species, q = 1)
rare_q2 <- rarefy(species, q = 2)
c(q0 = max(rare$expected), q1 = max(rare_q1$expected), q2 = max(rare_q2$expected))

## ----rarefy-anyq--------------------------------------------------------------
rare_q3 <- rarefy(species, q = 3)
rare_half <- rarefy(species, q = 0.5)
c(q0.5 = max(rare_half$expected), q3 = max(rare_q3$expected))

## ----analytical---------------------------------------------------------------
mt <- mao_tau(species)
cm <- coleman(species)
cl <- collector(species)
head(mt, 3)

## ----analytical-plot, eval = requireNamespace("ggplot2", quietly = TRUE)------
df <- rbind(
  data.frame(sites = mt$sites, S = mt$expected, kind = "Mao Tau"),
  data.frame(sites = cm$sites, S = cm$expected, kind = "Coleman"),
  data.frame(sites = cl$sites, S = cl$species,  kind = "Collector")
)
ggplot2::ggplot(df, ggplot2::aes(sites, S, color = kind)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(x = "Sites", y = "Species", color = NULL) +
  ggplot2::theme_minimal(base_size = 17) + transparent

## ----coverage, eval = requireNamespace("ggplot2", quietly = TRUE)-------------
cov_result <- spaccCoverage(species, coords, n_seeds = 30, progress = FALSE)
print(cov_result)

## ----coverage-plot, eval = requireNamespace("ggplot2", quietly = TRUE)--------
plot(cov_result) + transparent

## ----interpolate--------------------------------------------------------------
interp <- interpolateCoverage(cov_result, target = c(0.90, 0.95, 0.99))
round(colMeans(interp), 2)
round(apply(interp, 2, sd), 3)

## ----extrapolate--------------------------------------------------------------
extrap <- extrapolateCoverage(cov_result, target_coverage = c(0.95, 0.99), q = 0)
print(extrap)
summary(extrap)

## ----extrapolate-plot, eval = requireNamespace("ggplot2", quietly = TRUE)-----
plot(extrap) + transparent

## ----hill-coverage, eval = requireNamespace("ggplot2", quietly = TRUE)--------
hc <- spaccHillCoverage(species, coords, q = c(0, 1, 2),
                        target_coverage = 0.9,
                        n_seeds = 20, progress = FALSE)
print(hc)

## ----hill-coverage-std--------------------------------------------------------
sapply(hc$standardized, mean)

## ----hill-coverage-plot, eval = requireNamespace("ggplot2", quietly = TRUE)----
plot(hc, xaxis = "coverage") + transparent

## ----subsample----------------------------------------------------------------
set.seed(7)
keep <- subsample(coords, method = "grid", cell_size = 0.2, seed = 7)
length(keep)

## ----thinned-vs-full, eval = requireNamespace("ggplot2", quietly = TRUE)------
full  <- spacc(species, coords, n_seeds = 30, progress = FALSE)
thin  <- spacc(species[keep, ], coords[keep, ], n_seeds = 30, progress = FALSE)
df_full <- transform(as.data.frame(full), set = "Full")
df_thin <- transform(as.data.frame(thin), set = "Thinned (grid)")

## ----thinned-plot, eval = requireNamespace("ggplot2", quietly = TRUE)---------
ggplot2::ggplot(rbind(df_full, df_thin),
                ggplot2::aes(sites, mean, color = set)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(x = "Sites", y = "Species", color = NULL) +
  ggplot2::theme_minimal(base_size = 17) + transparent

## ----spatial-rare-------------------------------------------------------------
sr <- spatialRarefaction(species, coords, n_perm = 50)
tail(sr[, c("sites", "mean", "lower", "upper")], 2)

## ----tier-split---------------------------------------------------------------
lo <- species[lambdas == 1, ]
hi <- species[lambdas == 5, ]
raw <- c(low = sum(colSums(lo) > 0), high = sum(colSums(hi) > 0))
raw

## ----matched-effort-----------------------------------------------------------
n_match <- min(sum(lo), sum(hi))
S_lo <- rarefy(lo, n_individuals = n_match)$expected
S_hi <- rarefy(hi, n_individuals = n_match)$expected
c(low = round(S_lo, 1), high = round(S_hi, 1), n = n_match)

## ----matched-coverage---------------------------------------------------------
co_lo <- spaccCoverage(lo, coords[lambdas == 1, ], n_seeds = 20, progress = FALSE)
co_hi <- spaccCoverage(hi, coords[lambdas == 5, ], n_seeds = 20, progress = FALSE)
S_lo_c <- mean(interpolateCoverage(co_lo, target = 0.9)[, 1])
S_hi_c <- mean(interpolateCoverage(co_hi, target = 0.9)[, 1])
c(low = round(S_lo_c, 1), high = round(S_hi_c, 1))

## ----comparison-table---------------------------------------------------------
data.frame(
  standardization = c("Raw richness", "Matched individuals", "Matched coverage (0.90)"),
  low  = c(raw["low"],  round(S_lo, 1), round(S_lo_c, 1)),
  high = c(raw["high"], round(S_hi, 1), round(S_hi_c, 1)),
  row.names = NULL
)

