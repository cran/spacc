## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dev = "svglite",
  fig.ext = "svg"
)

## ----data---------------------------------------------------------------------
library(spacc)
set.seed(1)

S_true <- 120
rel <- exp(rnorm(S_true, mean = 0, sd = 1.5))
rel <- rel / sum(rel)

n_sites <- 80
species <- t(sapply(seq_len(n_sites),
                    function(i) rmultinom(1, size = 25, prob = rel)))
species <- species[, colSums(species) > 0, drop = FALSE]  # detected species
colnames(species) <- paste0("sp", seq_len(ncol(species)))

## ----data-summary-------------------------------------------------------------
S_obs <- ncol(species)
ab <- colSums(species)
data.frame(S_true = S_true, S_obs = S_obs,
           f1 = sum(ab == 1), f2 = sum(ab == 2),
           f3 = sum(ab == 3), f4 = sum(ab == 4))

## ----chao1--------------------------------------------------------------------
chao1(species)
spacc::ace(species)

## ----chao2--------------------------------------------------------------------
chao2(species)
jackknife(species, order = 1)
jackknife(species, order = 2)

## ----bootstrap----------------------------------------------------------------
bootstrap_richness(species, n_boot = 100)

## ----comparison---------------------------------------------------------------
estimators <- list(
  chao1(species), chao2(species), spacc::ace(species),
  jackknife(species, order = 1), jackknife(species, order = 2),
  bootstrap_richness(species, n_boot = 100)
)
tab <- do.call(rbind, lapply(estimators, as.data.frame))
tab$gap_to_truth <- S_true - tab$estimate
tab

## ----ichao--------------------------------------------------------------------
r_ichao1 <- iChao1(species)
r_ichao2 <- iChao2(species)
r_ichao1
r_ichao2

## ----ichao-compare------------------------------------------------------------
ich <- rbind(
  as.data.frame(chao1(species)), as.data.frame(r_ichao1),
  as.data.frame(chao2(species)), as.data.frame(r_ichao2)
)
ich$gap_to_truth <- S_true - ich$estimate
ich

## ----convergence--------------------------------------------------------------
effort <- c(10, 20, 30, 40, 50, 60, 70, 80)
conv <- do.call(rbind, lapply(effort, function(k) {
  reps <- t(vapply(1:30, function(r) {
    e <- chao1(species[sample(nrow(species), k), , drop = FALSE])
    c(e$S_obs, e$estimate, e$se)
  }, numeric(3)))
  data.frame(sites = k, S_obs = mean(reps[, 1]),
             estimate = mean(reps[, 2]), se = mean(reps[, 3]))
}))
round(conv, 1)

## ----convergence-plot, eval = requireNamespace("ggplot2", quietly = TRUE)-----
library(ggplot2)
ggplot(conv, aes(sites)) +
  geom_ribbon(aes(ymin = estimate - se, ymax = estimate + se),
              fill = "#4CAF50", alpha = 0.2) +
  geom_line(aes(y = estimate), color = "#4CAF50", linewidth = 1) +
  geom_line(aes(y = S_obs), color = "#78909C", linewidth = 1) +
  geom_hline(yintercept = S_true, linetype = "dashed") +
  labs(x = "Sampling units", y = "Richness",
       title = "Chao1 estimate (green) and observed count (grey) vs effort") +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----ace-threshold------------------------------------------------------------
do.call(rbind, lapply(c(5, 10, 15, 20), function(th) {
  e <- spacc::ace(species, threshold = th)
  data.frame(threshold = th, estimate = round(e$estimate, 1),
             se = round(e$se, 1))
}))

## ----boot-nboot---------------------------------------------------------------
set.seed(7)
do.call(rbind, lapply(c(50, 100, 500), function(nb) {
  e <- bootstrap_richness(species, n_boot = nb)
  data.frame(n_boot = nb, estimate = round(e$estimate, 1),
             se = round(e$se, 2))
}))

## ----inext, eval = requireNamespace("iNEXT", quietly = TRUE)------------------
cr <- iNEXT::ChaoRichness(as.numeric(colSums(species)), datatype = "abundance")
sp_c1 <- chao1(species)
data.frame(
  source = c("spacc::chao1", "iNEXT::ChaoRichness", "truth"),
  estimate = round(c(sp_c1$estimate, cr[["Estimator"]], S_true), 2),
  se = round(c(sp_c1$se, cr[["Est_s.e."]], NA), 2)
)

## ----completeness-------------------------------------------------------------
cp <- completenessProfile(species)
cp

## ----completeness-plot, eval = requireNamespace("ggplot2", quietly = TRUE)----
plot(cp) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

## ----guidance-table-----------------------------------------------------------
data.frame(
  data_type = c("abundance", "abundance", "abundance",
                "incidence", "incidence", "incidence"),
  estimator = c("chao1", "ace", "iChao1",
                "chao2", "jackknife", "iChao2"),
  use_when = c("counts reliable, f2 > 0",
               "many rare species, heterogeneous",
               "f3, f4 > 0, sample incomplete",
               "presence/absence, m >= 20",
               "incidence, simple bias correction",
               "Q3, Q4 > 0, sample incomplete"),
  unreliable_when = c("f2 = 0", "C_ace <= 0 (all singletons)", "f4 = 0",
                      "Q2 = 0 or m < 10", "few units", "Q4 = 0")
)

