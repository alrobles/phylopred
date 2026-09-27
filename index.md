# phylopred

Predicting species interactions from phylogenetic distance under
presence-only data.

`phylopred` models the probability that two species share a parasite (or
other interaction) as a function of their phylogenetic distance, and
predicts new hosts for a parasite from its known hosts:

- Pairwise host-sharing logistic models, parasite-cluster bootstrap, and
  formal slope comparison between clades (generalist vs specialist).
- Presence-only performance metrics: threshold efficiency
  `f_gamma = a / (gamma * b + c)` and a phylogenetic partial ROC
  (sensitivity vs fraction of the phylogeny predicted suitable).
- Positive-unlabeled (PU) host prediction with
  [xplus](https://CRAN.R-project.org/package=xplus) on phylogenetic
  features.
- Leave-hosts-out validation of predictive performance.

## Installation

``` r

remotes::install_github("alrobles/phylopred")
```

Documentation site: <https://alrobles.github.io/phylopred/>

## Example

``` r

library(phylopred)

pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
fit <- fit_phylopred_model(pairs)
summary(fit)

parasite <- names(sort(table(pairs$parasite), decreasing = TRUE))[1]
one <- pairs[pairs$parasite == parasite, ]
one$pred <- predict(fit$fit, newdata = one, type = "response")
m <- host_threshold_metric(one$pred, one$suscept)
plot(m)
roc <- phylo_partial_roc(one$pred, one$suscept, error_rate = 0.95)
plot(roc)
```

See the vignette `caligus-vs-lepeophtheirus` for a complete real-data
case study with the `cofid` copepod-fish interaction database and the
Fish Tree of Life.

## References

- Robles-Fernandez, A. L. & Lira-Noriega, A. (2017). Combining
  phylogenetic and occurrence information for risk assessment of pest
  and pathogen interactions with host plants. *Frontiers in Applied
  Mathematics and Statistics*, 3:17. <doi:10.3389/fams.2017.00017>
- Peterson, A. T., Papes, M. & Soberon, J. (2008). Rethinking receiver
  operating characteristic analysis applications in ecological niche
  modeling. *Ecological Modelling*, 213, 63-72.
- Rabosky, D. L. and colleagues (2018). An inverse latitudinal gradient
  in speciation rate for marine fishes. *Nature*, 559, 392-395.
- Morales-Serna, F. N. (2025). Global patterns of modularity and narrow
  host use in fish-parasitic copepods (Crustacea). *Biodiversity Data
  Journal*, 13, e163693. <https://doi.org/10.3897/BDJ.13.e163693> —
  basis of the cofid database (<https://github.com/alrobles/cofid>)
