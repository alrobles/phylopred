# Changelog

## phylopred 0.1.0

Initial release. Clean reimplementation of the modern `geotax` workflow
for predicting species interactions from phylogenetic distance:

- [`prepare_pair_data()`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md),
  [`fit_phylopred_model()`](https://alrobles.github.io/phylopred/reference/fit_phylopred_model.md),
  [`compare_phylopred_slopes()`](https://alrobles.github.io/phylopred/reference/compare_phylopred_slopes.md)
  — pairwise host-sharing logistic models and formal slope comparison
  between clades.
- [`bootstrap_phylopred_model()`](https://alrobles.github.io/phylopred/reference/bootstrap_phylopred_model.md),
  [`cluster_bootstrap_phylopred()`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md),
  [`phylopred_confidence_ellipse()`](https://alrobles.github.io/phylopred/reference/phylopred_confidence_ellipse.md)
  — parasite-cluster bootstrap and joint intercept-slope uncertainty.
- [`evaluate_phylopred_model()`](https://alrobles.github.io/phylopred/reference/evaluate_phylopred_model.md)
  — McFadden/Tjur pseudo-R2, AUC, and grouped cross-validated AUC.
- [`host_threshold_metric()`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)
  — presence-only threshold efficiency f_gamma = a / (gamma \* b + c)
  over the prediction threshold.
- [`phylo_partial_roc()`](https://alrobles.github.io/phylopred/reference/phylo_partial_roc.md)
  — phylogenetic partial ROC (sensitivity vs fraction of the phylogeny
  predicted suitable) with bootstrap ratio against a random-ranking
  null.
- [`fit_phylopred_pu()`](https://alrobles.github.io/phylopred/reference/fit_phylopred_pu.md)
  — positive-unlabeled host prediction with `xplus` on phylogenetic
  features.
- [`validate_host_prediction()`](https://alrobles.github.io/phylopred/reference/validate_host_prediction.md)
  — leave-hosts-out validation for both the GLM and PU routes.
- [`poulin_std()`](https://alrobles.github.io/phylopred/reference/poulin_std.md)
  — Poulin and Mouillot’s (2003) average taxonomic/ phylogenetic
  distinctness index S_TD and its sampling variance, for any distance
  matrix (taxonomic steps or phylogenetic distances).
- [`taxonomic_step_matrix()`](https://alrobles.github.io/phylopred/reference/taxonomic_step_matrix.md)
  — builds the classic 1-5 taxonomic step matrix from a
  genus/family/order/class hierarchy table.
- [`host_specificity_index()`](https://alrobles.github.io/phylopred/reference/host_specificity_index.md)
  — applies
  [`poulin_std()`](https://alrobles.github.io/phylopred/reference/poulin_std.md)
  to every parasite in an interaction table, for contrasting S_TD
  against
  [`host_threshold_metric()`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)’s
  model-based host range b\*.
- Data: `fish_tree` (Fish Tree of Life subtree, CC0),
  `beetleTreeInteractions`, `phy_dist`.
