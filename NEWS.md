# phylopred 0.1.0

Initial release. Clean reimplementation of the modern `geotax` workflow
for predicting species interactions from phylogenetic distance:

* `prepare_pair_data()`, `fit_phylopred_model()`,
  `compare_phylopred_slopes()` — pairwise host-sharing logistic models
  and formal slope comparison between clades.
* `bootstrap_phylopred_model()`, `cluster_bootstrap_phylopred()`,
  `phylopred_confidence_ellipse()` — parasite-cluster bootstrap and
  joint intercept-slope uncertainty.
* `evaluate_phylopred_model()` — McFadden/Tjur pseudo-R2, AUC, and
  grouped cross-validated AUC.
* `host_threshold_metric()` — presence-only threshold efficiency
  f_gamma = a / (gamma * b + c) over the prediction threshold.
* `phylo_partial_roc()` — phylogenetic partial ROC (sensitivity vs
  fraction of the phylogeny predicted suitable) with bootstrap
  ratio against a random-ranking null.
* `fit_phylopred_pu()` — positive-unlabeled host prediction with
  `xplus` on phylogenetic features.
* `validate_host_prediction()` — leave-hosts-out validation for both
  the GLM and PU routes.
* `poulin_std()` — Poulin and Mouillot's (2003) average taxonomic/
  phylogenetic distinctness index S_TD and its sampling variance, for
  any distance matrix (taxonomic steps or phylogenetic distances).
* `taxonomic_step_matrix()` — builds the classic 1-5 taxonomic step
  matrix from a genus/family/order/class hierarchy table.
* `host_specificity_index()` — applies `poulin_std()` to every
  parasite in an interaction table, for contrasting S_TD against
  `host_threshold_metric()`'s model-based host range b*.
* Data: `fish_tree` (Fish Tree of Life subtree, CC0),
  `beetleTreeInteractions`, `phy_dist`.
