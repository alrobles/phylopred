# Leave-hosts-out validation of host predictions

Presence-only predictive validation: a fraction of the parasite's known
hosts is hidden, the model is refitted on the remaining hosts, and the
hidden hosts are scored against the unlabeled species. Because true
absences do not exist, performance is measured as the AUC of the hidden
(held-out) hosts versus the unlabeled species — the probability that a
truly interacting species outranks a random unlabeled one. This measures
actual predictive ability, not in-sample fit.

## Usage

``` r
validate_host_prediction(
  interactions,
  phydist,
  parasite,
  method = c("glm", "pu"),
  prop_holdout = 0.25,
  n_rep = 10,
  seed = NULL,
  ...
)
```

## Arguments

- interactions:

  A data.frame whose first two columns are the parasite and host species
  names.

- phydist:

  Square numeric matrix of phylogenetic distances with species names as
  row and column names.

- parasite:

  Character scalar: the parasite whose hosts are predicted.

- method:

  `"glm"` for the pairwise logistic model (target species scored by
  their maximum predicted sharing probability across the training focal
  hosts) or `"pu"` for the xplus PU-learning model
  ([`fit_phylopred_pu`](https://alrobles.github.io/phylopred/reference/fit_phylopred_pu.md)).

- prop_holdout:

  Proportion of known hosts hidden in each replicate. Default 0.25.

- n_rep:

  Number of holdout replicates. Default 10.

- seed:

  Optional random seed.

- ...:

  Further arguments passed to
  [`fit_phylopred_pu`](https://alrobles.github.io/phylopred/reference/fit_phylopred_pu.md)
  when `method = "pu"`.

## Value

An object of class `"phylopred_holdout"`: a list with `auc` (vector of
holdout AUCs, one per replicate), `mean_auc`, `method`, `parasite`,
`prop_holdout`, and `n_rep`.

## See also

[`fit_phylopred_pu`](https://alrobles.github.io/phylopred/reference/fit_phylopred_pu.md),
[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md),
[`phylo_partial_roc`](https://alrobles.github.io/phylopred/reference/phylo_partial_roc.md)

## Examples

``` r
parasite <- names(sort(table(beetleTreeInteractions[[1]]),
                       decreasing = TRUE))[1]
validate_host_prediction(beetleTreeInteractions, phy_dist, parasite,
                         method = "glm", n_rep = 3, seed = 1)
#> Phylopred leave-hosts-out validation (glm)
#> Parasite: Xyleborus affinis
#> Holdout AUC (hidden hosts vs unlabeled): mean 0.563 over 3 reps [0.519, 0.595]
```
