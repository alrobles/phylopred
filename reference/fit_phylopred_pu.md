# Positive-unlabeled host prediction via xplus (PLUS algorithm)

Predicts the hosts of a parasite from the phylogeny using
positive-unlabeled (PU) learning, as a parallel alternative to the
pairwise logistic model. Known hosts are treated as positives and every
other species in the phylogeny as unlabeled (never as a true absence).
Each species is described by phylogenetic features and the PLUS
algorithm
([`xplus::xplus`](https://alrobles.github.io/xplus/reference/xplus.html))
iteratively re-labels the unlabeled species while fitting a penalized
logistic model.

## Usage

``` r
fit_phylopred_pu(
  interactions,
  phydist,
  parasite,
  n_eigen = 10,
  seed = NULL,
  ...
)
```

## Arguments

- interactions:

  A data.frame whose first two columns are the parasite and host species
  names (as in
  [`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)).

- phydist:

  Square numeric matrix of phylogenetic distances with species names as
  row and column names.

- parasite:

  Character scalar: the parasite whose hosts are predicted.

- n_eigen:

  Number of PCoA axes used as features. Default 10.

- seed:

  Optional random seed passed to
  [`xplus::xplus`](https://alrobles.github.io/xplus/reference/xplus.html).

- ...:

  Further arguments passed to
  [`xplus::xplus`](https://alrobles.github.io/xplus/reference/xplus.html)
  (e.g. `max_iter`, `learning_rate`, `qq`).

## Value

An object of class `"phylopred_pu"`: a list with

- prediction:

  data.frame with `species`, `pred` (predicted host suitability) and
  `is_host`.

- fit:

  The underlying `xplus` object.

- parasite:

  The parasite name.

The `prediction` columns can be passed directly to
[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)
and
[`phylo_partial_roc`](https://alrobles.github.io/phylopred/reference/phylo_partial_roc.md).

## Details

The feature matrix contains, for every species in `phydist`:

- the minimum and mean phylogenetic distance to the parasite's known
  hosts (excluding the species itself), and

- the first `n_eigen` principal coordinates (PCoA) of the phylogenetic
  distance matrix, describing the global position of the species in the
  phylogeny (the analogue of environmental layers in an SDM).

## See also

[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md),
[`phylo_partial_roc`](https://alrobles.github.io/phylopred/reference/phylo_partial_roc.md),
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)

## Examples

``` r
pu <- fit_phylopred_pu(beetleTreeInteractions, phy_dist,
                     parasite = beetleTreeInteractions[[1]][1],
                     max_iter = 20, seed = 1)
host_threshold_metric(pu$prediction$pred, pu$prediction$is_host)
#> Phylopred threshold metric f_gamma (gamma = 1)
#> Optimal threshold: 0.6020  f_gamma: 39.1667
#> Evaluated at 190 thresholds.
```
