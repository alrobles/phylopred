# Compare host-sharing slopes between two clades

Fits a single logistic regression with a group-by-phylogenetic-distance
interaction to a pairwise data set built with
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)
and tests whether the slopes of the two groups differ. A more negative
slope means the probability of sharing a parasite decays faster with
phylogenetic distance, i.e. the clade is more of a phylogenetic
specialist. Uncertainty for the slope difference comes from a
parasite-level cluster bootstrap.

## Usage

``` r
compare_phylopred_slopes(pairs, n_boot = 1000, seed = NULL, conf_level = 0.95)
```

## Arguments

- pairs:

  A data.frame from
  [`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)
  with a `group` column containing exactly two groups.

- n_boot:

  Integer. Number of cluster bootstrap replicates. Default 1000.

- seed:

  Integer or NULL. Random seed for reproducibility. Default NULL.

- conf_level:

  Confidence level. Default 0.95.

## Value

An object of class `"phylopred_slope_comparison"`, a list with:

- groups:

  The two group labels (reference first).

- slopes:

  Named vector with the slope of each group.

- slope_difference:

  Slope of the second group minus the reference.

- interaction_p_value:

  Wald p-value of the interaction term.

- bootstrap:

  The
  [`cluster_bootstrap_phylopred`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md)
  object.

- slope_cis:

  Bootstrap percentile CIs for each group's slope and for the
  difference.

- fitted_model:

  The underlying `glm` object.

## See also

[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md),
[`cluster_bootstrap_phylopred`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md)

## Examples

``` r
# Compare Xyleborus beetles against other genera (subset kept small
# so the bootstrap example runs quickly)
parasites <- unique(as.character(beetleTreeInteractions[[1]]))
keep <- c(grep("^Xyleborus", parasites, value = TRUE)[1:3],
          parasites[!grepl("^Xyleborus", parasites)][1:3])
sub <- beetleTreeInteractions[
  as.character(beetleTreeInteractions[[1]]) %in% keep, ]
pairs <- prepare_pair_data(
  sub, phy_dist,
  group = function(x) ifelse(grepl("^Xyleborus", x), "Xyleborus", "other")
)
#> Warning: 4 interaction record(s) dropped: host not in 'phydist'.
compare_phylopred_slopes(pairs, n_boot = 5, seed = 42)
#> Phylopred slope comparison: Xyleborus vs other 
#> Slope Xyleborus: -0.0005  [-0.0006, -0.0003]
#> Slope other: -0.0010  [-0.0015, -0.0002]
#> Difference: -0.0006  [-0.0009, 0.0004]  (Wald p = 0.553)
```
