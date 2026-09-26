# Cluster bootstrap for pairwise host-sharing models

Nonparametric bootstrap that resamples whole parasites (clusters) with
replacement and refits the logistic regression on each replicate. This
respects the non-independence of rows that share a parasite in a
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)
data set, unlike a naive row-level bootstrap. When a `group` column is
present, resampling is stratified within groups and the model includes a
group-by-distance interaction.

## Usage

``` r
cluster_bootstrap_phylopred(
  pairs,
  n_boot = 1000,
  seed = NULL,
  conf_level = 0.95
)
```

## Arguments

- pairs:

  A data.frame from
  [`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md).

- n_boot:

  Integer. Number of bootstrap replicates. Default 1000.

- seed:

  Integer or NULL. Random seed for reproducibility. Default NULL.

- conf_level:

  Confidence level for percentile intervals. Default 0.95.

## Value

An object of class `"phylopred_cluster_bootstrap"`, a list with:

- draws:

  Matrix of bootstrap coefficient draws (replicates x terms).

- estimates:

  Coefficients of the model fitted to the full data.

- confidence_intervals:

  Percentile confidence intervals per term.

- n_boot, seed, conf_level:

  Bootstrap settings.

## See also

[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md),
[`compare_phylopred_slopes`](https://alrobles.github.io/phylopred/reference/compare_phylopred_slopes.md)

## Examples

``` r
# Three parasites keep the example fast
top <- names(sort(table(beetleTreeInteractions[[1]]),
                  decreasing = TRUE))[1:3]
sub <- beetleTreeInteractions[
  as.character(beetleTreeInteractions[[1]]) %in% top, ]
pairs <- prepare_pair_data(sub, phy_dist)
#> Warning: 5 interaction record(s) dropped: host not in 'phydist'.
cluster_bootstrap_phylopred(pairs, n_boot = 5, seed = 42)
#> Phylopred cluster bootstrap (5 replicates)
#>             estimate    2.5%   97.5%
#> (Intercept)  -0.6507 -0.7935 -0.6153
#> phydist      -0.0007 -0.0013 -0.0003
```
