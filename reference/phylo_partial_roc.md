# Phylogenetic partial ROC (Peterson, Papes & Soberon 2008 analogue)

Partial ROC test for presence-only host-suitability predictions,
following Peterson, Papes & Soberon (2008, Ecological Modelling 213,
63-72). The ROC space is redefined for presence-only data: the y axis is
the sensitivity \\a(t)\\ (proportion of known hosts predicted at cutoff
\\t\\) and the x axis is the proportion of ALL species in the phylogeny
predicted as hosts, \\b(t)\\ (the analogue of proportional predicted
area). The curve is restricted to the low-omission region \\a(t) \ge 1 -
E\\ for a user-chosen omission tolerance `error_rate` \\E\\, and
performance is summarized as the ratio of the partial AUC of the
observed curve to the partial AUC of the null (random-ranking)
expectation, the diagonal \\a = b\\. A ratio of 1 equals random
performance; ratios above 1 indicate predictive signal. Significance is
assessed by bootstrap resampling of the known hosts.

## Usage

``` r
phylo_partial_roc(
  pred,
  is_host,
  error_rate = 0.05,
  n_boot = 500,
  boot_prop = 0.5,
  seed = NULL
)
```

## Arguments

- pred:

  Numeric vector of predicted suitability scores, one per species in the
  phylogeny.

- is_host:

  Logical or 0/1 vector of the same length: TRUE/1 for known hosts.

- error_rate:

  Allowed omission error \\E\\ in \\\[0, 1)\\. The partial AUC is
  computed over the region with sensitivity \\\ge 1 - E\\. Default 0.05.

- n_boot:

  Number of bootstrap replicates. Default 500.

- boot_prop:

  Proportion of known hosts resampled (with replacement) in each
  bootstrap replicate, as in Peterson et al. Default 0.5.

- seed:

  Integer or NULL. Random seed for the bootstrap. Default NULL.

## Value

An object of class `"phylopred_partial_roc"`, a list with:

- auc_ratio:

  Observed partial AUC ratio (observed / null).

- p_value:

  Bootstrap proportion of replicates with ratio \\\le 1\\.

- boot_ratios:

  Vector of bootstrap AUC ratios.

- curve:

  data.frame with `threshold`, `a` (sensitivity) and `b` (predicted
  breadth).

- error_rate, n_boot, boot_prop:

  The configuration used.

## See also

[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)

## Examples

``` r
pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
#> Warning: 22 interaction record(s) dropped: host not in 'phydist'.
fit <- stats::glm(suscept ~ phydist, data = pairs, family = stats::binomial())
one <- pairs[pairs$parasite == pairs$parasite[1], ]
one <- one[!duplicated(one$target), ]
pred <- stats::predict(fit, newdata = one, type = "response")
phylo_partial_roc(pred, one$suscept, n_boot = 50, seed = 42)
#> Phylopred phylogenetic partial ROC (E = 0.05)
#> Partial AUC ratio: 1.0284
#> Bootstrap p-value (ratio <= 1): 0.0000  [50 replicates]
```
