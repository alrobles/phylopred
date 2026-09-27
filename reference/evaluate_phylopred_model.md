# Evaluate the predictive capacity of a pairwise host-sharing model

Computes goodness-of-fit and discrimination metrics for a logistic
host-sharing model fitted to a
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md)
data set: McFadden and Tjur pseudo-R-squared, the area under the ROC
curve (AUC), and optionally a grouped cross-validated AUC where whole
parasites are held out (parasite-grouped k-fold), which measures the
ability to predict the host range of unseen parasites.

## Usage

``` r
evaluate_phylopred_model(pairs, cv = FALSE, k = 5, seed = NULL)
```

## Arguments

- pairs:

  A data.frame from
  [`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md).

- cv:

  Logical. If TRUE, also compute the parasite-grouped cross-validated
  AUC. Default FALSE (refitting can be slow on large data sets).

- k:

  Integer. Number of cross-validation folds. Default 5.

- seed:

  Integer or NULL. Random seed for fold assignment. Default NULL.

## Value

An object of class `"phylopred_evaluation"`, a list with:

- mcfadden_r2:

  McFadden pseudo-R-squared, \\1 - logLik(model)/logLik(null)\\.

- tjur_r2:

  Tjur coefficient of discrimination (mean fitted probability for 1s
  minus mean for 0s).

- auc:

  In-sample area under the ROC curve.

- cv_auc:

  Parasite-grouped k-fold cross-validated AUC (NA if `cv = FALSE`).

- fitted_model:

  The underlying `glm` object.

- n_pairs, n_parasites:

  Data set size.

## See also

[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md),
[`cluster_bootstrap_phylopred`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md)

## Examples

``` r
# Three parasites keep the example fast
top <- names(sort(table(beetleTreeInteractions[[1]]),
                  decreasing = TRUE))[1:3]
sub <- beetleTreeInteractions[
  as.character(beetleTreeInteractions[[1]]) %in% top, ]
pairs <- prepare_pair_data(sub, phy_dist)
#> Warning: 5 interaction record(s) dropped: host not in 'phydist'.
evaluate_phylopred_model(pairs)
#> Phylopred model evaluation
#> Pairs: 44460  Parasites: 3
#> McFadden pseudo-R2: 0.0010
#> Tjur R2:            0.0011
#> AUC (in-sample):    0.5223
evaluate_phylopred_model(pairs, cv = TRUE, k = 3, seed = 42)
#> Phylopred model evaluation
#> Pairs: 44460  Parasites: 3
#> McFadden pseudo-R2: 0.0010
#> Tjur R2:            0.0011
#> AUC (in-sample):    0.5223
#> AUC (grouped CV):   0.4684
```
