# Bootstrap a phylopred logistic regression model

Runs the phylopred logistic regression `n` times using random focal host
sampling, returning a distribution of coefficient estimates.

## Usage

``` r
bootstrap_phylopred_model(incidence, phydist, n = 1000, seed = NULL)
```

## Arguments

- incidence:

  A binary incidence matrix (rows = parasites, cols = hosts).

- phydist:

  A square, numeric phylogenetic distance matrix among hosts.

- n:

  Integer. Number of bootstrap iterations. Default 1000.

- seed:

  Integer or NULL. Random seed for reproducibility. Default NULL.

## Value

An S3 object of class `"phylopred_bootstrap"`, a list with:

- call:

  The matched call.

- n:

  Number of bootstrap iterations.

- seed:

  The seed used.

- coefficients:

  A matrix of `n` rows x 12 columns (coefficient info from each
  iteration).

- summary:

  A data.frame with the column-wise means (same 12 columns).

- metadata:

  List with `n_parasites` and `n_hosts`.

## See also

[`fit_phylopred_model`](https://alrobles.github.io/phylopred/reference/fit_phylopred_model.md),
[`predict_phylopred_probability`](https://alrobles.github.io/phylopred/reference/predict_phylopred_probability.md)

## Examples

``` r
incidence <- prepare_incidence_matrix(beetleTreeInteractions)
aligned   <- align_phylopred_inputs(incidence, phy_dist)
bootstrap_phylopred_model(aligned$incidence, aligned$phydist, n = 5, seed = 42)
#> Phylopred bootstrap (5 iterations)
#> Intercept (mean): -1.6305
#> Slope     (mean): -0.0016
```
