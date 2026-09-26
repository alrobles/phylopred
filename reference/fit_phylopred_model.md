# Fit a phylopred logistic regression model

Fits a logistic regression model relating host susceptibility to
phylogenetic distance.

## Usage

``` r
fit_phylopred_model(incidence, phydist, seed = NULL)
```

## Arguments

- incidence:

  A binary incidence matrix (rows = parasites/pathogens, cols = hosts).

- phydist:

  A square, numeric phylogenetic distance matrix among hosts. Row/col
  names must match the column names of `incidence`.

- seed:

  Integer or NULL. Random seed for reproducibility. Default NULL.

## Value

An S3 object of class `"phylopred_model"`, a list with:

- call:

  The matched call.

- coefficients:

  Named numeric vector `c(intercept, slope)`.

- standard_errors:

  Named numeric vector `c(intercept, slope)`.

- confidence_intervals:

  Matrix with rows `c("intercept","slope")` and cols
  `c("2.5%","97.5%")`.

- fitted_model:

  The underlying `glm` object.

- training_data:

  The data.frame used for model fitting.

- convergence:

  Logical; whether the GLM converged.

- metadata:

  List with `n_parasites`, `n_hosts`, and `seed`.

## See also

[`bootstrap_phylopred_model`](https://alrobles.github.io/phylopred/reference/bootstrap_phylopred_model.md),
[`predict_phylopred_probability`](https://alrobles.github.io/phylopred/reference/predict_phylopred_probability.md)

## Examples

``` r
incidence <- prepare_incidence_matrix(beetleTreeInteractions)
aligned   <- align_phylopred_inputs(incidence, phy_dist)
fit_phylopred_model(aligned$incidence, aligned$phydist, seed = 42)
#> Phylopred model
#> Intercept: -1.2397 (SE: 0.1791)
#> Slope:     -0.0025 (SE: 0.0004)
#> Convergence: TRUE
```
