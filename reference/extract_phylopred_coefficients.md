# Extract coefficients from a phylopred model or bootstrap object

Returns a named list of model coefficients and their associated
statistics from a `phylopred_model` or `phylopred_bootstrap` object.

## Usage

``` r
extract_phylopred_coefficients(model)
```

## Arguments

- model:

  A `phylopred_model` or `phylopred_bootstrap` object.

## Value

A named list with elements:

- intercept:

  Intercept coefficient (mean for bootstrap).

- slope:

  Slope coefficient (mean for bootstrap).

- intercept_se:

  Standard error of intercept.

- slope_se:

  Standard error of slope.

- intercept_ci:

  Two-element vector with 2.5% and 97.5% CI for intercept.

- slope_ci:

  Two-element vector with 2.5% and 97.5% CI for slope.

## Examples

``` r
inc <- prepare_incidence_matrix(beetleTreeInteractions)
al  <- align_phylopred_inputs(inc, phy_dist)
model <- fit_phylopred_model(al$incidence, al$phydist, seed = 42)
extract_phylopred_coefficients(model)
#> $intercept
#> [1] -1.239663
#> 
#> $slope
#> [1] -0.002466312
#> 
#> $intercept_se
#> [1] 0.1790548
#> 
#> $slope_se
#> [1] 0.000381976
#> 
#> $intercept_ci
#>       2.5%      97.5% 
#> -1.5906042 -0.8887224 
#> 
#> $slope_ci
#>         2.5%        97.5% 
#> -0.003214972 -0.001717653 
#> 
```
