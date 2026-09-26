# Predict host susceptibility probabilities from logistic regression coefficients

Applies a logit function using intercept and slope coefficients to a
vector of phylogenetic distances. Accepts a numeric vector, a
`phylopred_model`, or a `phylopred_bootstrap` object.

## Usage

``` r
predict_phylopred_probability(coef, dist)
```

## Arguments

- coef:

  A named numeric vector with elements `"intercept"` and `"slope"`, OR a
  `phylopred_model` object, OR a `phylopred_bootstrap` object (in which
  case the bootstrap summary means are used). An unnamed vector of
  length \>= 2 is interpreted as `c(intercept, slope)`.

- dist:

  A numeric vector of phylogenetic distances.

## Value

A named numeric vector of predicted probabilities of the same length as
`dist`.

## See also

[`fit_phylopred_model`](https://alrobles.github.io/phylopred/reference/fit_phylopred_model.md)

## Examples

``` r
inc <- prepare_incidence_matrix(beetleTreeInteractions)
al  <- align_phylopred_inputs(inc, phy_dist)
model <- fit_phylopred_model(al$incidence, al$phydist, seed = 42)
dist_vec <- seq(0, max(al$phydist), length.out = 50)
predict_phylopred_probability(model, dist_vec)
#>  [1] 0.22449459 0.21885016 0.21330860 0.20787002 0.20253441 0.19730164
#>  [7] 0.19217149 0.18714363 0.18221764 0.17739302 0.17266917 0.16804541
#> [13] 0.16352100 0.15909511 0.15476684 0.15053525 0.14639932 0.14235798
#> [19] 0.13841011 0.13455454 0.13079007 0.12711546 0.12352941 0.12003061
#> [25] 0.11661773 0.11328940 0.11004422 0.10688080 0.10379771 0.10079352
#> [31] 0.09786678 0.09501605 0.09223986 0.08953677 0.08690530 0.08434401
#> [37] 0.08185144 0.07942614 0.07706667 0.07477160 0.07253950 0.07036897
#> [43] 0.06825861 0.06620703 0.06421286 0.06227475 0.06039137 0.05856138
#> [49] 0.05678350 0.05505644
```
