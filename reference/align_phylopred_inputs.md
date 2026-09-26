# Align incidence matrix and phylogenetic distance matrix

Validates both inputs and reorders the columns of `incidence` to match
the column order of `phydist`.

## Usage

``` r
align_phylopred_inputs(incidence, phydist)
```

## Arguments

- incidence:

  A binary incidence matrix (rows = parasites, cols = hosts).

- phydist:

  A square phylogenetic distance matrix among hosts.

## Value

A named list with:

- incidence:

  The incidence matrix with columns reordered to match `phydist`.

- phydist:

  The (unchanged) phylogenetic distance matrix.

## Examples

``` r
incidence <- prepare_incidence_matrix(beetleTreeInteractions)
aligned <- align_phylopred_inputs(incidence, phy_dist)
```
