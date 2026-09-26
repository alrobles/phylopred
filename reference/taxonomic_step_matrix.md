# Build a taxonomic step-distance matrix (Poulin & Mouillot scale)

Constructs a symmetric taxonomic distance matrix among species from a
taxonomic hierarchy table, following the step scale of Poulin and
Mouillot (2003): 1 if two species share a genus, 2 if they share a
family but not a genus, 3 if they share an order but not a family, 4 if
they share a class but not an order, and 5 if they share none of the
above (different class / phylum). Any subset of `genus`, `family`,
`order`, `class` may be supplied; missing ranks are skipped and the
maximum step for that pair falls back to `length(ranks) + 1`.

## Usage

``` r
taxonomic_step_matrix(tax_table)
```

## Arguments

- tax_table:

  A data.frame with one row per species and columns `species` plus at
  least two of `genus`, `family`, `order`, `class`. No NA values are
  allowed in the rank columns (see
  [`host_specificity_index`](https://alrobles.github.io/phylopred/reference/host_specificity_index.md)
  for handling hosts with incomplete taxonomy).

## Value

A square numeric matrix with row/column names equal to
`tax_table$species`, diagonal 0, suitable as the `phydist` argument of
[`poulin_std`](https://alrobles.github.io/phylopred/reference/poulin_std.md)
or
[`host_specificity_index`](https://alrobles.github.io/phylopred/reference/host_specificity_index.md).

## References

Poulin, R. and Mouillot, D. (2003) Parasite specialization from a
phylogenetic perspective: a new index of host specificity.
*Parasitology* 126, 473-480.
[doi:10.1017/S0031182003002993](https://doi.org/10.1017/S0031182003002993)

## See also

[`poulin_std`](https://alrobles.github.io/phylopred/reference/poulin_std.md),
[`host_specificity_index`](https://alrobles.github.io/phylopred/reference/host_specificity_index.md)

## Examples

``` r
tax <- data.frame(
  species = c("Caligus elongatus", "Caligus curtus",
              "Lepeophtheirus salmonis", "Ergasilus sieboldi"),
  genus  = c("Caligus", "Caligus", "Lepeophtheirus", "Ergasilus"),
  family = c("Caligidae", "Caligidae", "Caligidae", "Ergasilidae"),
  order  = c("Siphonostomatoida", "Siphonostomatoida",
              "Siphonostomatoida", "Cyclopoida")
)
taxonomic_step_matrix(tax)
#>                         Caligus elongatus Caligus curtus
#> Caligus elongatus                       0              1
#> Caligus curtus                          1              0
#> Lepeophtheirus salmonis                 2              2
#> Ergasilus sieboldi                      4              4
#>                         Lepeophtheirus salmonis Ergasilus sieboldi
#> Caligus elongatus                             2                  4
#> Caligus curtus                                2                  4
#> Lepeophtheirus salmonis                       0                  4
#> Ergasilus sieboldi                            4                  0
```
