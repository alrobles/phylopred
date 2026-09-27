# Host specificity index (S_TD) for every parasite in an interaction table

Applies
[`poulin_std`](https://alrobles.github.io/phylopred/reference/poulin_std.md)
to every parasite in an interaction table, following the same expansion
logic as
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md).
Useful for contrasting the Poulin and Mouillot (2003) specificity index
against `phylopred`'s model-based host-range metric `b*`
([`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md))
across many parasites.

## Usage

``` r
host_specificity_index(interactions, phydist, group = NULL)
```

## Arguments

- interactions:

  A data.frame with at least two columns: the first the parasite (or
  symbiont) species, the second the host species. Extra columns are
  ignored. Duplicated rows are removed.

- phydist:

  A square, numeric distance matrix among hosts (taxonomic steps or
  phylogenetic distances), with matching row and column names.

- group:

  Optional. Either the name of a column in `interactions` giving a
  grouping factor for each parasite (e.g. genus or clade), or a function
  applied to the parasite name to derive the group.

## Value

A data.frame of class `"phylopred_specificity"` with columns `parasite`,
`n_hosts` (known hosts present in `phydist`), `std` (\\S\_{TD}\\),
`var_std` (\\Var(S\_{TD})\\), and `group` (only if `group` is supplied).

## References

Poulin, R. and Mouillot, D. (2003) Parasite specialization from a
phylogenetic perspective: a new index of host specificity.
*Parasitology* 126, 473-480.
[doi:10.1017/S0031182003002993](https://doi.org/10.1017/S0031182003002993)

## See also

[`poulin_std`](https://alrobles.github.io/phylopred/reference/poulin_std.md),
[`prepare_pair_data`](https://alrobles.github.io/phylopred/reference/prepare_pair_data.md),
[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)

## Examples

``` r
std_by_parasite <- host_specificity_index(beetleTreeInteractions, phy_dist)
head(std_by_parasite)
#>                    parasite n_hosts      std  var_std
#> 1     Ambrosiodmus obliquus       6 477.4053 11989.14
#> 2  Ambrosiodmus rubricollis       7 477.6867 22962.69
#> 3 Ambrosiodmus tachygraphus       4 467.0644 24030.53
#> 4   Coptoborus pseudotenuis      10 394.9607  5115.38
#> 5     Euwallacea fornicatus      31 457.2600 13711.36
#> 6        Euwallacea validus      22 511.3490 22369.51
```
