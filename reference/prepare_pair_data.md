# Build a pairwise host-sharing data set from an interaction table

Expands an interaction table into all (focal host, target host) pairs
per parasite. For every parasite, each of its known hosts is used as a
focal host, and every other host in `phydist` is a target host whose
response is 1 if it is also a known host of that parasite and 0
otherwise.

## Usage

``` r
prepare_pair_data(interactions, phydist, group = NULL)
```

## Arguments

- interactions:

  A data.frame with at least two columns: the first the parasite (or
  symbiont) species, the second the host species. Extra columns are
  ignored. Duplicated rows are removed.

- phydist:

  A square, numeric phylogenetic (or taxonomic) distance matrix among
  hosts, with matching row and column names. Hosts in `interactions`
  absent from `phydist` are dropped with a warning.

- group:

  Optional. Either the name of a column in `interactions` giving a
  grouping factor for each parasite (e.g. genus or clade), or a function
  applied to the parasite name to derive the group (e.g.
  `function(x) sub(" .*", "", x)` for the genus).

## Value

A data.frame of class `"phylopred_pairs"` with columns:

- parasite:

  Parasite species (cluster identifier).

- focal:

  Focal (known) host.

- target:

  Target host.

- phydist:

  Distance between focal and target host.

- suscept:

  1 if the target host is a known host of the parasite.

- group:

  Grouping factor (only if `group` is supplied).

## Details

Unlike the legacy workflow (which samples a single focal host per
parasite), this uses every known host as a focal host, so no information
is discarded and results are deterministic. Rows belonging to the same
parasite are not independent; account for this with
[`cluster_bootstrap_phylopred`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md)
when computing uncertainty.

## See also

[`fit_phylopred_model`](https://alrobles.github.io/phylopred/reference/fit_phylopred_model.md),
[`cluster_bootstrap_phylopred`](https://alrobles.github.io/phylopred/reference/cluster_bootstrap_phylopred.md),
[`compare_phylopred_slopes`](https://alrobles.github.io/phylopred/reference/compare_phylopred_slopes.md)

## Examples

``` r
pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
#> Warning: 22 interaction record(s) dropped: host not in 'phydist'.
head(pairs)
#>                parasite    focal        target  phydist suscept
#> 1 Ambrosiodmus obliquus Lysiloma     Ailanthus 443.1818       0
#> 2 Ambrosiodmus obliquus Lysiloma Phellodendron 443.1818       0
#> 3 Ambrosiodmus obliquus Lysiloma        Citrus 443.1818       0
#> 4 Ambrosiodmus obliquus Lysiloma   Zanthoxylum 443.1818       0
#> 5 Ambrosiodmus obliquus Lysiloma       Murraya 443.1818       0
#> 6 Ambrosiodmus obliquus Lysiloma   Chloroxylon 443.1818       0
```
