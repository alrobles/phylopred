# Poulin & Mouillot's taxonomic/phylogenetic distinctness index (S_TD)

Computes the average taxonomic (or phylogenetic) distinctness
\\S\_{TD}\\ of Poulin and Mouillot (2003) for a single parasite's known
host set: \$\$S\_{TD} = \frac{2 \sum\_{i\<j} \omega\_{ij}}{s(s-1)}\$\$
where \\\omega\_{ij}\\ is the distance between hosts \\i\\ and \\j\\
taken from `phydist`, and \\s\\ is the number of known hosts. When
`phydist` is a taxonomic step matrix (see
[`taxonomic_step_matrix`](https://alrobles.github.io/phylopred/reference/taxonomic_step_matrix.md)),
\\S\_{TD}\\ is the original 1-5 taxonomic index; when `phydist` is a
matrix of phylogenetic (patristic) distances, this is the continuous
phylogenetic generalization noted as possible by Poulin and Mouillot
(2003, after Clarke and Warwick, 2001), expressed directly in
branch-length units. The corresponding sampling variance \\Var(S\_{TD})
= \sum\_{i \ne j} (\omega\_{ij} - S\_{TD})^2 / \[s(s-1)\]\\ is also
returned.

## Usage

``` r
poulin_std(phydist, hosts)
```

## Arguments

- phydist:

  A square, numeric distance matrix among hosts with matching row and
  column names (taxonomic steps from
  [`taxonomic_step_matrix`](https://alrobles.github.io/phylopred/reference/taxonomic_step_matrix.md),
  or phylogenetic/patristic distances such as `phy_dist`).

- hosts:

  Character vector of known host names for one parasite. Hosts not
  present in `phydist` are dropped with a warning.

## Value

An object of class `"phylopred_std"`: a list with elements `n_hosts`,
`std` (\\S\_{TD}\\), and `var_std` (\\Var(S\_{TD})\\, `NA` when
`n_hosts < 3`). \\S\_{TD}\\ is undefined (`NA`, with a warning) when
`n_hosts < 2`, as noted by Poulin and Mouillot (2003) for single-host
parasites.

## References

Poulin, R. and Mouillot, D. (2003) Parasite specialization from a
phylogenetic perspective: a new index of host specificity.
*Parasitology* 126, 473-480.
[doi:10.1017/S0031182003002993](https://doi.org/10.1017/S0031182003002993)

Clarke, K.R. and Warwick, R.M. (2001) A further biodiversity index
applicable to species lists: variation in taxonomic distinctness.
*Marine Ecology Progress Series* 216, 265-278.

## See also

[`host_specificity_index`](https://alrobles.github.io/phylopred/reference/host_specificity_index.md),
[`taxonomic_step_matrix`](https://alrobles.github.io/phylopred/reference/taxonomic_step_matrix.md),
[`host_threshold_metric`](https://alrobles.github.io/phylopred/reference/host_threshold_metric.md)

## Examples

``` r
hosts <- unique(beetleTreeInteractions$genus[
  beetleTreeInteractions$beetle == beetleTreeInteractions$beetle[1]])
poulin_std(phy_dist, hosts)
#> Poulin & Mouillot S_TD
#>   n_hosts: 6
#>   S_TD:    477.4053
#>   Var(S_TD): 11989.1383
```
