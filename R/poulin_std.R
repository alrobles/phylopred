#' Build a taxonomic step-distance matrix (Poulin & Mouillot scale)
#'
#' Constructs a symmetric taxonomic distance matrix among species from a
#' taxonomic hierarchy table, following the step scale of
#' Poulin and Mouillot (2003): 1 if two species share a genus, 2 if they
#' share a family but not a genus, 3 if they share an order but not a
#' family, 4 if they share a class but not an order, and 5 if they share
#' none of the above (different class / phylum). Any subset of
#' \code{genus}, \code{family}, \code{order}, \code{class} may be supplied;
#' missing ranks are skipped and the maximum step for that pair falls back
#' to \code{length(ranks) + 1}.
#'
#' @param tax_table A data.frame with one row per species and columns
#'   \code{species} plus at least two of \code{genus}, \code{family},
#'   \code{order}, \code{class}. No NA values are allowed in the rank
#'   columns (see \code{\link{host_specificity_index}} for handling hosts
#'   with incomplete taxonomy).
#'
#' @return A square numeric matrix with row/column names equal to
#'   \code{tax_table$species}, diagonal 0, suitable as the \code{phydist}
#'   argument of \code{\link{poulin_std}} or \code{\link{host_specificity_index}}.
#'
#' @references Poulin, R. and Mouillot, D. (2003) Parasite specialization
#'   from a phylogenetic perspective: a new index of host specificity.
#'   \emph{Parasitology} 126, 473-480. \doi{10.1017/S0031182003002993}
#'
#' @seealso \code{\link{poulin_std}}, \code{\link{host_specificity_index}}
#'
#' @export
#'
#' @examples
#' tax <- data.frame(
#'   species = c("Caligus elongatus", "Caligus curtus",
#'               "Lepeophtheirus salmonis", "Ergasilus sieboldi"),
#'   genus  = c("Caligus", "Caligus", "Lepeophtheirus", "Ergasilus"),
#'   family = c("Caligidae", "Caligidae", "Caligidae", "Ergasilidae"),
#'   order  = c("Siphonostomatoida", "Siphonostomatoida",
#'               "Siphonostomatoida", "Cyclopoida")
#' )
#' taxonomic_step_matrix(tax)
taxonomic_step_matrix <- function(tax_table) {
  .validate_tax_table(tax_table)
  if (!"species" %in% colnames(tax_table)) {
    stop("'tax_table' must include a 'species' column of species names.")
  }
  sp <- as.character(tax_table$species)
  if (anyDuplicated(sp)) {
    stop("'tax_table' must have one row per unique species.")
  }
  ranks <- c("genus", "family", "order", "class")
  ranks <- ranks[ranks %in% colnames(tax_table)]
  if (length(ranks) < 2) {
    stop("'tax_table' must supply at least two of: genus, family, order, class.")
  }

  n <- length(sp)
  D <- matrix(0, n, n, dimnames = list(sp, sp))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      step <- .taxonomic_step(tax_table[i, ranks, drop = FALSE],
                               tax_table[j, ranks, drop = FALSE], ranks)
      D[i, j] <- step
      D[j, i] <- step
    }
  }
  D
}

# Step distance between two taxonomy rows on the Poulin & Mouillot scale.
.taxonomic_step <- function(row_i, row_j, ranks) {
  for (k in seq_along(ranks)) {
    if (identical(row_i[[ranks[k]]], row_j[[ranks[k]]])) {
      return(k)
    }
  }
  length(ranks) + 1
}

#' Poulin & Mouillot's taxonomic/phylogenetic distinctness index (S_TD)
#'
#' Computes the average taxonomic (or phylogenetic) distinctness
#' \eqn{S_{TD}} of Poulin and Mouillot (2003) for a single parasite's known
#' host set:
#' \deqn{S_{TD} = \frac{2 \sum_{i<j} \omega_{ij}}{s(s-1)}}
#' where \eqn{\omega_{ij}} is the distance between hosts \eqn{i} and
#' \eqn{j} taken from \code{phydist}, and \eqn{s} is the number of known
#' hosts. When \code{phydist} is a taxonomic step matrix (see
#' \code{\link{taxonomic_step_matrix}}), \eqn{S_{TD}} is the original
#' 1-5 taxonomic index; when \code{phydist} is a matrix of phylogenetic
#' (patristic) distances, this is the continuous phylogenetic
#' generalization noted as possible by Poulin and Mouillot (2003, after
#' Clarke and Warwick, 2001), expressed directly in branch-length units.
#' The corresponding sampling variance
#' \eqn{Var(S_{TD}) = \sum_{i \ne j} (\omega_{ij} - S_{TD})^2 / [s(s-1)]}
#' is also returned.
#'
#' @param phydist A square, numeric distance matrix among hosts with
#'   matching row and column names (taxonomic steps from
#'   \code{\link{taxonomic_step_matrix}}, or phylogenetic/patristic
#'   distances such as \code{phy_dist}).
#' @param hosts Character vector of known host names for one parasite.
#'   Hosts not present in \code{phydist} are dropped with a warning.
#'
#' @return An object of class \code{"phylopred_std"}: a list with elements
#'   \code{n_hosts}, \code{std} (\eqn{S_{TD}}), and \code{var_std}
#'   (\eqn{Var(S_{TD})}, \code{NA} when \code{n_hosts < 3}). \eqn{S_{TD}}
#'   is undefined (\code{NA}, with a warning) when \code{n_hosts < 2}, as
#'   noted by Poulin and Mouillot (2003) for single-host parasites.
#'
#' @references Poulin, R. and Mouillot, D. (2003) Parasite specialization
#'   from a phylogenetic perspective: a new index of host specificity.
#'   \emph{Parasitology} 126, 473-480. \doi{10.1017/S0031182003002993}
#'
#'   Clarke, K.R. and Warwick, R.M. (2001) A further biodiversity index
#'   applicable to species lists: variation in taxonomic distinctness.
#'   \emph{Marine Ecology Progress Series} 216, 265-278.
#'
#' @seealso \code{\link{host_specificity_index}},
#'   \code{\link{taxonomic_step_matrix}}, \code{\link{host_threshold_metric}}
#'
#' @export
#'
#' @examples
#' hosts <- unique(beetleTreeInteractions$genus[
#'   beetleTreeInteractions$beetle == beetleTreeInteractions$beetle[1]])
#' poulin_std(phy_dist, hosts)
poulin_std <- function(phydist, hosts) {
  .validate_phydist(phydist)
  hosts <- unique(as.character(hosts))
  in_mat <- hosts %in% rownames(phydist)
  if (!all(in_mat)) {
    warning(sprintf("%d host(s) not found in 'phydist' and were dropped.",
                    sum(!in_mat)))
  }
  hosts <- hosts[in_mat]
  s <- length(hosts)
  if (s == 0) {
    stop("None of 'hosts' are present in 'phydist'.")
  }
  if (s == 1) {
    warning("S_TD is undefined for a single-host parasite; returning NA.")
    return(structure(list(n_hosts = 1L, std = NA_real_, var_std = NA_real_),
                     class = "phylopred_std"))
  }

  sub <- phydist[hosts, hosts, drop = FALSE]
  upper <- sub[upper.tri(sub)]
  std <- mean(upper)
  var_std <- if (s >= 3) sum((upper - std)^2) * 2 / (s * (s - 1)) else NA_real_

  structure(list(n_hosts = s, std = std, var_std = var_std),
            class = "phylopred_std")
}

#' @export
print.phylopred_std <- function(x, ...) {
  cat("Poulin & Mouillot S_TD\n")
  cat(sprintf("  n_hosts: %d\n", x$n_hosts))
  cat(sprintf("  S_TD:    %s\n",
              if (is.na(x$std)) "NA" else sprintf("%.4f", x$std)))
  cat(sprintf("  Var(S_TD): %s\n",
              if (is.na(x$var_std)) "NA" else sprintf("%.4f", x$var_std)))
  invisible(x)
}

#' Host specificity index (S_TD) for every parasite in an interaction table
#'
#' Applies \code{\link{poulin_std}} to every parasite in an interaction
#' table, following the same expansion logic as
#' \code{\link{prepare_pair_data}}. Useful for contrasting Poulin and
#' Mouillot's (2003) specificity index against \code{phylopred}'s
#' model-based host-range metric \code{b*}
#' (\code{\link{host_threshold_metric}}) across many parasites.
#'
#' @param interactions A data.frame with at least two columns: the first
#'   the parasite (or symbiont) species, the second the host species.
#'   Extra columns are ignored. Duplicated rows are removed.
#' @param phydist A square, numeric distance matrix among hosts (taxonomic
#'   steps or phylogenetic distances), with matching row and column names.
#' @param group Optional. Either the name of a column in \code{interactions}
#'   giving a grouping factor for each parasite (e.g. genus or clade), or a
#'   function applied to the parasite name to derive the group.
#'
#' @return A data.frame of class \code{"phylopred_specificity"} with columns
#'   \code{parasite}, \code{n_hosts} (known hosts present in \code{phydist}),
#'   \code{std} (\eqn{S_{TD}}), \code{var_std} (\eqn{Var(S_{TD})}), and
#'   \code{group} (only if \code{group} is supplied).
#'
#' @references Poulin, R. and Mouillot, D. (2003) Parasite specialization
#'   from a phylogenetic perspective: a new index of host specificity.
#'   \emph{Parasitology} 126, 473-480. \doi{10.1017/S0031182003002993}
#'
#' @seealso \code{\link{poulin_std}}, \code{\link{prepare_pair_data}},
#'   \code{\link{host_threshold_metric}}
#'
#' @export
#'
#' @examples
#' std_by_parasite <- host_specificity_index(beetleTreeInteractions, phy_dist)
#' head(std_by_parasite)
host_specificity_index <- function(interactions, phydist, group = NULL) {
  if (!inherits(interactions, "data.frame")) {
    stop("'interactions' must be a data.frame.")
  }
  if (ncol(interactions) < 2) {
    stop("'interactions' must have at least two columns (parasite, host).")
  }
  .validate_phydist(phydist)

  group_col <- NULL
  if (!is.null(group)) {
    if (is.character(group) && length(group) == 1) {
      if (!group %in% colnames(interactions)) {
        stop(sprintf("Column '%s' not found in 'interactions'.", group))
      }
      group_col <- group
    } else if (!is.function(group)) {
      stop("'group' must be a column name or a function of the parasite name.")
    }
  }

  parasite <- as.character(interactions[[1]])
  host <- as.character(interactions[[2]])
  df0 <- unique(data.frame(parasite = parasite, host = host,
                           stringsAsFactors = FALSE))
  hosts_by_parasite <- split(df0$host, df0$parasite)

  out <- lapply(names(hosts_by_parasite), function(p) {
    hs <- hosts_by_parasite[[p]]
    hs <- hs[hs %in% rownames(phydist)]
    if (length(hs) < 2) {
      return(data.frame(parasite = p, n_hosts = length(hs),
                         std = NA_real_, var_std = NA_real_,
                         stringsAsFactors = FALSE))
    }
    m <- suppressWarnings(poulin_std(phydist, hs))
    data.frame(parasite = p, n_hosts = m$n_hosts, std = m$std,
               var_std = m$var_std, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, out)
  rownames(result) <- NULL

  if (!is.null(group_col)) {
    map <- unique(data.frame(parasite = parasite,
                             group = as.character(interactions[[group_col]]),
                             stringsAsFactors = FALSE))
    if (anyDuplicated(map$parasite)) {
      stop("Each parasite must map to a single group.")
    }
    group_of <- stats::setNames(map$group, map$parasite)
    result$group <- unname(group_of[result$parasite])
  } else if (is.function(group)) {
    result$group <- vapply(result$parasite, group, character(1))
  }

  class(result) <- c("phylopred_specificity", class(result))
  result
}
