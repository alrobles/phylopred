#' Build a pairwise host-sharing data set from an interaction table
#'
#' Expands an interaction table into all (focal host, target host) pairs per
#' parasite. For every parasite, each of its known hosts is used as a focal
#' host, and every other host in \code{phydist} is a target host whose
#' response is 1 if it is also a known host of that parasite and 0 otherwise.
#'
#' Unlike the legacy workflow (which samples a single focal host per parasite),
#' this uses every known host as a focal host, so no information is discarded
#' and results are deterministic. Rows belonging to the same parasite are not
#' independent; account for this with \code{\link{cluster_bootstrap_phylopred}}
#' when computing uncertainty.
#'
#' @param interactions A data.frame with at least two columns: the first the
#'   parasite (or symbiont) species, the second the host species. Extra
#'   columns are ignored. Duplicated rows are removed.
#' @param phydist A square, numeric phylogenetic (or taxonomic) distance
#'   matrix among hosts, with matching row and column names. Hosts in
#'   \code{interactions} absent from \code{phydist} are dropped with a warning.
#' @param group Optional. Either the name of a column in \code{interactions}
#'   giving a grouping factor for each parasite (e.g. genus or clade), or a
#'   function applied to the parasite name to derive the group (e.g.
#'   \code{function(x) sub(" .*", "", x)} for the genus).
#'
#' @return A data.frame of class \code{"phylopred_pairs"} with columns:
#' \describe{
#'   \item{parasite}{Parasite species (cluster identifier).}
#'   \item{focal}{Focal (known) host.}
#'   \item{target}{Target host.}
#'   \item{phydist}{Distance between focal and target host.}
#'   \item{suscept}{1 if the target host is a known host of the parasite.}
#'   \item{group}{Grouping factor (only if \code{group} is supplied).}
#' }
#'
#' @seealso \code{\link{fit_phylopred_model}}, \code{\link{cluster_bootstrap_phylopred}},
#'   \code{\link{compare_phylopred_slopes}}
#'
#' @export
#'
#' @examples
#' pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
#' head(pairs)
prepare_pair_data <- function(interactions, phydist, group = NULL) {
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
  host     <- as.character(interactions[[2]])
  hosts_all <- rownames(phydist)

  in_tree <- host %in% hosts_all
  if (!any(in_tree)) {
    stop("No hosts in 'interactions' match the rownames of 'phydist'.")
  }
  if (!all(in_tree)) {
    warning(sprintf("%d interaction record(s) dropped: host not in 'phydist'.",
                    sum(!in_tree)))
  }

  keep <- which(in_tree)
  df0 <- unique(data.frame(parasite = parasite[keep], host = host[keep],
                           stringsAsFactors = FALSE))

  group_of <- NULL
  if (!is.null(group_col)) {
    map <- unique(data.frame(parasite = parasite[keep],
                             group = as.character(interactions[[group_col]][keep]),
                             stringsAsFactors = FALSE))
    if (anyDuplicated(map$parasite)) {
      stop("Each parasite must map to a single group.")
    }
    group_of <- stats::setNames(map$group, map$parasite)
  } else if (is.function(group)) {
    group_of <- stats::setNames(vapply(unique(df0$parasite), group, character(1)),
                                unique(df0$parasite))
  }

  hosts_by_parasite <- split(df0$host, df0$parasite)
  out <- lapply(names(hosts_by_parasite), function(p) {
    H <- hosts_by_parasite[[p]]
    do.call(rbind, lapply(H, function(i) {
      targets <- setdiff(hosts_all, i)
      data.frame(parasite = p, focal = i, target = targets,
                 phydist = phydist[i, targets],
                 suscept = as.numeric(targets %in% H),
                 stringsAsFactors = FALSE)
    }))
  })
  pairs <- do.call(rbind, out)
  rownames(pairs) <- NULL
  if (!is.null(group_of)) {
    pairs$group <- unname(group_of[pairs$parasite])
  }
  class(pairs) <- c("phylopred_pairs", class(pairs))
  pairs
}
