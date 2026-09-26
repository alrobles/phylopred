#' Prepare an incidence matrix from an interaction table
#'
#' Builds an incidence (presence/absence) matrix from a two-column data.frame
#' of species interactions.
#'
#' @param db A data.frame with exactly two columns (e.g., parasite and host).
#'
#' @return An incidence matrix (rows = first column unique values, cols =
#'   second column unique values) filled with 0/1.
#'
#' @export
#'
#' @examples
#' prepare_incidence_matrix(beetleTreeInteractions)
prepare_incidence_matrix <- function(db) {
  if (!inherits(db, "data.frame")) {
    stop("'db' must be a data.frame.")
  }
  if (ncol(db) != 2) {
    stop("'db' must have exactly two columns.")
  }

  m <- as.matrix(db)
  interactors  <- unique(m[, 1])
  interactuans <- unique(m[, 2])
  incidence_mat <- sapply(interactors, function(x) {
    as.numeric(interactuans %in% m[(m[, 1] %in% x), 2])
  })
  rownames(incidence_mat) <- interactuans
  colnames(incidence_mat) <- interactors
  t(incidence_mat)
}

#' Align incidence matrix and phylogenetic distance matrix
#'
#' Validates both inputs and reorders the columns of \code{incidence} to match
#' the column order of \code{phydist}.
#'
#' @param incidence A binary incidence matrix (rows = parasites, cols = hosts).
#' @param phydist A square phylogenetic distance matrix among hosts.
#'
#' @return A named list with:
#' \describe{
#'   \item{incidence}{The incidence matrix with columns reordered to match \code{phydist}.}
#'   \item{phydist}{The (unchanged) phylogenetic distance matrix.}
#' }
#'
#' @export
#'
#' @examples
#' incidence <- prepare_incidence_matrix(beetleTreeInteractions)
#' aligned <- align_phylopred_inputs(incidence, phy_dist)
align_phylopred_inputs <- function(incidence, phydist) {
  .validate_incidence(incidence)
  .validate_phydist(phydist)

  common_cols <- intersect(colnames(phydist), colnames(incidence))
  if (length(common_cols) == 0) {
    stop("No shared column names between 'incidence' and 'phydist'.")
  }
  incidence_aligned <- incidence[, common_cols, drop = FALSE]
  phydist_aligned   <- phydist[common_cols, common_cols]

  list(incidence = incidence_aligned, phydist = phydist_aligned)
}
