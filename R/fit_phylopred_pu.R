#' Positive-unlabeled host prediction via xplus (PLUS algorithm)
#'
#' Predicts the hosts of a parasite from the phylogeny using
#' positive-unlabeled (PU) learning, as a parallel alternative to the
#' pairwise logistic model. Known hosts are treated as positives and every
#' other species in the phylogeny as unlabeled (never as a true absence).
#' Each species is described by phylogenetic features and the PLUS
#' algorithm (\code{xplus::xplus}) iteratively re-labels the unlabeled
#' species while fitting a penalized logistic model.
#'
#' The feature matrix contains, for every species in \code{phydist}:
#' \itemize{
#'   \item the minimum and mean phylogenetic distance to the parasite's
#'     known hosts (excluding the species itself), and
#'   \item the first \code{n_eigen} principal coordinates (PCoA) of the
#'     phylogenetic distance matrix, describing the global position of the
#'     species in the phylogeny (the analogue of environmental layers in
#'     an SDM).
#' }
#'
#' @param interactions A data.frame whose first two columns are the
#'   parasite and host species names (as in
#'   \code{\link{prepare_pair_data}}).
#' @param phydist Square numeric matrix of phylogenetic distances with
#'   species names as row and column names.
#' @param parasite Character scalar: the parasite whose hosts are
#'   predicted.
#' @param n_eigen Number of PCoA axes used as features. Default 10.
#' @param seed Optional random seed passed to \code{xplus::xplus}.
#' @param ... Further arguments passed to \code{xplus::xplus}
#'   (e.g. \code{max_iter}, \code{learning_rate}, \code{qq}).
#'
#' @return An object of class \code{"phylopred_pu"}: a list with
#' \describe{
#'   \item{prediction}{data.frame with \code{species}, \code{pred}
#'     (predicted host suitability) and \code{is_host}.}
#'   \item{fit}{The underlying \code{xplus} object.}
#'   \item{parasite}{The parasite name.}
#' }
#'   The \code{prediction} columns can be passed directly to
#'   \code{\link{host_threshold_metric}} and
#'   \code{\link{phylo_partial_roc}}.
#'
#' @seealso \code{\link{host_threshold_metric}},
#'   \code{\link{phylo_partial_roc}}, \code{\link{prepare_pair_data}}
#'
#' @export
#'
#' @examples
#' pu <- fit_phylopred_pu(beetleTreeInteractions, phy_dist,
#'                      parasite = beetleTreeInteractions[[1]][1],
#'                      max_iter = 20, seed = 1)
#' host_threshold_metric(pu$prediction$pred, pu$prediction$is_host)
fit_phylopred_pu <- function(interactions, phydist, parasite,
                          n_eigen = 10, seed = NULL, ...) {
  if (!is.data.frame(interactions) || ncol(interactions) < 2) {
    stop("'interactions' must be a data.frame with at least two columns.")
  }
  if (!is.matrix(phydist) || nrow(phydist) != ncol(phydist) ||
      is.null(rownames(phydist))) {
    stop("'phydist' must be a square named matrix of phylogenetic distances.")
  }
  parasite <- as.character(parasite)
  if (length(parasite) != 1 || is.na(parasite)) {
    stop("'parasite' must be a single parasite name.")
  }
  if (!is.numeric(n_eigen) || n_eigen < 1) stop("'n_eigen' must be >= 1.")

  species <- rownames(phydist)
  hosts <- unique(as.character(interactions[[2]])[as.character(interactions[[1]]) == parasite])
  hosts <- intersect(hosts, species)
  if (length(hosts) == 0) {
    stop(sprintf("No known hosts of '%s' are present in 'phydist'.", parasite))
  }
  if (length(hosts) == length(species)) {
    stop("All species in 'phydist' are known hosts; nothing to predict.")
  }

  dist_to_hosts <- phydist[, hosts, drop = FALSE]
  # Exclude self-distance (0) so known hosts are not trivially separable.
  min_dist <- vapply(seq_along(species), function(i) {
    d <- dist_to_hosts[i, setdiff(hosts, species[i]), drop = TRUE]
    if (length(d) == 0) 0 else min(d)
  }, numeric(1))
  mean_dist <- vapply(seq_along(species), function(i) {
    d <- dist_to_hosts[i, setdiff(hosts, species[i]), drop = TRUE]
    if (length(d) == 0) 0 else mean(d)
  }, numeric(1))

  n_eigen <- min(as.integer(n_eigen), length(species) - 1L)
  pcoa <- stats::cmdscale(stats::as.dist(phydist), k = n_eigen)

  x <- cbind(min_dist = min_dist, mean_dist = mean_dist, pcoa)
  colnames(x) <- c("min_dist", "mean_dist", paste0("pcoa", seq_len(n_eigen)))
  y <- as.integer(species %in% hosts)

  fit <- xplus::xplus(x = x, y = y, seed = seed, ...)
  pred <- as.numeric(stats::predict(fit, newx = x, type = "response"))

  structure(
    list(
      prediction = data.frame(species = species, pred = pred,
                              is_host = y, row.names = NULL),
      fit = fit,
      parasite = parasite
    ),
    class = "phylopred_pu"
  )
}

#' @export
print.phylopred_pu <- function(x, ...) {
  cat("Phylopred PU-learning host prediction (xplus/PLUS)\n")
  cat(sprintf("Parasite: %s\n", x$parasite))
  cat(sprintf("Known hosts: %d of %d species\n",
              sum(x$prediction$is_host), nrow(x$prediction)))
  cat(sprintf("Iterations: %d (stop: %s)\n",
              x$fit$n_iter, x$fit$stop_reason))
  invisible(x)
}
