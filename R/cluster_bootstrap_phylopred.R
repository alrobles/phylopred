#' Cluster bootstrap for pairwise host-sharing models
#'
#' Nonparametric bootstrap that resamples whole parasites (clusters) with
#' replacement and refits the logistic regression on each replicate. This
#' respects the non-independence of rows that share a parasite in a
#' \code{\link{prepare_pair_data}} data set, unlike a naive row-level
#' bootstrap. When a \code{group} column is present, resampling is stratified
#' within groups and the model includes a group-by-distance interaction.
#'
#' @param pairs A data.frame from \code{\link{prepare_pair_data}}.
#' @param n_boot Integer. Number of bootstrap replicates. Default 1000.
#' @param seed Integer or NULL. Random seed for reproducibility. Default NULL.
#' @param conf_level Confidence level for percentile intervals. Default 0.95.
#'
#' @return An object of class \code{"phylopred_cluster_bootstrap"}, a list with:
#' \describe{
#'   \item{draws}{Matrix of bootstrap coefficient draws (replicates x terms).}
#'   \item{estimates}{Coefficients of the model fitted to the full data.}
#'   \item{confidence_intervals}{Percentile confidence intervals per term.}
#'   \item{n_boot, seed, conf_level}{Bootstrap settings.}
#' }
#'
#' @seealso \code{\link{prepare_pair_data}}, \code{\link{compare_phylopred_slopes}}
#'
#' @export
#'
#' @examples
#' # Three parasites keep the example fast
#' top <- names(sort(table(beetleTreeInteractions[[1]]),
#'                   decreasing = TRUE))[1:3]
#' sub <- beetleTreeInteractions[
#'   as.character(beetleTreeInteractions[[1]]) %in% top, ]
#' pairs <- prepare_pair_data(sub, phy_dist)
#' cluster_bootstrap_phylopred(pairs, n_boot = 5, seed = 42)
cluster_bootstrap_phylopred <- function(pairs, n_boot = 1000, seed = NULL,
                                     conf_level = 0.95) {
  .validate_pairs(pairs)
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) {
    stop("'n_boot' must be a positive integer.")
  }
  n_boot <- as.integer(n_boot)
  if (!is.null(seed)) set.seed(seed)

  has_group <- "group" %in% colnames(pairs)
  fml <- if (has_group) suscept ~ phydist * group else suscept ~ phydist

  full_fit <- stats::glm(fml, data = pairs, family = stats::binomial())
  terms <- names(stats::coef(full_fit))

  idx_by_parasite <- split(seq_len(nrow(pairs)), pairs$parasite)
  parasites <- names(idx_by_parasite)
  strata <- if (has_group) {
    vapply(idx_by_parasite, function(i) as.character(pairs$group[i[1]]),
           character(1))
  } else {
    rep("all", length(parasites))
  }

  draws <- matrix(NA_real_, nrow = n_boot, ncol = length(terms),
                  dimnames = list(NULL, terms))
  for (b in seq_len(n_boot)) {
    samp <- unlist(lapply(unique(strata), function(s) {
      pool <- parasites[strata == s]
      sample(pool, length(pool), replace = TRUE)
    }), use.names = FALSE)
    rows <- unlist(idx_by_parasite[samp], use.names = FALSE)
    fit_b <- try(stats::glm(fml, data = pairs[rows, , drop = FALSE],
                            family = stats::binomial()), silent = TRUE)
    if (!inherits(fit_b, "try-error")) {
      cf <- stats::coef(fit_b)
      draws[b, names(cf)] <- cf
    }
  }

  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  ci <- t(apply(draws, 2, stats::quantile, probs = probs, na.rm = TRUE))

  structure(
    list(
      draws = draws,
      estimates = stats::coef(full_fit),
      confidence_intervals = ci,
      n_boot = n_boot,
      seed = seed,
      conf_level = conf_level
    ),
    class = "phylopred_cluster_bootstrap"
  )
}

#' @export
print.phylopred_cluster_bootstrap <- function(x, ...) {
  cat("Phylopred cluster bootstrap (", x$n_boot, " replicates)\n", sep = "")
  out <- cbind(estimate = x$estimates, x$confidence_intervals)
  print(round(out, 4))
  invisible(x)
}
