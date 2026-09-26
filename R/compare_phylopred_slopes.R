#' Compare host-sharing slopes between two clades
#'
#' Fits a single logistic regression with a group-by-phylogenetic-distance
#' interaction to a pairwise data set built with \code{\link{prepare_pair_data}}
#' and tests whether the slopes of the two groups differ. A more negative
#' slope means the probability of sharing a parasite decays faster with
#' phylogenetic distance, i.e. the clade is more of a phylogenetic specialist.
#' Uncertainty for the slope difference comes from a parasite-level cluster
#' bootstrap.
#'
#' @param pairs A data.frame from \code{\link{prepare_pair_data}} with a
#'   \code{group} column containing exactly two groups.
#' @param n_boot Integer. Number of cluster bootstrap replicates. Default 1000.
#' @param seed Integer or NULL. Random seed for reproducibility. Default NULL.
#' @param conf_level Confidence level. Default 0.95.
#'
#' @return An object of class \code{"phylopred_slope_comparison"}, a list with:
#' \describe{
#'   \item{groups}{The two group labels (reference first).}
#'   \item{slopes}{Named vector with the slope of each group.}
#'   \item{slope_difference}{Slope of the second group minus the reference.}
#'   \item{interaction_p_value}{Wald p-value of the interaction term.}
#'   \item{bootstrap}{The \code{\link{cluster_bootstrap_phylopred}} object.}
#'   \item{slope_cis}{Bootstrap percentile CIs for each group's slope and
#'     for the difference.}
#'   \item{fitted_model}{The underlying \code{glm} object.}
#' }
#'
#' @seealso \code{\link{prepare_pair_data}}, \code{\link{cluster_bootstrap_phylopred}}
#'
#' @export
#'
#' @examples
#' # Compare Xyleborus beetles against other genera (subset kept small
#' # so the bootstrap example runs quickly)
#' parasites <- unique(as.character(beetleTreeInteractions[[1]]))
#' keep <- c(grep("^Xyleborus", parasites, value = TRUE)[1:3],
#'           parasites[!grepl("^Xyleborus", parasites)][1:3])
#' sub <- beetleTreeInteractions[
#'   as.character(beetleTreeInteractions[[1]]) %in% keep, ]
#' pairs <- prepare_pair_data(
#'   sub, phy_dist,
#'   group = function(x) ifelse(grepl("^Xyleborus", x), "Xyleborus", "other")
#' )
#' compare_phylopred_slopes(pairs, n_boot = 5, seed = 42)
compare_phylopred_slopes <- function(pairs, n_boot = 1000, seed = NULL,
                                  conf_level = 0.95) {
  .validate_pairs(pairs)
  if (!"group" %in% colnames(pairs)) {
    stop("'pairs' must have a 'group' column; see prepare_pair_data(group = ...).")
  }
  groups <- sort(unique(pairs$group))
  if (length(groups) != 2) {
    stop("'pairs$group' must contain exactly two groups.")
  }

  pairs$group <- factor(pairs$group, levels = groups)
  fit <- stats::glm(suscept ~ phydist * group, data = pairs,
                    family = stats::binomial())
  ct <- stats::coef(summary(fit))
  int_term <- paste0("phydist:group", groups[2])

  boot <- cluster_bootstrap_phylopred(pairs, n_boot = n_boot, seed = seed,
                                   conf_level = conf_level)

  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  slope_ref  <- boot$draws[, "phydist"]
  slope_diff <- boot$draws[, int_term]
  slope_alt  <- slope_ref + slope_diff
  slope_cis <- rbind(
    stats::quantile(slope_ref, probs, na.rm = TRUE),
    stats::quantile(slope_alt, probs, na.rm = TRUE),
    stats::quantile(slope_diff, probs, na.rm = TRUE)
  )
  rownames(slope_cis) <- c(groups, "difference")

  cf <- stats::coef(fit)
  slopes <- c(cf["phydist"], cf["phydist"] + cf[int_term])
  names(slopes) <- groups

  structure(
    list(
      groups = groups,
      slopes = slopes,
      slope_difference = unname(cf[int_term]),
      interaction_p_value = unname(ct[int_term, "Pr(>|z|)"]),
      bootstrap = boot,
      slope_cis = slope_cis,
      fitted_model = fit
    ),
    class = "phylopred_slope_comparison"
  )
}

#' @export
print.phylopred_slope_comparison <- function(x, ...) {
  cat("Phylopred slope comparison:", x$groups[1], "vs", x$groups[2], "\n")
  cat(sprintf("Slope %s: %.4f  [%.4f, %.4f]\n", x$groups[1], x$slopes[1],
              x$slope_cis[1, 1], x$slope_cis[1, 2]))
  cat(sprintf("Slope %s: %.4f  [%.4f, %.4f]\n", x$groups[2], x$slopes[2],
              x$slope_cis[2, 1], x$slope_cis[2, 2]))
  cat(sprintf("Difference: %.4f  [%.4f, %.4f]  (Wald p = %.3g)\n",
              x$slope_difference, x$slope_cis[3, 1], x$slope_cis[3, 2],
              x$interaction_p_value))
  invisible(x)
}
