#' Leave-hosts-out validation of host predictions
#'
#' Presence-only predictive validation: a fraction of the parasite's known
#' hosts is hidden, the model is refitted on the remaining hosts, and the
#' hidden hosts are scored against the unlabeled species. Because true
#' absences do not exist, performance is measured as the AUC of the hidden
#' (held-out) hosts versus the unlabeled species — the probability that a
#' truly interacting species outranks a random unlabeled one. This measures
#' actual predictive ability, not in-sample fit.
#'
#' @param interactions A data.frame whose first two columns are the
#'   parasite and host species names.
#' @param phydist Square numeric matrix of phylogenetic distances with
#'   species names as row and column names.
#' @param parasite Character scalar: the parasite whose hosts are
#'   predicted.
#' @param method \code{"glm"} for the pairwise logistic model (target
#'   species scored by their maximum predicted sharing probability across
#'   the training focal hosts) or \code{"pu"} for the xplus PU-learning
#'   model (\code{\link{fit_phylopred_pu}}).
#' @param prop_holdout Proportion of known hosts hidden in each replicate.
#'   Default 0.25.
#' @param n_rep Number of holdout replicates. Default 10.
#' @param seed Optional random seed.
#' @param ... Further arguments passed to \code{\link{fit_phylopred_pu}}
#'   when \code{method = "pu"}.
#'
#' @return An object of class \code{"phylopred_holdout"}: a list with
#'   \code{auc} (vector of holdout AUCs, one per replicate),
#'   \code{mean_auc}, \code{method}, \code{parasite},
#'   \code{prop_holdout}, and \code{n_rep}.
#'
#' @seealso \code{\link{fit_phylopred_pu}}, \code{\link{host_threshold_metric}},
#'   \code{\link{phylo_partial_roc}}
#'
#' @export
#'
#' @examples
#' parasite <- names(sort(table(beetleTreeInteractions[[1]]),
#'                        decreasing = TRUE))[1]
#' validate_host_prediction(beetleTreeInteractions, phy_dist, parasite,
#'                          method = "glm", n_rep = 3, seed = 1)
validate_host_prediction <- function(interactions, phydist, parasite,
                                     method = c("glm", "pu"),
                                     prop_holdout = 0.25, n_rep = 10,
                                     seed = NULL, ...) {
  method <- match.arg(method)
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
  if (!is.numeric(prop_holdout) || prop_holdout <= 0 || prop_holdout >= 1) {
    stop("'prop_holdout' must be in (0, 1).")
  }
  if (!is.numeric(n_rep) || n_rep < 1) stop("'n_rep' must be >= 1.")
  if (!is.null(seed)) set.seed(seed)

  species <- rownames(phydist)
  hosts <- unique(as.character(interactions[[2]])[as.character(interactions[[1]]) == parasite])
  hosts <- intersect(hosts, species)
  if (length(hosts) < 3) {
    stop("Need at least 3 known hosts in 'phydist' for holdout validation.")
  }

  n_hide <- max(1L, floor(prop_holdout * length(hosts)))
  if (n_hide >= length(hosts)) n_hide <- length(hosts) - 1L

  auc <- vapply(seq_len(n_rep), function(r) {
    hidden <- sample(hosts, n_hide)
    train_hosts <- setdiff(hosts, hidden)

    scores <- if (method == "pu") {
      train_int <- data.frame(parasite = parasite, host = train_hosts)
      pu <- fit_phylopred_pu(train_int, phydist, parasite, ...)
      stats::setNames(pu$prediction$pred, pu$prediction$species)
    } else {
      train_int <- data.frame(parasite = parasite, host = train_hosts)
      pairs <- suppressWarnings(prepare_pair_data(train_int, phydist))
      fit <- stats::glm(suscept ~ phydist, data = pairs,
                        family = stats::binomial())
      # Max predicted sharing probability across training focal hosts.
      p <- stats::predict(fit,
                          newdata = data.frame(phydist = as.vector(
                            phydist[train_hosts, , drop = FALSE])),
                          type = "response")
      pm <- matrix(p, nrow = length(train_hosts))
      stats::setNames(apply(pm, 2, max), colnames(phydist))
    }

    unlabeled <- setdiff(species, hosts)
    s_pos <- scores[hidden]
    s_neg <- scores[unlabeled]
    r_all <- rank(c(s_pos, s_neg))
    n1 <- length(s_pos)
    (sum(r_all[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(s_neg))
  }, numeric(1))

  structure(
    list(auc = auc, mean_auc = mean(auc), method = method,
         parasite = parasite, prop_holdout = prop_holdout, n_rep = n_rep),
    class = "phylopred_holdout"
  )
}

#' @export
print.phylopred_holdout <- function(x, ...) {
  cat("Phylopred leave-hosts-out validation (", x$method, ")\n", sep = "")
  cat(sprintf("Parasite: %s\n", x$parasite))
  cat(sprintf("Holdout AUC (hidden hosts vs unlabeled): mean %.3f over %d reps [%.3f, %.3f]\n",
              x$mean_auc, x$n_rep, min(x$auc), max(x$auc)))
  invisible(x)
}
