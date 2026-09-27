#' Evaluate the predictive capacity of a pairwise host-sharing model
#'
#' Computes goodness-of-fit and discrimination metrics for a logistic
#' host-sharing model fitted to a \code{\link{prepare_pair_data}} data set:
#' McFadden and Tjur pseudo-R-squared, the area under the ROC curve (AUC),
#' and optionally a grouped cross-validated AUC where whole parasites are
#' held out (parasite-grouped k-fold), which measures the ability to predict
#' the host range of unseen parasites.
#'
#' @param pairs A data.frame from \code{\link{prepare_pair_data}}.
#' @param cv Logical. If TRUE, also compute the parasite-grouped
#'   cross-validated AUC. Default FALSE (refitting can be slow on large
#'   data sets).
#' @param k Integer. Number of cross-validation folds. Default 5.
#' @param seed Integer or NULL. Random seed for fold assignment. Default NULL.
#'
#' @return An object of class \code{"phylopred_evaluation"}, a list with:
#' \describe{
#'   \item{mcfadden_r2}{McFadden pseudo-R-squared, \eqn{1 - logLik(model)/logLik(null)}.}
#'   \item{tjur_r2}{Tjur coefficient of discrimination (mean fitted
#'     probability for 1s minus mean for 0s).}
#'   \item{auc}{In-sample area under the ROC curve.}
#'   \item{cv_auc}{Parasite-grouped k-fold cross-validated AUC (NA if
#'     \code{cv = FALSE}).}
#'   \item{fitted_model}{The underlying \code{glm} object.}
#'   \item{n_pairs, n_parasites}{Data set size.}
#' }
#'
#' @seealso \code{\link{prepare_pair_data}}, \code{\link{cluster_bootstrap_phylopred}}
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
#' evaluate_phylopred_model(pairs)
#' evaluate_phylopred_model(pairs, cv = TRUE, k = 3, seed = 42)
evaluate_phylopred_model <- function(pairs, cv = FALSE, k = 5, seed = NULL) {
  .validate_pairs(pairs)
  if (!is.numeric(k) || length(k) != 1 || k < 2) {
    stop("'k' must be an integer >= 2.")
  }
  k <- as.integer(k)

  has_group <- "group" %in% colnames(pairs)
  fml <- if (has_group) suscept ~ phydist * group else suscept ~ phydist

  fit <- stats::glm(fml, data = pairs, family = stats::binomial())
  null_fit <- stats::glm(suscept ~ 1, data = pairs, family = stats::binomial())

  mcfadden <- 1 - as.numeric(stats::logLik(fit)) / as.numeric(stats::logLik(null_fit))
  p_hat <- stats::fitted(fit)
  tjur <- mean(p_hat[pairs$suscept == 1]) - mean(p_hat[pairs$suscept == 0])
  auc <- .auc(pairs$suscept, p_hat)

  cv_auc <- NA_real_
  if (cv) {
    if (!is.null(seed)) set.seed(seed)
    parasites <- unique(pairs$parasite)
    if (length(parasites) < k) {
      stop("Number of parasites must be at least 'k' for grouped CV.")
    }
    folds <- sample(rep_len(seq_len(k), length(parasites)))
    names(folds) <- parasites
    pred <- rep(NA_real_, nrow(pairs))
    for (fold in seq_len(k)) {
      test <- pairs$parasite %in% parasites[folds == fold]
      fit_f <- stats::glm(fml, data = pairs[!test, , drop = FALSE],
                          family = stats::binomial())
      pred[test] <- stats::predict(fit_f, newdata = pairs[test, , drop = FALSE],
                                   type = "response")
    }
    cv_auc <- .auc(pairs$suscept, pred)
  }

  structure(
    list(
      mcfadden_r2 = mcfadden,
      tjur_r2 = tjur,
      auc = auc,
      cv_auc = cv_auc,
      fitted_model = fit,
      n_pairs = nrow(pairs),
      n_parasites = length(unique(pairs$parasite))
    ),
    class = "phylopred_evaluation"
  )
}

# Mann-Whitney AUC: P(score_1 > score_0) + 0.5 P(tie)
.auc <- function(y, score) {
  ok <- !is.na(score)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @export
print.phylopred_evaluation <- function(x, ...) {
  cat("Phylopred model evaluation\n")
  cat(sprintf("Pairs: %d  Parasites: %d\n", x$n_pairs, x$n_parasites))
  cat(sprintf("McFadden pseudo-R2: %.4f\n", x$mcfadden_r2))
  cat(sprintf("Tjur R2:            %.4f\n", x$tjur_r2))
  cat(sprintf("AUC (in-sample):    %.4f\n", x$auc))
  if (!is.na(x$cv_auc)) {
    cat(sprintf("AUC (grouped CV):   %.4f\n", x$cv_auc))
  }
  invisible(x)
}
