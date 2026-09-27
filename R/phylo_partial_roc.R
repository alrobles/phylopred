#' Phylogenetic partial ROC (Peterson, Papes & Soberon 2008 analogue)
#'
#' Partial ROC test for presence-only host-suitability predictions,
#' following Peterson, Papes & Soberon (2008, Ecological Modelling 213,
#' 63-72). The ROC space is redefined for presence-only data: the y axis is
#' the sensitivity \eqn{a(t)} (proportion of known hosts predicted at
#' cutoff \eqn{t}) and the x axis is the proportion of ALL species in the
#' phylogeny predicted as hosts, \eqn{b(t)} (the analogue of proportional
#' predicted area). The curve is restricted to the low-omission region
#' \eqn{a(t) \ge 1 - E} for a user-chosen omission tolerance \code{error_rate}
#' \eqn{E}, and performance is summarized as the ratio of the partial AUC of
#' the observed curve to the partial AUC of the null (random-ranking)
#' expectation, the diagonal \eqn{a = b}. A ratio of 1 equals random
#' performance; ratios above 1 indicate predictive signal. Significance is
#' assessed by bootstrap resampling of the known hosts.
#'
#' @param pred Numeric vector of predicted suitability scores, one per
#'   species in the phylogeny.
#' @param is_host Logical or 0/1 vector of the same length: TRUE/1 for
#'   known hosts.
#' @param error_rate Allowed omission error \eqn{E} in \eqn{[0, 1)}. The
#'   partial AUC is computed over the region with sensitivity
#'   \eqn{\ge 1 - E}. Default 0.05.
#' @param n_boot Number of bootstrap replicates. Default 500.
#' @param boot_prop Proportion of known hosts resampled (with replacement)
#'   in each bootstrap replicate, as in Peterson, Papes and Soberon
#'   (2008). Default 0.5.
#' @param seed Integer or NULL. Random seed for the bootstrap. Default NULL.
#'
#' @return An object of class \code{"phylopred_partial_roc"}, a list with:
#' \describe{
#'   \item{auc_ratio}{Observed partial AUC ratio (observed / null).}
#'   \item{p_value}{Bootstrap proportion of replicates with ratio \eqn{\le 1}.}
#'   \item{boot_ratios}{Vector of bootstrap AUC ratios.}
#'   \item{curve}{data.frame with \code{threshold}, \code{a} (sensitivity)
#'     and \code{b} (predicted breadth).}
#'   \item{error_rate, n_boot, boot_prop}{The configuration used.}
#' }
#'
#' @seealso \code{\link{host_threshold_metric}}
#'
#' @export
#'
#' @examples
#' pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
#' fit <- stats::glm(suscept ~ phydist, data = pairs, family = stats::binomial())
#' one <- pairs[pairs$parasite == pairs$parasite[1], ]
#' one <- one[!duplicated(one$target), ]
#' pred <- stats::predict(fit, newdata = one, type = "response")
#' phylo_partial_roc(pred, one$suscept, n_boot = 50, seed = 42)
phylo_partial_roc <- function(pred, is_host, error_rate = 0.05,
                              n_boot = 500, boot_prop = 0.5, seed = NULL) {
  curve <- .threshold_curve(pred, is_host)
  if (!is.numeric(error_rate) || error_rate < 0 || error_rate >= 1) {
    stop("'error_rate' must be in [0, 1).")
  }
  if (!is.numeric(n_boot) || n_boot < 1) stop("'n_boot' must be >= 1.")
  if (!is.numeric(boot_prop) || boot_prop <= 0 || boot_prop > 1) {
    stop("'boot_prop' must be in (0, 1].")
  }
  if (!is.null(seed)) set.seed(seed)

  is_host <- as.logical(is_host)
  observed <- .partial_auc_ratio(pred, is_host, error_rate)

  host_idx <- which(is_host)
  n_take <- max(1L, floor(boot_prop * length(host_idx)))
  boot_ratios <- vapply(seq_len(n_boot), function(b) {
    take <- sample(host_idx, n_take, replace = TRUE)
    keep <- c(take, which(!is_host))
    .partial_auc_ratio(pred[keep],
                       c(rep(TRUE, n_take), rep(FALSE, sum(!is_host))),
                       error_rate)
  }, numeric(1))
  boot_ratios <- boot_ratios[!is.na(boot_ratios)]

  structure(
    list(
      auc_ratio = observed,
      p_value = if (length(boot_ratios) > 0) mean(boot_ratios <= 1) else NA_real_,
      boot_ratios = boot_ratios,
      curve = curve[, c("threshold", "a", "b")],
      error_rate = error_rate,
      n_boot = n_boot,
      boot_prop = boot_prop
    ),
    class = "phylopred_partial_roc"
  )
}

# Partial AUC ratio over the region with sensitivity >= 1 - error_rate.
# Observed AUC integrates a db (trapezoid); null AUC is the diagonal a = b
# over the same b range: (b_hi^2 - b_lo^2) / 2.
.partial_auc_ratio <- function(pred, is_host, error_rate) {
  curve <- .threshold_curve(pred, is_host)
  keep <- curve$a >= 1 - error_rate
  a <- curve$a[keep]
  b <- curve$b[keep]
  if (length(b) < 2) return(NA_real_)
  ord <- order(b)
  a <- a[ord]
  b <- b[ord]
  obs <- sum(diff(b) * (utils::head(a, -1) + utils::tail(a, -1)) / 2)
  null <- (max(b)^2 - min(b)^2) / 2
  if (null <= 0) return(NA_real_)
  obs / null
}

#' @export
plot.phylopred_partial_roc <- function(x, ...) {
  ord <- order(x$curve$b)
  graphics::plot(x$curve$b[ord], x$curve$a[ord], type = "s",
                 xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "proportion of phylogeny predicted (b)",
                 ylab = "sensitivity (a)",
                 main = sprintf("Partial ROC (E = %g, AUC ratio = %.2f)",
                                x$error_rate, x$auc_ratio),
                 ...)
  graphics::abline(0, 1, lty = 2)
  graphics::abline(h = 1 - x$error_rate, lty = 3)
  invisible(x)
}

#' @export
print.phylopred_partial_roc <- function(x, ...) {
  cat("Phylopred phylogenetic partial ROC (E = ", x$error_rate, ")\n", sep = "")
  cat(sprintf("Partial AUC ratio: %.4f\n", x$auc_ratio))
  cat(sprintf("Bootstrap p-value (ratio <= 1): %.4f  [%d replicates]\n",
              x$p_value, length(x$boot_ratios)))
  invisible(x)
}
