#' Threshold metric f_gamma for predicted host suitability
#'
#' Computes the threshold-dependent performance metric
#' \deqn{f_\gamma(t) = a(t) / (\gamma b(t) + c(t))}
#' where, for a cutoff \eqn{t} on the predicted host-suitability scores:
#' \describe{
#'   \item{\code{a(t)}}{proportion of known hosts with prediction \eqn{\ge t}
#'     (sensitivity / recall);}
#'   \item{\code{b(t)}}{proportion of ALL species in the phylogeny with
#'     prediction \eqn{\ge t} (predicted host breadth, the phylogenetic
#'     analogue of predicted area in SDM threshold metrics);}
#'   \item{\code{c(t) = 1 - a(t)}}{proportion of known hosts lost (omission).}
#' }
#' The metric is designed for presence-only / positive-unlabeled settings:
#' it never uses absences, and \code{gamma} tunes the trade-off between
#' retaining known hosts and predicting a compact host set. The maximum of
#' \eqn{f_\gamma(t)} gives an optimal cutoff \eqn{t^*}.
#'
#' @param pred Numeric vector of predicted suitability scores, one per
#'   species in the phylogeny (any real values; typically probabilities).
#' @param is_host Logical or 0/1 vector of the same length: TRUE/1 for
#'   species that are known hosts (presences).
#' @param gamma Positive scalar weight on the predicted breadth term.
#'   Default 1.
#'
#' @return An object of class \code{"phylopred_threshold_metric"}: a data.frame
#'   with columns \code{threshold}, \code{a}, \code{b}, \code{c}, and
#'   \code{f_gamma}, evaluated at every distinct prediction value. The
#'   optimal cutoff is stored in attributes \code{"optimal_threshold"} and
#'   \code{"optimal_f"}, and \code{gamma} in attribute \code{"gamma"}.
#'
#' @seealso \code{\link{phylo_partial_roc}}, \code{\link{evaluate_phylopred_model}}
#'
#' @export
#'
#' @examples
#' pairs <- prepare_pair_data(beetleTreeInteractions, phy_dist)
#' fit <- stats::glm(suscept ~ phydist, data = pairs, family = stats::binomial())
#' # Suitability of every host species for one focal parasite
#' one <- pairs[pairs$parasite == pairs$parasite[1], ]
#' one <- one[!duplicated(one$target), ]
#' pred <- stats::predict(fit, newdata = one, type = "response")
#' m <- host_threshold_metric(pred, one$suscept, gamma = 1)
#' plot(m)
host_threshold_metric <- function(pred, is_host, gamma = 1) {
  metrics <- .threshold_curve(pred, is_host)
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("'gamma' must be a positive scalar.")
  }

  metrics$f_gamma <- metrics$a / (gamma * metrics$b + metrics$c)
  best <- which.max(metrics$f_gamma)

  structure(
    metrics,
    optimal_threshold = metrics$threshold[best],
    optimal_f = metrics$f_gamma[best],
    gamma = gamma,
    class = c("phylopred_threshold_metric", class(metrics))
  )
}

# Step curve of sensitivity a(t) and predicted breadth b(t) at every
# distinct prediction value.
.threshold_curve <- function(pred, is_host) {
  if (!is.numeric(pred)) stop("'pred' must be a numeric vector.")
  is_host <- as.logical(is_host)
  if (length(pred) != length(is_host)) {
    stop("'pred' and 'is_host' must have the same length.")
  }
  if (anyNA(pred) || anyNA(is_host)) {
    stop("'pred' and 'is_host' must not contain NA values.")
  }
  n_pos <- sum(is_host)
  if (n_pos == 0) stop("'is_host' must contain at least one known host.")
  if (n_pos == length(pred)) {
    stop("'is_host' must contain at least one non-host (unlabeled) species.")
  }

  thresholds <- sort(unique(pred))
  a <- vapply(thresholds, function(t) mean(pred[is_host] >= t), numeric(1))
  b <- vapply(thresholds, function(t) mean(pred >= t), numeric(1))
  data.frame(threshold = thresholds, a = a, b = b, c = 1 - a)
}

#' @export
plot.phylopred_threshold_metric <- function(x, ...) {
  graphics::plot(x$threshold, x$f_gamma, type = "s",
                 xlab = "threshold", ylab = expression(f[gamma]),
                 main = sprintf("f_gamma curve (gamma = %g)", attr(x, "gamma")),
                 ...)
  graphics::abline(v = attr(x, "optimal_threshold"), lty = 2)
  invisible(x)
}

#' @export
print.phylopred_threshold_metric <- function(x, ...) {
  cat("Phylopred threshold metric f_gamma (gamma = ",
      attr(x, "gamma"), ")\n", sep = "")
  cat(sprintf("Optimal threshold: %.4f  f_gamma: %.4f\n",
              attr(x, "optimal_threshold"), attr(x, "optimal_f")))
  cat(sprintf("Evaluated at %d thresholds.\n", nrow(x)))
  invisible(x)
}
