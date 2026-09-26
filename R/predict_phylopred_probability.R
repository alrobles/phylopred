#' Predict host susceptibility probabilities from logistic regression coefficients
#'
#' Applies a logit function using intercept and slope coefficients to a vector
#' of phylogenetic distances. Accepts a numeric vector, a \code{phylopred_model},
#' or a \code{phylopred_bootstrap} object.
#'
#' @param coef A named numeric vector with elements \code{"intercept"} and
#'   \code{"slope"}, OR a \code{phylopred_model} object, OR a \code{phylopred_bootstrap}
#'   object (in which case the bootstrap summary means are used).
#'   An unnamed vector of length >= 2 is interpreted as \code{c(intercept, slope)}.
#' @param dist A numeric vector of phylogenetic distances.
#'
#' @return A named numeric vector of predicted probabilities of the same length
#'   as \code{dist}.
#'
#' @seealso \code{\link{fit_phylopred_model}}
#'
#' @export
#'
#' @examples
#' inc <- prepare_incidence_matrix(beetleTreeInteractions)
#' al  <- align_phylopred_inputs(inc, phy_dist)
#' model <- fit_phylopred_model(al$incidence, al$phydist, seed = 42)
#' dist_vec <- seq(0, max(al$phydist), length.out = 50)
#' predict_phylopred_probability(model, dist_vec)
predict_phylopred_probability <- function(coef, dist) {
  if (inherits(coef, "phylopred_model")) {
    intercept <- coef$coefficients["intercept"]
    slope     <- coef$coefficients["slope"]
  } else if (inherits(coef, "phylopred_bootstrap")) {
    intercept <- coef$summary$intercept
    slope     <- coef$summary$slope
  } else {
    if (!is.numeric(coef) || length(coef) < 2) {
      stop("'coef' must be a numeric vector of length >= 2, a phylopred_model, or a phylopred_bootstrap.")
    }
    if (!is.null(names(coef)) && "intercept" %in% names(coef) && "slope" %in% names(coef)) {
      intercept <- coef["intercept"]
      slope     <- coef["slope"]
    } else {
      intercept <- coef[1]
      slope     <- coef[2]
    }
  }

  if (!is.numeric(dist)) stop("'dist' must be a numeric vector.")

  logit <- intercept + slope * dist
  probs <- exp(logit) / (1 + exp(logit))
  names(probs) <- names(dist)
  probs
}

#' Extract coefficients from a phylopred model or bootstrap object
#'
#' Returns a named list of model coefficients and their associated statistics
#' from a \code{phylopred_model} or \code{phylopred_bootstrap} object.
#'
#' @param model A \code{phylopred_model} or \code{phylopred_bootstrap} object.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{intercept}{Intercept coefficient (mean for bootstrap).}
#'   \item{slope}{Slope coefficient (mean for bootstrap).}
#'   \item{intercept_se}{Standard error of intercept.}
#'   \item{slope_se}{Standard error of slope.}
#'   \item{intercept_ci}{Two-element vector with 2.5\% and 97.5\% CI for intercept.}
#'   \item{slope_ci}{Two-element vector with 2.5\% and 97.5\% CI for slope.}
#' }
#'
#' @export
#'
#' @examples
#' inc <- prepare_incidence_matrix(beetleTreeInteractions)
#' al  <- align_phylopred_inputs(inc, phy_dist)
#' model <- fit_phylopred_model(al$incidence, al$phydist, seed = 42)
#' extract_phylopred_coefficients(model)
extract_phylopred_coefficients <- function(model) {
  if (inherits(model, "phylopred_model")) {
    list(
      intercept    = unname(model$coefficients["intercept"]),
      slope        = unname(model$coefficients["slope"]),
      intercept_se = unname(model$standard_errors["intercept"]),
      slope_se     = unname(model$standard_errors["slope"]),
      intercept_ci = model$confidence_intervals["intercept", ],
      slope_ci     = model$confidence_intervals["slope", ]
    )
  } else if (inherits(model, "phylopred_bootstrap")) {
    s <- model$summary
    list(
      intercept    = s$intercept,
      slope        = s$slope,
      intercept_se = s$intercept_se,
      slope_se     = s$slope_se,
      intercept_ci = c("2.5%" = s$intercept_ci_lower, "97.5%" = s$intercept_ci_upper),
      slope_ci     = c("2.5%" = s$slope_ci_lower,     "97.5%" = s$slope_ci_upper)
    )
  } else {
    stop("'model' must be a phylopred_model or phylopred_bootstrap object.")
  }
}
