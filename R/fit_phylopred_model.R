#' Fit a phylopred logistic regression model
#'
#' Fits a logistic regression model relating host susceptibility to phylogenetic
#' distance.
#'
#' @param incidence A binary incidence matrix (rows = parasites/pathogens,
#'   cols = hosts).
#' @param phydist A square, numeric phylogenetic distance matrix among hosts.
#'   Row/col names must match the column names of \code{incidence}.
#' @param seed Integer or NULL. Random seed for reproducibility. Default NULL.
#'
#' @return An S3 object of class \code{"phylopred_model"}, a list with:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{coefficients}{Named numeric vector \code{c(intercept, slope)}.}
#'   \item{standard_errors}{Named numeric vector \code{c(intercept, slope)}.}
#'   \item{confidence_intervals}{Matrix with rows \code{c("intercept","slope")} and
#'     cols \code{c("2.5\%","97.5\%")}.}
#'   \item{fitted_model}{The underlying \code{glm} object.}
#'   \item{training_data}{The data.frame used for model fitting.}
#'   \item{convergence}{Logical; whether the GLM converged.}
#'   \item{metadata}{List with \code{n_parasites}, \code{n_hosts}, and \code{seed}.}
#' }
#'
#' @seealso \code{\link{bootstrap_phylopred_model}}, \code{\link{predict_phylopred_probability}}
#'
#' @export
#'
#' @examples
#' incidence <- prepare_incidence_matrix(beetleTreeInteractions)
#' aligned   <- align_phylopred_inputs(incidence, phy_dist)
#' fit_phylopred_model(aligned$incidence, aligned$phydist, seed = 42)
fit_phylopred_model <- function(incidence, phydist, seed = NULL) {
  .validate_phylopred_inputs(incidence, phydist)

  if (!is.null(seed)) set.seed(seed)

  df      <- .build_glm_dataframe(incidence, phydist)
  glm_fit <- .fit_glm(df)
  ct      <- .extract_coef_table(glm_fit)

  ci <- matrix(
    c(ct$intercept_ci_lower, ct$intercept_ci_upper,
      ct$slope_ci_lower,     ct$slope_ci_upper),
    nrow = 2, ncol = 2, byrow = TRUE,
    dimnames = list(c("intercept", "slope"), c("2.5%", "97.5%"))
  )

  structure(
    list(
      call               = match.call(),
      coefficients       = c(intercept = ct$intercept, slope = ct$slope),
      standard_errors    = c(intercept = ct$intercept_se, slope = ct$slope_se),
      confidence_intervals = ci,
      fitted_model       = glm_fit,
      training_data      = df,
      convergence        = glm_fit$converged,
      metadata           = list(
        n_parasites = nrow(incidence),
        n_hosts     = ncol(incidence),
        seed        = seed
      )
    ),
    class = "phylopred_model"
  )
}

#' @export
print.phylopred_model <- function(x, ...) {
  cat("Phylopred model\n")
  cat(sprintf("Intercept: %.4f (SE: %.4f)\n",
              x$coefficients["intercept"], x$standard_errors["intercept"]))
  cat(sprintf("Slope:     %.4f (SE: %.4f)\n",
              x$coefficients["slope"], x$standard_errors["slope"]))
  cat(sprintf("Convergence: %s\n", x$convergence))
  invisible(x)
}

#' @export
summary.phylopred_model <- function(object, ...) {
  ct <- .extract_coef_table(object$fitted_model)
  result <- list(call = object$call, coefficients = ct)
  cat("Phylopred model summary\n")
  cat("Call: "); print(object$call)
  cat("\nCoefficients:\n")
  print(ct)
  invisible(result)
}
