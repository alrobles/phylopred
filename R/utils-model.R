# Internal model building helpers - not exported

.build_glm_dataframe <- function(incidence, phydist) {
  bool_mat <- apply(incidence, 2, function(x) x == 1)
  if (!is.matrix(bool_mat)) {
    bool_mat <- matrix(bool_mat, nrow = 1,
                       dimnames = list(rownames(incidence), colnames(incidence)))
  }
  hosts_by_parasite <- lapply(seq_len(nrow(bool_mat)), function(i) {
    names(which(bool_mat[i, ]))
  })
  names(hosts_by_parasite) <- rownames(bool_mat)
  # Remove parasites with no hosts
  hosts_by_parasite <- hosts_by_parasite[sapply(hosts_by_parasite, length) > 0]
  if (length(hosts_by_parasite) == 0) {
    stop("No parasites with at least one host found in the incidence matrix.")
  }
  focal <- lapply(hosts_by_parasite, function(x) sample(x, 1))
  df <- do.call(rbind, lapply(names(focal), function(x) {
    data.frame(
      tohost   = colnames(incidence),
      fromhost = focal[[x]],
      phydist  = phydist[, focal[[x]]],
      suscept  = as.numeric(incidence[x, ]),
      focal    = x,
      stringsAsFactors = FALSE
    )
  }))
  rownames(df) <- NULL
  df
}

.fit_glm <- function(df) {
  stats::glm(suscept ~ phydist, data = df, family = stats::binomial(link = "logit"))
}

.extract_coef_table <- function(glm_fit) {
  stats_out <- stats::coef(summary(glm_fit))
  conf_int  <- stats::confint.default(glm_fit)
  data.frame(
    intercept          = glm_fit$coefficients[1],
    intercept_se       = stats_out[1, 2],
    intercept_z        = stats_out[1, 3],
    intercept_p        = stats_out[1, 4],
    intercept_ci_lower = conf_int[1, 1],
    intercept_ci_upper = conf_int[1, 2],
    slope              = glm_fit$coefficients[2],
    slope_se           = stats_out[2, 2],
    slope_z            = stats_out[2, 3],
    slope_p            = stats_out[2, 4],
    slope_ci_lower     = conf_int[2, 1],
    slope_ci_upper     = conf_int[2, 2],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
