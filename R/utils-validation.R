# Internal validation helpers - not exported

.validate_tax_table <- function(tax_table) {
  if (!inherits(tax_table, "data.frame")) {
    stop("'tax_table' must be a data.frame.")
  }
  fullNames <- c("class", "order", "family", "genus", "species")
  matchNames <- match(colnames(tax_table), fullNames)
  matchNames <- matchNames[!is.na(matchNames)]
  if (length(matchNames) < 3) {
    stop("'tax_table' must have at least 3 columns matching: class, order, family, genus, species.")
  }
  if (nrow(unique(tax_table)) < 2) {
    stop("'tax_table' must have at least 2 unique rows.")
  }
  matched_cols <- colnames(tax_table)[colnames(tax_table) %in% fullNames]
  if (length(matched_cols) > 0 && any(is.na(tax_table[matched_cols]))) {
    stop("'tax_table' must not contain NA values in taxonomic rank columns.")
  }
  invisible(NULL)
}

.validate_incidence <- function(incidence) {
  if (!is.matrix(incidence) && !is.data.frame(incidence)) {
    stop("'incidence' must be a matrix or data.frame.")
  }
  if (is.null(rownames(incidence)) || is.null(colnames(incidence))) {
    stop("'incidence' must have both rownames and colnames.")
  }
  dup_rows <- duplicated(rownames(incidence))
  if (any(dup_rows)) {
    stop(paste0("'incidence' must not have duplicate rownames: ",
                paste(unique(rownames(incidence)[dup_rows]), collapse = ", ")))
  }
  dup_cols <- duplicated(colnames(incidence))
  if (any(dup_cols)) {
    stop(paste0("'incidence' must not have duplicate colnames: ",
                paste(unique(colnames(incidence)[dup_cols]), collapse = ", ")))
  }
  vals <- as.numeric(as.matrix(incidence))
  if (!all(vals %in% c(0, 1, NA))) {
    stop("'incidence' must be a binary (0/1) matrix.")
  }
  if (any(apply(incidence, 1, function(r) all(is.na(r))))) {
    stop("'incidence' must not contain all-NA rows.")
  }
  invisible(NULL)
}

.validate_phydist <- function(phydist) {
  if (!is.matrix(phydist)) {
    stop("'phydist' must be a matrix.")
  }
  if (nrow(phydist) != ncol(phydist)) {
    stop("'phydist' must be a square matrix.")
  }
  if (!identical(rownames(phydist), colnames(phydist))) {
    stop("'phydist' rownames must equal colnames.")
  }
  if (!is.numeric(phydist)) {
    stop("'phydist' must be numeric.")
  }
  if (anyNA(phydist)) {
    stop("'phydist' must not contain NA values.")
  }
  if (!all(diag(phydist) == 0)) {
    stop("'phydist' diagonal should be 0.")
  }
  off_diag <- phydist[row(phydist) != col(phydist)]
  if (any(off_diag < 0)) {
    stop("'phydist' off-diagonal values must be non-negative.")
  }
  invisible(NULL)
}

.validate_spatial_inputs <- function(shp, points = NULL, res = NULL) {
  if (is.null(shp)) {
    stop("'shp' must be provided.")
  }
  if (!is.null(points)) {
    if (!is.data.frame(points) && !is.matrix(points)) {
      stop("'points' must be a data.frame or matrix.")
    }
    if (ncol(points) != 3) {
      stop("'points' must have exactly 3 columns (species, longitude, latitude).")
    }
    if (anyNA(points[, 2]) || anyNA(points[, 3])) {
      stop("'points' must not contain NA values in coordinate columns.")
    }
  }
  if (!is.null(res)) {
    if (!is.numeric(res) || any(res <= 0)) {
      stop("'res' must be a positive numeric value.")
    }
  }
  invisible(NULL)
}

.check_incidence_phydist_alignment <- function(incidence, phydist) {
  inc_cols <- colnames(incidence)
  phy_names <- colnames(phydist)
  missing <- setdiff(inc_cols, phy_names)
  if (length(missing) > 0) {
    stop(paste0("The following incidence columns are not in phydist: ",
                paste(missing, collapse = ", ")))
  }
  invisible(NULL)
}

.validate_pairs <- function(pairs) {
  if (!inherits(pairs, "data.frame")) {
    stop("'pairs' must be a data.frame (see prepare_pair_data()).")
  }
  needed <- c("parasite", "focal", "target", "phydist", "suscept")
  missing <- setdiff(needed, colnames(pairs))
  if (length(missing) > 0) {
    stop(paste0("'pairs' is missing required column(s): ",
                paste(missing, collapse = ", "),
                ". Build it with prepare_pair_data()."))
  }
  if (nrow(pairs) == 0) {
    stop("'pairs' must have at least one row.")
  }
  if (!all(pairs$suscept %in% c(0, 1))) {
    stop("'pairs$suscept' must be binary (0/1).")
  }
  invisible(NULL)
}

# Combined validator for incidence + phydist inputs used by modeling functions.
# Validates each input individually and then checks alignment between them.
.validate_phylopred_inputs <- function(incidence, phydist) {
  .validate_incidence(incidence)
  .validate_phydist(phydist)
  .check_incidence_phydist_alignment(incidence, phydist)
  invisible(NULL)
}
