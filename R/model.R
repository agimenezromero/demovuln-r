#' Matrix population model
#'
#' Create a matrix population model and define the demographic targets used by
#' perturbation functions.
#'
#' @param A Numeric square projection matrix. Columns are source stages at time
#'   `t` and rows are destination stages at time `t + 1`.
#' @param fecundity_mask Optional logical matrix with the same dimensions as `A`,
#'   identifying fecundity entries.
#' @param fecundity_rows Integer vector identifying rows interpreted as newborn
#'   or reproductive-output rows. Defaults to the first row.
#' @param adult_stages Optional integer vector with source-stage columns
#'   interpreted as adult or reproductive stages.
#' @param juvenile_stages Optional integer vector with source-stage columns
#'   interpreted as juvenile or pre-reproductive stages.
#' @param name Optional model or species label.
#'
#' @return An object of class `demovuln_model`.
#' @export
matrix_population_model <- function(A,
                                    fecundity_mask = NULL,
                                    fecundity_rows = 1L,
                                    adult_stages = NULL,
                                    juvenile_stages = NULL,
                                    name = NULL) {
  A <- as.matrix(A)

  if (!is.numeric(A)) {
    stop("A must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop("A must be a square projection matrix.", call. = FALSE)
  }
  if (any(!is.finite(A))) {
    stop("A must contain only finite values.", call. = FALSE)
  }
  if (any(A < 0)) {
    stop("A must be non-negative.", call. = FALSE)
  }

  n_stages <- nrow(A)

  if (is.null(fecundity_mask)) {
    fecundity_rows <- as.integer(fecundity_rows)
    validate_indices(fecundity_rows, n_stages, "fecundity_rows")

    fecundity_mask <- matrix(FALSE, nrow = n_stages, ncol = n_stages)
    fecundity_mask[fecundity_rows, ] <- A[fecundity_rows, , drop = FALSE] > 0
  } else {
    fecundity_mask <- as.matrix(fecundity_mask)
    if (!is.logical(fecundity_mask)) {
      stop("fecundity_mask must be a logical matrix.", call. = FALSE)
    }
    if (!identical(dim(fecundity_mask), dim(A))) {
      stop("fecundity_mask must have the same dimensions as A.", call. = FALSE)
    }
  }

  if (is.null(adult_stages)) {
    adult_stages <- which(colSums(fecundity_mask) > 0)
  } else {
    adult_stages <- as.integer(adult_stages)
    validate_indices(adult_stages, n_stages, "adult_stages")
  }

  if (is.null(juvenile_stages)) {
    juvenile_stages <- setdiff(seq_len(n_stages), adult_stages)
  } else {
    juvenile_stages <- as.integer(juvenile_stages)
    validate_indices(juvenile_stages, n_stages, "juvenile_stages")
  }

  adult_stage_mask <- rep(FALSE, n_stages)
  adult_stage_mask[adult_stages] <- TRUE

  juvenile_stage_mask <- rep(FALSE, n_stages)
  juvenile_stage_mask[juvenile_stages] <- TRUE

  out <- list(
    matrix = A,
    n_stages = n_stages,
    fecundity_element_mask = fecundity_mask,
    adult_stage_mask = adult_stage_mask,
    juvenile_stage_mask = juvenile_stage_mask,
    lambda_ = dominant_eigenvalue(A),
    name = name
  )

  class(out) <- "demovuln_model"
  out
}

#' Dominant eigenvalue of a projection matrix
#'
#' @param A Numeric square projection matrix.
#' @return The real part of the eigenvalue with largest modulus.
#' @export
dominant_eigenvalue <- function(A) {
  A <- as.matrix(A)
  eig <- eigen(A)
  idx <- which.max(Mod(eig$values))
  Re(eig$values[idx])
}

#' Stable stage distribution
#'
#' @param model A `demovuln_model` object or a numeric square projection matrix.
#' @return Numeric vector normalized to sum to one.
#' @export
stable_stage_distribution <- function(model) {
  A <- if (inherits(model, "demovuln_model")) model$matrix else as.matrix(model)
  eig <- eigen(A)
  idx <- which.max(Mod(eig$values))
  v <- Re(eig$vectors[, idx])

  if (sum(v) < 0) {
    v <- -v
  }
  if (any(v < 0)) {
    v <- abs(v)
  }
  if (sum(v) <= 0 || any(!is.finite(v))) {
    v <- rep(1 / length(v), length(v))
  } else {
    v <- v / sum(v)
  }

  as.numeric(v)
}

as_demovuln_model <- function(model) {
  if (inherits(model, "demovuln_model")) {
    return(model)
  }
  matrix_population_model(model)
}

validate_indices <- function(x, n, name) {
  if (length(x) == 0) {
    return(invisible(NULL))
  }
  if (any(is.na(x)) || any(x < 1) || any(x > n)) {
    stop(sprintf("%s must contain valid 1-based stage indices.", name), call. = FALSE)
  }
  invisible(NULL)
}
