#' Build a perturbation target mask
#'
#' Return a logical matrix selecting the entries affected by a perturbation.
#'
#' @param model A `demovuln_model` object or a numeric square projection matrix.
#' @param target One of `"adult_survival"`, `"juvenile_survival"`, `"fecundity"`,
#'   `"all"`, or `"custom"`.
#' @param survival_affects_fecundity Logical. If `TRUE`, survival perturbations
#'   scale whole source-stage columns, including fecundity entries.
#' @param custom_mask Optional logical matrix used when `target = "custom"`.
#'
#' @return Logical matrix with the same dimensions as the projection matrix.
#' @export
build_target_mask <- function(model,
                              target,
                              survival_affects_fecundity = TRUE,
                              custom_mask = NULL) {
  model <- as_demovuln_model(model)
  target <- match.arg(
    target,
    choices = c("adult_survival", "juvenile_survival", "fecundity", "all", "custom")
  )

  n <- model$n_stages

  if (target == "custom") {
    if (is.null(custom_mask)) {
      stop("custom_mask must be supplied when target = 'custom'.", call. = FALSE)
    }
    custom_mask <- as.matrix(custom_mask)
    if (!is.logical(custom_mask) || !identical(dim(custom_mask), dim(model$matrix))) {
      stop("custom_mask must be a logical matrix with the same dimensions as A.", call. = FALSE)
    }
    return(custom_mask)
  }

  if (target == "fecundity") {
    return(model$fecundity_element_mask)
  }

  mask <- matrix(FALSE, nrow = n, ncol = n)

  if (target == "adult_survival") {
    mask[, model$adult_stage_mask] <- TRUE
  } else if (target == "juvenile_survival") {
    mask[, model$juvenile_stage_mask] <- TRUE
  } else if (target == "all") {
    mask[,] <- TRUE
  }

  if (!survival_affects_fecundity && target %in% c("adult_survival", "juvenile_survival")) {
    mask <- mask & !model$fecundity_element_mask
  }

  mask
}

#' Apply a proportional perturbation to a projection matrix
#'
#' @param model A `demovuln_model` object or a numeric square projection matrix.
#' @param target One of `"adult_survival"`, `"juvenile_survival"`, `"fecundity"`,
#'   `"all"`, or `"custom"`.
#' @param magnitude Proportional reduction applied to the selected entries.
#'   Must lie in `[0, 1]`.
#' @param survival_affects_fecundity Logical. If `TRUE`, survival perturbations
#'   scale whole source-stage columns, including fecundity entries.
#' @param custom_mask Optional logical matrix used when `target = "custom"`.
#'
#' @return Perturbed projection matrix.
#' @export
apply_perturbation <- function(model,
                               target,
                               magnitude,
                               survival_affects_fecundity = TRUE,
                               custom_mask = NULL) {
  model <- as_demovuln_model(model)
  target <- match.arg(
    target,
    choices = c("adult_survival", "juvenile_survival", "fecundity", "all", "custom")
  )
  validate_magnitude(magnitude)

  factor <- 1 - magnitude
  B <- model$matrix

  if (target == "all" && isTRUE(survival_affects_fecundity)) {
    B <- B * factor
    B[model$fecundity_element_mask] <- B[model$fecundity_element_mask] * factor
    return(B)
  }

  mask <- build_target_mask(
    model = model,
    target = target,
    survival_affects_fecundity = survival_affects_fecundity,
    custom_mask = custom_mask
  )

  B[mask] <- B[mask] * factor
  B
}

validate_magnitude <- function(magnitude) {
  if (!is.numeric(magnitude) || length(magnitude) != 1 || !is.finite(magnitude)) {
    stop("magnitude must be a single finite numeric value.", call. = FALSE)
  }
  if (magnitude < 0 || magnitude > 1) {
    stop("magnitude must lie in the interval [0, 1].", call. = FALSE)
  }
  invisible(NULL)
}
