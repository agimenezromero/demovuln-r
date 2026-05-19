#' Compute percent population reduction
#'
#' @param final_population Final population size under perturbed dynamics.
#' @param baseline_final_population Final population size under unperturbed
#'   baseline dynamics.
#'
#' @return Percent population reduction relative to the baseline.
#' @export
population_reduction <- function(final_population, baseline_final_population) {
  if (!is.numeric(baseline_final_population) || baseline_final_population <= 0) {
    stop("baseline_final_population must be positive.", call. = FALSE)
  }
  100 * (1 - final_population / baseline_final_population)
}

#' Simulate dynamics under a temporally structured perturbation
#'
#' @param model A `demovuln_model` object or a numeric square projection matrix.
#' @param target One of `"adult_survival"`, `"juvenile_survival"`, `"fecundity"`,
#'   `"all"`, or `"custom"`.
#' @param magnitude Proportional reduction applied to the selected entries.
#' @param duration Number of consecutive projection intervals during which each
#'   perturbation event is active.
#' @param period Number of projection intervals between perturbation onsets.
#' @param t_max Number of projection intervals in the perturbation-forcing window.
#' @param recovery_steps Number of additional unperturbed projection intervals
#'   after the forcing window.
#' @param start Projection interval at which the first perturbation event starts.
#'   The default `0` means that forcing can begin at the first projection step.
#' @param initial_state Optional initial population vector. If omitted, the
#'   stable stage distribution of the unperturbed model is used.
#' @param normalize_by_lambda Logical. If `TRUE`, baseline and perturbed matrices
#'   are divided by the dominant eigenvalue of the unperturbed projection matrix.
#' @param survival_affects_fecundity Logical. If `TRUE`, survival perturbations
#'   scale whole source-stage columns, including fecundity entries.
#' @param custom_mask Optional logical matrix used when `target = "custom"`.
#' @param return_stage_vectors Logical. If `TRUE`, return the full stage-vector
#'   trajectories.
#' @param force_during_recovery Logical. If `TRUE`, scheduled perturbations
#'   continue during the recovery window. The default is `FALSE`.
#'
#' @return An object of class `demovuln_simulation`.
#' @export
simulate_dynamics <- function(model,
                              target = "adult_survival",
                              magnitude,
                              duration,
                              period,
                              t_max,
                              recovery_steps = 0L,
                              start = 0L,
                              initial_state = NULL,
                              normalize_by_lambda = TRUE,
                              survival_affects_fecundity = TRUE,
                              custom_mask = NULL,
                              return_stage_vectors = FALSE,
                              force_during_recovery = FALSE) {
  model <- as_demovuln_model(model)
  target <- match.arg(
    target,
    choices = c("adult_survival", "juvenile_survival", "fecundity", "all", "custom")
  )

  validate_magnitude(magnitude)
  duration <- validate_nonnegative_integer(duration, "duration")
  period <- validate_positive_integer(period, "period")
  t_max <- validate_nonnegative_integer(t_max, "t_max")
  recovery_steps <- validate_nonnegative_integer(recovery_steps, "recovery_steps")
  start <- validate_nonnegative_integer(start, "start")

  if (duration > period) {
    stop("duration must be less than or equal to period.", call. = FALSE)
  }

  total_steps <- t_max + recovery_steps

  A_base <- model$matrix
  A_pert <- apply_perturbation(
    model = model,
    target = target,
    magnitude = magnitude,
    survival_affects_fecundity = survival_affects_fecundity,
    custom_mask = custom_mask
  )

  if (isTRUE(normalize_by_lambda)) {
    if (!is.finite(model$lambda_) || model$lambda_ == 0) {
      stop("The dominant eigenvalue must be finite and non-zero for normalization.", call. = FALSE)
    }
    A_base <- A_base / model$lambda_
    A_pert <- A_pert / model$lambda_
  }

  n <- model$n_stages
  if (is.null(initial_state)) {
    state <- stable_stage_distribution(model)
  } else {
    state <- as.numeric(initial_state)
    if (length(state) != n || any(!is.finite(state)) || any(state < 0)) {
      stop("initial_state must be a non-negative finite vector with one entry per stage.", call. = FALSE)
    }
    if (sum(state) <= 0) {
      stop("initial_state must have positive total abundance.", call. = FALSE)
    }
    state <- state / sum(state)
  }

  baseline_state <- state

  abundance <- numeric(total_steps + 1L)
  baseline_abundance <- numeric(total_steps + 1L)
  abundance[1L] <- sum(state)
  baseline_abundance[1L] <- sum(baseline_state)

  if (isTRUE(return_stage_vectors)) {
    stage_vectors <- matrix(NA_real_, nrow = total_steps + 1L, ncol = n)
    baseline_stage_vectors <- matrix(NA_real_, nrow = total_steps + 1L, ncol = n)
    stage_vectors[1L, ] <- state
    baseline_stage_vectors[1L, ] <- baseline_state
  } else {
    stage_vectors <- NULL
    baseline_stage_vectors <- NULL
  }

  if (total_steps > 0) {
    for (step in seq_len(total_steps)) {
      t <- step - 1L
      active <- is_perturbation_active(
        t = t,
        duration = duration,
        period = period,
        start = start,
        t_max = t_max,
        force_during_recovery = force_during_recovery
      )

      state <- if (active) A_pert %*% state else A_base %*% state
      baseline_state <- A_base %*% baseline_state

      abundance[step + 1L] <- sum(state)
      baseline_abundance[step + 1L] <- sum(baseline_state)

      if (isTRUE(return_stage_vectors)) {
        stage_vectors[step + 1L, ] <- as.numeric(state)
        baseline_stage_vectors[step + 1L, ] <- as.numeric(baseline_state)
      }
    }
  }

  final_population <- abundance[total_steps + 1L]
  baseline_final_population <- baseline_abundance[total_steps + 1L]

  out <- list(
    abundance = abundance,
    baseline_abundance = baseline_abundance,
    stage_vectors = stage_vectors,
    baseline_stage_vectors = baseline_stage_vectors,
    reduction = population_reduction(final_population, baseline_final_population),
    final_population = final_population,
    baseline_final_population = baseline_final_population,
    magnitude = magnitude,
    duration = duration,
    period = period,
    target = target
  )

  class(out) <- "demovuln_simulation"
  out
}

is_perturbation_active <- function(t,
                                   duration,
                                   period,
                                   start,
                                   t_max,
                                   force_during_recovery) {
  if (duration <= 0) {
    return(FALSE)
  }
  if (t < start) {
    return(FALSE)
  }
  if (!force_during_recovery && t >= t_max) {
    return(FALSE)
  }
  ((t - start) %% period) < duration
}

validate_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x != as.integer(x)) {
    stop(sprintf("%s must be a non-negative integer.", name), call. = FALSE)
  }
  as.integer(x)
}

validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0 || x != as.integer(x)) {
    stop(sprintf("%s must be a positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}
