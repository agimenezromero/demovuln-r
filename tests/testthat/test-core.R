test_that("matrix_population_model accepts valid projection matrices", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  expect_s3_class(model, "demovuln_model")
  expect_equal(model$n_stages, 2)
  expect_true(model$lambda_ > 0)
  expect_equal(sum(stable_stage_distribution(model)), 1, tolerance = 1e-10)
})

test_that("matrix_population_model rejects invalid matrices", {
  expect_error(matrix_population_model(matrix(1, nrow = 2, ncol = 3)))
  expect_error(matrix_population_model(matrix(c(-1, 0, 0, 1), nrow = 2)))
})

test_that("zero magnitude gives zero reduction", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  sim <- simulate_dynamics(
    model,
    target = "adult_survival",
    magnitude = 0,
    duration = 1,
    period = 2,
    t_max = 20,
    recovery_steps = 5
  )

  expect_lt(abs(sim$reduction), 1e-10)
})

test_that("positive magnitude can reduce population size", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  sim <- simulate_dynamics(
    model,
    target = "adult_survival",
    magnitude = 0.5,
    duration = 1,
    period = 2,
    t_max = 20,
    recovery_steps = 5
  )

  expect_gt(sim$reduction, 0)
})

test_that("adult survival perturbation scales adult source columns", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  B <- apply_perturbation(model, target = "adult_survival", magnitude = 0.5)

  expected <- A
  expected[, 2] <- expected[, 2] * 0.5

  expect_equal(B, expected)
})

test_that("grid vulnerability is the mean reduction", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  grid <- perturbation_grid(
    magnitudes = c(0, 0.5),
    durations = 1,
    periods = c(2, 3)
  )

  result <- run_grid(
    model,
    target = "adult_survival",
    grid = grid,
    t_max = 20,
    recovery_steps = 5
  )

  expect_equal(nrow(result$table), 4)
  expect_equal(result$vulnerability, mean(result$table$population_reduction), tolerance = 1e-10)
})

test_that("infeasible scenarios can be retained", {
  A <- matrix(c(0.0, 0.4, 2.0, 0.7), nrow = 2, byrow = FALSE)
  model <- matrix_population_model(A)

  grid <- perturbation_grid(
    magnitudes = 0.5,
    durations = 3,
    periods = 2
  )

  result <- run_grid(
    model,
    target = "adult_survival",
    grid = grid,
    t_max = 20,
    skip_infeasible = FALSE
  )

  expect_equal(nrow(result$table), 1)
  expect_false(result$table$feasible[1])
  expect_true(is.na(result$table$population_reduction[1]))
})
