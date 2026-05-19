# demovuln for R

`demovuln` is a lightweight R implementation of the demographic vulnerability framework for matrix population models.

The main Python package is available at:

- PyPI: https://pypi.org/project/demovuln/
- Documentation: https://demovuln.readthedocs.io/en/latest/
- GitHub: https://github.com/agimenezromero/demovuln

This R package is intended as a companion implementation that can be installed directly from GitHub.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("agimenezromero/demovuln-r")
```

## Basic usage

```r
library(demovuln)

A <- matrix(
  c(0.0, 0.4,
    2.0, 0.7),
  nrow = 2,
  byrow = FALSE
)

model <- matrix_population_model(A)

sim <- simulate_dynamics(
  model,
  target = "adult_survival",
  magnitude = 0.25,
  duration = 1,
  period = 3,
  t_max = 50,
  recovery_steps = 10
)

sim$reduction
sim$abundance
```

The projection matrix follows the standard matrix-population-model convention: columns are source stages at time `t`, and rows are destination stages at time `t + 1`.

The example above produces the matrix:

```text
     [,1] [,2]
[1,]  0.0  2.0
[2,]  0.4  0.7
```

## Perturbation-grid analysis

```r
grid <- perturbation_grid(
  magnitudes = seq(0, 1, length.out = 11),
  durations = c(0, 1, 2, 3),
  periods = c(1, 2, 3, 5, 10)
)

out <- run_grid(
  model,
  target = "adult_survival",
  grid = grid,
  t_max = 50,
  recovery_steps = 10
)

out$vulnerability
head(out$table)
```

## Demographic targets

The package supports perturbations to:

- `adult_survival`
- `juvenile_survival`
- `fecundity`
- `all`
- `custom`

By default, adult stages are inferred as source-stage columns with at least one fecundity entry, and juvenile stages are inferred as the remaining source-stage columns. These definitions can be specified explicitly:

```r
model <- matrix_population_model(
  A,
  adult_stages = 2,
  juvenile_stages = 1
)
```

Custom perturbation targets can be defined with logical masks:

```r
custom_mask <- matrix(
  c(FALSE, TRUE,
    FALSE, FALSE),
  nrow = 2,
  byrow = FALSE
)

sim <- simulate_dynamics(
  model,
  target = "custom",
  custom_mask = custom_mask,
  magnitude = 0.5,
  duration = 1,
  period = 3,
  t_max = 50
)
```

## Conceptual summary

For a given perturbation regime, population reduction is computed as:

```text
rho = 100 * (1 - N_perturbed(T) / N_baseline(T))
```

where `N_perturbed(T)` is the final population size under perturbed dynamics and `N_baseline(T)` is the final population size under the unperturbed baseline.

Integrated vulnerability is the mean population reduction across the simulated perturbation space:

```text
Phi = mean(rho)
```

## Example vignettes

Example vignettes are available in the `vignettes/` directory.

To install the package and build the vignettes locally:

```r
remotes::install_github(
  "agimenezromero/demovuln-r",
  build_vignettes = TRUE
)

browseVignettes("demovuln")

## Development checks

After cloning the repository, run:

```r
devtools::test()
devtools::check()
```

## License

This package is distributed under the MIT License.
