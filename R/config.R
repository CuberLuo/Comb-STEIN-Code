if (!requireNamespace("Iso", quietly = TRUE)) {
  stop("Package 'Iso' is required. Install it with install.packages('Iso').")
}
if (!requireNamespace("rjags", quietly = TRUE)) {
  stop(
    "Package 'rjags' is required. Install JAGS first, then install.packages('rjags')."
  )
}

DEFAULT_SEED <- 2026
DEFAULT_NTRIAL <- 1000L

## Time-to-event settings used by TITE-STEIN (Sun et al., 2025).
## Times are expressed in days.  Delayed outcomes are opt-in so the original
## complete-outcome simulation remains the default.
DEFAULT_DELAY_CONTROL <- list(
  accrual = 10,
  tox_window = 30,
  eff_window = 90,
  suspension = 0.5,
  use_suspension = TRUE,
  accrual_random = FALSE,
  tox_distribution = "uniform",
  eff_distribution = "uniform"
)

## Design parameters (Section 3 of the manuscript)
phi0 <- 0.3
phi1 <- 0.5
psi1 <- 0.3
psi2 <- 0.8
Tthresh <- 0.8
Ethresh <- 0.9
w1 <- 0.33
w2 <- 1.09
cohortsize <- 2
COHORTS_PER_COMBINATION <- 2
RUNIN_COHORTS <- 4
EFF_STOP_N <- 6
EFF_STOP_CUT <- 0.7
EFF_STOP_PROB <- 0.9

## Bayesian quadratic-logistic efficacy model (Cai, Yuan, and Ji, 2014).
EFF_MCMC_BURNIN <- 1000L
EFF_MCMC_SAMPLES <- 2000L
EFF_MCMC_CHAINS <- 2L
EFF_MCMC_THIN <- 1L

## Derive the toxicity decision boundary from the acceptable (`phi0`) and
## unacceptable (`phi1`) toxicity probabilities.
calc.phi <- function(phi0, phi1) {
  log((1 - phi0) / (1 - phi1)) /
    log(phi1 * (1 - phi0) / (phi0 * (1 - phi1)))
}

## Derive the efficacy decision boundary from the unacceptable (`psi1`) and
## desirable (`psi2`) efficacy probabilities.
calc.psi <- function(psi1, psi2) {
  log((1 - psi1) / (1 - psi2)) /
    log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))
}

phi <- calc.phi(phi0, phi1)
psi <- calc.psi(psi1, psi2)
