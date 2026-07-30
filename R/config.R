if (!requireNamespace("Iso", quietly = TRUE)) {
  stop("Package 'Iso' is required. Install it with install.packages('Iso').")
}

DEFAULT_SEED <- 2026

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

calc.phi <- function(phi0, phi1) {
  log((1 - phi0) / (1 - phi1)) /
    log(phi1 * (1 - phi0) / (phi0 * (1 - phi1)))
}

calc.psi <- function(psi1, psi2) {
  log((1 - psi1) / (1 - psi2)) /
    log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))
}

phi <- calc.phi(phi0, phi1)
psi <- calc.psi(psi1, psi2)
