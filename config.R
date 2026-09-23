## Global design configuration for Comb-STEIN / RunIn-STEIN.
## Everything here is read at source() time by the simulation and plotting scripts.

DEFAULT_SEED <- 2026
DEFAULT_NTRIAL <- 1000L  # replicated trials per scenario

# NULL -> auto: available cores minus one, never fewer than one worker;
# overridable via run.scenario(workers=) or --workers=
DEFAULT_WORKERS <- NULL

## All intermediate diagnostics go through these two wrappers, so --quiet
## (option comb.stein.quiet) silences them at once.
debug.cat <- function(...) {
  if (!isTRUE(getOption("comb.stein.quiet", FALSE))) cat(...)
}

debug.print <- function(x, ...) {
  if (!isTRUE(getOption("comb.stein.quiet", FALSE))) print(x, ...)
}

## TITE-STEIN time-to-event setup (Sun et al., 2025). All times in days.
## Delayed outcomes are opt-in, so complete-data simulation stays the default.
DEFAULT_DELAY_CONTROL <- list(
  accrual = 10,          # mean days between successive cohort arrivals
  tox_window = 30,       # days of follow-up before toxicity is scored
  eff_window = 90,       # days of follow-up before efficacy is scored
  suspension = 0.5,      # fraction of follow-up suspending new accrual
  use_suspension = TRUE, # FALSE disables accrual suspension entirely
  accrual_random = FALSE,# FALSE -> fixed spacing; TRUE -> random arrival times
  tox_distribution = "uniform",
  eff_distribution = "uniform"
)

## Toxicity anchors (phi) and efficacy anchors (psi) defining the boundaries
phi1 <- 0.3  # highest toxicity rate still deemed acceptable
phi2 <- 0.5  # toxicity rate deemed unacceptably high
psi1 <- 0.3  # efficacy rate deemed unacceptably low
psi2 <- 0.8  # efficacy rate targeted as desirable

## Posterior probability cut-offs for the safety / futility stopping rules
Tthresh <- 0.8  # Pr(toxicity > phi) above which a cell is over-toxic
Ethresh <- 0.9  # Pr(efficacy < psi1) above which a cell is futile

## Utility: q - w1*p - w2*p*I(p > phi), so phi is also the toxicity penalty knot
w1 <- 0.33
w2 <- 1.09
# Utility values within this distance of the global maximum are OBDCs.
OBDC_UTILITY_TOLERANCE <- 0.03

cohortsize <- 3                  # patients enrolled per cohort
COHORTS_PER_COMBINATION <- 2     # cohorts allowed at each combination before moving on
RUNIN_COHORTS <- 4               # number of cohorts following the run-in path

## Early stopping for efficacy: needs >= EFF_STOP_N treated patients showing
## Pr(efficacy > EFF_STOP_CUT) > EFF_STOP_PROB
EFF_STOP_N <- 6
EFF_STOP_CUT <- 0.7
EFF_STOP_PROB <- 0.9

## Bayesian quadratic logistic efficacy model (Cai, Yuan and Ji, 2014).
EFF_MCMC_BURNIN <- 1000L
EFF_MCMC_SAMPLES <- 2000L  # posterior draws kept in total, summed over chains
EFF_MCMC_CHAINS <- 2L
EFF_MCMC_THIN <- 1L

## Toxicity decision boundary implied by the acceptable (phi1) and
## unacceptable (phi2) toxicity rates.
calc.phi <- function(phi1, phi2) {
  log((1 - phi1) / (1 - phi2)) /
    log(phi2 * (1 - phi1) / (phi1 * (1 - phi2)))
}

## Efficacy decision boundary implied by the unacceptably low (psi1) and
## ideal (psi2) efficacy rates.
calc.psi <- function(psi1, psi2) {
  log((1 - psi1) / (1 - psi2)) /
    log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))
}

phi <- calc.phi(phi1, phi2)
psi <- calc.psi(psi1, psi2)
