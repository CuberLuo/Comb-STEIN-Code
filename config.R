if (!requireNamespace("Iso", quietly = TRUE)) {
  stop("Package 'Iso' is required. Install it with install.packages('Iso').")
}
if (!requireNamespace("rjags", quietly = TRUE)) {
  stop(
    "Package 'rjags' is required. Install JAGS first, then install.packages('rjags')."
  )
}

DEFAULT_SEED <- 2026
DEFAULT_NTRIAL <- 1L#试验次数

## TITE-STEIN 使用的时间到事件设置（Sun 等，2025）。
## 时间单位为天。延迟结局需要显式开启，因此原有的完整结局模拟仍是默认行为。
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

## 设计参数（论文第 3 节）。
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

## 贝叶斯二次逻辑回归疗效模型（Cai、Yuan 和 Ji，2014）。
EFF_MCMC_BURNIN <- 1000L
EFF_MCMC_SAMPLES <- 2000L
EFF_MCMC_CHAINS <- 2L
EFF_MCMC_THIN <- 1L

## 根据可接受（`phi0`）和不可接受（`phi1`）的毒性概率推导毒性决策边界。
calc.phi <- function(phi0, phi1) {
  log((1 - phi0) / (1 - phi1)) /
    log(phi1 * (1 - phi0) / (phi0 * (1 - phi1)))
}

## 根据不可接受（`psi1`）和理想（`psi2`）的疗效概率推导疗效决策边界。
calc.psi <- function(psi1, psi2) {
  log((1 - psi1) / (1 - psi2)) /
    log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))
}

phi <- calc.phi(phi0, phi1)
psi <- calc.psi(psi1, psi2)
