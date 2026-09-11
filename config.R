DEFAULT_SEED <- 2026
DEFAULT_NTRIAL <- 1000L#试验次数

## 并行 worker 数。NULL 表示自动使用可用 CPU 核数减一，
## 至少保留一个 worker；可由 run.scenario(workers=) 或 --workers= 覆盖。
DEFAULT_WORKERS <- NULL

## 中间调试输出统一通过这两个函数控制。
## 默认保留调试信息；命令行添加 --quiet 后关闭。
debug.cat <- function(...) {
  if (!isTRUE(getOption("comb.stein.quiet", FALSE))) cat(...)
}

debug.print <- function(x, ...) {
  if (!isTRUE(getOption("comb.stein.quiet", FALSE))) print(x, ...)
}

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

## 设计参数
phi1 <- 0.3
phi2 <- 0.5
psi1 <- 0.3
psi2 <- 0.8
Tthresh <- 0.8
Ethresh <- 0.9
w1 <- 0.33
w2 <- 1.09
cohortsize <- 3
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

## 根据可接受（`phi1`）和不可接受（`phi2`）的毒性概率推导毒性决策边界。
calc.phi <- function(phi1, phi2) {
  log((1 - phi1) / (1 - phi2)) /
    log(phi2 * (1 - phi1) / (phi1 * (1 - phi2)))
}

## 根据不可接受（`psi1`）和理想（`psi2`）的疗效概率推导疗效决策边界。
calc.psi <- function(psi1, psi2) {
  log((1 - psi1) / (1 - psi2)) /
    log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))
}

phi <- calc.phi(phi1, phi2)
psi <- calc.psi(psi1, psi2)
