## 简化的逐行 Waterfall 对照方法。

## 使用完整或延迟结局模拟一次逐行 Waterfall 试验。
## 返回分配计数、结局、选定组合、停止指标和试验时长。
simulate.trial.waterfall <- function(
  pT,
  pE,
  method = c("Waterfall-L", "Waterfall-H"),
  delayed = FALSE,
  delay.control = NULL
) {
  if (isTRUE(delayed)) {
    return(simulate.trial.waterfall.tite(
      pT,
      pE,
      method,
      normalize.delay.control(delay.control)
    ))
  }
  delay.control <- normalize.delay.control(delay.control)
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  tox_monitor_triggered <- FALSE
  time_current <- 0
  duration <- 0

  row.order <- if (method == "Waterfall-L") seq_len(I) else I:1
  start.col <- if (method == "Waterfall-L") 1 else J
  cohorts.per.row <- rep(ncohort %/% I, I)
  remainder <- ncohort %% I
  if (remainder > 0) {
    cohorts.per.row[seq_len(remainder)] <-
      cohorts.per.row[seq_len(remainder)] + 1
  }

  cohort.count <- 0
  current.col <- start.col
  outer.stop <- FALSE

  for (r in seq_along(row.order)) {
    i <- row.order[r]
    j <- current.col
    n.this.row <- cohorts.per.row[r]
    if (n.this.row == 0) next

    for (k in seq_len(n.this.row)) {
      cohort.count <- cohort.count + 1
      tox_outcome <- as.integer(runif(cohortsize) < pT[i, j])
      eff_outcome <- as.integer(runif(cohortsize) < pE[i, j])
      wT <- sum(tox_outcome)
      wE <- sum(eff_outcome)
      x_mat[i, j] <- x_mat[i, j] + wT
      y_mat[i, j] <- y_mat[i, j] + wE
      n_mat[i, j] <- n_mat[i, j] + cohortsize
      timing <- complete.cohort.timing(
        time_current, tox_outcome, eff_outcome, delay.control
      )
      time_current <- timing$next_time
      duration <- timing$duration

      flags <- monitor.flags(x_mat, y_mat, n_mat)
      tox_monitor_triggered <- tox_monitor_triggered || any(flags$overtox)
      admissible <- update.admissible(x_mat, y_mat, n_mat, admissible)

      if (!any(admissible)) {
        earlystop_tox <- TRUE
        outer.stop <- TRUE
        break
      }
      if (!any(admissible[i, ])) break

      if (
        n_mat[i, j] >= EFF_STOP_N &&
          pr.eff.gt(y_mat[i, j], n_mat[i, j], EFF_STOP_CUT) > EFF_STOP_PROB
      ) {
        earlystop_eff <- TRUE
        outer.stop <- TRUE
        break
      }
      if (cohort.count >= ncohort || k == n.this.row) break

      p_ij_hat <- x_mat[i, j] / n_mat[i, j]
      candidates <- list()
      if (p_ij_hat > phi) {
        if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
      } else {
        if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
        candidates <- c(candidates, list(c(i, j)))
        if (j < J) candidates <- c(candidates, list(c(i, j + 1)))
      }
      candidates <- Filter(
        function(ij) admissible[ij[1], ij[2]],
        candidates
      )

      if (length(candidates) == 0) {
        row.admissible <- which(admissible[i, ])
        nearest <- row.admissible[
          abs(row.admissible - j) == min(abs(row.admissible - j))
        ]
        candidates <- lapply(nearest, function(jj) c(i, jj))
      }
      j <- pick.best.by.eff(candidates, y_mat, n_mat)[2]
    }

    current.col <- j
    if (outer.stop || cohort.count >= ncohort) break
  }

  selected <- if (earlystop_tox) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(x_mat, y_mat, n_mat, admissible)
  }

  list(
    sel = selected,
    x = x_mat,
    y = y_mat,
    n = n_mat,
    earlystop_tox = earlystop_tox,
    earlystop_eff = earlystop_eff,
    tox_monitor_triggered = tox_monitor_triggered,
    delayed = FALSE,
    duration = duration
  )
}

## 运行一次 Waterfall 试验的时间到事件版本。
## `delay.control` 必须已经由 `normalize.delay.control()` 规范化。
simulate.trial.waterfall.tite <- function(pT, pE, method, delay.control) {
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  state <- new.tite.state()
  admissible <- matrix(TRUE, I, J)
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  tox_monitor_triggered <- FALSE
  time_current <- 0

  row.order <- if (method == "Waterfall-L") seq_len(I) else I:1
  start.col <- if (method == "Waterfall-L") 1 else J
  cohorts.per.row <- rep(ncohort %/% I, I)
  remainder <- ncohort %% I
  if (remainder > 0L) {
    cohorts.per.row[seq_len(remainder)] <-
      cohorts.per.row[seq_len(remainder)] + 1
  }

  cohort.count <- 0L
  current.col <- start.col
  outer.stop <- FALSE

  for (r in seq_along(row.order)) {
    i <- row.order[r]
    j <- current.col
    n.this.row <- cohorts.per.row[r]
    if (n.this.row == 0L) next

    for (k in seq_len(n.this.row)) {
      cohort.count <- cohort.count + 1L
      current <- c(i, j)
      state <- tite.enroll.cohort(
        state,
        cohort.count,
        current,
        pT[i, j],
        pE[i, j],
        time_current,
        cohortsize,
        delay.control
      )
      decision_time <- tite.next.decision.time(
        state,
        current,
        time_current,
        cohortsize,
        delay.control
      )
      interim <- tite.interim(state, decision_time, I, J, delay.control)
      time_current <- decision_time

      flags <- monitor.flags(
        interim$x,
        interim$y,
        interim$ess_t,
        interim$ess_e
      )
      tox_monitor_triggered <- tox_monitor_triggered || any(flags$overtox)
      admissible <- update.admissible(
        interim$x,
        interim$y,
        interim$ess_t,
        admissible,
        n_eff = interim$ess_e
      )
      if (!any(admissible)) {
        earlystop_tox <- TRUE
        outer.stop <- TRUE
        break
      }
      if (!any(admissible[i, ])) break

      if (
        interim$n[i, j] >= EFF_STOP_N &&
          pr.eff.gt(
            interim$y[i, j],
            interim$ess_e[i, j],
            EFF_STOP_CUT
          ) > EFF_STOP_PROB
      ) {
        earlystop_eff <- TRUE
        outer.stop <- TRUE
        break
      }
      if (cohort.count >= ncohort || k == n.this.row) break

      p_ij_hat <- interim$x[i, j] / interim$ess_t[i, j]
      candidates <- list()
      if (p_ij_hat > phi) {
        if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
      } else {
        if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
        candidates <- c(candidates, list(c(i, j)))
        if (j < J) candidates <- c(candidates, list(c(i, j + 1)))
      }
      candidates <- Filter(
        function(ij) admissible[ij[1], ij[2]],
        candidates
      )
      if (length(candidates) == 0L) {
        row.admissible <- which(admissible[i, ])
        nearest <- row.admissible[
          abs(row.admissible - j) == min(abs(row.admissible - j))
        ]
        candidates <- lapply(nearest, function(jj) c(i, jj))
      }
      j <- pick.best.by.eff(candidates, interim$y, interim$ess_e)[2]
    }

    current.col <- j
    if (outer.stop || cohort.count >= ncohort) break
  }

  final <- tite.final.counts(state, I, J)
  selected <- if (earlystop_tox) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(final$x, final$y, final$n, admissible)
  }

  list(
    sel = selected,
    x = final$x,
    y = final$y,
    n = final$n,
    earlystop_tox = earlystop_tox,
    earlystop_eff = earlystop_eff,
    tox_monitor_triggered = tox_monitor_triggered,
    delayed = TRUE,
    duration = final$duration
  )
}
