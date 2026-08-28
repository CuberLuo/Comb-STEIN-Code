## 使用完整或延迟数据模拟一次 Comb-STEIN 试验。
## 返回分配计数、结局、选定组合、停止指标和试验时长。
simulate.trial.comb.stein <- function(
  pT,
  pE,
  method = c("CombSTEIN-L", "CombSTEIN-H"),
  delayed = FALSE,
  delay.control = NULL
) {
  method <- match.arg(method)
  ## 延迟结局模式使用单独的 TITE 实现；这里先完成分流。
  if (isTRUE(delayed)) {
    return(simulate.trial.comb.stein.tite(
      pT,
      pE,
      method,
      normalize.delay.control(delay.control)
    ))
  }
  delay.control <- normalize.delay.control(delay.control)
  ## I、J 是剂量网格的行数和列数。
  ## 每个剂量组合预计入组 COHORTS_PER_COMBINATION 个队列，
  ## 因此总队列数由剂量组合数量决定。
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  ## x_mat：毒性事件数；y_mat：疗效事件数；n_mat：入组人数。
  ## admissible 表示当前仍然允许使用的剂量组合。
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)

  ## CombSTEIN-H 从最高剂量开始，CombSTEIN-L 从最低剂量开始。
  ## current 保存当前剂量组合的行、列坐标。
  current <- if (method == "CombSTEIN-H") c(I, J) else c(1, 1)
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  earlystop_no_candidate <- FALSE
  time_current <- 0
  duration <- 0

  for (cohort in seq_len(ncohort)) {
    ## 每次循环模拟一个队列，并在队列结束后进行一次剂量决策。
    debug.cat("current:", current, "\n")
    i <- current[1]
    j <- current[2]

    ## 根据当前剂量组合真实的毒性率和疗效率，
    ## 为队列中的每位患者随机生成毒性和疗效结局。
    tox_outcome <- as.integer(runif(cohortsize) < pT[i, j])
    eff_outcome <- as.integer(runif(cohortsize) < pE[i, j])
    # debug.cat("tox_outcome:",tox_outcome,"\n")
    # debug.cat("eff_outcome:",eff_outcome,"\n")
    wT <- sum(tox_outcome)#当前cohort出现毒性的人数
    wE <- sum(eff_outcome)#当前cohort出现疗效的人数
    debug.cat("wT:", wT, "\n")
    debug.cat("wE:", wE, "\n")
    x_mat[i, j] <- x_mat[i, j] + wT
    y_mat[i, j] <- y_mat[i, j] + wE
    n_mat[i, j] <- n_mat[i, j] + cohortsize
    # debug.cat("x_mat:","\n")
    # debug.print(x_mat)
    # debug.cat("y_mat:","\n")
    # debug.print(y_mat)
    # debug.cat("n_mat:","\n")
    # debug.print(n_mat)
    timing <- complete.cohort.timing(
      time_current, tox_outcome, eff_outcome, delay.control
    )
    time_current <- timing$next_time#完成评估且下一队列可以开始入组的最早时间
    duration <- timing$duration#当前已入组患者完成全部评估的时间

    ## 更新可接受剂量集合：超毒性组合及其右上方向可能被排除，
    ## 疗效不足的组合会被排除。
    admissible <- update.admissible(x_mat, y_mat, n_mat, admissible)
    # cat("admissible:\n")
    # print(admissible)
    if (!any(admissible)) {
      ## 如果没有任何可接受剂量，试验因毒性规则停止。
      earlystop_tox <- TRUE
      break
    }

    ## 当前剂量组合达到最少疗效观察人数，且疗效超过停止阈值时，
    ## 试验因疗效充分而提前停止。
    if (
      n_mat[i, j] >= EFF_STOP_N &&
        pr.eff.gt(y_mat[i, j], n_mat[i, j], EFF_STOP_CUT) > EFF_STOP_PROB
    ) {
      earlystop_eff <- TRUE
      break
    }

    if (cohort == ncohort) break

    ## Comb-STEIN 每次都在当前剂量的局部候选集合中选择下一剂量。
    p_ij_hat <- x_mat[i, j] / n_mat[i, j]
    # 筛选1：根据当前毒性概率得到 2剂量/5剂量 局部候选集
    candidates <- alt.set.C(i, j, I, J, p_ij_hat)
    # 筛选2：局部候选集∩可接受集
    candidates <- Filter(
      function(ij) admissible[ij[1], ij[2]],
      candidates
    )
    debug.cat("candidates:\n")
    debug.print(candidates)
    if (length(candidates) == 0) {
      if (admissible[i, j]) {
        candidates <- list(c(i, j))
      } else {
        earlystop_no_candidate <- TRUE
        break
      }
    }
    current <- pick.best.by.eff(candidates, y_mat, n_mat)
  }

  ## 只有正常完成试验时才进行最终剂量选择；
  ## 因毒性或无候选组合停止时，不报告最终选择。
  selected <- if (earlystop_tox || earlystop_no_candidate) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(x_mat, y_mat, n_mat, admissible)
  }

  ## 返回一次试验的全部信息，供 run.scenario() 汇总多个重复试验。
  list(
    sel = selected,
    x = x_mat,
    y = y_mat,
    n = n_mat,
    earlystop_tox = earlystop_tox,
    earlystop_eff = earlystop_eff,
    earlystop_no_candidate = earlystop_no_candidate,
    delayed = FALSE,
    duration = duration
  )
}

## 运行一次 STEIN 试验的时间到事件版本。
## `delay.control` 必须已经由 `normalize.delay.control()` 规范化。
simulate.trial.comb.stein.tite <- function(pT, pE, method, delay.control) {
  ## TITE 版本与完整结局版本使用相同的剂量移动规则，
  ## 但每次决策前先根据随访进度计算有效样本量。
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  ## state 保存每位患者的入组时间、结局观察时间和当前状态。
  state <- new.tite.state()
  admissible <- matrix(TRUE, I, J)
  ## H 版本从最高剂量开始，L 版本从最低剂量开始。
  current <- if (method == "CombSTEIN-H") c(I, J) else c(1, 1)
  time_current <- 0
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  earlystop_no_candidate <- FALSE

  for (cohort in seq_len(ncohort)) {
    ## 入组当前队列，并按照延迟结局机制推进试验时间。
    i <- current[1]
    j <- current[2]
    state <- tite.enroll.cohort(
      state,
      cohort,
      current,
      pT[i, j],
      pE[i, j],
      time_current,
      cohortsize,
      delay.control
    )
    ## 找到当前队列完成后可以进行下一次决策的时间点。
    decision_time <- tite.next.decision.time(
      state,
      current,
      time_current,
      cohortsize,
      delay.control
    )
    ## 在决策时刻计算每个剂量组合的有效毒性/疗效样本量和事件数。
    interim <- tite.interim(state, decision_time, I, J, delay.control)
    time_current <- decision_time

    ## 用有效样本量执行与完整结局版本相同的监测和可接受性更新。
    admissible <- update.admissible(
      interim$x,
      interim$y,
      interim$ess_t,
      admissible,
      n_eff = interim$ess_e
    )
    if (!any(admissible)) {
      ## 所有组合都不可接受，因毒性停止。
      earlystop_tox <- TRUE
      break
    }

    ## 如果当前组合的有效疗效样本量和疗效证据达到停止条件，
    ## 则因疗效充分而停止。
    if (
      interim$n[i, j] >= EFF_STOP_N &&
        pr.eff.gt(
          interim$y[i, j],
          interim$ess_e[i, j],
          EFF_STOP_CUT
        ) > EFF_STOP_PROB
    ) {
      earlystop_eff <- TRUE
      break
    }
    if (cohort == ncohort) break

    ## 用有效毒性样本量估计当前组合的毒性率，
    ## 决定下一步允许探索的方向。
    p_ij_hat <- interim$x[i, j] / interim$ess_t[i, j]
    ## 从当前组合的局部邻域中选择后验疗效最优的组合。
    candidates <- alt.set.C(i, j, I, J, p_ij_hat)
    candidates <- Filter(
      function(ij) admissible[ij[1], ij[2]],
      candidates
    )
    if (length(candidates) == 0L) {
      if (admissible[i, j]) {
        candidates <- list(c(i, j))
      } else {
        earlystop_no_candidate <- TRUE
        break
      }
    }
    current <- pick.best.by.eff(candidates, interim$y, interim$ess_e)
  }

  ## 试验结束后，将所有患者的最终随访结果整理成完整计数矩阵。
  final <- tite.final.counts(state, I, J)
  ## 与完整结局版本一致：异常停止时不报告最终选择，
  ## 正常结束时根据最终完整数据选择 utility 最大的可接受组合。
  selected <- if (earlystop_tox || earlystop_no_candidate) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(final$x, final$y, final$n, admissible)
  }

  ## 返回最终计数、选择结果、停止原因和试验持续时间。
  list(
    sel = selected,
    x = final$x,
    y = final$y,
    n = final$n,
    earlystop_tox = earlystop_tox,
    earlystop_eff = earlystop_eff,
    earlystop_no_candidate = earlystop_no_candidate,
    delayed = TRUE,
    duration = final$duration
  )
}

