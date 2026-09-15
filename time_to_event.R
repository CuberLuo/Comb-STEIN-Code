## 时间到事件数据生成和有效样本量计算。
##
## 本模块遵循 TITE-STEIN 论文及其参考实现：
##   * 事件时间在各自的评估窗口内服从均匀分布；
##   * 尚未发生事件的待定结局按随访时间/评估窗口计入有效样本量；
##   * 已观测事件或已完成的无事件评估计为 1；
##   * 当前组合的毒性和疗效结局达到要求的确认比例后，才能继续入组。

## 将用户延迟设置与默认值合并，返回完整的延迟控制列表。
normalize.delay.control <- function(control = NULL) {
  if (is.null(control)) control <- list()
  out <- utils::modifyList(DEFAULT_DELAY_CONTROL, control)
  for (field in c("use_suspension", "accrual_random")) {
    out[[field]] <- as.logical(out[[field]])
  }
  for (field in c("tox_distribution", "eff_distribution")) {
    out[[field]] <- tolower(as.character(out[[field]])[1])
  }
  out
}

## 推进完整结局（non-TITE）设计的一个队列时间。
## 参考 TITE-STEIN 实现中 suspension = 1 的定义：必须等当前队列
## 所有毒性和疗效结局均已确认，才允许对下一队列作出分配决策。
complete.cohort.timing <- function(
  time_current,
  tox_outcome,
  eff_outcome,
  control
) {
  cohort_size <- length(tox_outcome)
  if (length(eff_outcome) != cohort_size) {
    stop("tox_outcome and eff_outcome must have the same length")
  }

  if (control$accrual_random) {
    enroll <- time_current + cumsum(
      stats::runif(cohort_size, 0, 2 * control$accrual)
    ) + 1
  } else {
    enroll <- time_current + seq.int(0, cohort_size - 1L) *
      control$accrual + 1
  }

  ## Complete 模式只使用已确认的二元结局，不额外抽取事件发生时刻。
  ## 这样不会改变原有完整结局模拟的随机数序列；持续时间按评估窗结束定义。
  tox_end <- enroll + control$tox_window
  eff_end <- enroll + control$eff_window

  next_accrual <- if (control$accrual_random) {
    max(enroll) + stats::runif(1L, 0, 2 * control$accrual) + 1
  } else {
    time_current + cohort_size * control$accrual + 1
  }

  list(
    next_time = max(next_accrual, tox_end, eff_end),
    duration = max(tox_end, eff_end)
  )
}

## 创建时间到事件试验的空患者级状态表。
new.tite.state <- function() {
  data.frame(
    cohort = integer(),
    dose_i = integer(),
    dose_j = integer(),
    enroll = numeric(),
    tox_end = numeric(),
    eff_end = numeric(),
    tox_outcome = integer(),
    eff_outcome = integer(),
    tox_event = numeric(),
    eff_event = numeric()
  )
}

## 入组一个队列，并生成其毒性/疗效结局及事件时间。
## 返回追加了每位新入组患者一行记录的 `state`。
tite.enroll.cohort <- function(
  state,
  cohort,
  dose,
  p_tox,
  p_eff,
  time_current,
  cohort_size,
  control
) {
  if (control$accrual_random) {
    enroll <- time_current + cumsum(
      stats::runif(cohort_size, 0, 2 * control$accrual)
    ) + 1
  } else {
    enroll <- time_current + seq.int(0, cohort_size - 1L) *
      control$accrual + 1
  }

  tox_outcome <- as.integer(stats::runif(cohort_size) < p_tox)
  eff_outcome <- as.integer(stats::runif(cohort_size) < p_eff)
  tox_event <- rep(Inf, cohort_size)
  eff_event <- rep(Inf, cohort_size)
  has_tox <- tox_outcome == 1L
  has_eff <- eff_outcome == 1L
  tox_event[has_tox] <- enroll[has_tox] + stats::runif(
    sum(has_tox), 0, control$tox_window
  )
  eff_event[has_eff] <- enroll[has_eff] + stats::runif(
    sum(has_eff), 0, control$eff_window
  )

  new_rows <- data.frame(
    cohort = rep.int(as.integer(cohort), cohort_size),
    dose_i = rep.int(as.integer(dose[1]), cohort_size),
    dose_j = rep.int(as.integer(dose[2]), cohort_size),
    enroll = enroll,
    tox_end = enroll + control$tox_window,
    eff_end = enroll + control$eff_window,
    tox_outcome = tox_outcome,
    eff_outcome = eff_outcome,
    tox_event = tox_event,
    eff_event = eff_event
  )
  rbind(state, new_rows)
}

## 确定最近一个队列之后允许决策的最早时间。
## 启用暂停时，当前组合的两个终点都必须达到要求的确认比例后才能继续入组。
tite.next.decision.time <- function(
  state,
  current,
  time_current,
  cohort_size,
  control
) {
  if (control$accrual_random) {
    base_time <- max(state$enroll) +
      stats::runif(1L, 0, 2 * control$accrual) + 1
  } else {
    base_time <- time_current + cohort_size * control$accrual + 1
  }
  if (!control$use_suspension) return(base_time)

  at_current <- state$dose_i == current[1] & state$dose_j == current[2]
  current_data <- state[at_current, , drop = FALSE]
  n_current <- nrow(current_data)
  required <- min(
    floor(n_current * control$suspension + 1),
    n_current
  )
  tox_confirm <- pmin(current_data$tox_end, current_data$tox_event)
  eff_confirm <- pmin(current_data$eff_end, current_data$eff_event)
  max(
    base_time,
    sort(tox_confirm)[required],
    sort(eff_confirm)[required]
  )
}

## 汇总中期分析时点的所有患者记录。
## 每名患者在决策时点的状态按参考实现（fun_TITE_STEIN_update.R 第 3.1 步）
## 分为三类：已观察到事件（delta = 1）、评估窗结束且无事件（delta = 0）
## 以及尚未判定（delta = -1）。
## 返回每个剂量组合的观测计数、各终点有效样本量、疗效结局已判定的例数
## 和待定结局计数。
tite.interim <- function(state, time, I, J, control) {
  dims <- c(I, J)
  n_mat <- x_mat <- y_mat <- ess_t <- ess_e <- asc_e <- matrix(0, I, J)

  tox_event_seen <- state$tox_event <= time
  eff_event_seen <- state$eff_event <= time
  tox_no_event_complete <- state$tox_end <= time & !tox_event_seen
  eff_no_event_complete <- state$eff_end <= time & !eff_event_seen
  tox_pending <- !tox_event_seen & !tox_no_event_complete
  eff_pending <- !eff_event_seen & !eff_no_event_complete
  follow_up <- pmax(0, time - state$enroll)
  tox_weight <- ifelse(
    tox_pending,
    pmin(follow_up / control$tox_window, 1),
    0
  )
  eff_weight <- ifelse(
    eff_pending,
    pmin(follow_up / control$eff_window, 1),
    0
  )

  cell <- cbind(state$dose_i, state$dose_j)
  for (k in seq_len(nrow(state))) {
    i <- cell[k, 1]
    j <- cell[k, 2]
    n_mat[i, j] <- n_mat[i, j] + 1
    x_mat[i, j] <- x_mat[i, j] + tox_event_seen[k]
    y_mat[i, j] <- y_mat[i, j] + eff_event_seen[k]
    ess_t[i, j] <- ess_t[i, j] + tox_event_seen[k] +
      tox_no_event_complete[k] + tox_weight[k]
    ess_e[i, j] <- ess_e[i, j] + eff_event_seen[k] +
      eff_no_event_complete[k] + eff_weight[k]
    ## ESS_e 中权重为 1 的部分：已观察到应答或评估窗结束且无应答；
    ## 尚未判定的患者不计入。
    asc_e[i, j] <- asc_e[i, j] + eff_event_seen[k] +
      eff_no_event_complete[k]
  }

  list(
    n = n_mat,
    x = x_mat,
    y = y_mat,
    ess_t = ess_t,
    ess_e = ess_e,
    asc_e = asc_e,
    pending_t = matrix(
      tabulate(
        (state$dose_j[tox_pending] - 1L) * I +
          state$dose_i[tox_pending],
        nbins = prod(dims)
      ),
      I,
      J
    ),
    pending_e = matrix(
      tabulate(
        (state$dose_j[eff_pending] - 1L) * I +
          state$dose_i[eff_pending],
        nbins = prod(dims)
      ),
      I,
      J
    )
  )
}

## 将患者级状态汇总为完整观测的最终计数矩阵。
## 返回列表还包含以天为单位的总随访时长。
tite.final.counts <- function(state, I, J) {
  n_mat <- x_mat <- y_mat <- matrix(0, I, J)
  for (k in seq_len(nrow(state))) {
    i <- state$dose_i[k]
    j <- state$dose_j[k]
    n_mat[i, j] <- n_mat[i, j] + 1
    x_mat[i, j] <- x_mat[i, j] + state$tox_outcome[k]
    y_mat[i, j] <- y_mat[i, j] + state$eff_outcome[k]
  }
  list(
    n = n_mat,
    x = x_mat,
    y = y_mat,
    duration = if (nrow(state) == 0L) {
      0
    } else {
      max(c(state$tox_end, state$eff_end))
    }
  )
}
