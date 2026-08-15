## Time-to-event data generation and effective-sample-size calculations.
##
## This follows the TITE-STEIN paper and its reference implementation:
##   * event times are uniform within their assessment windows;
##   * a pending no-event outcome contributes follow-up/window to the ESS;
##   * an observed event or a completed no-event assessment contributes 1;
##   * accrual can pause until the requested fraction of toxicity AND efficacy
##     outcomes at the current combination has been ascertained.

## Merge user delay settings with defaults and validate every supported field.
## Returns a complete, normalized delay-control list.
normalize.delay.control <- function(control = NULL) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("delay.control must be a list")

  unknown <- setdiff(names(control), names(DEFAULT_DELAY_CONTROL))
  if (length(unknown) > 0L) {
    stop("Unknown delay.control field(s): ", paste(unknown, collapse = ", "))
  }
  out <- utils::modifyList(DEFAULT_DELAY_CONTROL, control)

  positive <- c("accrual", "tox_window", "eff_window")
  for (field in positive) {
    value <- out[[field]]
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(field, " must be a positive finite number")
    }
  }
  if (
    length(out$suspension) != 1L ||
      !is.finite(out$suspension) ||
      out$suspension <= 0 ||
      out$suspension > 1
  ) {
    stop("suspension must lie in (0, 1]")
  }
  for (field in c("use_suspension", "accrual_random")) {
    if (length(out[[field]]) != 1L || is.na(out[[field]])) {
      stop(field, " must be TRUE or FALSE")
    }
    out[[field]] <- as.logical(out[[field]])
  }
  for (field in c("tox_distribution", "eff_distribution")) {
    out[[field]] <- tolower(as.character(out[[field]])[1])
    if (out[[field]] != "uniform") {
      stop(field, " currently supports only 'uniform'")
    }
  }
  out
}

## Create an empty patient-level state table for a time-to-event trial.
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

## Enroll one cohort and generate its toxicity/efficacy outcomes and event times.
## Returns `state` with one row appended for each newly enrolled patient.
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

## Determine the earliest permitted decision time after the latest cohort.
## When suspension is enabled, the required fraction of both endpoints at the
## current combination must be observed before accrual can continue.
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

## Summarize all patient records at an interim analysis time.
## Returns observed counts, endpoint-specific effective sample sizes, and
## pending-outcome counts for every dose combination.
tite.interim <- function(state, time, I, J, control) {
  dims <- c(I, J)
  n_mat <- x_mat <- y_mat <- ess_t <- ess_e <- matrix(0, I, J)

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
  }

  list(
    n = n_mat,
    x = x_mat,
    y = y_mat,
    ess_t = ess_t,
    ess_e = ess_e,
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

## Collapse patient-level state into fully observed final count matrices.
## The returned list also contains total follow-up duration in days.
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
