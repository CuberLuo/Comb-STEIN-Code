## Time-to-event outcome generation and effective sample size (ESS) computation.
## Follows TITE-STEIN:
##   * event times are uniform within the assessment window;
##   * pending outcomes contribute follow_up / window to ESS;
##   * an observed event or a completed event-free assessment counts as 1;
##   * with suspension, accrual waits for the required confirmation fraction.

# merges user delay settings with the defaults
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

## Advances cohort timing for the complete-outcome (non-TITE) design.
## suspension = 1: the next cohort is assigned only once all toxicity and
## efficacy outcomes of the current cohort are confirmed.
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

  # Complete mode uses confirmed binary outcomes only, no event times are drawn,
  # so the RNG stream of the original complete-outcome simulation is unchanged;
  # duration is defined by the end of the assessment window.
  tox_end <- enroll + control$tox_window
  eff_end <- enroll + control$eff_window

  next_accrual <- if (control$accrual_random) {
    max(enroll) + stats::runif(1L, 0, 2 * control$accrual) + 1
  } else {
    time_current + cohort_size * control$accrual + 1
  }

  list(
    enroll = enroll,
    next_time = max(next_accrual, tox_end, eff_end),
    duration = max(tox_end, eff_end)
  )
}

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

## Earliest time at which the next decision may be made after the latest cohort.
## With suspension, accrual is held until the required confirmation fraction of the
## current combination is reached for both endpoints.
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

## Summarizes all patient records at an interim decision time.
## Each patient is classified as event observed (delta = 1), window closed with
## no event (delta = 0), or not yet adjudicated (delta = -1); pending patients
## contribute the follow-up fraction follow_up / window to ESS.
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
    # asc_e: weight-1 part of ESS_e (response observed or window closed without
    # response); pending patients are excluded.
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

## Collapses patient-level state into fully observed final count matrices.
## `duration` is total follow-up time in days.
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

## Per-patient trajectory table for single-trial trajectory plots: one row per
## patient with the adjudication status of both endpoints at the decision time
## of that patient's cohort.
##   tox_state: dlt (DLT observed) | pending (window still open) | none (window closed, no DLT)
##   eff_state: benefit | pending | none
tite.cohort.trace <- function(state, decision_times) {
  if (nrow(state) == 0L) {
    return(NULL)
  }
  t_dec <- decision_times[state$cohort]
  tox_seen <- state$tox_event <= t_dec
  eff_seen <- state$eff_event <= t_dec
  tox_done <- state$tox_end <= t_dec
  eff_done <- state$eff_end <= t_dec
  data.frame(
    cohort = state$cohort,
    dose_i = state$dose_i,
    dose_j = state$dose_j,
    enroll = state$enroll,
    decision_time = t_dec,
    tox_state = ifelse(tox_seen, "dlt", ifelse(tox_done, "none", "pending")),
    eff_state = ifelse(eff_seen, "benefit", ifelse(eff_done, "none", "pending")),
    stringsAsFactors = FALSE
  )
}
