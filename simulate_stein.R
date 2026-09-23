## Simulates one Comb-STEIN trial under complete outcomes, dispatching to the
## TITE implementation when delayed outcomes are requested.
## Returns allocation counts, outcomes, selected combination, stopping flags and duration.
simulate.trial.comb.stein <- function(
  pT,
  pE,
  method = c("CombSTEIN-L", "CombSTEIN-H"),
  delayed = FALSE,
  delay.control = NULL,
  trace = FALSE
) {
  method <- match.arg(method)
  if (isTRUE(delayed)) {
    return(simulate.trial.comb.stein.tite(
      pT,
      pE,
      method,
      normalize.delay.control(delay.control),
      trace = trace
    ))
  }
  delay.control <- normalize.delay.control(delay.control)
  # COHORTS_PER_COMBINATION cohorts are planned per dose combination.
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  # x_mat / y_mat / n_mat: toxicity events, efficacy events, enrolled patients.
  # admissible: combinations still allowed for use.
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)

  # CombSTEIN-H starts at the highest dose, CombSTEIN-L at the lowest.
  current <- if (method == "CombSTEIN-H") {
    c(I, J)
  } else {
    c(1, 1)
  }
  earlystop_tox <- FALSE
  earlystop_no_candidate <- FALSE
  time_current <- 0
  duration <- 0
  # trace: per-cohort records for single-trial trajectory plots; NULL otherwise.
  cohort.trace <- if (trace) vector("list", ncohort) else NULL

  for (cohort in seq_len(ncohort)) {
    debug.cat("current:", current, "\n")
    i <- current[1]
    j <- current[2]

    tox_outcome <- as.integer(runif(cohortsize) < pT[i, j])
    eff_outcome <- as.integer(runif(cohortsize) < pE[i, j])
    wT <- sum(tox_outcome)
    wE <- sum(eff_outcome)
    debug.cat("wT:", wT, "\n")
    debug.cat("wE:", wE, "\n")
    x_mat[i, j] <- x_mat[i, j] + wT
    y_mat[i, j] <- y_mat[i, j] + wE
    n_mat[i, j] <- n_mat[i, j] + cohortsize
    timing <- complete.cohort.timing(
      time_current, tox_outcome, eff_outcome, delay.control
    )
    # next_time: earliest time the next cohort may start enrolling (days).
    time_current <- timing$next_time
    # duration: time until all enrolled patients complete assessment (days).
    duration <- timing$duration

    if (trace) {
      cohort.trace[[cohort]] <- data.frame(
        cohort = cohort,
        dose_i = i,
        dose_j = j,
        enroll_start = min(timing$enroll),
        enroll_end = max(timing$enroll),
        decision_time = time_current,
        tox = wT,
        eff = wE,
        pend_tox = 0L,
        pend_eff = 0L
      )
    }

    # overly toxic combinations (and those up/right of them) and inefficacious ones are dropped
    admissible <- update.admissible(x_mat, y_mat, n_mat, admissible)
    if (!any(admissible)) {
      # no admissible combination left => stop for toxicity
      earlystop_tox <- TRUE
      break
    }

    # early stop for sufficient efficacy: needs EFF_STOP_N patients and
    # P(efficacy rate > EFF_STOP_CUT) > EFF_STOP_PROB
    if (
      n_mat[i, j] >= EFF_STOP_N &&
        pr.eff.gt(y_mat[i, j], n_mat[i, j], EFF_STOP_CUT) > EFF_STOP_PROB
    ) {
      break
    }

    if (cohort == ncohort) break

    # local candidate set from the current toxicity estimate, then intersect with admissible;
    # empty => stay at current dose, or stop if it is itself inadmissible
    p_ij_hat <- x_mat[i, j] / n_mat[i, j]
    candidates <- alt.set.C(i, j, I, J, p_ij_hat)
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

  selected <- if (earlystop_tox || earlystop_no_candidate) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(x_mat, y_mat, n_mat, admissible)
  }

  list(
    sel = selected,
    x = x_mat,
    y = y_mat,
    n = n_mat,
    earlystop_no_candidate = earlystop_no_candidate,
    delayed = FALSE,
    duration = duration,
    trace = if (trace) do.call(rbind, Filter(Negate(is.null), cohort.trace))
  )
}

## Time-to-event (TITE) version of the STEIN trial.
## `delay.control` must already be normalized by `normalize.delay.control()`.
simulate.trial.comb.stein.tite <- function(
  pT,
  pE,
  method,
  delay.control,
  trace = FALSE
) {
  # same dose-moving rules as the complete-outcome version, but each decision is
  # based on effective sample sizes computed from current follow-up
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  state <- new.tite.state()
  admissible <- matrix(TRUE, I, J)
  current <- if (method == "CombSTEIN-H") {
    c(I, J)
  } else {
    c(1, 1)
  }
  time_current <- 0
  earlystop_tox <- FALSE
  earlystop_no_candidate <- FALSE
  # decision_times[k]: time at which cohort k's dose decision is made
  decision_times <- numeric(0)

  for (cohort in seq_len(ncohort)) {
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
    decision_time <- tite.next.decision.time(
      state,
      current,
      time_current,
      cohortsize,
      delay.control
    )
    interim <- tite.interim(state, decision_time, I, J, delay.control)
    time_current <- decision_time
    decision_times[cohort] <- decision_time

    admissible <- update.admissible(
      interim$x,
      interim$y,
      interim$ess_t,
      admissible,
      n_eff = interim$ess_e
    )
    if (!any(admissible)) {
      # no admissible combination left => stop for toxicity
      earlystop_tox <- TRUE
      break
    }

    # Efficacy stopping uses only adjudicated outcomes (delta_e of ref impl step 3.1):
    # pending outcomes enter neither the numerator nor the ESS-weighted denominator,
    # otherwise efficacy evidence is deflated and early stopping almost never fires.
    if (
      interim$asc_e[i, j] >= EFF_STOP_N &&
        pr.eff.gt(
          interim$y[i, j],
          interim$asc_e[i, j],
          EFF_STOP_CUT
        ) > EFF_STOP_PROB
    ) {
      break
    }
    if (cohort == ncohort) break

    # toxicity rate from ESS_t decides which directions may be explored next
    p_ij_hat <- interim$x[i, j] / interim$ess_t[i, j]
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

  final <- tite.final.counts(state, I, J)
  selected <- if (earlystop_tox || earlystop_no_candidate) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(final$x, final$y, final$n, admissible)
  }

  list(
    sel = selected,
    x = final$x,
    y = final$y,
    n = final$n,
    earlystop_no_candidate = earlystop_no_candidate,
    delayed = TRUE,
    duration = final$duration,
    trace = if (trace) tite.cohort.trace(state, decision_times)
  )
}

