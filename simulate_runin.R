## RunIn-STEIN simulation: complete-outcome and time-to-event (TITE) variants.
## Both follow the run-in path for the first RUNIN_COHORTS cohorts, then switch
## to the local candidate set, stopping early if nothing stays admissible.

## Complete-outcome version; `delayed = TRUE` delegates to the TITE simulator.
simulate.trial.runin.stein <- function(
  pT,
  pE,
  delayed = FALSE,
  delay.control = NULL
) {
  if (isTRUE(delayed)) {
    return(simulate.trial.runin.stein.tite(
      pT, pE, normalize.delay.control(delay.control)
    ))
  }
  delay.control <- normalize.delay.control(delay.control)
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)
  current <- c(1, 1)
  runin_path <- make.runin.path(I, J)
  debug.print(runin_path)
  path_pos <- 1
  earlystop_tox <- earlystop_no_candidate <- FALSE
  time_current <- 0
  duration <- 0

  for (cohort in seq_len(ncohort)) {
    i <- current[1]; j <- current[2]
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
    admissible <- update.admissible(x_mat, y_mat, n_mat, admissible)
    if (!any(admissible)) { earlystop_tox <- TRUE; break }
    if (n_mat[i, j] >= EFF_STOP_N &&
        pr.eff.gt(y_mat[i, j], n_mat[i, j], EFF_STOP_CUT) > EFF_STOP_PROB) {
      break
    }
    if (cohort == ncohort) break

    p_ij_hat <- x_mat[i, j] / n_mat[i, j]
    # first RUNIN_COHORTS cohorts follow the run-in path
    if (cohort < RUNIN_COHORTS) {
      if (path_pos < length(runin_path)) path_pos <- path_pos + 1L
      next_dose <- runin_path[[path_pos]]
      if (p_ij_hat > phi || !admissible[next_dose[1], next_dose[2]]) {
        candidates <- Filter(
          function(ij) admissible[ij[1], ij[2]],
          alt.set.C(i, j, I, J, p_ij_hat)
        )
        if (length(candidates) == 0 && admissible[i, j]) {
          candidates <- list(c(i, j))
        }
        if (length(candidates) == 0) {
          earlystop_no_candidate <- TRUE
          break
        }
        next_dose <- pick.best.by.eff(candidates, y_mat, n_mat)
      }
      current <- next_dose
    } else {
      candidates <- Filter(
        function(ij) admissible[ij[1], ij[2]],
        alt.set.C(i, j, I, J, p_ij_hat)
      )
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
  }

  selected <- if (earlystop_tox || earlystop_no_candidate) {
    c(NA_integer_, NA_integer_)
  } else {
    final.selection(x_mat, y_mat, n_mat, admissible)
  }
  list(
    sel = selected, x = x_mat, y = y_mat, n = n_mat,
    earlystop_no_candidate = earlystop_no_candidate,
    delayed = FALSE, duration = duration
  )
}

## TITE version. Monitoring uses effective sample sizes rather than enrolled
## counts, because pending outcomes contribute partial weight.
simulate.trial.runin.stein.tite <- function(pT, pE, delay.control) {
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  state <- new.tite.state()
  admissible <- matrix(TRUE, I, J)
  current <- c(1, 1)
  runin_path <- make.runin.path(I, J)
  path_pos <- 1L
  time_current <- 0
  earlystop_tox <- FALSE
  earlystop_no_candidate <- FALSE

  for (cohort in seq_len(ncohort)) {
    i <- current[1]
    j <- current[2]
    state <- tite.enroll.cohort(
      state, cohort, current, pT[i, j], pE[i, j],
      time_current, cohortsize, delay.control
    )
    decision_time <- tite.next.decision.time(
      state, current, time_current, cohortsize, delay.control
    )
    interim <- tite.interim(state, decision_time, I, J, delay.control)
    time_current <- decision_time

    admissible <- update.admissible(
      interim$x, interim$y, interim$ess_t,
      admissible, n_eff = interim$ess_e
    )
    if (!any(admissible)) {
      earlystop_tox <- TRUE
      break
    }
    ## As in Comb-STEIN: the efficacy stop uses only resolved efficacy outcomes.
    if (
      interim$asc_e[i, j] >= EFF_STOP_N &&
        pr.eff.gt(
          interim$y[i, j], interim$asc_e[i, j], EFF_STOP_CUT
        ) > EFF_STOP_PROB
    ) {
      break
    }
    if (cohort == ncohort) break

    p_ij_hat <- interim$x[i, j] / interim$ess_t[i, j]
    if (cohort < RUNIN_COHORTS) {
      if (path_pos < length(runin_path)) path_pos <- path_pos + 1L
      next_dose <- runin_path[[path_pos]]
      if (p_ij_hat > phi || !admissible[next_dose[1], next_dose[2]]) {
        candidates <- Filter(
          function(ij) admissible[ij[1], ij[2]],
          alt.set.C(i, j, I, J, p_ij_hat)
        )
        if (length(candidates) == 0L && admissible[i, j]) {
          candidates <- list(c(i, j))
        }
        if (length(candidates) == 0L) {
          earlystop_no_candidate <- TRUE
          break
        }
        next_dose <- pick.best.by.eff(
          candidates, interim$y, interim$ess_e
        )
      }
      current <- next_dose
    } else {
      candidates <- Filter(
        function(ij) admissible[ij[1], ij[2]],
        alt.set.C(i, j, I, J, p_ij_hat)
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
    duration = final$duration
  )
}
