## Simulate one trial under Comb-STEIN (Rule A) or RunIn-STEIN (Rule B).

simulate.trial.stein <- function(
  pT,
  pE,
  method = c("CombSTEIN-L", "CombSTEIN-H", "RunIn-STEIN")
) {
  method <- match.arg(method)
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)

  current <- if (method == "CombSTEIN-H") c(I, J) else c(1, 1)
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  earlystop_no_candidate <- FALSE
  tox_monitor_triggered <- FALSE

  if (method == "RunIn-STEIN") {
    runin_path <- make.runin.path(I, J)
    path_pos <- 1
  }

  for (cohort in seq_len(ncohort)) {
    i <- current[1]
    j <- current[2]

    wT <- sum(runif(cohortsize) < pT[i, j])
    wE <- sum(runif(cohortsize) < pE[i, j])
    x_mat[i, j] <- x_mat[i, j] + wT
    y_mat[i, j] <- y_mat[i, j] + wE
    n_mat[i, j] <- n_mat[i, j] + cohortsize

    flags <- monitor.flags(x_mat, y_mat, n_mat)
    tox_monitor_triggered <- tox_monitor_triggered || any(flags$overtox)
    admissible <- update.admissible(x_mat, y_mat, n_mat, admissible)
    if (!any(admissible)) {
      earlystop_tox <- TRUE
      break
    }

    if (
      n_mat[i, j] >= EFF_STOP_N &&
        pr.eff.gt(y_mat[i, j], n_mat[i, j], EFF_STOP_CUT) > EFF_STOP_PROB
    ) {
      earlystop_eff <- TRUE
      break
    }

    if (cohort == ncohort) break

    if (method == "RunIn-STEIN" && cohort < RUNIN_COHORTS) {
      p_ij_hat <- x_mat[i, j] / n_mat[i, j]
      if (path_pos < length(runin_path)) path_pos <- path_pos + 1L
      next_dose <- runin_path[[path_pos]]

      if (p_ij_hat > phi || !admissible[next_dose[1], next_dose[2]]) {
        candidates <- alt.set.C(i, j, I, J, p_ij_hat)
        candidates <- Filter(
          function(ij) admissible[ij[1], ij[2]],
          candidates
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
      p_ij_hat <- x_mat[i, j] / n_mat[i, j]
      candidates <- alt.set.C(i, j, I, J, p_ij_hat)
      candidates <- Filter(
        function(ij) admissible[ij[1], ij[2]],
        candidates
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
    sel = selected,
    x = x_mat,
    y = y_mat,
    n = n_mat,
    earlystop_tox = earlystop_tox,
    earlystop_eff = earlystop_eff,
    earlystop_no_candidate = earlystop_no_candidate,
    tox_monitor_triggered = tox_monitor_triggered
  )
}
