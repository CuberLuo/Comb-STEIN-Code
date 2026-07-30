## Simplified row-by-row waterfall comparator.

simulate.trial.waterfall <- function(
  pT,
  pE,
  method = c("Waterfall-L", "Waterfall-H")
) {
  method <- match.arg(method)
  I <- nrow(pT)
  J <- ncol(pT)
  ncohort <- I * J * COHORTS_PER_COMBINATION
  x_mat <- y_mat <- n_mat <- matrix(0, I, J)
  admissible <- matrix(TRUE, I, J)
  earlystop_tox <- FALSE
  earlystop_eff <- FALSE
  tox_monitor_triggered <- FALSE

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
    tox_monitor_triggered = tox_monitor_triggered
  )
}
