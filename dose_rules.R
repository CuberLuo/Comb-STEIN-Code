## Posterior monitoring, admissibility updates and dose-selection rules.

## All three assume the conjugate Beta(1, 1) prior; arguments may be scalars,
## vectors or conformable matrices.
## Pr(toxicity probability > cut)
pr.overtox <- function(x, n, cut) {
  1 - pbeta(cut, x + 1, n - x + 1)
}

## Pr(efficacy probability < cut)
pr.ineff <- function(y, n, cut) {
  pbeta(cut, y + 1, n - y + 1)
}

## Pr(efficacy probability > cut)
pr.eff.gt <- function(y, n, cut) {
  1 - pbeta(cut, y + 1, n - y + 1)
}

## Unobserved cells (n == 0) are never flagged: monitoring relies on data, not
## on the prior, so escalation is not blocked before a combination is tried.
monitor.flags <- function(x_mat, y_mat, n_tox, n_eff = n_tox) {
  observed_tox <- n_tox > 0
  observed_eff <- n_eff > 0
  overtox <- observed_tox &
    (pr.overtox(x_mat, n_tox, phi) > Tthresh)
  ineff <- observed_eff &
    (pr.ineff(y_mat, n_eff, psi1) > Ethresh)
  list(overtox = overtox, ineff = ineff)
}

get.true.utility <- function(pT, pE) {
  utility.fn(pT, pE, w1, w2, phi)
}

## True (oracle) OBDCs: every combination whose utility is within
## OBDC_UTILITY_TOLERANCE of the maximum, so ties are all counted as correct.
get.true.obdc <- function(pT, pE) {
  utility <- get.true.utility(pT, pE)
  which(
    max(utility) - utility <= OBDC_UTILITY_TOLERANCE,
    arr.ind = TRUE
  )
}

## An over-toxic cell also eliminates its north-east quadrant, by the monotonic
## toxicity assumption; a futile cell only eliminates itself.
update.admissible <- function(
  x_mat,
  y_mat,
  n_tox,
  admissible,
  n_eff = n_tox
) {
  I <- nrow(n_tox)
  J <- ncol(n_tox)
  flags <- monitor.flags(x_mat, y_mat, n_tox, n_eff)

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      if (flags$overtox[i, j]) admissible[i:I, j:J] <- FALSE
      if (flags$ineff[i, j]) admissible[i, j] <- FALSE
    }
  }
  admissible
}

## Local candidate set; no escalation is allowed once the observed toxicity
## rate at (i, j) exceeds phi, which leaves only the two de-escalations.
alt.set.C <- function(i, j, I, J, p_ij_hat) {
  candidates <- list()
  if (p_ij_hat > phi) {
    if (i > 1) candidates <- c(candidates, list(c(i - 1, j)))
    if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
  } else {
    if (i > 1) candidates <- c(candidates, list(c(i - 1, j)))
    if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
    candidates <- c(candidates, list(c(i, j)))
    if (j < J) candidates <- c(candidates, list(c(i, j + 1)))
    if (i < I) candidates <- c(candidates, list(c(i + 1, j)))
  }
  candidates
}

## The 1e-6 * (i + j) term breaks ties deterministically by dose index, so
## equal posteriors are resolved reproducibly rather than by Vapply order.
pick.best.by.eff <- function(cand_list, y_mat, n_mat) {
  posterior <- sapply(cand_list, function(ij) {
    score <- pr.eff.gt(
      y_mat[ij[1], ij[2]],
      n_mat[ij[1], ij[2]],
      psi
    )
    score + 1e-6 * (ij[1] + ij[2])
  })
  # debug.cat("y_mat:\n")
  # debug.print(y_mat)
  # debug.cat("n_mat:\n")
  # debug.print(n_mat)
  # debug.cat("posterior:\n")
  # debug.print(posterior)
  cand_list[[which.max(posterior)]]
}

## Final pick among combinations that were both treated and admissible.
## Pseudo-counts (x + 0.05) / (n + 0.1) are also passed to Iso::biviso as weights,
## which enforces the partial ordering of both agents before scoring utility.
## Returns c(row, column), or two NA integers when nothing is eligible.
final.selection <- function(
  x_mat,
  y_mat,
  n_mat,
  admissible
) {
  eligible <- admissible & n_mat > 0
  if (!any(eligible)) return(c(NA_integer_, NA_integer_))

  phat <- (x_mat + 0.05) / (n_mat + 0.1)
  p_star <- Iso::biviso(phat, w = n_mat + 0.1)
  q_star <- eff.estimate(y_mat, n_mat)
  utility <- utility.fn(p_star, q_star, w1, w2, phi)
  utility[!eligible] <- -Inf

  if (all(!is.finite(utility))) {
    return(c(NA_integer_, NA_integer_))
  }
  which(utility == max(utility), arr.ind = TRUE)[1, ]
}

## Run-in path: row-first, alternating between the two agents so both axes are
## exercised early. Hitting a grid edge stops that direction, so the returned
## list can be shorter than ncohort.
make.runin.path <- function(I, J, ncohort = RUNIN_COHORTS) {
  path <- list(c(1L, 1L))
  current <- c(1L, 1L)
  prefer_row <- TRUE

  while (length(path) < ncohort && (current[1] < I || current[2] < J)) {
    if (prefer_row && current[1] < I) {
      current[1] <- current[1] + 1L
    } else if (!prefer_row && current[2] < J) {
      current[2] <- current[2] + 1L
    } else if (current[2] < J) {
      current[2] <- current[2] + 1L
    } else {
      current[1] <- current[1] + 1L
    }
    path[[length(path) + 1L]] <- current
    prefer_row <- !prefer_row
  }
  path
}
