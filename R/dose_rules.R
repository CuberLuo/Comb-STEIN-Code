## Posterior monitoring, admissibility updates, and dose-selection rules.

## Compute Pr(toxicity probability > `cut`) under a Beta(1, 1) prior.
## `x` and `n` may be scalars, vectors, or conformable matrices.
pr.overtox <- function(x, n, cut) {
  1 - pbeta(cut, x + 1, n - x + 1)
}

## Compute Pr(efficacy probability < `cut`) under a Beta(1, 1) prior.
## `y` and `n` may be scalars, vectors, or conformable matrices.
pr.ineff <- function(y, n, cut) {
  pbeta(cut, y + 1, n - y + 1)
}

## Compute Pr(efficacy probability > `cut`) under a Beta(1, 1) prior.
## `y` and `n` may be scalars, vectors, or conformable matrices.
pr.eff.gt <- function(y, n, cut) {
  1 - pbeta(cut, y + 1, n - y + 1)
}

## Flag observed dose combinations that cross toxicity or futility thresholds.
## Returns logical matrices named `overtox` and `ineff`.
monitor.flags <- function(x_mat, y_mat, n_tox, n_eff = n_tox) {
  observed_tox <- n_tox > 0
  observed_eff <- n_eff > 0
  overtox <- observed_tox &
    (pr.overtox(x_mat, n_tox, phi) > Tthresh)
  ineff <- observed_eff &
    (pr.ineff(y_mat, n_eff, psi1) > Ethresh)
  list(overtox = overtox, ineff = ineff)
}

## Identify all true optimal biological dose combinations by maximum utility.
## Returns a two-column matrix of row and column indices, including exact ties.
get.true.obdc <- function(pT, pE) {
  utility <- utility.fn(pT, pE, w1, w2, phi)
  which(
    abs(utility - max(utility)) <= sqrt(.Machine$double.eps),
    arr.ind = TRUE
  )
}

## Apply posterior safety and futility rules to the admissibility matrix.
## An overtoxic cell also eliminates combinations northeast of that cell;
## an inefficacious cell eliminates only itself.
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

## Build the local candidate set for the next dose-combination assignment.
## Escalation is excluded when the current observed toxicity rate exceeds phi.
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

## Select the candidate with the largest posterior probability of efficacy.
## A deterministic dose-index term breaks otherwise equal posterior scores.
pick.best.by.eff <- function(cand_list, y_mat, n_mat) {
  posterior <- sapply(cand_list, function(ij) {
    score <- pr.eff.gt(
      y_mat[ij[1], ij[2]],
      n_mat[ij[1], ij[2]],
      psi
    )
    score + 1e-6 * (ij[1] + ij[2])
  })
  cand_list[[which.max(posterior)]]
}

## Select the final treated, admissible combination with maximum fitted utility.
## Returns c(row, column), or two NA integers when no combination is eligible.
final.selection <- function(x_mat, y_mat, n_mat, admissible) {
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
  selected <- which(utility == max(utility), arr.ind = TRUE)
  selected[1, ]
}

## Construct the alternating row/column escalation path used during run-in.
## The returned list contains at most `ncohort` integer index pairs.
make.runin.path <- function(I, J, ncohort = RUNIN_COHORTS) {
  path <- list(c(1L, 1L))
  current <- c(1L, 1L)
  prefer_col <- TRUE

  while (length(path) < ncohort && (current[1] < I || current[2] < J)) {
    if (prefer_col && current[2] < J) {
      current[2] <- current[2] + 1L
    } else if (!prefer_col && current[1] < I) {
      current[1] <- current[1] + 1L
    } else if (current[2] < J) {
      current[2] <- current[2] + 1L
    } else {
      current[1] <- current[1] + 1L
    }
    path[[length(path) + 1L]] <- current
    prefer_col <- !prefer_col
  }
  path
}
