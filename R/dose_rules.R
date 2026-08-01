## Posterior monitoring, admissibility updates, and dose-selection rules.

pr.overtox <- function(x, n, cut) {
  1 - pbeta(cut, x + 1, n - x + 1)
}

pr.ineff <- function(y, n, cut) {
  pbeta(cut, y + 1, n - y + 1)
}

pr.eff.gt <- function(y, n, cut) {
  1 - pbeta(cut, y + 1, n - y + 1)
}

monitor.flags <- function(x_mat, y_mat, n_mat) {
  observed <- n_mat > 0
  overtox <- observed & (pr.overtox(x_mat, n_mat, phi) > Tthresh)
  ineff <- observed & (pr.ineff(y_mat, n_mat, psi1) > Ethresh)
  list(overtox = overtox, ineff = ineff)
}

get.true.obdc <- function(pT, pE) {
  utility <- utility.fn(pT, pE, w1, w2, phi)
  which(
    abs(utility - max(utility)) <= sqrt(.Machine$double.eps),
    arr.ind = TRUE
  )
}

update.admissible <- function(x_mat, y_mat, n_mat, admissible) {
  I <- nrow(n_mat)
  J <- ncol(n_mat)
  flags <- monitor.flags(x_mat, y_mat, n_mat)

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      if (n_mat[i, j] > 0) {
        if (flags$overtox[i, j]) admissible[i:I, j:J] <- FALSE
        if (flags$ineff[i, j]) admissible[i, j] <- FALSE
      }
    }
  }
  admissible
}

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
