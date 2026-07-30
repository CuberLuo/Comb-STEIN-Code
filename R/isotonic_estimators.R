## Unimodal efficacy estimators (Section 2.2.4).

unimodal.1d <- function(y, w, peak) {
  if (peak < 1L || peak > length(y)) stop("peak is outside y")
  as.numeric(Iso::ufit(
    y,
    imode = peak,
    x = seq_along(y),
    w = w
  )[[2]])
}

unimodal.2d <- function(
  qhat,
  n_mat,
  peak_row,
  peak_col,
  max.iter = 100,
  tol = 1e-8
) {
  I <- nrow(qhat)
  J <- ncol(qhat)
  q <- qhat
  w <- n_mat + 0.5

  for (iter in seq_len(max.iter)) {
    q.old <- q
    for (i in seq_len(I)) {
      q[i, ] <- unimodal.1d(q[i, ], w[i, ], peak_col)
    }
    for (j in seq_len(J)) {
      q[, j] <- unimodal.1d(q[, j], w[, j], peak_row)
    }
    if (max(abs(q - q.old)) < tol) break
  }
  q
}

eff.estimate <- function(y_mat, n_mat) {
  I <- nrow(y_mat)
  J <- ncol(y_mat)
  qhat <- (y_mat + 0.05) / (n_mat + 0.1)
  loglik <- matrix(-Inf, I, J)
  fits <- vector("list", I * J)
  idx <- 1

  for (m in seq_len(I)) {
    for (nn in seq_len(J)) {
      fit <- unimodal.2d(qhat, n_mat, m, nn)
      fit <- pmin(pmax(fit, 1e-4), 1 - 1e-4)
      loglik[m, nn] <- sum(dbinom(
        as.vector(y_mat),
        as.vector(n_mat),
        as.vector(fit),
        log = TRUE
      ))
      fits[[idx]] <- fit
      idx <- idx + 1
    }
  }

  loglik <- loglik - max(loglik)
  model.weights <- exp(loglik)
  model.weights <- model.weights / sum(model.weights)
  q.final <- matrix(0, I, J)
  idx <- 1

  for (m in seq_len(I)) {
    for (nn in seq_len(J)) {
      q.final <- q.final + model.weights[m, nn] * fits[[idx]]
      idx <- idx + 1
    }
  }
  q.final
}

utility.fn <- function(p, q, w1, w2, phi) {
  q - w1 * p - w2 * p * (p > phi)
}
