## Bayesian quadratic logistic efficacy estimator of Cai, Yuan and Ji (2014),
## as used by Comb-BOIN12:
##
##   logit(q[i,j]) = gamma[1] + gamma[2]*a[i] + gamma[3]*b[j]
##                 + gamma[4]*a[i]^2 + gamma[5]*b[j]^2.
##
## Dose scores are standardised to mean 0, SD 0.5. Only main effects and squared
## terms are kept: no agent-by-agent interaction and no per-cell random effect,
## so the fit borrows strength across the whole grid.

.efficacy_fit_cache <- new.env(parent = emptyenv())

## A single-level agent gets the neutral score 0, because its SD is undefined.
.standardize.dose.scores <- function(n.levels) {
  scores <- seq_len(n.levels)
  if (n.levels == 1L) return(0)
  (scores - mean(scores)) / (2 * stats::sd(scores))
}

## Cache key over counts + MCMC settings: identical data refits only once.
.efficacy.cache.key <- function(
  y_mat,
  n_mat,
  n_burnin,
  n_samples,
  n_chains,
  thin
) {
  paste(
    nrow(y_mat),
    ncol(y_mat),
    paste(as.integer(y_mat), collapse = ","),
    paste(as.integer(n_mat), collapse = ","),
    n_burnin,
    n_samples,
    n_chains,
    thin,
    sep = "|"
  )
}

## Seed derived from the data (not the clock), so fits are reproducible.
.efficacy.jags.seed <- function(y_mat, n_mat) {
  values <- c(
    nrow(y_mat),
    ncol(y_mat),
    as.integer(y_mat),
    as.integer(n_mat)
  )
  positions <- seq_along(values)
  weighted <- (values + 1) * ((positions * 104729) %% 2147483629)
  as.integer((sum(weighted) %% 2147483646) + 1)
}

## Posterior mean efficacy per cell, from the hierarchical quadratic logistic
## model above (weakly informative t priors on gamma).
eff.estimate <- function(
  y_mat,
  n_mat,
  n_burnin = EFF_MCMC_BURNIN,
  n_samples = EFF_MCMC_SAMPLES,
  n_chains = EFF_MCMC_CHAINS,
  thin = EFF_MCMC_THIN,
  use_cache = TRUE
) {
  I <- nrow(y_mat)
  J <- ncol(y_mat)
  cache.key <- .efficacy.cache.key(
    y_mat,
    n_mat,
    n_burnin,
    n_samples,
    n_chains,
    thin
  )
  if (use_cache && exists(
    cache.key,
    envir = .efficacy_fit_cache,
    inherits = FALSE
  )) {
    return(get(cache.key, envir = .efficacy_fit_cache, inherits = FALSE))
  }

  a <- .standardize.dose.scores(I)
  b <- .standardize.dose.scores(J)
  design <- cbind(
    intercept = 1,
    agent_a = rep(a, times = J),
    agent_b = rep(b, each = I),
    agent_a_sq = rep(a^2, times = J),
    agent_b_sq = rep(b^2, each = I)
  )

  model.text <- "
    model {
      for (d in 1:N) {
        y[d] ~ dbin(q[d], n[d])
        logit(q[d]) <- inprod(gamma[1:5], X[d, 1:5])
      }
      gamma[1] ~ dt(0, 0.01, 1)
      for (k in 2:5) {
        gamma[k] ~ dt(0, 0.16, 1)
      }
    }
  "

  base.seed <- .efficacy.jags.seed(y_mat, n_mat)
  initial.values <- lapply(seq_len(n_chains), function(chain) {
    chain.seed <- as.integer(
      ((base.seed - 1 + (chain - 1) * 104729) %% 2147483646) + 1
    )
    list(
      gamma = rep(0, 5),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = chain.seed
    )
  })

  model.connection <- textConnection(model.text)
  model <- tryCatch(
    rjags::jags.model(
      model.connection,
      data = list(
        N = length(y_mat),
        y = as.integer(y_mat),
        n = as.integer(n_mat),
        X = unname(design)
      ),
      inits = initial.values,
      n.chains = n_chains,
      quiet = TRUE
    ),
    finally = close(model.connection)
  )
  stats::update(
    model,
    n.iter = n_burnin,
    progress.bar = "none"
  )

  draws.per.chain <- ceiling(n_samples / n_chains)
  ## Read the mcarray from rjags directly: for large fits coda::as.matrix()
  ## can retain expression objects as rownames, which then breaks
  ## base::as.data.frame(). The mcarray is already cell x iteration x chain.
  samples <- rjags::jags.samples(
    model,
    variable.names = "q",
    n.iter = draws.per.chain * thin,
    thin = thin,
    type = "trace",
    progress.bar = "none"
  )
  q.samples <- samples[["q"]]
  expected.cells <- length(y_mat)
  if (
    is.null(q.samples) ||
      length(dim(q.samples)) != 3L ||
      dim(q.samples)[1] != expected.cells
  ) {
    stop("JAGS output has unexpected dimensions for efficacy probabilities.")
  }
  q.draws <- matrix(
    as.numeric(q.samples),
    nrow = expected.cells,
    dimnames = NULL
  )
  if (ncol(q.draws) > n_samples) {
    q.draws <- q.draws[, seq_len(n_samples), drop = FALSE]
  }
  if (ncol(q.draws) < n_samples) {
    stop("JAGS returned fewer efficacy draws than requested.")
  }

  q.final <- matrix(
    rowMeans(q.draws),
    nrow = I,
    ncol = J
  )

  if (use_cache) {
    assign(cache.key, q.final, envir = .efficacy_fit_cache)
  }
  q.final
}

## Design utility. The extra w2-weighted penalty applies only above phi, making
## phi both the toxicity decision boundary and the kink of the utility.
utility.fn <- function(p, q, w1, w2, phi) {
  q - w1 * p - w2 * p * (p > phi)
}
