## Cai、Yuan 和 Ji（2014）及 Comb-BOIN12 使用的贝叶斯二次逻辑回归疗效估计器：
##
##   logit(q[i,j]) = gamma[1] + gamma[2] * a[i] + gamma[3] * b[j]
##                   + gamma[4] * a[i]^2 + gamma[5] * b[j]^2
##                   + gamma[6] * a[i] * b[j].
##
## 剂量得分 a 和 b 标准化为均值 0、标准差 0.5。
## 交互项补全二阶响应面，使疗效峰值可以随另一个药物的剂量变化而移动。
## 层级化的剂量单元偏差可以避免全局二次曲面对局部组合效应过度平滑。
## 数据稀疏时这些偏差会收缩到 0；数据充足时则可以恢复局部峰值。

.efficacy_fit_cache <- new.env(parent = emptyenv())

## 将有序剂量水平标准化为均值 0、标准差 0.5。
## 只有一个剂量水平的药物使用中性得分 0。
.standardize.dose.scores <- function(n.levels) {
  scores <- seq_len(n.levels)
  if (n.levels == 1L) return(0)
  (scores - mean(scores)) / (2 * stats::sd(scores))
}

## 将疗效计数和 MCMC 设置序列化为确定性的缓存键。
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

## 根据观测计数矩阵生成可复现且为正的 JAGS 随机种子。
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

## 使用层级二次逻辑回归 JAGS 模型估计疗效概率。
## 输入为维度一致的疗效计数矩阵和样本量矩阵；返回后验均值概率矩阵。
eff.estimate <- function(
  y_mat,
  n_mat,
  n_burnin = EFF_MCMC_BURNIN,
  n_samples = EFF_MCMC_SAMPLES,
  n_chains = EFF_MCMC_CHAINS,
  thin = EFF_MCMC_THIN,
  use_cache = TRUE
) {
  if (!is.matrix(y_mat) || !is.matrix(n_mat)) {
    stop("y_mat and n_mat must be matrices.")
  }
  if (!identical(dim(y_mat), dim(n_mat))) {
    stop("y_mat and n_mat must have identical dimensions.")
  }
  if (any(!is.finite(y_mat)) || any(!is.finite(n_mat))) {
    stop("y_mat and n_mat must contain only finite values.")
  }
  if (any(y_mat < 0) || any(n_mat < 0) || any(y_mat > n_mat)) {
    stop("Counts must satisfy 0 <= y_mat <= n_mat.")
  }
  if (any(y_mat != floor(y_mat)) || any(n_mat != floor(n_mat))) {
    stop("y_mat and n_mat must contain integer counts.")
  }

  mcmc.settings <- c(n_burnin, n_samples, n_chains, thin)
  if (any(!is.finite(mcmc.settings)) ||
      any(mcmc.settings < 1) ||
      any(mcmc.settings != floor(mcmc.settings))) {
    stop("MCMC settings must be positive integers.")
  }

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
    agent_b_sq = rep(b^2, each = I),
    agent_a_by_b = rep(a, times = J) * rep(b, each = I)
  )

  model.text <- "
    model {
      for (d in 1:N) {
        y[d] ~ dbin(q[d], n[d])
        logit(q[d]) <- inprod(gamma[1:6], X[d, 1:6]) + delta[d]
        delta[d] ~ dnorm(0, tau_delta)
      }
      gamma[1] ~ dt(0, 0.01, 1)
      for (k in 2:6) {
        gamma[k] ~ dt(0, 0.16, 1)
      }
      sigma_delta ~ dunif(0.05, 2.5)
      tau_delta <- pow(sigma_delta, -2)
    }
  "

  base.seed <- .efficacy.jags.seed(y_mat, n_mat)
  initial.values <- lapply(seq_len(n_chains), function(chain) {
    chain.seed <- as.integer(
      ((base.seed - 1 + (chain - 1) * 104729) %% 2147483646) + 1
    )
    list(
      gamma = rep(0, 6),
      delta = rep(0, length(y_mat)),
      sigma_delta = 0.5,
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
  ## 直接从 rjags 读取监测数组。对大量拟合结果调用 coda::as.matrix()
  ## 偶尔会把 expression 对象保留为行名，导致 base::as.data.frame() 失败。
  ## 原始 mcarray 的维度为 q 单元格 x 迭代 x 链，不需要行名。
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

## 根据毒性（`p`）和疗效（`q`）计算设计效用。
## 毒性高于 `phi` 时，额外施加由 `w2` 加权的惩罚。
utility.fn <- function(p, q, w1, w2, phi) {
  q - w1 * p - w2 * p * (p > phi)
}
