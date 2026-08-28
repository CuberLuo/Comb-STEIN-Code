## 后验监测、可接受性更新和剂量选择规则。

## 在 Beta(1, 1) 先验下计算 Pr（毒性概率 > `cut`）。
## `x` 和 `n` 可以是标量、向量或维度相容的矩阵。
pr.overtox <- function(x, n, cut) {
  1 - pbeta(cut, x + 1, n - x + 1)
}

## 在 Beta(1, 1) 先验下计算 Pr（疗效概率 < `cut`）。
## `y` 和 `n` 可以是标量、向量或维度相容的矩阵。
pr.ineff <- function(y, n, cut) {
  pbeta(cut, y + 1, n - y + 1)
}

## 在 Beta(1, 1) 先验下计算 Pr（疗效概率 > `cut`）。
## `y` 和 `n` 可以是标量、向量或维度相容的矩阵。
pr.eff.gt <- function(y, n, cut) {
  1 - pbeta(cut, y + 1, n - y + 1)
}

## 标记超过毒性或无效性阈值的已观测剂量组合。
## 返回名为 `overtox` 和 `ineff` 的逻辑矩阵。
monitor.flags <- function(x_mat, y_mat, n_tox, n_eff = n_tox) {
  observed_tox <- n_tox > 0
  observed_eff <- n_eff > 0
  overtox <- observed_tox &
    (pr.overtox(x_mat, n_tox, phi) > Tthresh)
  ineff <- observed_eff &
    (pr.ineff(y_mat, n_eff, psi1) > Ethresh)
  list(overtox = overtox, ineff = ineff)
}

## 根据场景的真实毒性率和疗效率计算真实效用矩阵。
get.true.utility <- function(pT, pE) {
  utility.fn(pT, pE, w1, w2, phi)
}

## 按效用最大值识别所有真实最优生物剂量组合。
## 返回包含行、列索引的两列矩阵，并保留完全并列的组合。
get.true.obdc <- function(pT, pE) {
  utility <- get.true.utility(pT, pE)
  which(
    abs(utility - max(utility)) <= sqrt(.Machine$double.eps),
    arr.ind = TRUE
  )
}

## 将后验安全性和无效性规则应用于可接受性矩阵。
## 过度毒性的单元格还会排除其东北方向的组合；
## 疗效不足的单元格只排除自身。
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

## 构造下一次剂量组合分配的局部候选集合。
## 当前观测毒性率超过 phi 时，不允许升级。
alt.set.C <- function(i, j, I, J, p_ij_hat) {
  candidates <- list()
  if (p_ij_hat > phi) {# 高毒性则候选剂量集合为2个剂量
    if (i > 1) candidates <- c(candidates, list(c(i - 1, j)))
    if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
  } else {#低毒性则候选剂量集合为5个剂量
    if (i > 1) candidates <- c(candidates, list(c(i - 1, j)))
    if (j > 1) candidates <- c(candidates, list(c(i, j - 1)))
    candidates <- c(candidates, list(c(i, j)))
    if (j < J) candidates <- c(candidates, list(c(i, j + 1)))
    if (i < I) candidates <- c(candidates, list(c(i + 1, j)))
  }
  candidates
}

## 选择疗效后验概率最大的候选组合。
## 使用确定性的剂量索引项打破后验得分相同的情况。
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

## 在已治疗且可接受的组合中选择拟合效用最大的最终组合。
## 返回 c(row, column)；没有符合条件的组合时返回两个 NA 整数。
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

## 构造入组阶段使用的行列交替升级路径。
## 返回列表中最多包含 `ncohort` 个整数索引对。
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
