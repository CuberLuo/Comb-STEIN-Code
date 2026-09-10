## 根据论文设置构造的多样化模拟场景。
## 所有毒性和疗效概率最多保留两位小数。毒性随两个药物剂量均不下降。
## 按照效用规则，S1-S4 有 1 个 OBDC，S5-S8 有 2 个 OBDC，S9-S12 有 3 个 OBDC，
## S13-S14 有 1 个 OBDC；每个场景的完整 OBDC 位置集合都不同。
##
## 每个场景的毒性真值（pT）与疗效真值（pE）分开定义为独立矩阵，
## 便于单独阅读、修改和复用。行对应药物 A 的剂量水平 i，
## 列对应药物 B 的剂量水平 j。

## 场景 S1（3 x 5）
pT.true1 <- matrix(c(
  0.05, 0.10, 0.15, 0.20, 0.30,
  0.10, 0.15, 0.20, 0.30, 0.40,
  0.15, 0.20, 0.30, 0.40, 0.50
), nrow = 3, byrow = TRUE)
pE.true1 <- matrix(c(
  0.35, 0.80, 0.65, 0.50, 0.40,
  0.40, 0.60, 0.55, 0.45, 0.35,
  0.45, 0.50, 0.40, 0.35, 0.30
), nrow = 3, byrow = TRUE)

## 场景 S2（3 x 5）
pT.true2 <- matrix(c(
  0.05, 0.10, 0.15, 0.20, 0.30,
  0.10, 0.15, 0.20, 0.25, 0.35,
  0.15, 0.20, 0.25, 0.35, 0.45
), nrow = 3, byrow = TRUE)
pE.true2 <- matrix(c(
  0.25, 0.35, 0.50, 0.65, 0.70,
  0.30, 0.45, 0.60, 0.85, 0.65,
  0.35, 0.50, 0.55, 0.60, 0.50
), nrow = 3, byrow = TRUE)

## 场景 S3（3 x 5）
pT.true3 <- matrix(c(
  0.05, 0.10, 0.15, 0.25, 0.35,
  0.10, 0.15, 0.20, 0.30, 0.40,
  0.15, 0.25, 0.30, 0.40, 0.50
), nrow = 3, byrow = TRUE)
pE.true3 <- matrix(c(
  0.20, 0.30, 0.40, 0.50, 0.45,
  0.30, 0.45, 0.60, 0.65, 0.55,
  0.40, 0.85, 0.70, 0.55, 0.40
), nrow = 3, byrow = TRUE)

## 场景 S4（3 x 5）
pT.true4 <- matrix(c(
  0.05, 0.10, 0.15, 0.20, 0.30,
  0.10, 0.15, 0.20, 0.30, 0.40,
  0.15, 0.20, 0.30, 0.40, 0.50
), nrow = 3, byrow = TRUE)
pE.true4 <- matrix(c(
  0.20, 0.35, 0.50, 0.70, 0.90,
  0.25, 0.40, 0.55, 0.65, 0.75,
  0.30, 0.45, 0.60, 0.70, 0.80
), nrow = 3, byrow = TRUE)

## 场景 S5（3 x 5）
pT.true5 <- matrix(c(
  0.05, 0.10, 0.20, 0.30, 0.40,
  0.10, 0.20, 0.30, 0.40, 0.50,
  0.20, 0.30, 0.40, 0.50, 0.60
), nrow = 3, byrow = TRUE)
pE.true5 <- matrix(c(
  0.25, 0.45, 0.80, 0.65, 0.50,
  0.40, 0.80, 0.70, 0.55, 0.40,
  0.55, 0.65, 0.60, 0.45, 0.30
), nrow = 3, byrow = TRUE)

## 场景 S6（3 x 4）
pT.true6 <- matrix(c(
  0.05, 0.10, 0.20, 0.30,
  0.10, 0.20, 0.30, 0.40,
  0.20, 0.30, 0.40, 0.50
), nrow = 3, byrow = TRUE)
pE.true6 <- matrix(c(
  0.25, 0.45, 0.65, 0.85,
  0.40, 0.60, 0.85, 0.70,
  0.50, 0.65, 0.70, 0.55
), nrow = 3, byrow = TRUE)

## 场景 S7（3 x 4）
pT.true7 <- matrix(c(
  0.05, 0.10, 0.20, 0.30,
  0.10, 0.20, 0.30, 0.40,
  0.20, 0.30, 0.40, 0.50
), nrow = 3, byrow = TRUE)
pE.true7 <- matrix(c(
  0.25, 0.45, 0.60, 0.55,
  0.40, 0.80, 0.70, 0.60,
  0.80, 0.65, 0.55, 0.45
), nrow = 3, byrow = TRUE)

## 场景 S8（3 x 4）
pT.true8 <- matrix(c(
  0.05, 0.10, 0.15, 0.25,
  0.10, 0.15, 0.25, 0.35,
  0.15, 0.25, 0.35, 0.45
), nrow = 3, byrow = TRUE)
pE.true8 <- matrix(c(
  0.20, 0.35, 0.50, 0.60,
  0.30, 0.50, 0.65, 0.85,
  0.40, 0.60, 0.85, 0.70
), nrow = 3, byrow = TRUE)

## 场景 S9（3 x 3）
pT.true9 <- matrix(c(
  0.05, 0.20, 0.30,
  0.20, 0.20, 0.35,
  0.30, 0.35, 0.45
), nrow = 3, byrow = TRUE)
pE.true9 <- matrix(c(
  0.30, 0.80, 0.60,
  0.80, 0.80, 0.65,
  0.50, 0.60, 0.55
), nrow = 3, byrow = TRUE)

## 场景 S10（3 x 3）
pT.true10 <- matrix(c(
  0.05, 0.10, 0.25,
  0.10, 0.25, 0.35,
  0.25, 0.35, 0.45
), nrow = 3, byrow = TRUE)
pE.true10 <- matrix(c(
  0.30, 0.55, 0.85,
  0.55, 0.85, 0.70,
  0.85, 0.70, 0.60
), nrow = 3, byrow = TRUE)

## 场景 S11（2 x 3）
pT.true11 <- matrix(c(
  0.10, 0.10, 0.25,
  0.10, 0.25, 0.40
), nrow = 2, byrow = TRUE)
pE.true11 <- matrix(c(
  0.75, 0.75, 0.60,
  0.75, 0.60, 0.50
), nrow = 2, byrow = TRUE)

## 场景 S12（2 x 3）
pT.true12 <- matrix(c(
  0.05, 0.15, 0.30,
  0.15, 0.30, 0.30
), nrow = 2, byrow = TRUE)
pE.true12 <- matrix(c(
  0.25, 0.50, 0.85,
  0.45, 0.85, 0.85
), nrow = 2, byrow = TRUE)

## 场景 S13（3 x 4，OBDC 位于 (i=2, j=2)）
pT.true13 <- matrix(c(
  0.05, 0.10, 0.15, 0.19,
  0.11, 0.15, 0.20, 0.25,
  0.16, 0.21, 0.25, 0.30
), nrow = 3, byrow = TRUE)
pE.true13 <- matrix(c(
  0.18, 0.28, 0.38, 0.48,
  0.25, 0.77, 0.45, 0.55,
  0.32, 0.42, 0.52, 0.62
), nrow = 3, byrow = TRUE)

## 场景 S14（3 x 4，OBDC 位于 (i=1, j=3)）
pT.true14 <- matrix(c(
  0.05, 0.10, 0.15, 0.19,
  0.11, 0.15, 0.20, 0.25,
  0.16, 0.21, 0.25, 0.30
), nrow = 3, byrow = TRUE)
pE.true14 <- matrix(c(
  0.20, 0.30, 0.75, 0.50,
  0.26, 0.36, 0.46, 0.56,
  0.33, 0.43, 0.53, 0.63
), nrow = 3, byrow = TRUE)

## 将各场景的毒性/疗效真值汇总为按场景名（S1-S14）索引的列表，
## 下游仿真代码（main.R、run_scenarios.R、Cai、SFU、BOINETC 等）继续
## 通过 names(scenarios) 与 build.pTpE(name) 统一访问。
scenarios <- list(
  S1 = list(pT = pT.true1, pE = pE.true1),
  S2 = list(pT = pT.true2, pE = pE.true2),
  S3 = list(pT = pT.true3, pE = pE.true3),
  S4 = list(pT = pT.true4, pE = pE.true4),
  S5 = list(pT = pT.true5, pE = pE.true5),
  S6 = list(pT = pT.true6, pE = pE.true6),
  S7 = list(pT = pT.true7, pE = pE.true7),
  S8 = list(pT = pT.true8, pE = pE.true8),
  S9 = list(pT = pT.true9, pE = pE.true9),
  S10 = list(pT = pT.true10, pE = pE.true10),
  S11 = list(pT = pT.true11, pE = pE.true11),
  S12 = list(pT = pT.true12, pE = pE.true12),
  S13 = list(pT = pT.true13, pE = pE.true13),
  S14 = list(pT = pT.true14, pE = pE.true14)
)

## 每个场景的网格维度，从真值矩阵直接导出，避免重复维护。
nrow_map <- vapply(scenarios, function(s) nrow(s$pT), integer(1))
ncol_map <- vapply(scenarios, function(s) ncol(s$pT), integer(1))

## 返回命名场景的毒性和疗效真值矩阵。
build.pTpE <- function(name) {
  key <- as.character(name)
  list(pT = scenarios[[key]]$pT, pE = scenarios[[key]]$pE)
}
