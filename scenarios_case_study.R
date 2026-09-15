## 真实案例情景（论文表 4 与表 5）。
##
## 来源：乳腺癌患者中诱导 CD8+ T 细胞的疫苗联合曲妥珠单抗的 Ib 期试验
## （NCT03014076，Clifton 等，2017）。论文第 6 节依据该试验的实际结果构造了
## 三个情景，用于说明所提设计在真实临床场景下的可行性。
##
## 剂量网格为 3×2：疫苗 A 三个剂量水平，药物 B 两个剂量水平。
## 论文表 4 的排版是“行 = 药物 B、列 = 疫苗 A”；本文件与 code/scenarios.R
## 保持一致，改为“行 = 疫苗 A、列 = 药物 B”，即 pT[i, j] 是疫苗 A 第 i 个
## 剂量与药物 B 第 j 个剂量组合的真实毒性概率，pE[i, j] 为真实疗效概率。
## 因此表 4 中第 j 行第 i 列的 (毒性, 疗效) 对应这里的 pT[i, j]、pE[i, j]。
##
## 表 4 中红色粗体单元格是各情景的真实最优生物学剂量组合（OBDC）。OBDC 不在
## 这里写死：它由设计效用函数决定，用 case.study.obdc() 调用 get.true.obdc()
## 计算得到，与 run.scenario() 计算选择率时的口径完全一致。对这三个情景算出
## 的 OBDC 都是最高剂量组合 (3, 2)，与论文“每种情景的 OBDC 均设为最高剂量”
## 的说明相符。
##
## 这些情景单独保存、不并入 scenarios.R 的 S1-S12，避免编号冲突。
## 运行方式：Rscript code/main.R C1 C2 C3 1000

## 与 code/scenarios.R 中的同名辅助函数定义一致：按行填写矩阵。
mk <- function(values, nrow) {
  matrix(values, nrow = nrow, byrow = TRUE)
}

case.study.scenarios <- list(
  ## 表 4 情景 1：疫苗组合毒性极低、疗效中等，红色粗体单元格为 (0.10, 0.60)。
  C1 = list(
    label = "Case study 1 (Table 4, S1)",
    pT = mk(c(0.01, 0.05,
              0.03, 0.07,
              0.05, 0.10), nrow = 3),
    pE = mk(c(0.10, 0.12,
              0.15, 0.25,
              0.30, 0.60), nrow = 3)
  ),
  ## 表 4 情景 2：毒性上升、疗效显著增强，红色粗体单元格为 (0.25, 0.80)。
  C2 = list(
    label = "Case study 2 (Table 4, S2)",
    pT = mk(c(0.10, 0.15,
              0.15, 0.20,
              0.20, 0.25), nrow = 3),
    pE = mk(c(0.01, 0.15,
              0.15, 0.30,
              0.25, 0.80), nrow = 3)
  ),
  ## 表 4 情景 3：中间情形，红色粗体单元格为 (0.20, 0.75)。
  C3 = list(
    label = "Case study 3 (Table 4, S3)",
    pT = mk(c(0.01, 0.10,
              0.05, 0.15,
              0.15, 0.20), nrow = 3),
    pE = mk(c(0.05, 0.15,
              0.15, 0.25,
              0.35, 0.75), nrow = 3)
  )
)

## 按名称取回真实案例情景，与 scenarios.R 的 build.pTpE() 用法一致。
build.case.study.pTpE <- function(name) {
  key <- as.character(name)
  if (!key %in% names(case.study.scenarios)) {
    stop(sprintf("Unknown case study scenario: %s", key))
  }
  case.study.scenarios[[key]]
}

## 把三个真实案例情景并入 scenarios，使它们在本次会话中可直接按名字运行。
## 默认不合并：main.R 的默认情景列表仍是 scenarios.R 中的 S1-S12。
register.case.study.scenarios <- function(envir = globalenv()) {
  scenarios <- get("scenarios", envir = envir)
  scenarios[names(case.study.scenarios)] <- case.study.scenarios
  assign("scenarios", scenarios, envir = envir)
  assign(
    "nrow_map",
    vapply(scenarios, function(s) nrow(s$pT), integer(1)),
    envir = envir
  )
  assign(
    "ncol_map",
    vapply(scenarios, function(s) ncol(s$pT), integer(1)),
    envir = envir
  )
  invisible(names(scenarios))
}

## 计算真实案例情景的真实 OBDC。
## 与 run.scenario() 判定选择率时使用同一套规则：先由 get.true.utility() 按
## 设计效用 utility.fn(p, q, w1, w2, phi) 算出真实效用矩阵，再把效用与最大值
## 相差不超过 OBDC_UTILITY_TOLERANCE 的组合全部判为 OBDC。因此这里不需要、
## 也不应该在情景里预先指定 OBDC。
##
## name 缺省时返回所有情景的命名列表；给出单个名字时直接返回该情景的
## 两列索引矩阵（第 1 列为疫苗 A 剂量，第 2 列为药物 B 剂量）。
case.study.obdc <- function(name = NULL, verbose = TRUE) {
  keys <- if (is.null(name)) {
    names(case.study.scenarios)
  } else {
    as.character(name)
  }
  unknown <- setdiff(keys, names(case.study.scenarios))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown case study scenario: %s", paste(unknown, collapse = ", ")))
  }
  obdcs <- lapply(keys, function(key) {
    scenario <- case.study.scenarios[[key]]
    get.true.obdc(scenario$pT, scenario$pE)
  })
  names(obdcs) <- keys
  if (verbose) {
    for (key in keys) {
      found <- obdcs[[key]]
      cat(sprintf(
        "%s (%s): 真实 OBDC = %s\n",
        key, case.study.scenarios[[key]]$label,
        paste(sprintf("(%d, %d)", found[, 1], found[, 2]), collapse = ", ")
      ))
    }
  }
  if (length(keys) == 1L) obdcs[[1]] else obdcs
}
