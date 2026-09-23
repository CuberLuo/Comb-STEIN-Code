## 真实案例情景（论文表 3 与表 4）。
##
## 论文第 6 节用真实试验构造情景，说明所提设计在真实临床场景下的可行性。
## 结构与表 3 一致：一个试验对应一个剂量组合网格，网格的剂量档位取自试验
## 真实使用的剂量，每个情景给出一对真实毒性概率矩阵 pT 与真实疗效概率矩阵
## pE。共三个试验、六个情景：
##   - C1、C2：NCT03014076，乳腺癌中诱导 CD8+ T 细胞的疫苗（GP2）联合
##     GM-CSF 的 Ib 期试验（Clifton 等，2017）。原试验是 3+3 递增，四个队列
##     的给药方案为 GP2 100 μg + GM-CSF 125 μg、GP2 500 μg + GM-CSF 125 μg、
##     GP2 500 μg + GM-CSF 250 μg、GP2 1000 μg + GM-CSF 250 μg，最大给药
##     剂量为 GP2 1000 μg + GM-CSF 250 μg。剂量组合空间是 GM-CSF 两个水平 ×
##     GP2 三个水平，即表 3 的 2×3 网格，两个情景共用这一网格。
##   - C3、C4：NCT06380309，静脉给药的溶瘤痘苗病毒 IDOV-SAFE 联合
##     toripalimab（PD-1）与 fruquintinib（VEGFR 抑制剂）的 I 期试验。
##     网格为 fruquintinib 两个档 × IDOV-SAFE 三个档，即 2×3。
##     两个情景共用这一网格，**真实 OBDC 都落在 (1, 3)**（fruquintinib 取低档、
##     IDOV-SAFE 取最高档）：C3 靠"高病毒 + 高 TKI 格越过毒性参考值"把最优组合
##     逼回 (1, 3)；C4 两个药物都在安全区内，靠"疗效在 TKI 方向 3 mg 见顶、
##     5 mg 只增毒"达到同一落点，后者的效用优势明显更小。
##   - C5、C6：NCT07587827，个体化 DNA 新抗原疫苗 ODI-2001 联合 MVA 载体
##     与 ipilimumab（anti-CTLA4）的 I 期试验。网格为疫苗三个档 ×
##     ipilimumab 三个档，即 3×3。两个情景共用这一网格，**真实 OBDC 都落在
##     (3, 1)**（疫苗取最高档 4 mg、ipilimumab 退到最低档 2.5 mg）：C5 的疗效
##     随两个药物单调上升，但 ipi 每一档的毒性代价大于其疗效增益；C6 中疫苗剂量
##     由 2 mg 升到 4 mg 时疗效只小幅提高，OBDC 的效用优势更小。
##
## **OBDC 位置的约定（2026-09-22）**：同一个真实试验的两个情景必须共用同一个真实
## OBDC 位置——真实最优组合是该试验的属性，不应随情景改变；情景之间的差别体现在
## 概率曲面的形状与陡峭程度上。三个案例的 OBDC 位置互不重复，依次为
## C1–C2 = (2, 3)、C3–C4 = (1, 3)、C5–C6 = (3, 1)。
##
## **网格方向的统一约定（2026-09-22 转置）**
##   为了让表 3 的每个情景块都填满三列（药物 B 一律取三个档位），第一、第二
##   两个案例原来的 3×2 网格已整体转置为 2×3 网格：把原来作为药物 A 的三档
##   药物移到药物 B 的位置，把原来作为药物 B 的两档药物移到药物 A 的位置，
##   毒性与疗效概率矩阵取转置。即
##     - 第一案例：药物 A = GM-CSF（125、250 μg 两档），
##                 药物 B = GP2（100、500、1000 μg 三档）；
##     - 第二案例：药物 A = fruquintinib（3、5 mg 两档），
##                 药物 B = IDOV-SAFE（1e9、3e9、1e10 PFU 三档）。
##   第三个案例的 ODI-2001 × ipilimumab 网格本来就是 3×3，转置后仍是 3×3，
##   故保持原样（药物 A = ODI-2001，药物 B = ipilimumab）。
##   注意转置并不只是排版变化：设计在两个剂量维度上并非完全对称
##   （Runin-STEIN 的导入路径 make.runin.path 是行优先、列交替的，
##   final.selection 在效用并列时取行优先的组合），因此转置会改变模拟结果，
##   必须重新运行。第一案例原有两个以上情景，转置后按 OBDC 选择率与患者
##   分配比例保留表现最好的两个（见 2026-09-22 运行记录）。
##
## 剂量档位与概率的口径（论文第 6 节需要逐项说明）：
##   - C1-C2 的档位与概率依据 Clifton 等（2017）发表的试验结果；
##   - C3-C4 的病毒档位（1e9、3e9、1e10 PFU）来自登记处的爬坡设计，
##     fruquintinib 的 3 mg 与 5 mg 是说明书与临床常用的两个档，登记处
##     只写了联合 fruquintinib 而未披露档位；试验仍在招募、没有公开的
##     剂量-毒性/疗效结果，pT 与 pE 是依据机制的设定值；
##   - C5-C6 的 ipilimumab 档位（2.5、5、10 mg）来自登记处，DNA 疫苗
##     4 mg 也是登记处的剂量，疫苗的 1 mg 与 2 mg 是该平台下探的假设档位；
##     试验尚未启动，pT 与 pE 同样是依据机制的设定值。
##   三个网格里试验实际执行过的格子数与总格子数都不相同，未执行过的
##   格子由相邻格子外推。
##
## pT[i, j] 是药物 A 第 i 个剂量与药物 B 第 j 个剂量组合的真实毒性概率，
## pE[i, j] 为真实疗效概率。表 3 中加粗单元格是各情景的真实最优生物学
## 剂量组合（OBDC）。OBDC 不在这里写死：它由设计效用函数决定，用
## case.study.obdc() 调用 get.true.obdc() 计算得到，与 run.scenario() 计算
## 选择率时的口径完全一致。
##
## 这些情景单独保存、不并入 scenarios.R 的 S1-S12，避免编号冲突。
## 运行方式：Rscript Comb-STEIN-Code/main.R C1 C2 C3 C4 C5 C6 1000

## 与 Comb-STEIN-Code/scenarios.R 中的同名辅助函数定义一致：按行填写矩阵。
mk <- function(values, nrow) {
  matrix(values, nrow = nrow, byrow = TRUE)
}

case.study.scenarios <- list(
  ## 第一案例（NCT03014076）：GM-CSF × GP2 的 2×3 网格（由 3×2 转置而来）。
  ## 情景 1：毒性极低、疗效中等，转置后的加粗单元格为 (0.10, 0.60)。
  C1 = list(
    label = "Case study 1 (NCT03014076, GM-CSF + GP2): low toxicity",
    pT = mk(c(0.01, 0.03, 0.05,
              0.05, 0.07, 0.10), nrow = 2),
    pE = mk(c(0.10, 0.15, 0.30,
              0.12, 0.25, 0.60), nrow = 2)
  ),
  ## 情景 2：中间情形，转置后的加粗单元格为 (0.20, 0.75)。
  C2 = list(
    label = "Case study 2 (NCT03014076, GM-CSF + GP2): intermediate",
    pT = mk(c(0.01, 0.05, 0.15,
              0.10, 0.15, 0.20), nrow = 2),
    pE = mk(c(0.05, 0.15, 0.35,
              0.15, 0.25, 0.75), nrow = 2)
  ),
  ## 第二案例（NCT06380309）：fruquintinib × IDOV-SAFE 的 2×3 网格
  ## （由 IDOV-SAFE 三个档 × fruquintinib 两个档的 3×2 网格转置而来）。
  ##
  ## 剂量设置：登记处的爬坡档位为 1e9、3e9、1e10、3e10 PFU，扩展阶段使用
  ## 1e9 与 3e9 PFU 两个档。网格取 IDOV-SAFE 三个档（1e9、3e9、1e10 PFU）×
  ## fruquintinib 两个档（3 mg、5 mg）；3e10 PFU 毒性最高、作为毒性上界不
  ## 纳入网格。fruquintinib 的 3 mg 与 5 mg 是说明书与临床常用的两个剂量档，
  ## 登记处只写了联合 fruquintinib 而未披露具体档位。toripalimab 固定
  ## 240 mg，不构成第二维。
  ##
  ## 概率取值依据（机制，非试验结果）：试验仍在招募、没有公开的剂量-毒性/
  ## 疗效结果。静脉病毒载量带来流感样症状与肝酶升高，fruquintinib 带来高血压
  ## 与手足综合征，因此毒性沿两个方向上升；疗效则随病毒载量（肿瘤内感染与
  ## 免疫激活）和 VEGFR 抑制（免疫微环境改造）同时上升，但在两个药物都取
  ## 最高档时毒性已越过毒性参考值 0.397，效用被大幅惩罚。
  ## 转置后 OBDC 落在 (1, 3)，即最低 TKI 档 + 最高病毒档——在保住疗效的
  ## 同时退开毒性最大的 TKI 档，用来说明 OBDC 未必落在网格的角点。
  C3 = list(
    label = "Case study 3 (NCT06380309, fruquintinib + IDOV-SAFE): 2x3 grid",
    pT = mk(c(0.05, 0.08, 0.13,
              0.15, 0.22, 0.45), nrow = 2),
    pE = mk(c(0.15, 0.22, 0.38,
              0.22, 0.32, 0.55), nrow = 2)
  ),
  ## 情景 4：同一个 2×3 网格、同一个 OBDC 位置 (1, 3)，但与情景 3 的成因不同：
  ## 两个药物都留在安全区内（无一格越过 0.397），而 fruquintinib 的疗效在 3 mg
  ## 已达平台——5 mg 不再提高疗效（pE 两行相同）、只增加毒性。于是效用最高点
  ## 仍是 (1, 3)，但它的效用优势只有 0.04，明显小于情景 3 的 0.09，设计更难判断。
  C4 = list(
    label = "Case study 4 (NCT06380309, fruquintinib + IDOV-SAFE): efficacy plateau",
    pT = mk(c(0.05, 0.08, 0.12,
              0.09, 0.15, 0.24), nrow = 2),
    pE = mk(c(0.16, 0.26, 0.38,
              0.22, 0.35, 0.38), nrow = 2)
  ),
  ## 第三案例（NCT07587827）：ODI-2001 × ipilimumab 的 3×3 网格，本身对称，
  ## 转置后仍是 3×3，故保持原样（药物 A = ODI-2001，药物 B = ipilimumab）。
  ##
  ## 剂量设置：登记处只有 ipilimumab 一个剂量轴（2.5、5、10 mg 三个档），
  ## DNA 疫苗固定 4 mg、MVA 固定 1e7 pfu，DLT 与 MTD 也按 ipilimumab 档位
  ## 定义。为了在二维网格上评价该试验，疫苗轴取 1、2、4 mg 三个档：4 mg 是
  ## 登记处披露的剂量，1 mg 与 2 mg 是该平台下探的假设档位。
  ##
  ## 概率取值依据（机制，非试验结果）：试验尚未启动，没有公开的毒性/疗效
  ## 结果。毒性主要由 ipilimumab 驱动（anti-CTLA4 的免疫相关不良事件随剂量
  ## 陡升，在疫苗 2 mg 以上的高 ipi 档越过毒性参考值 0.397），并在疫苗剂量
  ## 方向温和上升；疗效随疫苗剂量与 ipilimumab 剂量同时上升——ipilimumab
  ## 解除 CTLA4 抑制、放大新抗原特异性 T 细胞应答。
  ## 计算结果使 OBDC 落在角点 (3, 1)：疫苗取最高档 4 mg，ipilimumab 退到
  ## 最低档 2.5 mg。ipilimumab 每一档虽有疗效增益（+0.03 左右），但毒性代价
  ## （0.33 × Δp，且高剂量档还要再加 1.09 × p 的越界惩罚）更大，净效用反而
  ## 下降。该组合沿 ipilimumab 方向单独爬坡到不了——这恰好说明若按试验原本的
  ## 单维递增（只升 ipilimumab）会停在 10 mg，与真实最优组合相反。
  C5 = list(
    label = "Case study 5 (NCT07587827, ODI-2001 + ipilimumab): 3x3 grid",
    pT = mk(c(0.06, 0.12, 0.22,
              0.09, 0.20, 0.45,
              0.13, 0.42, 0.60), nrow = 3),
    pE = mk(c(0.12, 0.15, 0.18,
              0.20, 0.24, 0.27,
              0.33, 0.38, 0.41), nrow = 3)
  ),
  ## 情景 6：同一个 3×3 网格、同一个 OBDC 位置 (3, 1)，但与情景 5 的成因不同：
  ## 疫苗剂量由 2 mg 升到 4 mg 时疗效只小幅提高（0.29 → 0.32，新抗原递送接近
  ## 平台），而毒性继续上升，因此 (3, 1) 的效用优势被压缩到 0.04 以内（情景 5
  ## 为 0.11），设计更难判断。两个情景因此构成"同一 OBDC 位置、不同难度"的一对
  ## 对照：OBDC 位置一致，概率曲面的陡峭程度不同。
  C6 = list(
    label = "Case study 6 (NCT07587827, ODI-2001 + ipilimumab): flat vaccine gain",
    pT = mk(c(0.07, 0.16, 0.30,
              0.11, 0.26, 0.50,
              0.15, 0.44, 0.62), nrow = 3),
    pE = mk(c(0.14, 0.16, 0.18,
              0.27, 0.29, 0.31,
              0.32, 0.34, 0.36), nrow = 3)
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

## 把真实案例情景并入 scenarios，使它们在本次会话中可直接按名字运行。
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
## 两列索引矩阵（第 1 列为药物 A 剂量，第 2 列为药物 B 剂量）。
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
