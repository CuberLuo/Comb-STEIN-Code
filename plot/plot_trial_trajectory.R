## 单次组合用药试验的过程记录。
##
## 组合试验的剂量是二维的（药物 A 的剂量 + 药物 B 的剂量），过程记录不适合
## 画成单药那样的剂量阶梯图，因此逐队列的过程用表格记录，图形只保留两张：
##   - 表格：每个队列分配到的剂量组合、决策时间、该队列的毒性和疗效结局，
##     输出 output/tables/ 下的 CSV、Markdown 与抬头说明三个文件；
##   - 图 A：两个药物各自随队列入组的剂量阶梯；
##   - 图 B：延迟结局（TITE）与等待完整结局两种运行方式的日历时间对比。
##
## 面向期刊投稿的几点约定：
##   * 图形尺寸与字号独立：文字字号固定，用 --width/--height 按期刊栏宽出图，
##     因此缩到目标栏宽后字号仍然可读；
##   * PDF 走 cairo 设备（字体嵌入完整），PNG 默认 600 dpi；
##   * 线条颜色取色盲友好的蓝／朱红，并同时用实线／虚线区分，黑白打印也能分辨；
##   * 剂量单位用 μg 而不是 ug；
##   * 图内的解释性文字尽量短，完整的图注（含随机种子、软件版本等可重复性
##     信息）单独写到 caption 文件，便于直接放进论文。
##
## 用法：Rscript Comb-STEIN-Code/plot/plot_trial_trajectory.R [项目根目录]
##       [--scenario=C1] [--method=CombSTEIN-L] [--seed=2026]
##       [--width=10.5] [--height=5.2] [--dpi=600]
##       [--figure-dir=output/figures] [--table-dir=output/tables]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

## 取出 --name=value 形式的选项，其余位置参数留给调用方。
.pull.option <- function(args, name) {
  pattern <- paste0("^--", name, "=")
  hit <- grep(pattern, args, value = TRUE)
  if (length(hit) == 0) {
    return(list(value = NULL, args = args))
  }
  list(value = sub(pattern, "", hit[[1]]), args = args[!grepl(pattern, args)])
}

option.value <- function(args, name, default) {
  pulled <- .pull.option(args, name)
  list(
    value = if (is.null(pulled$value)) default else pulled$value,
    args = pulled$args
  )
}

scenario.opt <- option.value(args, "scenario", "C1")
args <- scenario.opt$args
method.opt <- option.value(args, "method", "CombSTEIN-L")
args <- method.opt$args
seed.opt <- option.value(args, "seed", NULL)
args <- seed.opt$args
width.opt <- option.value(args, "width", "10.5")
args <- width.opt$args
height.opt <- option.value(args, "height", "5.2")
args <- height.opt$args
dpi.opt <- option.value(args, "dpi", "600")
args <- dpi.opt$args
figure.opt <- option.value(args, "figure-dir", "output/figures")
args <- figure.opt$args
table.opt <- option.value(args, "table-dir", "output/tables")
args <- table.opt$args
## 延迟结局的时间设置可覆盖，默认沿用 config.R 的 DEFAULT_DELAY_CONTROL，
## 这样论文里的过程图能与主分析使用同一套评估窗与入组速度。
accrual.opt <- option.value(args, "accrual", NULL)
args <- accrual.opt$args
tox.opt <- option.value(args, "tox-window", NULL)
args <- tox.opt$args
eff.opt <- option.value(args, "eff-window", NULL)
args <- eff.opt$args
susp.opt <- option.value(args, "suspension", NULL)
args <- susp.opt$args

scenario <- scenario.opt$value
method <- method.opt$value
figure.width <- as.numeric(width.opt$value)
figure.height <- as.numeric(height.opt$value)
figure.dpi <- as.numeric(dpi.opt$value)

project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

## 相对路径按项目根目录展开，绝对路径原样使用。
resolve.dir <- function(value) {
  if (grepl("^(/|[A-Za-z]:)", value)) value else file.path(project_root, value)
}
figure_dir <- resolve.dir(figure.opt$value)
table_dir <- resolve.dir(table.opt$value)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

## 加载模拟代码。模块顺序与 main.R 保持一致。
code_dir <- file.path(project_root, "Comb-STEIN-Code")
for (module in c(
  "config.R", "time_to_event.R", "isotonic_estimators.R", "dose_rules.R",
  "scenarios.R", "scenarios_case_study.R", "simulate_stein.R"
)) {
  sys.source(file.path(code_dir, module), envir = globalenv())
}

seed <- if (is.null(seed.opt$value)) DEFAULT_SEED else as.integer(seed.opt$value)

## 用命令行覆盖（若给出）组装本次模拟使用的延迟设置。
delay.control <- DEFAULT_DELAY_CONTROL
if (!is.null(accrual.opt$value)) delay.control$accrual <- as.numeric(accrual.opt$value)
if (!is.null(tox.opt$value)) delay.control$tox_window <- as.numeric(tox.opt$value)
if (!is.null(eff.opt$value)) delay.control$eff_window <- as.numeric(eff.opt$value)
if (!is.null(susp.opt$value)) delay.control$suspension <- as.numeric(susp.opt$value)
days_per_month <- 30.4375

## R 在中文 Windows 上通常运行在 C 语言环境，此时直接用 writeLines/write.csv
## 写非 ASCII 字符会被转义或丢失（μg 会写成 <U+03BC>g）。这里统一用二进制
## 连接按 UTF-8 字节写出，保证表格与图注文件里的单位正确。
micro.sign <- intToUtf8(0x3bc)
write.utf8 <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con, useBytes = TRUE)
}
## write.csv 会把 UTF-8 字符串翻译成当前本地编码（这里是 C / ASCII），μ 同样
## 会变成转义串，因此 CSV 也逐行拼好再按字节写出。
csv.field <- function(x) {
  x <- as.character(x)
  ifelse(grepl('["\r\n]', x), paste0('"', gsub('"', '""', x), '"'), x)
}
to.csv.lines <- function(data) {
  ## 先整体转成字符矩阵：直接对 data.frame 用 apply() 会对数值列做对齐补空格。
  chars <- as.matrix(data.frame(
    lapply(data, as.character),
    stringsAsFactors = FALSE, check.names = FALSE
  ))
  c(
    paste(csv.field(colnames(chars)), collapse = ","),
    apply(chars, 1, function(row) paste(csv.field(row), collapse = ","))
  )
}

## 各真实案例情景的剂量档位标签，与 Comb-STEIN-Code/scenarios_case_study.R 的网格一一
## 对应。这里的文字只用于作图与制表，不参与模拟。
## 注意（2026-09-22 转置）：第一、第二案例的网格已把两个药物对调，
## 药物 A 现在是两档的药物（GM-CSF、fruquintinib），药物 B 是三档的药物
## （GP2、IDOV-SAFE）。这里的顺序必须与 Comb-STEIN-Code/scenarios_case_study.R 一致。
dose.labels <- list(
  ## 第一案例（GM-CSF × GP2，2×3）：药物 A = GM-CSF 两档，药物 B = GP2 三档。
  C1 = list(drugA = "GM-CSF", a = c(125, 250), a.unit = paste0(micro.sign, "g"),
            drugB = "GP2", b = c(100, 500, 1000), b.unit = paste0(micro.sign, "g")),
  C2 = list(drugA = "GM-CSF", a = c(125, 250), a.unit = paste0(micro.sign, "g"),
            drugB = "GP2", b = c(100, 500, 1000), b.unit = paste0(micro.sign, "g")),
  ## 第二案例（fruquintinib × IDOV-SAFE，2×3）：药物 A = fruquintinib 两档，
  ## 药物 B = IDOV-SAFE 三个档（1e9、3e9、1e10 PFU）。作图时只能有一个连续
  ## 刻度，因此刻度值取 1 / 3 / 10，同时用 b.text 给出完整的字符串标签，
  ## 表格与图注里的剂量写法就不会出现 1e+09 这种科学计数法。
  C3 = list(drugA = "Fruquintinib", a = c(3, 5), a.unit = "mg",
            drugB = "IDOV-SAFE", b = c(1, 3, 10),
            b.unit = paste0("\u00d710", "\u2079 PFU"),
            b.text = c("IDOV-SAFE 1\u00d710\u2079 PFU",
                       "IDOV-SAFE 3\u00d710\u2079 PFU",
                       "IDOV-SAFE 1\u00d710\u00b9\u2070 PFU")),
  C4 = list(drugA = "Fruquintinib", a = c(3, 5), a.unit = "mg",
            drugB = "IDOV-SAFE", b = c(1, 3, 10),
            b.unit = paste0("\u00d710", "\u2079 PFU"),
            b.text = c("IDOV-SAFE 1\u00d710\u2079 PFU",
                       "IDOV-SAFE 3\u00d710\u2079 PFU",
                       "IDOV-SAFE 1\u00d710\u00b9\u2070 PFU")),
  ## 第三案例（ODI-2001 × ipilimumab，3×3）本身对称，转置后不变。
  C5 = list(drugA = "ODI-2001", a = c(1, 2, 4), a.unit = "mg",
            drugB = "Ipilimumab", b = c(2.5, 5, 10), b.unit = "mg"),
  C6 = list(drugA = "ODI-2001", a = c(1, 2, 4), a.unit = "mg",
            drugB = "Ipilimumab", b = c(2.5, 5, 10), b.unit = "mg")
)
if (!scenario %in% names(dose.labels)) {
  stop(sprintf("情景 %s 缺少剂量档位标签，请在 dose.labels 中补充", scenario))
}
doses <- dose.labels[[scenario]]

## 用固定的随机种子分别模拟一次延迟结局试验和一次完整结局试验。
truth <- build.case.study.pTpE(scenario)
set.seed(seed)
tite.trial <- simulate.trial.comb.stein(
  truth$pT, truth$pE, method, delayed = TRUE,
  delay.control = delay.control, trace = TRUE
)
set.seed(seed)
complete.trial <- simulate.trial.comb.stein(
  truth$pT, truth$pE, method, delayed = FALSE,
  delay.control = delay.control, trace = TRUE
)

## 真实最优生物学剂量组合（OBDC）由设计效用函数算出，不在情景里写死。
obdc <- get.true.obdc(truth$pT, truth$pE)

dose.text <- function(dose, drug, values, unit, texts = NULL) {
  ## texts 给出该轴每一档的完整写法（用于 PFU 这类不便于用连续刻度表示的
  ## 剂量）；没有给出时按“药物名 + 数值 + 单位”拼接。
  if (!is.null(texts)) return(texts[dose])
  sprintf("%s %g %s", drug, values[dose], unit)
}

## ---------------------------------------------------------------------------
## 逐队列的过程表。
## 每位患者一行（来自 TITE 轨迹），先按队列汇总成每个队列的结局计数。
trace <- tite.trial$trace
process <- do.call(rbind, lapply(split(trace, trace$cohort), function(chunk) {
  data.frame(
    cohort = chunk$cohort[[1]],
    dose_a_level = chunk$dose_i[[1]],
    dose_b_level = chunk$dose_j[[1]],
    dose_a = dose.text(chunk$dose_i[[1]], doses$drugA, doses$a, doses$a.unit,
                       texts = doses$a.text),
    dose_b = dose.text(chunk$dose_j[[1]], doses$drugB, doses$b, doses$b.unit,
                       texts = doses$b.text),
    enroll_start = min(chunk$enroll),
    enroll_end = max(chunk$enroll),
    decision_day = chunk$decision_time[[1]],
    decision_month = chunk$decision_time[[1]] / days_per_month,
    n = nrow(chunk),
    dlt = sum(chunk$tox_state == "dlt"),
    pending_dlt = sum(chunk$tox_state == "pending"),
    benefit = sum(chunk$eff_state == "benefit"),
    pending_eff = sum(chunk$eff_state == "pending"),
    stringsAsFactors = FALSE
  )
}))
process <- process[order(process$cohort), ]
process$no_dlt <- process$n - process$dlt - process$pending_dlt
process$no_benefit <- process$n - process$benefit - process$pending_eff
## 天数是时间到事件模拟里的连续值，表里保留一位小数即可。
process$enroll_start <- round(process$enroll_start, 1)
process$enroll_end <- round(process$enroll_end, 1)
process$decision_day <- round(process$decision_day, 1)
process$decision_month <- round(process$decision_month, 1)
process <- process[, c(
  "cohort", "dose_a_level", "dose_b_level", "dose_a", "dose_b",
  "enroll_start", "enroll_end",
  "decision_day", "decision_month", "n",
  "dlt", "pending_dlt", "no_dlt",
  "benefit", "pending_eff", "no_benefit"
)]
rownames(process) <- NULL

## 试验层面的汇总信息，写在表格文件抬头。
tite.months <- tite.trial$duration / days_per_month
complete.months <- complete.trial$duration / days_per_month
saved.months <- complete.months - tite.months
if (any(is.na(tite.trial$sel))) {
  selected.text <- "none (early termination without a recommendation)"
} else {
  selected.text <- sprintf(
    "%s + %s",
    dose.text(tite.trial$sel[1], doses$drugA, doses$a, doses$a.unit,
              texts = doses$a.text),
    dose.text(tite.trial$sel[2], doses$drugB, doses$b, doses$b.unit,
              texts = doses$b.text)
  )
}
obdc.text <- paste(
  sprintf(
    "%s + %s",
    dose.text(obdc[, 1], doses$drugA, doses$a, doses$a.unit,
              texts = doses$a.text),
    dose.text(obdc[, 2], doses$drugB, doses$b, doses$b.unit,
              texts = doses$b.text)
  ),
  collapse = "; "
)
header <- c(
  sprintf("Scenario: %s (%s)", scenario, method),
  sprintf("Random seed: %d", seed),
  sprintf("True OBDC: %s", obdc.text),
  sprintf("Combination selected at the end: %s", selected.text),
  sprintf(
    "Calendar time: %.1f months with delayed outcomes (TITE), %.1f months waiting for complete outcomes, %.1f months saved",
    tite.months, complete.months, saved.months
  ),
  sprintf(
    "Enrolment %g days per patient, toxicity window %g days, efficacy window %g days, suspension %.0f%%",
    delay.control$accrual,
    delay.control$tox_window,
    delay.control$eff_window,
    100 * delay.control$suspension
  )
)

file.stem <- sprintf("trial_cohort_process_%s_%s", scenario, method)
write.utf8(header, file.path(table_dir, paste0(file.stem, ".txt")))
write.utf8(
  to.csv.lines(process),
  file.path(table_dir, paste0(file.stem, ".csv"))
)

## Markdown 版本：便于直接粘进论文或文档。
write.utf8(
  c(
    sprintf("### %s - %s (seed %d)", scenario, method, seed),
    "",
    paste0("- ", header),
    "",
    paste(
      "| cohort | drug A dose | drug B dose | cohort entry (day) |",
      "dose decision (day) | dose decision (month) | n | DLT | pending DLT |",
      "no DLT | benefit | pending benefit | no benefit |"
    ),
    paste(
      "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
      "--- | --- | --- |"
    ),
    vapply(seq_len(nrow(process)), function(r) {
      sprintf(
        "| %d | %s | %s | %g-%g | %g | %.1f | %d | %d | %d | %d | %d | %d | %d |",
        process$cohort[r], process$dose_a[r], process$dose_b[r],
        process$enroll_start[r], process$enroll_end[r],
        process$decision_day[r], process$decision_month[r],
        process$n[r], process$dlt[r], process$pending_dlt[r],
        process$no_dlt[r], process$benefit[r],
        process$pending_eff[r], process$no_benefit[r]
      )
    }, character(1))
  ),
  file.path(table_dir, paste0(file.stem, ".md"))
)

## ---------------------------------------------------------------------------
## 统一的期刊用主题：白底、细边框、极浅网格，字号固定不随图幅缩放。
journal.theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.border = element_rect(colour = "grey40", linewidth = 0.4),
    axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.title = element_text(size = 9.5, colour = "black"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines"),
    plot.margin = margin(5, 8, 5, 5)
  )

## ---------------------------------------------------------------------------
## 面板 A：两个药物各自随队列入组的剂量阶梯。
## 两个药物共用横轴但剂量量纲不同，因此拆成上下两张子图，各自把纵轴刻度
## 落在该药物真实的剂量档位上。
assignments <- process[, c("cohort", "dose_a_level", "dose_b_level")]
cohort.max <- max(assignments$cohort)
## 两个药物的单位不同（例如 IDOV-SAFE 用 ×10⁹ PFU、fruquintinib 用 mg），
## 因此单位必须随各自的轨道传入，不能共用 doses$a.unit。
track.list <- list(
  list(
    drug = doses$drugA,
    unit = doses$a.unit,
    dose = doses$a[assignments$dose_a_level],
    breaks = doses$a
  ),
  list(
    drug = doses$drugB,
    unit = doses$b.unit,
    dose = doses$b[assignments$dose_b_level],
    breaks = doses$b
  )
)

build.track <- function(track, show.x) {
  frame <- data.frame(cohort = assignments$cohort, dose = track$dose)
  ## 只在剂量发生变化（或首次出现）的队列上标注具体剂量，避免同一档位重复标。
  frame$change <- c(TRUE, frame$dose[-1] != frame$dose[-nrow(frame)])
  pad <- 0.14 * diff(range(track$breaks))
  plot <- ggplot(frame, aes(x = cohort, y = dose)) +
    geom_line(colour = "grey60", linewidth = 0.4) +
    geom_point(size = 1.7, colour = "black") +
    geom_text(
      data = subset(frame, change),
      aes(label = dose), vjust = -0.9, size = 2.5, colour = "grey25"
    ) +
    scale_x_continuous(breaks = seq_len(cohort.max)) +
    scale_y_continuous(
      breaks = track$breaks,
      limits = range(track$breaks) + c(-pad, pad),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    labs(y = sprintf("%s dose (%s)", track$drug, track$unit)) +
    journal.theme
  if (show.x) {
    plot + labs(x = "Cohort entry")
  } else {
    plot +
      labs(x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
}

track.a <- build.track(track.list[[1]], show.x = FALSE) +
  labs(title = sprintf("A: Dose escalation of each agent (%s, %s)", scenario, method))
track.b <- build.track(track.list[[2]], show.x = TRUE)

panel.a <- track.a / track.b +
  plot_layout(heights = c(1, 1), axes = "collect_x")

## ---------------------------------------------------------------------------
## 面板 B：两种运行方式的日历时间对比。
## 颜色取色盲友好的蓝／朱红，同时用实线／虚线区分，黑白打印同样可辨。
timeline <- rbind(
  data.frame(
    mode = "Delayed outcomes (TITE)",
    months = sort(unique(trace$decision_time)) / days_per_month,
    finish = tite.months
  ),
  data.frame(
    mode = "Complete outcomes",
    months = sort(unique(complete.trial$trace$decision_time)) / days_per_month,
    finish = complete.months
  )
)
timeline$mode <- factor(
  timeline$mode,
  levels = c("Delayed outcomes (TITE)", "Complete outcomes")
)
mode.style <- c(
  "Delayed outcomes (TITE)" = "#0072B2",
  "Complete outcomes" = "#D55E00"
)
mode.linetype <- c(
  "Delayed outcomes (TITE)" = "solid",
  "Complete outcomes" = "dashed"
)
finish.labels <- unique(timeline[, c("mode", "finish")])
saved.arrow <- data.frame(
  x = tite.months,
  xend = complete.months,
  y = length(levels(timeline$mode)) + 0.4
)

panel.b <- ggplot(timeline, aes(y = mode, colour = mode, linetype = mode)) +
  geom_segment(
    data = finish.labels,
    aes(x = 0, xend = finish, y = mode, yend = mode),
    linewidth = 0.7
  ) +
  geom_point(aes(x = months), size = 1.7, show.legend = FALSE) +
  geom_point(
    data = finish.labels, aes(x = finish, y = mode),
    size = 3, shape = 21, fill = "white", stroke = 0.9, show.legend = FALSE
  ) +
  geom_text(
    data = finish.labels,
    aes(x = finish, y = mode, label = sprintf("%.1f", finish)),
    hjust = -0.25, vjust = -0.9, size = 2.6, colour = "black",
    show.legend = FALSE
  ) +
  geom_segment(
    data = saved.arrow,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE, colour = "grey20", linewidth = 0.4,
    arrow = arrow(ends = "both", length = unit(0.15, "cm"), type = "closed")
  ) +
  geom_text(
    data = saved.arrow,
    aes(x = (x + xend) / 2, y = y, label = sprintf("%.1f months saved", xend - x)),
    inherit.aes = FALSE, vjust = -0.7, size = 2.8, colour = "grey20"
  ) +
  scale_colour_manual(values = mode.style, guide = "none") +
  scale_linetype_manual(values = mode.linetype, guide = "none") +
  scale_x_continuous(
    limits = c(0, max(finish.labels$finish) * 1.30),
    breaks = scales::pretty_breaks(5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.9))) +
  labs(
    title = "B: Calendar time saved",
    x = "Calendar time (months)",
    y = NULL
  ) +
  journal.theme

figure <- panel.a | panel.b

figure.file <- sprintf("fig_trial_trajectory_%s_%s", scenario, method)
ggsave(
  file.path(figure_dir, paste0(figure.file, ".pdf")),
  figure, width = figure.width, height = figure.height,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(figure_dir, paste0(figure.file, ".png")),
  figure, width = figure.width, height = figure.height, dpi = figure.dpi
)

## ---------------------------------------------------------------------------
## 图注：期刊要求图片本身不承载解释性文字，完整说明放在图注里，并写明
## 可重复性信息（种子、软件版本、模拟设置）。
caption <- c(
  sprintf(
    "Figure.  Process of a single simulated combination trial under the %s design, scenario %s (%s, %s).",
    method, scenario, doses$drugA, doses$drugB
  ),
  "(A) Dose of each agent at every cohort entry. One point per cohort; the numbers give the dose actually used.",
  "(B) Calendar time to complete the trial under delayed outcomes (TITE, blue solid) and under the alternative of waiting for complete outcomes (orange dashed).",
  sprintf(
    "The trial is one simulated replicate with random seed %d; the operating characteristics of the design for this scenario are summarised over 1000 replicates in the accompanying table.",
    seed
  ),
  sprintf(
    "The true OBDC of this scenario is %s; the combination recommended at the end of this trial is %s.",
    obdc.text, selected.text
  ),
  sprintf(
    "Simulation settings: accrual %g days per patient, toxicity assessment window %g days, efficacy assessment window %g days, suspension rule active with proportion %.1f.",
    delay.control$accrual,
    delay.control$tox_window,
    delay.control$eff_window,
    delay.control$suspension
  ),
  sprintf(
    "Reproducibility: %s with ggplot2 %s. Code: Comb-STEIN-Code/plot/plot_trial_trajectory.R; command: Rscript Comb-STEIN-Code/plot/plot_trial_trajectory.R --scenario=%s --method=%s --seed=%d --accrual=%g --tox-window=%g --eff-window=%g --suspension=%g.",
    R.version.string,
    as.character(utils::packageVersion("ggplot2")),
    scenario, method, seed,
    delay.control$accrual, delay.control$tox_window,
    delay.control$eff_window, delay.control$suspension
  )
)
write.utf8(caption, file.path(table_dir, paste0(file.stem, "_caption.txt")))

## 控制台在 C 语言环境下会把 μg 显示成转义串，因此打印时换成 ASCII 写法。
cat(gsub(micro.sign, "u", header, fixed = TRUE), sep = "\n")
cat("\n队列过程表已保存到：\n")
cat(" ", file.path(table_dir, paste0(file.stem, ".csv")), "\n")
cat(" ", file.path(table_dir, paste0(file.stem, ".md")), "\n")
cat("图注已保存到：\n")
cat(" ", file.path(table_dir, paste0(file.stem, "_caption.txt")), "\n")
cat("图已保存到", file.path(figure_dir, paste0(figure.file, ".pdf")), "\n")
