## 表3：不同评估窗（AW）下的敏感性分析结果。
##
## 读取 output/sensitivity/runs.csv 登记的四次主分析运行：
##   main_c12  完整结局，AW=(1,2)，入组 1.5 例/月
##   main_c23  完整结局，AW=(2,3)
##   main_t12  TITE，AW=(1,2)
##   main_t23  TITE，AW=(2,3)
## 四次运行使用同一种子与同一批情景，因此逐试验一一配对。
##
## 输出两个文件：
##   output/tables/table3_windows_full.csv  逐方法明细（12 情景 x 3 方法 x 4 设置）
##   output/tables/table3_windows.md        论文排版：每个情景 3 行
##   output/tables/table3_windows.csv       与 md 同内容的 CSV
##
## 说明：完整结局（Comb-STEIN）的决策不依赖评估窗，因此其选择率、患者、
## 样本量与 DLT 数在 AW=(1,2) 与 AW=(2,3) 下相同，仅持续时间不同，
## 表中以 "a/b" 形式给出。表中的 3 行分别对三种 STEIN 类方法取平均，
## 逐方法结果见明细 CSV。
##
## 用法：Rscript code/tmp/make_table3.R [项目根目录]

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

manifest_path <- file.path(project_root, "output", "sensitivity", "runs.csv")
if (!file.exists(manifest_path)) {
  stop("找不到运行清单：", manifest_path, call. = FALSE)
}
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)

table_dir <- file.path(project_root, "output", "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

DAYS_PER_MONTH <- 30.4375

## 一次运行的路径信息。
run.info <- function(id) {
  hit <- manifest$id == id
  if (!any(hit)) stop("运行清单缺少配置：", id, call. = FALSE)
  row <- manifest[which(hit)[1], ]
  list(
    id = id,
    dir = file.path(project_root, row$run_dir),
    csv = row$csv
  )
}

## 读取一次运行的合并结果表（12 情景 x 3 方法）。
read.run <- function(info) {
  path <- file.path(info$dir, info$csv)
  if (!file.exists(path)) stop("找不到运行结果：", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

## 从场景报告 TXT 中读取每种方法的平均 DLT 例数。
## 报告里 "Number of toxicity observed at each dose (average):" 之后
## 是一张剂量网格，逐元素求和即为该方法的平均 DLT 例数。
read.dlt <- function(info) {
  result <- list()
  for (scenario in paste0("S", seq_len(12))) {
    txt <- file.path(
      info$dir,
      sub("_all\\.csv$", paste0("_", scenario, ".txt"), info$csv)
    )
    if (!file.exists(txt)) next
    lines <- readLines(txt, warn = FALSE)
    method <- NA_character_
    collecting <- FALSE
    total <- 0
    for (line in lines) {
      if (grepl("^=+\\s.*\\s=+$", line)) {
        method <- trimws(gsub("=", "", line))
        next
      }
      if (grepl("^Number of toxicity observed at each dose", line)) {
        collecting <- TRUE
        total <- 0
        next
      }
      if (collecting) {
        numbers <- suppressWarnings(
          as.numeric(strsplit(trimws(line), "\\s+")[[1]])
        )
        if (length(numbers) == 0 || all(is.na(numbers))) {
          collecting <- FALSE
          if (!is.na(method)) {
            result[[paste(scenario, method, sep = "/")]] <- total
          }
        } else {
          total <- total + sum(numbers, na.rm = TRUE)
        }
      }
    }
  }
  frame <- data.frame(
    key = names(result),
    DLT = round(unlist(result, use.names = FALSE), 2),
    stringsAsFactors = FALSE
  )
  frame$Scenario <- sub("/.*", "", frame$key)
  frame$Method <- sub(".*/", "", frame$key)
  frame[, c("Scenario", "Method", "DLT")]
}

RUNS <- list(
  "Complete AW=(1,2)" = run.info("main_c12"),
  "Complete AW=(2,3)" = run.info("main_c23"),
  "TITE (1, 2)"       = run.info("main_t12"),
  "TITE (2, 3)"       = run.info("main_t23")
)

rows <- list()
for (setting in names(RUNS)) {
  info <- RUNS[[setting]]
  frame <- read.run(info)
  frame$Setting <- setting
  frame$Run <- info$id
  frame <- merge(frame, read.dlt(info), by = c("Scenario", "Method"),
                 all.x = TRUE, sort = FALSE)
  rows[[setting]] <- frame
}
full <- do.call(rbind, rows)
full$DurationMonths <- round(full$Duration / DAYS_PER_MONTH, 1)
full <- full[order(
  as.integer(sub("^S", "", full$Scenario)),
  full$Run,
  full$Method
), ]
full <- full[, c("Scenario", "Method", "Run", "Setting", "Selection", "ET",
                 "Patient", "Sample", "Duration", "DurationMonths", "DLT")]

utils::write.csv(
  full, file.path(table_dir, "table3_windows_full.csv"), row.names = FALSE
)

## 论文排版的紧凑表：每个情景 3 行，分别为完整结局与两种评估窗下的 TITE，
## 每行对三种 STEIN 类方法（Comb-STEIN-L、Comb-STEIN-H、Runin-STEIN）取平均。
scenarios <- sort(unique(full$Scenario), method = "radix")
scenarios <- paste0("S", sort(as.integer(sub("^S", "", scenarios))))
compact <- list()
for (scenario in scenarios) {
  block <- full[full$Scenario == scenario, ]
  pick <- function(setting, columns) {
    sub <- block[block$Setting == setting, ]
    vapply(columns, function(column) round(mean(sub[[column]]), 1), numeric(1))
  }
  durations <- c(
    pick("Complete AW=(1,2)", "DurationMonths")[["DurationMonths"]],
    pick("Complete AW=(2,3)", "DurationMonths")[["DurationMonths"]]
  )
  complete <- pick("Complete AW=(1,2)",
                   c("Selection", "Patient", "Sample", "DLT"))
  t12 <- pick("TITE (1, 2)",
              c("Selection", "Patient", "Sample", "DurationMonths", "DLT"))
  t23 <- pick("TITE (2, 3)",
              c("Selection", "Patient", "Sample", "DurationMonths", "DLT"))
  compact[[length(compact) + 1]] <- data.frame(
    Scenario = scenario,
    Method = c("Comb-STEIN", "TITE (1, 2)", "TITE (2, 3)"),
    Selection = c(complete[["Selection"]], t12[["Selection"]], t23[["Selection"]]),
    Patient = c(complete[["Patient"]], t12[["Patient"]], t23[["Patient"]]),
    Sample = c(complete[["Sample"]], t12[["Sample"]], t23[["Sample"]]),
    Duration = c(
      sprintf("%.1f/%.1f", durations[1], durations[2]),
      sprintf("%.1f", t12[["DurationMonths"]]),
      sprintf("%.1f", t23[["DurationMonths"]])
    ),
    DLT = c(complete[["DLT"]], t12[["DLT"]], t23[["DLT"]]),
    stringsAsFactors = FALSE
  )
}
compact <- do.call(rbind, compact)
names(compact) <- c("Scenario", "Method", "SelectionRate(%)", "Patient(%)",
                    "SampleSize", "Duration(months)", "DLTcount")

utils::write.csv(
  compact, file.path(table_dir, "table3_windows.csv"), row.names = FALSE
)

md <- c(
  "# 表3 不同评估窗下的敏感性分析结果",
  "",
  "| 情景 | 方法 | 选择率% | 患者% | 样本量 | 持续时间 | DLT数 |",
  "| --- | --- | --- | --- | --- | --- | --- |"
)
for (index in seq_len(nrow(compact))) {
  row <- compact[index, ]
  prefix <- if (index %% 3 == 1) row$Scenario else ""
  md <- c(md, sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    prefix, row$Method, row$`SelectionRate(%)`, row$`Patient(%)`,
    row$`SampleSize`, row$`Duration(months)`, row$`DLTcount`
  ))
}
md <- c(
  md, "",
  "注：持续时间单位为月，由模拟输出的天数除以 30.4375 换算。",
  "Comb-STEIN 为完整结局设计，其决策不依赖评估窗，故其选择率、患者、样本量",
  "与 DLT 数在 AW=(1, 2) 与 AW=(2, 3) 下相同，表中持续时间以“a/b”形式给出，",
  "a 对应 AW=(1, 2)，b 对应 AW=(2, 3)。括号内数值依次为毒性与疗效评估窗，",
  "单位为月。表中三行分别对 Comb-STEIN-L、Comb-STEIN-H 与 Runin-STEIN 取平均，",
  "逐方法结果见 table3_windows_full.csv。"
)
## 源码为 UTF-8，而 Windows 上 R 的本地编码可能是 C/GBK，
## enc2utf8() 会把已经是 UTF-8 的字面量按本地编码二次转换而损坏，
## 因此这里直接按原始字节写盘。
writeLines(md, file.path(table_dir, "table3_windows.md"), useBytes = TRUE)

cat("表3 已生成：\n")
cat("  ", file.path(table_dir, "table3_windows.md"), "\n")
cat("  ", file.path(table_dir, "table3_windows.csv"), "\n")
cat("  ", file.path(table_dir, "table3_windows_full.csv"), "\n")
print(compact, row.names = FALSE)
