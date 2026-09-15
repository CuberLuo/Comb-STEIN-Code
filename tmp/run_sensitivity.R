## 主分析与敏感性分析的批量模拟驱动。
##
## 依次运行 code/main.R 的全部情景，并把每次运行的输出目录登记到
## output/sensitivity/runs.csv，供画图与制表脚本按“设置”取用。
## 所有配置都使用同一个随机种子和同一批情景顺序，因此不同设置之间
## 逐试验一一配对，可以直接比较。
##
## 用法：Rscript code/tmp/run_sensitivity.R [项目根目录]
## 重复运行时会跳过已经登记且输出仍存在的配置，便于中断后续跑。

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
main_script <- file.path(project_root, "code", "main.R")
if (!file.exists(main_script)) {
  stop("找不到 code/main.R：", main_script, call. = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")

output_dir <- file.path(project_root, "output")
sensitivity_dir <- file.path(output_dir, "sensitivity")
log_dir <- file.path(sensitivity_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(sensitivity_dir, "runs.csv")

NTRIAL <- 1000L
WORKERS <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

## 主分析：AW=(1,2) 与 AW=(2,3)，入组 1.5 例/月（20 天/例）。
## 敏感性分析：暂停比例、入组速度、队列大小，基准均为 AW=(1,2) + 1.5 例/月。
configs <- list(
  list(id = "main_c12", group = "main", label = "Complete AW=(1,2)",
       delay = FALSE, accrual = 20, tox = 30, eff = 60, susp = 0.5, csize = 3),
  list(id = "main_t12", group = "main", label = "TITE AW=(1,2)",
       delay = TRUE, accrual = 20, tox = 30, eff = 60, susp = 0.5, csize = 3),
  list(id = "main_c23", group = "main", label = "Complete AW=(2,3)",
       delay = FALSE, accrual = 20, tox = 60, eff = 90, susp = 0.5, csize = 3),
  list(id = "main_t23", group = "main", label = "TITE AW=(2,3)",
       delay = TRUE, accrual = 20, tox = 60, eff = 90, susp = 0.5, csize = 3),
  list(id = "sus_r1", group = "suspension", label = "r=1.0",
       delay = TRUE, accrual = 20, tox = 30, eff = 60, susp = 1, csize = 3),
  list(id = "acc_a30", group = "accrual", label = "1 patient/month",
       delay = TRUE, accrual = 30, tox = 30, eff = 60, susp = 0.5, csize = 3),
  list(id = "acc_a15", group = "accrual", label = "2 patients/month",
       delay = TRUE, accrual = 15, tox = 30, eff = 60, susp = 0.5, csize = 3),
  list(id = "cs_1", group = "cohortsize", label = "cohort size 1",
       delay = TRUE, accrual = 20, tox = 30, eff = 60, susp = 0.5, csize = 1),
  list(id = "cs_2", group = "cohortsize", label = "cohort size 2",
       delay = TRUE, accrual = 20, tox = 30, eff = 60, susp = 0.5, csize = 2)
)

manifest <- if (file.exists(manifest_path)) {
  utils::read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  NULL
}

already.done <- function(cfg) {
  if (is.null(manifest) || !any(manifest$id == cfg$id)) return(FALSE)
  row <- manifest[manifest$id == cfg$id, ][1, ]
  csv <- if (isTRUE(row$delay)) {
    "comb_stein_results_utility_tite_all.csv"
  } else {
    "comb_stein_results_utility_all.csv"
  }
  file.exists(file.path(project_root, row$run_dir, csv))
}

for (cfg in configs) {
  if (already.done(cfg)) {
    cat(sprintf("[skip] %-10s 已完成\n", cfg$id))
    next
  }
  cli <- c(
    as.character(NTRIAL),
    "--quiet",
    sprintf("--workers=%d", WORKERS),
    sprintf("--accrual=%g", cfg$accrual),
    sprintf("--tox-window=%g", cfg$tox),
    sprintf("--eff-window=%g", cfg$eff),
    sprintf("--suspension=%g", cfg$susp),
    sprintf("--cohort-size=%d", cfg$csize)
  )
  if (isTRUE(cfg$delay)) cli <- c(cli, "--tite")

  log_path <- file.path(log_dir, paste0(cfg$id, ".log"))
  before <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  started <- Sys.time()
  cat(sprintf(
    "[run ] %-10s %s accrual=%g AW=(%g,%g) susp=%g csize=%d ...\n",
    cfg$id, cfg$label, cfg$accrual, cfg$tox, cfg$eff, cfg$susp, cfg$csize
  ))
  status <- system2(
    rscript,
    args = c(shQuote(main_script), shQuote(cli)),
    stdout = log_path, stderr = log_path
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  if (!identical(as.integer(status), 0L)) {
    cat(sprintf("[fail] %-10s 退出码 %s，见 %s\n", cfg$id, status, log_path))
    stop("模拟运行失败：", cfg$id, call. = FALSE)
  }
  new_dirs <- setdiff(
    list.dirs(output_dir, recursive = FALSE, full.names = TRUE),
    before
  )
  new_dirs <- new_dirs[grepl("run_", basename(new_dirs))]
  if (length(new_dirs) == 0) {
    stop("未找到 ", cfg$id, " 的输出目录", call. = FALSE)
  }
  info <- file.info(new_dirs)
  run_dir <- basename(new_dirs[which.max(info$mtime)])
  csv <- if (isTRUE(cfg$delay)) {
    "comb_stein_results_utility_tite_all.csv"
  } else {
    "comb_stein_results_utility_all.csv"
  }
  if (!file.exists(file.path(output_dir, run_dir, csv))) {
    stop("输出目录缺少 ", csv, "：", run_dir, call. = FALSE)
  }

  row <- data.frame(
    id = cfg$id, group = cfg$group, label = cfg$label,
    delay = isTRUE(cfg$delay), accrual = cfg$accrual,
    tox_window = cfg$tox, eff_window = cfg$eff, suspension = cfg$susp,
    cohortsize = cfg$csize, ntrials = NTRIAL, workers = WORKERS,
    run_dir = file.path("output", run_dir), csv = csv,
    minutes = round(elapsed, 2), stringsAsFactors = FALSE
  )
  manifest <- if (is.null(manifest)) row else {
    rbind(manifest[manifest$id != cfg$id, , drop = FALSE], row)
  }
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
  cat(sprintf("[done] %-10s %s  用时 %.1f 分钟\n", cfg$id, run_dir, elapsed))
}

cat("\n全部配置完成，清单：", manifest_path, "\n")
print(manifest[, c("id", "group", "label", "run_dir", "minutes")], row.names = FALSE)
