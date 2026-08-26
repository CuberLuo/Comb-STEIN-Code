## Comb-STEIN 命令行入口和公共函数加载器。
##
## 用法：
##   Rscript main.R
##   Rscript main.R S1 S2 S3 50
##   Rscript main.R S1 50 --delayed=true
##   Rscript main.R S1 50 --tite --accrual=10 --tox-window=30
##     --eff-window=90 --suspension=0.5
##
## 每次命令行运行都会在 output 下创建新的时间戳文件夹，

## 定位包含入口脚本的目录。优先使用 source 框架中的路径，
## 再依次回退到 Rscript 的 --file 参数和当前工作目录。
## 返回规范化的绝对路径。
.comb_stein_script_dir <- function() {
  ## RStudio 的 debugSource() 可能在调用栈中加入一个工作目录路径，
  ## 因此不能盲目使用最后一个 sys.frames()$ofile。
  frame.files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1]
    },
    character(1)
  )
  frame.files <- frame.files[nzchar(frame.files)]
  if (length(frame.files) > 0) {
    frame.dirs <- unique(vapply(
      frame.files,
      function(path) dirname(normalizePath(path, mustWork = TRUE)),
      character(1)
    ))
    ## 选择真正包含本脚本模块的目录，例如 code/，而不是项目根目录。
    valid.dir <- frame.dirs[file.exists(file.path(frame.dirs, "config.R"))]
    if (length(valid.dir) > 0) return(valid.dir[length(valid.dir)])
  }

  file.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file.arg) > 0) {
    file.dir <- dirname(normalizePath(
      sub("^--file=", "", file.arg[1]),
      mustWork = TRUE
    ))
    if (file.exists(file.path(file.dir, "config.R"))) return(file.dir)
  }
  working.dir <- normalizePath(getwd(), mustWork = TRUE)
  if (file.exists(file.path(working.dir, "config.R"))) return(working.dir)
  if (file.exists(file.path(working.dir, "code", "config.R"))) {
    return(normalizePath(file.path(working.dir, "code"), mustWork = TRUE))
  }
  stop("无法定位 main.R 所在目录：未找到 config.R")
}

.comb_stein_root <- .comb_stein_script_dir()
.comb_stein_project_root <- normalizePath(
  file.path(.comb_stein_root, ".."),
  mustWork = TRUE
)
.comb_stein_modules <- c(
  "config.R",
  "time_to_event.R",
  "isotonic_estimators.R",
  "dose_rules.R",
  "simulate_stein.R",
  "simulate_waterfall.R",
  "scenarios.R",
  "run_scenarios.R"
)

for (.comb_stein_module in .comb_stein_modules) {
  sys.source(
    file.path(.comb_stein_root, .comb_stein_module),
    envir = environment()
  )
}

## 将命令行中的布尔值写法转换为单个逻辑值。
## `option` 仅用于生成明确的校验错误信息。
.parse.cli.boolean <- function(value, option) {
  value <- tolower(as.character(value)[1])
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  FALSE
}

## 从 `args` 中移除一个 --option=value 参数。
## 同时返回解析后的值（或 `default`）以及未消费的参数。
.pull.cli.option <- function(args, option, default = NULL) {
  pattern <- paste0("^--", option, "=")
  hits <- grep(pattern, args)
  if (length(hits) == 0L) return(list(value = default, args = args))
  list(
    value = sub(pattern, "", args[hits]),
    args = args[-hits]
  )
}

## 使用开始时间在 `output.root` 下创建唯一的运行目录。
## 如果同一时间戳已被占用，则追加数字后缀。
.create.timestamped.results.dir <- function(
  output.root,
  started.at = Sys.time()
) {
  dir.create(output.root, showWarnings = FALSE, recursive = TRUE)
  timestamp <- format(started.at, "%Y%m%d_%H%M%S")
  base.name <- paste0("run_", timestamp)
  attempt <- 0L

  repeat {
    folder.name <- if (attempt == 0L) {
      base.name
    } else {
      paste0(base.name, "_", sprintf("%02d", attempt))
    }
    candidate <- file.path(output.root, folder.name)
    if (dir.create(candidate, showWarnings = FALSE, recursive = FALSE)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    attempt <- attempt + 1L
  }
}

## 解析命令行参数，运行指定的模拟场景并保存 CSV 文件。
## 不可见地返回合并结果表；程序化调用者可用 `results.dir` 覆盖
## 默认的时间戳输出目录。
run.comb.stein.cli <- function(
  args = commandArgs(trailingOnly = TRUE),
  results.dir = NULL
) {
  ## 记录本次运行的开始时间。这个时间既用于屏幕输出，
  ## 也用于生成本次运行专属的结果文件夹名称。
  run.started.at <- Sys.time()
  if (is.null(results.dir)) {
    ## 如果没有指定输出目录，就在项目根目录的 output/ 下
    ## 自动创建一个带时间戳的文件夹，避免覆盖以前的结果。
    results.dir <- .create.timestamped.results.dir(
      file.path(.comb_stein_project_root, "output"),
      run.started.at
    )
  } else {
    ## 程序化调用时也可以手动指定输出目录。
    ## 如果目录不存在，则递归创建；normalizePath() 用于得到绝对路径。
    dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
    results.dir <- normalizePath(results.dir, mustWork = TRUE)
  }

  ## 固定随机种子，使同样的参数能够得到可重复的模拟结果。
  set.seed(DEFAULT_SEED)
  ## 输出本次运行的基本信息和由设计参数推导出的两个阈值：
  ## phi 是毒性阈值，psi 是疗效阈值。
  cat("Run started:", format(run.started.at, "%Y-%m-%d %H:%M:%S %Z"), "\n")
  cat("Output folder:", results.dir, "\n")
  cat(sprintf("Derived toxicity cutoff  phi = %.4f\n", phi))
  cat(sprintf("Derived efficacy  cutoff psi = %.4f\n", psi))

  ## 默认每个场景重复的试验次数。
  ## 如果命令行中提供了数字，后面会用该数字覆盖默认值。
  ntrial <- DEFAULT_NTRIAL

  ## delayed = FALSE 表示使用完整结局；
  ## delayed = TRUE 表示启用 TITE/延迟结局模拟。
  ## 默认采用完整结局，只有明确指定参数时才启用延迟结局。
  delayed <- FALSE
  delayed.flag <- .pull.cli.option(args, "delayed", NULL)
  args <- delayed.flag$args
  if (!is.null(delayed.flag$value)) {
    ## 处理 --delayed=true 或 --delayed=false 形式的参数。
    delayed <- .parse.cli.boolean(delayed.flag$value, "--delayed")
  }
  ## 也支持更简洁的开关形式 --tite 和 --no-tite。
  has.tite <- "--tite" %in% args
  has.no.tite <- "--no-tite" %in% args
  ## 读取开关的值后，从参数列表中删除它们，避免被误判为场景名称。
  if (has.tite) delayed <- TRUE
  if (has.no.tite) delayed <- FALSE
  args <- args[!args %in% c("--tite", "--no-tite")]

  ## delay.control 保存延迟结局相关的设置。
  ## 先建立“命令行参数名 -> R 内部名称”的对应关系。
  delay.control <- list()
  numeric.options <- c(
    accrual = "accrual",
    `tox-window` = "tox_window",
    `eff-window` = "eff_window",
    suspension = "suspension"
  )
  for (option in names(numeric.options)) {
    ## 读取数值型延迟参数，例如：
    ## --accrual=10、--tox-window=30、--eff-window=90。
    pulled <- .pull.cli.option(args, option, NULL)
    args <- pulled$args
    if (!is.null(pulled$value)) {
      delay.control[[numeric.options[[option]]]] <-
        as.numeric(pulled$value)
    }
  }
  logical.options <- c(
    `use-suspension` = "use_suspension",
    `accrual-random` = "accrual_random"
  )
  for (option in names(logical.options)) {
    ## 读取逻辑型延迟参数，例如：
    ## --use-suspension=true 或 --accrual-random=false。
    pulled <- .pull.cli.option(args, option, NULL)
    args <- pulled$args
    if (!is.null(pulled$value)) {
      delay.control[[logical.options[[option]]]] <-
        .parse.cli.boolean(pulled$value, paste0("--", option))
    }
  }
  ## Complete 和 TITE 都使用同一套入组速度与评估窗计算试验持续时间。
  delay.control <- normalize.delay.control(delay.control)

  ## 读取要运行的方法。多个方法使用逗号分隔。
  methods.flag <- .pull.cli.option(args, "methods", NULL)
  args <- methods.flag$args
  methods <- c(
    "Comb-STEIN-L",
    "Comb-STEIN-H",
    "Runin-STEIN",
    "Waterfall-L",
    "Waterfall-H"
  )
  if (!is.null(methods.flag$value)) {
    methods <- trimws(
      strsplit(methods.flag$value, ",", fixed = TRUE)[[1]]
    )
  }

  ## 默认运行 scenarios.R 中定义的全部场景。
  ## 用户也可以在命令行中指定一个或多个场景名称。
  scenarios.to.run <- names(scenarios)
  if (length(args) > 0) {
    ## 数字参数解释为每个场景的重复次数，例如：main.R S1 S2 50。
    ## 非数字参数解释为场景名称，例如 S1、S2。
    is.numeric.arg <- suppressWarnings(!is.na(as.numeric(args)))
    if (any(is.numeric.arg)) {
      ntrial <- as.numeric(args[is.numeric.arg])[1]
    }
    scenario.args <- args[!is.numeric.arg]
    if (length(scenario.args) > 0) scenarios.to.run <- scenario.args
  }

  ## 在正式开始模拟前，把实际采用的延迟设置打印出来，方便核对。
  if (delayed) {
    cat(sprintf(
      paste0(
        "Delayed outcomes: ON (accrual %.1f days, toxicity window %.1f ",
        "days, efficacy window %.1f days, suspension %.0f%%, rule %s)\n"
      ),
      delay.control$accrual,
      delay.control$tox_window,
      delay.control$eff_window,
      100 * delay.control$suspension,
      if (delay.control$use_suspension) "ON" else "OFF"
    ))
  } else {
    cat(sprintf(
      paste0(
        "Delayed outcomes: OFF (complete outcomes used at each decision; ",
        "accrual %.1f days, toxicity window %.1f days, efficacy window %.1f days)\n"
      ),
      delay.control$accrual,
      delay.control$tox_window,
      delay.control$eff_window
    ))
  }

  ## 用于保存所有场景的汇总结果。
  all.results <- NULL

  for (name in scenarios.to.run) {
    ## 每个场景都会调用 run.scenario()，运行 methods 中指定的方法。
    ## OBDC 始终定义为 utility 最大的组合
    cat("Running scenario", name, "with ntrial =", ntrial, "...\n")
    # *****执行scenario*****
    result <- run.scenario(
      name,
      ntrial = ntrial,
      delayed = delayed,
      delay.control = delay.control,
      methods = methods
    )
    ## run.scenario() 将当前场景的真实 OBDC 保存为结果属性。
    ## true.ij 是一个两列矩阵，每一行表示一个并列的最优组合 (行, 列)。
    true.ij <- attr(result, "true.obdcs")
    labels <- apply(
      true.ij,
      1,
      function(ij) sprintf("(%d,%d)", ij[1], ij[2])
    )
    cat(
      "  true OBDC(s) [utility maximum]: ",
      paste(labels, collapse = ", "),
      "\n",
      sep = ""
    )

    ## 每个场景单独保存一个 CSV 文件，便于后续单独查看。
    scenario.outfile <- file.path(
      results.dir,
      paste0(
        "comb_stein_results_",
        "utility",
        if (delayed) "_tite" else "",
        "_",
        name,
        ".csv"
      )
    )
    write.csv(result, scenario.outfile, row.names = FALSE)
    cat("  scenario results saved to", scenario.outfile, "\n")
    ## 同时把当前场景追加到总结果表中。
    all.results <- rbind(all.results, result)
  }

  ## 打印所有场景的汇总表。
  cat("\n\n===================== Table 2 simulation =====================\n")
  print(all.results, row.names = FALSE)

  ## 将所有场景合并后的结果保存为一个总 CSV 文件。
  outfile <- file.path(
    results.dir,
    paste0(
      "comb_stein_results_",
      "utility",
      if (delayed) "_tite" else "",
      "_all.csv"
    )
  )
  write.csv(all.results, outfile, row.names = FALSE)
  cat("\nResults saved to", outfile, "\n")
  ## 返回结果但不自动打印一遍，方便其他 R 函数接收并继续处理。
  invisible(all.results)
}

run.comb.stein.cli()
