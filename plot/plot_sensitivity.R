## 敏感性分析图：暂停比例、入组速度、队列大小。
##
## 三项敏感性分析共用同一批情景、同一随机种子，基准设置均为
## AW=(1,2)、入组 1.5 例/月、暂停比例 0.5、队列大小 3。
## 每项分析输出一张四面板图（选择率 / 患者分配 / 样本量 / 持续时间）
## 和一份对应的数据表。
##
## 柱形与图例的排列与正文图3/图4一致：同一方法的各个设置**相邻**（同一色系的
## 深浅色成对出现），方法之间再依次排列；图例按列填充，每列是一种方法的深浅色。
##
## 运行配置来自 output/sensitivity/runs.csv（由 Comb-STEIN-Code/tmp/run_sensitivity.R 生成）。
##
## 用法：Rscript Comb-STEIN-Code/plot/plot_sensitivity.R [项目根目录]
##       [--manifest=output/sensitivity/runs.csv] [--figure-dir=output/figures]

suppressPackageStartupMessages({
  library(ggplot2)
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

manifest.opt <- .pull.option(args, "manifest")
args <- manifest.opt$args
figure.opt <- .pull.option(args, "figure-dir")
args <- figure.opt$args

project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

figure_dir <- if (is.null(figure.opt$value)) {
  file.path(project_root, "output", "figures")
} else if (grepl("^(/|[A-Za-z]:)", figure.opt$value)) {
  figure.opt$value
} else {
  file.path(project_root, figure.opt$value)
}
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- if (is.null(manifest.opt$value)) {
  file.path(project_root, "output", "sensitivity", "runs.csv")
} else if (grepl("^(/|[A-Za-z]:)", manifest.opt$value)) {
  manifest.opt$value
} else {
  file.path(project_root, manifest.opt$value)
}
if (!file.exists(manifest_path)) {
  stop("找不到运行清单：", manifest_path, call. = FALSE)
}
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)

method_levels <- c("Comb-STEIN-L", "Comb-STEIN-H", "Runin-STEIN")
method_colors <- c(
  "Comb-STEIN-L" = "#1973B0",
  "Comb-STEIN-H" = "#769658",
  "Runin-STEIN"  = "#E3A472"
)
panel_labels <- c(
  Selection = "A: OBDC selection rate (%)",
  Patient   = "B: Patients allocated to OBDC (%)",
  Sample    = "C: Average sample size (patients)",
  Duration  = "D: Trial duration (days)"
)
metric_levels <- names(panel_labels)

## 把基准色按给定比例与白色混合，得到同一色系的浅色变体。
.mix.white <- function(color, prop) {
  rgb <- grDevices::col2rgb(color)[, 1]
  mixed <- round(rgb + (255 - rgb) * prop)
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

## 第 k 个设置（共 n 个）的浅化比例：第一个设置最深，最后一个最浅。
.shade_of <- function(k, n) {
  if (n <= 1) 0 else 0.55 * (k - 1) / (n - 1)
}

.load.run <- function(id, label) {
  hit <- manifest$id == id
  if (!any(hit)) stop("运行清单缺少配置：", id, call. = FALSE)
  row <- manifest[which(hit)[1], ]
  path <- file.path(project_root, row$run_dir, row$csv)
  if (!file.exists(path)) stop("找不到运行结果：", path, call. = FALSE)
  frame <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  frame$Run <- row$id
  frame$Setting <- label
  frame
}

## 三项敏感性分析：entries 的顺序决定颜色深浅（第一项最深）。
analyses <- list(
  list(
    name = "suspension",
    out = "fig_sensitivity_suspension",
    x_title = "Suspension ratio r",
    entries = list(
      list(id = "main_t12", label = "r = 0.5"),
      list(id = "sus_r1", label = "r = 1.0")
    )
  ),
  list(
    name = "accrual",
    out = "fig_sensitivity_accrual",
    entries = list(
      list(id = "acc_a30", label = "1 / month"),
      list(id = "main_t12", label = "1.5 / month"),
      list(id = "acc_a15", label = "2 / month")
    )
  ),
  list(
    name = "cohortsize",
    out = "fig_sensitivity_cohortsize",
    x_title = "Cohort size",
    entries = list(
      list(id = "cs_1", label = "size 1"),
      list(id = "cs_2", label = "size 2"),
      list(id = "main_t12", label = "size 3")
    )
  )
)

for (analysis in analyses) {
  ## 组装长表：每个方法一组颜色，方法内按设置展开。
  ## 系列顺序与正文图3/图4保持一致：**同一方法的各个设置相邻**（深浅成对），
  ## 方法之间再依次排列，即 [Comb-STEIN-L 各设置] [Comb-STEIN-H 各设置] [Runin-STEIN 各设置]。
  ## 注意：展开顺序必须是"方法在外、设置在内"，写成 lapply(entries, ...) 会把
  ## 同一设置的三根柱形排在一起、深浅色系被拆散，与正文图不一致。
  ## 这里不能用 outer() 生成：它会把设置下标扩成长向量，而 [[ 收到
  ## 向量下标时会做递归索引，导致报错。
  entry_labels <- vapply(analysis$entries, function(entry) entry$label, character(1))
  series_levels <- unlist(lapply(
    method_levels,
    function(method) paste(method, entry_labels, sep = "\n")
  ))
  long <- list()
  csv_rows <- list()
  for (entry_index in seq_along(analysis$entries)) {
    entry <- analysis$entries[[entry_index]]
    frame <- .load.run(entry$id, entry$label)
    frame$Series <- factor(
      paste(frame$Method, entry$label, sep = "\n"),
      levels = series_levels
    )
    csv_rows[[length(csv_rows) + 1]] <- frame
    long[[length(long) + 1]] <- do.call(rbind, lapply(metric_levels, function(metric) {
      data.frame(
        Scenario = frame$Scenario,
        Series = frame$Series,
        Metric = factor(metric, levels = metric_levels),
        Value = frame[[metric]]
      )
    }))
  }
  long <- do.call(rbind, long)
  csv <- do.call(rbind, csv_rows)

  scenario_number <- as.integer(sub("^S", "", csv$Scenario))
  scenario_levels <- unique(csv$Scenario[order(scenario_number)])
  long$Scenario <- factor(long$Scenario, levels = scenario_levels)

  ## 配色必须与 series_levels 的顺序**同一口径**（方法在外、设置在内），
  ## 否则 setNames 会把深浅色错配到别的系列上。
  series_colors <- setNames(
    unlist(lapply(
      method_levels,
      function(method) vapply(
        seq_along(analysis$entries),
        function(entry_index) {
          .mix.white(method_colors[[method]], .shade_of(entry_index, length(analysis$entries)))
        },
        character(1)
      )
    )),
    series_levels
  )

  plot <- ggplot(long, aes(x = Scenario, y = Value, fill = Series)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.78, colour = NA) +
    facet_wrap(~Metric, ncol = 2, scales = "free_y",
               labeller = labeller(Metric = panel_labels)) +
    scale_fill_manual(
      values = series_colors, breaks = series_levels, name = "Method",
      ## byrow = FALSE 让图例按**列**填充：每一列是一种方法的深浅色，
      ## 与正文图3/图4的图例排法一致（每列上下两根是同色系的深浅两档）。
      guide = guide_legend(nrow = length(analysis$entries), byrow = FALSE)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08)),
      labels = function(x) format(x, trim = TRUE, scientific = FALSE)
    ) +
    labs(x = "Scenario", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.border = element_rect(colour = "grey40", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = 10.5),
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.5),
      axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "lines"),
      legend.margin = margin(2, 2, 2, 2),
      panel.spacing = unit(0.9, "lines"),
      plot.margin = margin(6, 10, 6, 6)
    )

  ggsave(file.path(figure_dir, paste0(analysis$out, ".png")), plot,
         width = 10.5, height = 9, dpi = 300)
  ggsave(file.path(figure_dir, paste0(analysis$out, ".pdf")), plot,
         width = 10.5, height = 9, device = grDevices::cairo_pdf)
  table_columns <- c(
    "Scenario", "Setting", "Run", "Method", "Delay",
    "Selection", "ET", "Patient", "Sample", "Duration"
  )
  utils::write.csv(
    csv[order(csv$Scenario, csv$Setting, csv$Method),
        table_columns[table_columns %in% names(csv)]],
    file.path(figure_dir, paste0(sub("^fig_", "", analysis$out), "_metrics.csv")),
    row.names = FALSE
  )
  cat("saved:", analysis$out, "(", length(analysis$entries), "settings )\n")
}
