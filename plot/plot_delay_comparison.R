## 延迟（TITE）与完整结局（non-TITE）操作特征对比图。
##
## 读取 output/ 下最新一次 TITE 与完整结局运行的合并结果，绘制
## Selection、Patient、Sample、Duration 四个指标的对比图。
##
## 图形风格与论文图3/图4保持一致：
##   * 横轴为情景（S1-S12），每个情景内三种方法各两根柱形；
##   * 深色柱形为完整结局（non-TITE），浅色柱形为 TITE 扩展；
##   * 图例为两行三列，每列对应一种方法，上行为非 TITE、下行为 TITE。
##
## 用法：Rscript code/plot/plot_delay_comparison.R [项目根目录]

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!dir.exists(file.path(project_root, "output"))) {
  stop("No output/ directory under ", project_root, call. = FALSE)
}

output_root <- file.path(project_root, "output")
figure_dir <- file.path(output_root, "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

method_levels <- c("Comb-STEIN-L", "Comb-STEIN-H", "Runin-STEIN")
delay_levels <- c("Complete", "TITE")

## 与论文图3/图4完全相同的配色：深色为完整结局，浅色为 TITE。
method_colors <- c(
  "Comb-STEIN-L" = "#1973B0",
  "Comb-STEIN-H" = "#769658",
  "Runin-STEIN"  = "#E3A472"
)
tite_colors <- c(
  "Comb-STEIN-L" = "#A2CCE1",
  "Comb-STEIN-H" = "#A7C698",
  "Runin-STEIN"  = "#DEC2AF"
)

## 图例标签沿用论文中的命名习惯（TITE 版本加 TITE- 前缀）。
method_labels <- c(
  "Comb-STEIN-L" = "Comb-STEIN-L",
  "Comb-STEIN-H" = "Comb-STEIN-H",
  "Runin-STEIN"  = "Runin-STEIN"
)
tite_labels <- c(
  "Comb-STEIN-L" = "TITE-Comb-STEIN-L",
  "Comb-STEIN-H" = "TITE-Comb-STEIN-H",
  "Runin-STEIN"  = "TITE-Runin-STEIN"
)

## 每个情景内共 3 种方法 x 2 种结局 = 6 根柱形。
## ggplot2 会把 position_dodge(width=) 视为整组总宽度并按系列数均分，
## 因此这里把单根柱宽表示为组宽的比例。
group_width <- 0.88
bar_fraction <- 0.85

## 找出 output/ 下所有合并结果文件，并按延迟设置保留最新一次运行。
discover_runs <- function(output_root) {
  files <- list.files(
    output_root,
    pattern = "^comb_stein_results_utility.*_all\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No combined result CSV found under ", output_root, call. = FALSE)
  }
  ## 时间戳包含日期和时间（run_YYYYMMDD_HHMMSS），必须整体捕获，
  ## 否则同一天内先后两次运行会被当成同一个运行。
  run_stamp <- sub(".*run_([0-9]+(_[0-9]+)?).*", "\\1", files)
  run_stamp[!grepl("run_[0-9]+", files)] <- "0"
  frames <- lapply(files, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  ## 时间戳按文件逐一写回，不能先 rbind 再整体赋值：那样 R 会按行循环
  ## 复用长度等于文件数的时间戳向量，把每个文件的标记写错。
  stamped <- Map(
    function(frame, stamp) {
      frame$.run_stamp <- stamp
      frame
    },
    frames,
    run_stamp
  )
  do.call(rbind, stamped)
}

raw <- discover_runs(output_root)

required <- c("Scenario", "Method", "Delay", "Selection", "Patient", "Sample", "Duration")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

## 每个延迟设置只保留时间戳最新的运行，避免历史目录重复计入。
raw <- raw[order(raw$.run_stamp, decreasing = TRUE), ]
raw <- raw[!duplicated(raw[c("Scenario", "Method", "Delay")]), ]
raw$.run_stamp <- NULL

## 情景按编号排序（S1、S2、...、S12）。
scenario_number <- as.integer(sub("^S", "", raw$Scenario))
raw$Scenario <- factor(
  raw$Scenario,
  levels = unique(raw$Scenario[order(scenario_number)])
)
raw$Method <- factor(raw$Method, levels = method_levels)
raw$Delay <- factor(raw$Delay, levels = delay_levels)
if (any(is.na(raw$Delay))) {
  stop("Unexpected Delay values; expected only: ",
       paste(delay_levels, collapse = ", "), call. = FALSE)
}

expected_rows <- length(levels(raw$Scenario)) * length(method_levels) * length(delay_levels)
if (nrow(raw) != expected_rows) {
  stop(
    "Expected ", expected_rows, " rows but found ", nrow(raw),
    "; result files are incomplete or inconsistent.",
    call. = FALSE
  )
}

write.csv(
  raw[order(raw$Scenario, raw$Method, raw$Delay),
      c("Scenario", "Method", "Delay", "Selection", "Patient", "Sample", "Duration")],
  file.path(figure_dir, "delay_comparison_metrics.csv"),
  row.names = FALSE
)

## 系列顺序决定柱形与图例的排列：同一方法的两根柱形相邻，
## 且完整结局在前、TITE 在后。
series_levels <- as.vector(t(outer(
  method_levels, delay_levels,
  function(method, delay) paste(method, delay, sep = "\n")
)))
series_colors <- setNames(
  c(
    method_colors[["Comb-STEIN-L"]], tite_colors[["Comb-STEIN-L"]],
    method_colors[["Comb-STEIN-H"]], tite_colors[["Comb-STEIN-H"]],
    method_colors[["Runin-STEIN"]],  tite_colors[["Runin-STEIN"]]
  ),
  series_levels
)
series_labels <- setNames(
  c(
    method_labels[["Comb-STEIN-L"]], tite_labels[["Comb-STEIN-L"]],
    method_labels[["Comb-STEIN-H"]], tite_labels[["Comb-STEIN-H"]],
    method_labels[["Runin-STEIN"]],  tite_labels[["Runin-STEIN"]]
  ),
  series_levels
)

panel_labels <- c(
  Selection = "A: OBDC selection rate (%)",
  Patient   = "B: Patients allocated to OBDC (%)",
  Sample    = "C: Average sample size (patients)",
  Duration  = "D: Trial duration (days)"
)
metric_levels <- names(panel_labels)

long <- do.call(rbind, lapply(metric_levels, function(metric) {
  data.frame(
    Scenario = raw$Scenario,
    Series = factor(
      paste(raw$Method, raw$Delay, sep = "\n"),
      levels = series_levels
    ),
    Metric = factor(metric, levels = metric_levels),
    Value = raw[[metric]]
  )
}))

plot <- ggplot(long, aes(x = Scenario, y = Value, fill = Series)) +
  geom_col(
    position = position_dodge(width = group_width),
    width = group_width * bar_fraction,
    colour = NA
  ) +
  facet_wrap(~Metric, ncol = 2, scales = "free_y",
             labeller = labeller(Metric = panel_labels)) +
  scale_fill_manual(
    values = series_colors,
    breaks = series_levels,
    labels = series_labels,
    name = "Method",
    guide = guide_legend(nrow = 2, byrow = FALSE)
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
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(0.85, "lines"),
    legend.margin = margin(2, 2, 2, 2),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(6, 10, 6, 6)
  )

figure_width <- 10
figure_height <- 9

ggsave(file.path(figure_dir, "fig_delay_comparison_metrics.png"), plot,
       width = figure_width, height = figure_height, dpi = 300)
ggsave(file.path(figure_dir, "fig_delay_comparison_metrics.pdf"), plot,
       width = figure_width, height = figure_height, device = grDevices::cairo_pdf)

cat("Saved figure and table to:", figure_dir, "\n")
cat("Scenarios:", nlevels(raw$Scenario),
    " methods:", length(method_levels),
    " delays:", length(delay_levels), "\n")
