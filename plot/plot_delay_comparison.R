## Operating characteristics of delayed (TITE) vs complete (non-TITE) outcomes for
## Selection / Patient / Sample / Duration, from the latest combined results in output/.
## Two bars per method: dark = complete outcomes, light = TITE.
##
## Usage: Rscript Comb-STEIN-Code/plot/plot_delay_comparison.R [project root]
##       [--complete-run=output/run_...] [--tite-run=output/run_...]
##       [--out=fig_delay_comparison_metrics]
##
## Default: newest run of each delay setting under output/. When one batch holds
## several settings (e.g. AW=(1,2) vs AW=(2,3)), pass --complete-run/--tite-run.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

.pull.option <- function(args, name) {
  pattern <- paste0("^--", name, "=")
  hit <- grep(pattern, args, value = TRUE)
  if (length(hit) == 0) {
    return(list(value = NULL, args = args))
  }
  list(value = sub(pattern, "", hit[[1]]), args = args[!grepl(pattern, args)])
}

complete.opt <- .pull.option(args, "complete-run")
args <- complete.opt$args
tite.opt <- .pull.option(args, "tite-run")
args <- tite.opt$args
out.opt <- .pull.option(args, "out")
args <- out.opt$args

project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!dir.exists(file.path(project_root, "output"))) {
  stop("No output/ directory under ", project_root, call. = FALSE)
}

out_name <- if (is.null(out.opt$value)) {
  "fig_delay_comparison_metrics"
} else {
  out.opt$value
}

.run.csv <- function(run_dir, delay) {
  path <- if (grepl("^(/|[A-Za-z]:)", run_dir)) {
    run_dir
  } else {
    file.path(project_root, run_dir)
  }
  file <- file.path(
    path,
    if (delay == "TITE") {
      "comb_stein_results_utility_tite_all.csv"
    } else {
      "comb_stein_results_utility_all.csv"
    }
  )
  if (!file.exists(file)) {
    stop("Run results not found for ", delay, ": ", file, call. = FALSE)
  }
  file
}

output_root <- file.path(project_root, "output")
figure_dir <- file.path(output_root, "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

method_levels <- c("Comb-STEIN-L", "Comb-STEIN-H", "Runin-STEIN")
delay_levels <- c("Complete", "TITE")

## Same palette as manuscript Figs 3/4: dark = complete outcomes, light = TITE.
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

## Legend labels follow the manuscript naming; the TITE variants get a TITE- prefix.
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

## Six bars per scenario (3 methods x 2 outcome modes). ggplot2 treats
## position_dodge(width=) as the width of the whole group and splits it across the
## series, so the bar width is expressed as a fraction of the group width.
group_width <- 0.88
bar_fraction <- 0.85

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
  ## run_YYYYMMDD_HHMMSS includes the time; capture it as a whole, otherwise two
  ## runs on the same day collapse into one.
  run_stamp <- sub(".*run_([0-9]+(_[0-9]+)?).*", "\\1", files)
  run_stamp[!grepl("run_[0-9]+", files)] <- "0"
  frames <- lapply(files, function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  ## Stamps are written per file before rbind: assigning after rbind recycles the
  ## length-n_files vector along the rows and mislabels each frame.
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

## Keep only the newest run per delay setting, so stale directories are not counted twice.
raw <- raw[order(raw$.run_stamp, decreasing = TRUE), ]
raw <- raw[!duplicated(raw[c("Scenario", "Method", "Delay")]), ]

## An explicit run replaces the whole auto-discovered block of that Delay, so runs
## from different settings (e.g. AW=(1,2) vs AW=(2,3)) never mix.
.replace.delay <- function(raw, delay, run_dir) {
  frame <- utils::read.csv(
    .run.csv(run_dir, delay),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  frame$Delay <- delay
  frame$.run_stamp <- "explicit"
  rbind(raw[raw$Delay != delay, , drop = FALSE], frame)
}
if (!is.null(complete.opt$value)) {
  raw <- .replace.delay(raw, "Complete", complete.opt$value)
}
if (!is.null(tite.opt$value)) {
  raw <- .replace.delay(raw, "TITE", tite.opt$value)
}
raw$.run_stamp <- NULL

## Scenarios are ordered by number (S1...S12, or C1...C6 for case studies). The
## prefix must accept both S and C: stripping only S gives NA for the case
## studies, order() falls back to the original order, and the scenario order in
## the figure then depends on directory scan order.
scenario_number <- as.integer(sub("^[SC]", "", raw$Scenario))
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

## Figure is <out_name>.png; the data table drops the fig_ prefix and appends
## _metrics.csv, which keeps the existing file names. A trailing _metrics in
## out_name is stripped first to avoid ..._metrics_metrics.csv.
data_base <- sub("_metrics$", "", sub("^fig_", "", out_name))
data_path <- file.path(figure_dir, paste0(data_base, "_metrics.csv"))

write.csv(
  raw[order(raw$Scenario, raw$Method, raw$Delay),
      c("Scenario", "Method", "Delay", "Selection", "Patient", "Sample", "Duration")],
  data_path,
  row.names = FALSE
)

## Series order fixes both the bar and the legend order: the two bars of a method
## sit adjacent, complete outcomes first, then TITE.
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

ggsave(file.path(figure_dir, paste0(out_name, ".png")), plot,
       width = figure_width, height = figure_height, dpi = 300)
ggsave(file.path(figure_dir, paste0(out_name, ".pdf")), plot,
       width = figure_width, height = figure_height, device = grDevices::cairo_pdf)

cat("Saved figure and table to:", figure_dir, "\n")
cat("Scenarios:", nlevels(raw$Scenario),
    " methods:", length(method_levels),
    " delays:", length(delay_levels), "\n")
