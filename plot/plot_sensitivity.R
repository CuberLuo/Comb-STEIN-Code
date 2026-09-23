## Supplementary sensitivity figures: suspension ratio, accrual rate, cohort size.
## Same scenarios and seed throughout; baseline AW=(1,2), accrual 1.5/month,
## suspension ratio 0.5, cohort size 3. Each analysis writes a four-panel figure
## (selection / allocation / sample size / duration) plus a metrics CSV.
## Runs are listed in output/sensitivity/runs.csv (Comb-STEIN-Code/tmp/run_sensitivity.R).
## Usage: Rscript Comb-STEIN-Code/plot/plot_sensitivity.R [project root]
##        [--manifest=output/sensitivity/runs.csv] [--figure-dir=output/figures]

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
  stop("Run manifest not found: ", manifest_path, call. = FALSE)
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

## Lighten towards white so one method keeps one hue across its settings.
.mix.white <- function(color, prop) {
  rgb <- grDevices::col2rgb(color)[, 1]
  mixed <- round(rgb + (255 - rgb) * prop)
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

## Whitening fraction of setting k of n: first setting darkest, last lightest.
.shade_of <- function(k, n) {
  if (n <= 1) 0 else 0.55 * (k - 1) / (n - 1)
}

.load.run <- function(id, label) {
  hit <- manifest$id == id
  if (!any(hit)) stop("Run manifest is missing configuration: ", id, call. = FALSE)
  row <- manifest[which(hit)[1], ]
  path <- file.path(project_root, row$run_dir, row$csv)
  if (!file.exists(path)) stop("Run results not found: ", path, call. = FALSE)
  frame <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  frame$Run <- row$id
  frame$Setting <- label
  frame
}

## The three sensitivity analyses; order of `entries` fixes the shade order.
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
  ## Series order must match manuscript Figs 3-4: settings of one method stay
  ## adjacent (dark/light pair), methods in sequence. So the expansion has to be
  ## method-outer / setting-inner; lapply(entries, ...) would group the same
  ## setting together and break the colour pairing.
  ## outer() cannot be used here: it expands indices into a vector subscript,
  ## and `[[` does recursive indexing on vector subscripts.
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

  ## Colours must be built in the same method-outer / setting-inner order as
  ## series_levels, otherwise setNames mismatches shades across series.
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
      ## byrow = FALSE fills the legend by column, so each column is the
      ## dark/light pair of one method, as in manuscript Figs 3-4.
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
