## Trajectory of one simulated combination trial: cohort process table, dose
## escalation of each agent (panel A), calendar time with delayed (TITE) vs
## complete outcomes (panel B). Font sizes are fixed, so set --width/--height to
## the journal column width to keep text legible after scaling.
##
## Usage: Rscript Comb-STEIN-Code/plot/plot_trial_trajectory.R [project root]
##       [--scenario=C1] [--method=CombSTEIN-L] [--seed=2026]
##       [--width=10.5] [--height=5.2] [--dpi=600]
##       [--figure-dir=output/figures] [--table-dir=output/tables]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

## Returns the --name=value option together with the remaining args.
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
## Delay settings may be overridden from the CLI; defaults follow DEFAULT_DELAY_CONTROL
## so the figure uses the same windows and accrual rate as the main analysis.
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

## Relative paths expand against the project root; absolute paths are kept as-is.
resolve.dir <- function(value) {
  if (grepl("^(/|[A-Za-z]:)", value)) value else file.path(project_root, value)
}
figure_dir <- resolve.dir(figure.opt$value)
table_dir <- resolve.dir(table.opt$value)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

## Module order matches main.R.
code_dir <- file.path(project_root, "Comb-STEIN-Code")
for (module in c(
  "config.R", "time_to_event.R", "isotonic_estimators.R", "dose_rules.R",
  "scenarios.R", "scenarios_case_study.R", "simulate_stein.R"
)) {
  sys.source(file.path(code_dir, module), envir = globalenv())
}

seed <- if (is.null(seed.opt$value)) DEFAULT_SEED else as.integer(seed.opt$value)

delay.control <- DEFAULT_DELAY_CONTROL
if (!is.null(accrual.opt$value)) delay.control$accrual <- as.numeric(accrual.opt$value)
if (!is.null(tox.opt$value)) delay.control$tox_window <- as.numeric(tox.opt$value)
if (!is.null(eff.opt$value)) delay.control$eff_window <- as.numeric(eff.opt$value)
if (!is.null(susp.opt$value)) delay.control$suspension <- as.numeric(susp.opt$value)
days_per_month <- 30.4375

## R on Chinese Windows runs in the C locale, where writeLines escapes non-ASCII
## (μg becomes <U+03BC>g); write UTF-8 bytes through a binary connection instead.
micro.sign <- intToUtf8(0x3bc)
write.utf8 <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con, useBytes = TRUE)
}
## write.csv translates to the current locale, so CSV is assembled by hand too.
csv.field <- function(x) {
  x <- as.character(x)
  ifelse(grepl('["\r\n]', x), paste0('"', gsub('"', '""', x), '"'), x)
}
to.csv.lines <- function(data) {
  ## Cast to a character matrix first: apply() on a data.frame pads numeric columns.
  chars <- as.matrix(data.frame(
    lapply(data, as.character),
    stringsAsFactors = FALSE, check.names = FALSE
  ))
  c(
    paste(csv.field(colnames(chars)), collapse = ","),
    apply(chars, 1, function(row) paste(csv.field(row), collapse = ","))
  )
}

## Dose-level labels for the case-study scenarios, one entry per grid in
## Comb-STEIN-Code/scenarios_case_study.R; for plotting and tables only.
## After the 2026-09-22 transpose the agents are swapped: drug A is the two-level
## agent, drug B the three-level one. Order must match scenarios_case_study.R.
dose.labels <- list(
  C1 = list(drugA = "GM-CSF", a = c(125, 250), a.unit = paste0(micro.sign, "g"),
            drugB = "GP2", b = c(100, 500, 1000), b.unit = paste0(micro.sign, "g")),
  C2 = list(drugA = "GM-CSF", a = c(125, 250), a.unit = paste0(micro.sign, "g"),
            drugB = "GP2", b = c(100, 500, 1000), b.unit = paste0(micro.sign, "g")),
  ## IDOV-SAFE levels are 1e9 / 3e9 / 1e10 PFU. A plot axis needs one continuous
  ## scale, so ticks use 1 / 3 / 10 and b.text carries the full labels, keeping
  ## 1e+09-style notation out of tables and captions.
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
  ## The 3x3 ODI-2001 x ipilimumab grid is symmetric, so the transpose leaves it unchanged.
  C5 = list(drugA = "ODI-2001", a = c(1, 2, 4), a.unit = "mg",
            drugB = "Ipilimumab", b = c(2.5, 5, 10), b.unit = "mg"),
  C6 = list(drugA = "ODI-2001", a = c(1, 2, 4), a.unit = "mg",
            drugB = "Ipilimumab", b = c(2.5, 5, 10), b.unit = "mg")
)
if (!scenario %in% names(dose.labels)) {
  stop(sprintf("Scenario %s has no dose-level labels; add them to dose.labels", scenario))
}
doses <- dose.labels[[scenario]]

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

## True OBDC comes from the design utility function, not hard-coded in the scenario.
obdc <- get.true.obdc(truth$pT, truth$pE)

dose.text <- function(dose, drug, values, unit, texts = NULL) {
  ## texts gives the full spelling of each level (needed for PFU); otherwise
  ## "drug + value + unit" is pasted together.
  if (!is.null(texts)) return(texts[dose])
  sprintf("%s %g %s", drug, values[dose], unit)
}

## Per-cohort process table
## One row per cohort, aggregated from the per-patient TITE trace.
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
## Days are continuous in the time-to-event simulation; one decimal suffices.
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

## Trial-level summary, written as the table header.
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

## Shared journal theme
## Font sizes are fixed and do not scale with the figure dimensions.
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

## Panel A: dose escalation of each agent
## Both agents share the cohort axis but have different dose scales, so they are
## split into two stacked subplots, each with ticks at its own real dose levels.
assignments <- process[, c("cohort", "dose_a_level", "dose_b_level")]
cohort.max <- max(assignments$cohort)
## Units differ between agents (e.g. x10^9 PFU vs mg), so each track carries its own.
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
  ## Label only cohorts where the dose changes (plus the first), to avoid repeats.
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

## Panel B: calendar time of the two run modes
## Colour-blind blue/vermilion plus solid/dashed linetypes, so the two modes stay
## distinguishable in black-and-white printing.
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

## Figure caption
## The figure itself carries no explanatory text; the caption holds the full
## description plus reproducibility details (seed, R version, simulation settings).
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

## The C-locale console escapes μg, so print an ASCII spelling instead.
cat(gsub(micro.sign, "u", header, fixed = TRUE), sep = "\n")
cat("\nCohort-by-cohort table saved to:\n")
cat(" ", file.path(table_dir, paste0(file.stem, ".csv")), "\n")
cat(" ", file.path(table_dir, paste0(file.stem, ".md")), "\n")
cat("Figure caption saved to:\n")
cat(" ", file.path(table_dir, paste0(file.stem, "_caption.txt")), "\n")
cat("Figure saved to", file.path(figure_dir, paste0(figure.file, ".pdf")), "\n")
