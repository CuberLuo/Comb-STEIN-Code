## Comb-STEIN command-line entry point and public loader.
##
## Usage:
##   Rscript Comb-STEIN.R
##   Rscript Comb-STEIN.R S1 S2 S3 50
##   Rscript Comb-STEIN.R S2 50 --obdc=table1
##   Rscript Comb-STEIN.R S1 50 --delayed=true
##   Rscript Comb-STEIN.R S1 50 --tite --accrual=10 --tox-window=30
##     --eff-window=90 --suspension=0.5
##
## The file may also be sourced to load the public functions without starting
## a simulation:
##   source("Comb-STEIN.R")

.comb_stein_script_dir <- function() {
  frame.files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1]
    },
    character(1)
  )
  frame.files <- frame.files[nzchar(frame.files)]
  if (length(frame.files) > 0) {
    return(dirname(normalizePath(
      frame.files[length(frame.files)],
      mustWork = TRUE
    )))
  }

  file.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file.arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file.arg[1]),
      mustWork = TRUE
    )))
  }
  normalizePath(getwd(), mustWork = TRUE)
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
    file.path(.comb_stein_root, "R", .comb_stein_module),
    envir = environment()
  )
}

.parse.cli.boolean <- function(value, option) {
  value <- tolower(as.character(value)[1])
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(option, " must be true or false")
}

.pull.cli.option <- function(args, option, default = NULL) {
  pattern <- paste0("^--", option, "=")
  hits <- grep(pattern, args)
  if (length(hits) > 1L) stop("--", option, " may be specified only once")
  if (length(hits) == 0L) return(list(value = default, args = args))
  list(
    value = sub(pattern, "", args[hits]),
    args = args[-hits]
  )
}

run.comb.stein.cli <- function(
  args = commandArgs(trailingOnly = TRUE),
  results.dir = file.path(.comb_stein_project_root, "output")
) {
  set.seed(DEFAULT_SEED)
  cat(sprintf("Derived toxicity cutoff  phi = %.4f\n", phi))
  cat(sprintf("Derived efficacy  cutoff psi = %.4f\n", psi))

  ntrial <- DEFAULT_NTRIAL

  pulled <- .pull.cli.option(args, "obdc", "utility")
  obdc.definition <- pulled$value
  args <- pulled$args
  if (!obdc.definition %in% c("utility", "table1")) {
    stop("--obdc must be 'utility' or 'table1'")
  }

  delayed <- FALSE
  delayed.flag <- .pull.cli.option(args, "delayed", NULL)
  args <- delayed.flag$args
  if (!is.null(delayed.flag$value)) {
    delayed <- .parse.cli.boolean(delayed.flag$value, "--delayed")
  }
  has.tite <- "--tite" %in% args
  has.no.tite <- "--no-tite" %in% args
  if (has.tite && has.no.tite) {
    stop("Use only one of --tite and --no-tite")
  }
  if (!is.null(delayed.flag$value) && (has.tite || has.no.tite)) {
    stop("Do not combine --delayed=... with --tite or --no-tite")
  }
  if (has.tite) delayed <- TRUE
  if (has.no.tite) delayed <- FALSE
  args <- args[!args %in% c("--tite", "--no-tite")]

  delay.control <- list()
  numeric.options <- c(
    accrual = "accrual",
    `tox-window` = "tox_window",
    `eff-window` = "eff_window",
    suspension = "suspension"
  )
  for (option in names(numeric.options)) {
    pulled <- .pull.cli.option(args, option, NULL)
    args <- pulled$args
    if (!is.null(pulled$value)) {
      value <- suppressWarnings(as.numeric(pulled$value))
      if (!is.finite(value)) stop("--", option, " must be numeric")
      delay.control[[numeric.options[[option]]]] <- value
    }
  }
  logical.options <- c(
    `use-suspension` = "use_suspension",
    `accrual-random` = "accrual_random"
  )
  for (option in names(logical.options)) {
    pulled <- .pull.cli.option(args, option, NULL)
    args <- pulled$args
    if (!is.null(pulled$value)) {
      delay.control[[logical.options[[option]]]] <-
        .parse.cli.boolean(pulled$value, paste0("--", option))
    }
  }
  if (length(delay.control) > 0L && !delayed) {
    stop("Delay settings require --tite or --delayed=true")
  }
  if (delayed) delay.control <- normalize.delay.control(delay.control)

  unknown.options <- args[grepl("^--", args)]
  if (length(unknown.options) > 0L) {
    stop("Unknown option(s): ", paste(unknown.options, collapse = ", "))
  }

  scenarios.to.run <- names(scenarios)
  if (length(args) > 0) {
    is.numeric.arg <- suppressWarnings(!is.na(as.numeric(args)))
    if (any(is.numeric.arg)) {
      ntrial <- as.numeric(args[is.numeric.arg])[1]
    }
    scenario.args <- args[!is.numeric.arg]
    if (length(scenario.args) > 0) scenarios.to.run <- scenario.args
  }

  unknown <- setdiff(scenarios.to.run, names(scenarios))
  if (length(unknown) > 0) {
    stop("Unknown scenario(s): ", paste(unknown, collapse = ", "))
  }

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
    cat("Delayed outcomes: OFF (complete outcomes used at each decision)\n")
  }

  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  all.results <- NULL

  for (name in scenarios.to.run) {
    cat("Running scenario", name, "with ntrial =", ntrial, "...\n")
    result <- run.scenario(
      name,
      ntrial = ntrial,
      obdc.definition = obdc.definition,
      delayed = delayed,
      delay.control = delay.control
    )
    true.ij <- attr(result, "true.obdcs")
    labels <- apply(
      true.ij,
      1,
      function(ij) sprintf("(%d,%d)", ij[1], ij[2])
    )
    cat(
      "  true OBDC(s) [",
      obdc.definition,
      "]: ",
      paste(labels, collapse = ", "),
      "\n",
      sep = ""
    )

    scenario.outfile <- file.path(
      results.dir,
      paste0(
        "comb_stein_results_",
        obdc.definition,
        if (delayed) "_tite" else "",
        "_",
        name,
        ".csv"
      )
    )
    write.csv(result, scenario.outfile, row.names = FALSE)
    cat("  scenario results saved to", scenario.outfile, "\n")
    all.results <- rbind(all.results, result)
  }

  cat("\n\n===================== Table 2 simulation =====================\n")
  print(all.results, row.names = FALSE)

  outfile <- file.path(
    results.dir,
    paste0(
      "comb_stein_results_",
      obdc.definition,
      if (delayed) "_tite" else "",
      "_all.csv"
    )
  )
  write.csv(all.results, outfile, row.names = FALSE)
  cat("\nResults saved to", outfile, "\n")
  invisible(all.results)
}

rm(.comb_stein_module)
rm(.comb_stein_modules)
rm(.comb_stein_root)

if (sys.nframe() == 0L) {
  run.comb.stein.cli()
}
