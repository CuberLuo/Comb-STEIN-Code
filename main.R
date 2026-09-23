## Comb-STEIN command-line entry point and module loader.
##
## Usage:
##   Rscript Comb-STEIN-Code/main.R
##   Rscript Comb-STEIN-Code/main.R S1 S2 S3 50
##   Rscript Comb-STEIN-Code/main.R S1 50 --delayed=true
##   Rscript Comb-STEIN-Code/main.R S1 50 --tite --accrual=10 --tox-window=30
##     --eff-window=90 --suspension=0.5
##   Rscript Comb-STEIN-Code/main.R S1 50 --workers=4
##   Rscript Comb-STEIN-Code/main.R S1 50 --quiet
##
## Each run writes into a new timestamped folder under output/.

## Locate the directory holding this script: source() frame, then Rscript
## --file, then the working directory.
.comb_stein_script_dir <- function() {
  ## RStudio's debugSource() can push a working-directory path onto the stack,
  ## so don't blindly take the last sys.frames()$ofile.
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
    ## Prefer a dir that really holds the modules (Comb-STEIN-Code/), not the project root.
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
  if (file.exists(file.path(working.dir, "Comb-STEIN-Code", "config.R"))) {
    return(normalizePath(file.path(working.dir, "Comb-STEIN-Code"), mustWork = TRUE))
  }
  stop("Cannot locate the directory of main.R: config.R not found")
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
  "simulate_runin.R",
  "simulate_stein.R",
  "scenarios.R",
  "scenarios_case_study.R",
  "run_scenarios.R"
)

for (.comb_stein_module in .comb_stein_modules) {
  sys.source(
    file.path(.comb_stein_root, .comb_stein_module),
    envir = environment()
  )
}

.format.dose.matrix <- function(x, digits = 2) {
  if (is.null(dim(x))) return(format(round(x, digits), nsmall = digits))
  rows <- apply(x, 1, function(row) {
    paste(sprintf(paste0("%8.", digits, "f"), row), collapse = " ")
  })
  paste(rows, collapse = "\n")
}

.write.scenario.txt <- function(path, name, result, pT, pE, ntrial) {
  details <- attr(result, "dose.details")
  true.ij <- attr(result, "true.obdcs")
  labels <- apply(true.ij, 1, function(ij) sprintf("(%d,%d)", ij[1], ij[2]))
  lines <- c(
    "CombSTEIN simulation results",
    paste0("Scenario: ", name),
    paste0("Number of trials: ", ntrial),
    paste0(
      "Outcome mode: ",
      if (isTRUE(attr(result, "delayed"))) "TITE" else "Complete"
    ),
    paste0("True OBDC(s): ", paste(labels, collapse = ", ")),
    "Rows: Drug A dose levels i=1..I; columns: Drug B dose levels j=1..J",
    "",
    "True toxicity:",
    .format.dose.matrix(pT, 2),
    "",
    "True efficacy:",
    .format.dose.matrix(pE, 2),
    ""
  )

  for (method in names(details)) {
    d <- details[[method]]
    lines <- c(
      lines,
      paste0("================ ", method, " ================"),
      "Number of patients treated at each dose (average):",
      .format.dose.matrix(d$average_patients, 2),
      "",
      "Patient allocation at each dose (%):",
      .format.dose.matrix(d$allocation_percent, 2),
      "",
      "Number of toxicity observed at each dose (average):",
      .format.dose.matrix(d$average_toxicity, 2),
      "",
      "Number of response observed at each dose (average):",
      .format.dose.matrix(d$average_efficacy, 2),
      "",
      "Selection probabilities (%):",
      .format.dose.matrix(d$selection_percent, 2),
      ""
    )
  }
  writeLines(lines, path, useBytes = TRUE)
}

## `option` is kept for error reporting only.
.parse.cli.boolean <- function(value, option) {
  value <- tolower(as.character(value)[1])
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  FALSE
}

.pull.cli.option <- function(args, option, default = NULL) {
  pattern <- paste0("^--", option, "=")
  hits <- grep(pattern, args)
  if (length(hits) == 0L) return(list(value = default, args = args))
  list(
    value = sub(pattern, "", args[hits]),
    args = args[-hits]
  )
}

## Appends a numeric suffix when the timestamp folder is already taken.
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

## Runs the requested scenarios and writes CSV/TXT reports; returns the
## combined table invisibly. `results.dir` overrides the timestamped default.
run.comb.stein.cli <- function(
  args = commandArgs(trailingOnly = TRUE),
  results.dir = NULL
) {
  run.started.at <- Sys.time()
  if (is.null(results.dir)) {
    ## Timestamped folder under output/ so earlier runs are not overwritten.
    results.dir <- .create.timestamped.results.dir(
      file.path(.comb_stein_project_root, "output"),
      run.started.at
    )
  } else {
    dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
    results.dir <- normalizePath(results.dir, mustWork = TRUE)
  }

  ## --quiet silences interim simulation output only; summary, utility matrix
  ## and result tables still print.
  quiet.flag <- .pull.cli.option(args, "quiet", NULL)
  args <- quiet.flag$args
  quiet <- !is.null(quiet.flag$value) &&
    .parse.cli.boolean(quiet.flag$value, "--quiet")
  if ("--quiet" %in% args) {
    quiet <- TRUE
    args <- args[args != "--quiet"]
  }
  options(comb.stein.quiet = quiet)

  ## Fixed seed so identical arguments give reproducible results.
  set.seed(DEFAULT_SEED)
  ## phi: toxicity cutoff, psi: efficacy cutoff, both derived from config.
  cat("Run started:", format(run.started.at, "%Y-%m-%d %H:%M:%S %Z"), "\n")
  cat("Output folder:", results.dir, "\n")
  cat(sprintf("Derived toxicity cutoff  phi = %.4f\n", phi))
  cat(sprintf("Derived efficacy  cutoff psi = %.4f\n", psi))

  ntrial <- DEFAULT_NTRIAL
  workers <- DEFAULT_WORKERS
  workers.flag <- .pull.cli.option(args, "workers", NULL)
  args <- workers.flag$args
  if (!is.null(workers.flag$value)) {
    workers <- as.integer(workers.flag$value)
  }

  ## delayed = TRUE switches to TITE/delayed outcomes; default is complete outcomes.
  delayed <- FALSE
  delayed.flag <- .pull.cli.option(args, "delayed", NULL)
  args <- delayed.flag$args
  if (!is.null(delayed.flag$value)) {
    delayed <- .parse.cli.boolean(delayed.flag$value, "--delayed")
  }
  ## --tite / --no-tite shorthand; dropped from args so they are not read as scenario names.
  has.tite <- "--tite" %in% args
  has.no.tite <- "--no-tite" %in% args
  if (has.tite) delayed <- TRUE
  if (has.no.tite) delayed <- FALSE
  args <- args[!args %in% c("--tite", "--no-tite")]

  ## Maps CLI option names to internal delay.control names.
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
      delay.control[[numeric.options[[option]]]] <-
        as.numeric(pulled$value)
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
  ## Accrual rate and assessment windows apply to both complete and TITE modes.
  delay.control <- normalize.delay.control(delay.control)

  ## Methods may be given as a comma-separated list.
  methods.flag <- .pull.cli.option(args, "methods", NULL)
  args <- methods.flag$args
  methods <- c(
    "Comb-STEIN-L",
    "Comb-STEIN-H",
    "Runin-STEIN"
  )
  if (!is.null(methods.flag$value)) {
    methods <- trimws(
      strsplit(methods.flag$value, ",", fixed = TRUE)[[1]]
    )
  }

  ## Simulation code reads `cohortsize` from the global env, so the override
  ## must be assigned there.
  cohort.size.flag <- .pull.cli.option(args, "cohort-size", NULL)
  args <- cohort.size.flag$args
  if (!is.null(cohort.size.flag$value)) {
    cohort.size.value <- suppressWarnings(
      as.integer(cohort.size.flag$value)
    )
    if (is.na(cohort.size.value) || cohort.size.value < 1L) {
      stop("--cohort-size must be a positive integer", call. = FALSE)
    }
    assign("cohortsize", cohort.size.value, envir = globalenv())
    cat(sprintf("Cohort size overridden to %d\n", cohortsize))
  }

  ## Defaults to all scenarios; numeric args set ntrial, others name scenarios.
  scenarios.to.run <- names(scenarios)
  if (length(args) > 0) {
    is.numeric.arg <- suppressWarnings(!is.na(as.numeric(args)))
    if (any(is.numeric.arg)) {
      ntrial <- as.numeric(args[is.numeric.arg])[1]
    }
    scenario.args <- args[!is.numeric.arg]
    if (length(scenario.args) > 0) scenarios.to.run <- scenario.args
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

  all.results <- NULL

  for (name in scenarios.to.run) {
    ## OBDC = combinations within the configured utility tolerance of the maximum.
    cat("Running scenario", name, "with ntrial =", ntrial, "...\n")

    scenario.truth <- build.pTpE(name)
    true.utility <- get.true.utility(scenario.truth$pT, scenario.truth$pE)
    dimnames(true.utility) <- list(
      paste0("i=", seq_len(nrow(true.utility))),
      paste0("j=", seq_len(ncol(true.utility)))
    )
    cat("  true utility matrix (rows i, columns j):\n")
    utility.display <- matrix(
      sprintf("%.4f", true.utility),
      nrow = nrow(true.utility),
      ncol = ncol(true.utility),
      dimnames = dimnames(true.utility)
    )
    print(utility.display, quote = FALSE)

    result <- run.scenario(
      name,
      ntrial = ntrial,
      delayed = delayed,
      delay.control = delay.control,
      workers = workers,
      methods = methods
    )
    ## true.ij: two-column matrix, one row per tied optimal combination (row, column).
    true.ij <- attr(result, "true.obdcs")
    labels <- apply(
      true.ij,
      1,
      function(ij) sprintf("(%d,%d)", ij[1], ij[2])
    )
    cat(
      "  true OBDC(s) [within utility tolerance]: ",
      paste(labels, collapse = ", "),
      "\n",
      sep = ""
    )

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
    scenario.txtfile <- file.path(
      results.dir,
      paste0(
        "comb_stein_results_utility",
        if (delayed) "_tite" else "",
        "_",
        name,
        ".txt"
      )
    )
    .write.scenario.txt(
      scenario.txtfile,
      name,
      result,
      scenario.truth$pT,
      scenario.truth$pE,
      ntrial
    )
    cat("  scenario details saved to", scenario.txtfile, "\n")
    all.results <- rbind(all.results, result)
  }

  cat("\n\n===================== Table 2 simulation =====================\n")
  print(all.results, row.names = FALSE)

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
  invisible(all.results)
}

run.comb.stein.cli()
