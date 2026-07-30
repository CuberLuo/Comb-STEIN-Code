## Comb-STEIN command-line entry point and public loader.
##
## Usage:
##   Rscript Comb-STEIN.R
##   Rscript Comb-STEIN.R S1 S2 S3 200
##   Rscript Comb-STEIN.R S2 200 --obdc=table1
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

run.comb.stein.cli <- function(
  args = commandArgs(trailingOnly = TRUE),
  results.dir = file.path(.comb_stein_project_root, "output")
) {
  set.seed(DEFAULT_SEED)
  cat(sprintf("Derived toxicity cutoff  phi = %.4f\n", phi))
  cat(sprintf("Derived efficacy  cutoff psi = %.4f\n", psi))

  ntrial <- 1000
  obdc.definition <- "utility"
  obdc.arg <- grep("^--obdc=", args, value = TRUE)
  if (length(obdc.arg) > 0) {
    obdc.definition <- sub("^--obdc=", "", obdc.arg[1])
    args <- args[!grepl("^--obdc=", args)]
  }
  if (!obdc.definition %in% c("utility", "table1")) {
    stop("--obdc must be 'utility' or 'table1'")
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

  dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
  all.results <- NULL

  for (name in scenarios.to.run) {
    cat("Running scenario", name, "with ntrial =", ntrial, "...\n")
    result <- run.scenario(
      name,
      ntrial = ntrial,
      obdc.definition = obdc.definition
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
    paste0("comb_stein_results_", obdc.definition, "_all.csv")
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
