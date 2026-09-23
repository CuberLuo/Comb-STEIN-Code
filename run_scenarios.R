## Replicated trial simulation and aggregation of operating characteristics.

## Replicates every design method under one named scenario.
## Returns a method-level table; the true OBDC and delay settings are attached as attributes.
run.scenario <- function(
  name,
  ntrial = DEFAULT_NTRIAL,
  delayed = FALSE,
  delay.control = NULL,
  workers = DEFAULT_WORKERS,
  methods = c(
    "Comb-STEIN-L",
    "Comb-STEIN-H",
    "Runin-STEIN"
  )
) {
  delayed <- isTRUE(delayed)
  delay.control <- normalize.delay.control(delay.control)

  scenario <- build.pTpE(name)
  pT <- scenario$pT
  pE <- scenario$pE
  true.ij <- get.true.obdc(pT, pE)

  out <- data.frame(
    Scenario = name,
    Method = methods,
    Delay = if (delayed) "TITE" else "Complete",
    Selection = NA,
    ET = NA,
    Patient = NA,
    Sample = NA,
    Duration = NA
  )

  # method x replicate are expanded into independent parallel tasks; each worker keeps
  # its own RNG stream and JAGS/cache state in a separate R process, so the main
  # process aggregates only at the end and no shared counter can race.
  tasks <- expand.grid(
    method = methods,
    trial = seq_len(ntrial),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  quiet <- isTRUE(getOption("comb.stein.quiet", FALSE))
  task.fun <- function(task) {
    # multisession workers do not inherit options, sync the debug switch explicitly
    options(comb.stein.quiet = quiet)
    method <- task$method[[1]]
    trial <- as.integer(task$trial[[1]])
    result <- tryCatch(
      switch(
        method,
        "Comb-STEIN-L" = simulate.trial.comb.stein(
          pT, pE, "CombSTEIN-L", delayed, delay.control
        ),
        "Comb-STEIN-H" = simulate.trial.comb.stein(
          pT, pE, "CombSTEIN-H", delayed, delay.control
        ),
        "Runin-STEIN" = simulate.trial.runin.stein(
          pT, pE, delayed, delay.control
        ),
        stop(sprintf("Unknown method: %s", method), call. = FALSE)
      ),
      error = function(err) {
        stop(
          sprintf(
            "Scenario %s, method %s, trial %d failed: %s",
            name,
            method,
            trial,
            conditionMessage(err)
          ),
          call. = FALSE
        )
      }
    )
    list(method = method, trial = trial, result = result)
  }

  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  }
  previous.plan <- future::plan()
  on.exit(future::plan(previous.plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  task.results <- progressr::with_progress(
    {
      progress <- progressr::progressor(steps = nrow(tasks))
      future.apply::future_lapply(
        seq_len(nrow(tasks)),
        function(index) {
          task <- tasks[index, , drop = FALSE]
          result <- task.fun(task)
          progress(message = sprintf(
            "%s, trial %d/%d",
            task$method[[1]],
            task$trial[[1]],
            ntrial
          ))
          result
        },
        future.seed = TRUE,
        future.packages = "rjags"
      )
    },
    handlers = progressr::handler_txtprogressbar,
    # --quiet only suppresses simulation debug output; the progress bar stays on
    enable = TRUE
  )

  # per-combination counts for the scenario TXT report
  dose.details <- list()
  for (method in methods) {
    selection.count <- 0
    sample.sum <- 0
    duration.sum <- 0
    selection.matrix <- matrix(0, nrow = nrow(pT), ncol = ncol(pT))
    patient.matrix <- matrix(0, nrow = nrow(pT), ncol = ncol(pT))
    tox.matrix <- matrix(0, nrow = nrow(pT), ncol = ncol(pT))
    eff.matrix <- matrix(0, nrow = nrow(pT), ncol = ncol(pT))

    method.results <- Filter(
      function(task) identical(task$method, method),
      task.results
    )
    for (task in method.results) {
      result <- task$result

      # result$n: enrolled per combination; x / y: observed toxicity / efficacy events;
      # sel: final selection; duration: trial duration (days)
      total.n <- sum(result$n)
      patient.matrix <- patient.matrix + result$n
      tox.matrix <- tox.matrix + result$x
      eff.matrix <- eff.matrix + result$y
      sample.sum <- sample.sum + total.n
      duration.sum <- duration.sum + result$duration

      # true.ij may contain several tied optimal combinations, hence the row-wise compare
      if (!any(is.na(result$sel))) {
        selection.matrix[result$sel[1], result$sel[2]] <-
          selection.matrix[result$sel[1], result$sel[2]] + 1
        correct <- any(apply(
          true.ij,
          1,
          function(ij) all(result$sel == ij)
        ))
        if (correct) selection.count <- selection.count + 1
      }

    }

    ## Selection: % of trials correctly selecting an OBDC.
    ## ET: % of trials with no recommended dose = 1 - sum(selection.matrix) / ntrial.
    ## Patient: % of patients allocated to a true OBDC, pooled over replicates.
    ## Sample: mean enrollment. Duration: mean trial duration (days).
    row <- out$Method == method
    out[row, "Selection"] <- round(100 * selection.count / ntrial, 1)
    out[row, "ET"] <- round(100 * (1 - sum(selection.matrix) / ntrial), 1)
    total.patient.count <- sum(patient.matrix)
    true.index <- cbind(true.ij[, 1], true.ij[, 2])
    true.patient.count <- sum(patient.matrix[true.index])
    out[row, "Patient"] <- if (total.patient.count > 0) {
      round(100 * true.patient.count / total.patient.count, 1)
    } else {
      0
    }
    out[row, "Sample"] <- round(sample.sum / ntrial, 1)
    out[row, "Duration"] <- round(duration.sum / ntrial, 1)

    dose.details[[method]] <- list(
      average_patients = patient.matrix / ntrial,
      allocation_percent = if (total.patient.count > 0) {
        100 * patient.matrix / total.patient.count
      } else {
        matrix(0, nrow = nrow(pT), ncol = ncol(pT))
      },
      average_toxicity = tox.matrix / ntrial,
      average_efficacy = eff.matrix / ntrial,
      selection_percent = 100 * selection.matrix / ntrial
    )
  }

  # attached as attributes so main.R can report them without adding data columns
  attr(out, "true.obdcs") <- true.ij
  attr(out, "delayed") <- delayed
  attr(out, "delay.control") <- delay.control
  attr(out, "dose.details") <- dose.details
  out
}
