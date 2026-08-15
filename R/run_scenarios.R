## Repeated-trial simulation and summary aggregation.

## Simulate every design method repeatedly under one named scenario.
## Returns a method-level operating-characteristic table and stores the true
## OBDC definition and delay settings as attributes.
run.scenario <- function(
  name,
  ntrial = DEFAULT_NTRIAL,
  obdc.definition = c("utility", "table1"),
  delayed = FALSE,
  delay.control = NULL
) {
  obdc.definition <- match.arg(obdc.definition)
  if (length(delayed) != 1L || is.na(delayed)) {
    stop("delayed must be TRUE or FALSE")
  }
  delayed <- isTRUE(delayed)
  if (delayed) delay.control <- normalize.delay.control(delay.control)
  if (!name %in% names(scenarios)) stop("Unknown scenario: ", name)
  if (
    length(ntrial) != 1L ||
      !is.finite(ntrial) ||
      ntrial < 1 ||
      ntrial != as.integer(ntrial)
  ) {
    stop("ntrial must be a positive integer")
  }

  scenario <- build.pTpE(name)
  pT <- scenario$pT
  pE <- scenario$pE
  true.ij <- if (obdc.definition == "utility") {
    get.true.obdc(pT, pE)
  } else {
    table1.obdcs[[name]]
  }

  methods <- c(
    "Comb-STEIN-L",
    "Comb-STEIN-H",
    "Runin-STEIN",
    "Waterfall-L",
    "Waterfall-H"
  )
  out <- data.frame(
    Scenario = name,
    Method = methods,
    Delay = if (delayed) "TITE" else "Complete",
    Selection = NA,
    Patient = NA,
    Sample = NA,
    DLTs = NA,
    Tox = NA,
    Eff = NA,
    Duration = NA
  )

  for (method in methods) {
    selection.count <- 0
    patient.percent.sum <- 0
    sample.sum <- 0
    dlt.sum <- 0
    toxicity.stop.count <- 0
    efficacy.stop.count <- 0
    duration.sum <- 0

    for (trial in seq_len(ntrial)) {
      result <- tryCatch(
        switch(
          method,
          "Comb-STEIN-L" = simulate.trial.stein(
            pT, pE, "CombSTEIN-L", delayed, delay.control
          ),
          "Comb-STEIN-H" = simulate.trial.stein(
            pT, pE, "CombSTEIN-H", delayed, delay.control
          ),
          "Runin-STEIN" = simulate.trial.stein(
            pT, pE, "RunIn-STEIN", delayed, delay.control
          ),
          "Waterfall-L" = simulate.trial.waterfall(
            pT, pE, "Waterfall-L", delayed, delay.control
          ),
          "Waterfall-H" = simulate.trial.waterfall(
            pT, pE, "Waterfall-H", delayed, delay.control
          )
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

      total.n <- sum(result$n)
      sample.sum <- sample.sum + total.n
      dlt.sum <- dlt.sum + sum(result$x)
      if (result$earlystop_tox) {
        toxicity.stop.count <- toxicity.stop.count + 1
      }
      if (result$earlystop_eff) {
        efficacy.stop.count <- efficacy.stop.count + 1
      }
      if (delayed) duration.sum <- duration.sum + result$duration

      if (!any(is.na(result$sel))) {
        correct <- any(apply(
          true.ij,
          1,
          function(ij) all(result$sel == ij)
        ))
        if (correct) selection.count <- selection.count + 1
      }
      if (total.n > 0) {
        true.index <- cbind(true.ij[, 1], true.ij[, 2])
        n.at.true <- sum(result$n[true.index])
        patient.percent.sum <-
          patient.percent.sum + 100 * n.at.true / total.n
      }
    }

    row <- out$Method == method
    out[row, "Selection"] <- round(100 * selection.count / ntrial, 1)
    out[row, "Patient"] <- round(patient.percent.sum / ntrial, 1)
    out[row, "Sample"] <- round(sample.sum / ntrial, 1)
    out[row, "DLTs"] <- round(dlt.sum / ntrial, 1)
    out[row, "Tox"] <- round(100 * toxicity.stop.count / ntrial, 1)
    out[row, "Eff"] <- round(100 * efficacy.stop.count / ntrial, 1)
    if (delayed) {
      out[row, "Duration"] <- round(duration.sum / ntrial, 1)
    }
  }

  attr(out, "true.obdcs") <- true.ij
  attr(out, "obdc.definition") <- obdc.definition
  attr(out, "delayed") <- delayed
  if (delayed) attr(out, "delay.control") <- delay.control
  out
}
