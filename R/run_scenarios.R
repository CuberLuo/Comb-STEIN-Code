## Repeated-trial simulation and summary aggregation.

run.scenario <- function(
  name,
  ntrial = 1000,
  obdc.definition = c("utility", "table1")
) {
  obdc.definition <- match.arg(obdc.definition)
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
    Selection = NA,
    Patient = NA,
    Sample = NA,
    DLTs = NA,
    Tox = NA,
    Eff = NA
  )

  for (method in methods) {
    selection.count <- 0
    patient.percent.sum <- 0
    sample.sum <- 0
    dlt.sum <- 0
    toxicity.stop.count <- 0
    efficacy.stop.count <- 0

    for (trial in seq_len(ntrial)) {
      result <- switch(
        method,
        "Comb-STEIN-L" = simulate.trial.stein(
          pT, pE, "CombSTEIN-L"
        ),
        "Comb-STEIN-H" = simulate.trial.stein(
          pT, pE, "CombSTEIN-H"
        ),
        "Runin-STEIN" = simulate.trial.stein(
          pT, pE, "RunIn-STEIN"
        ),
        "Waterfall-L" = simulate.trial.waterfall(
          pT, pE, "Waterfall-L"
        ),
        "Waterfall-H" = simulate.trial.waterfall(
          pT, pE, "Waterfall-H"
        )
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
  }

  attr(out, "true.obdcs") <- true.ij
  attr(out, "obdc.definition") <- obdc.definition
  out
}
