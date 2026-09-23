## Real-trial case-study scenarios (paper Tables 3-4).
## pT/pE are constructed, not trial-reported: plausible surfaces built from each
## agent's mechanism plus the registered (or, where flagged, assumed) dose levels.
## pT[i, j], pE[i, j] = true toxicity / efficacy for drug A level i x drug B level j.
## OBDC is never hard-coded: case.study.obdc() recomputes it with the same rule
## run.scenario() uses for selection percentages. Axis order matters (row-first
## escalation, row-priority ties), so a transposed grid is a different design.
## Run: Rscript Comb-STEIN-Code/main.R C1 C2 C3 C4 C5 C6 1000

# row-wise fill; duplicates the helper of the same name in scenarios.R
mk <- function(values, nrow) {
  matrix(values, nrow = nrow, byrow = TRUE)
}

case.study.scenarios <- list(
  ## Case 1 (NCT03014076): GM-CSF 125/250 ug x GP2 100/500/1000 ug, 2x3.
  ## All four levels are the published 3+3 cohorts (Clifton et al. 2017), and C1/C2
  ## are the only scenarios backed by reported outcomes - hence the very low
  ## toxicity. GP2 drives the CD8+ T-cell response, GM-CSF is only the adjuvant, so
  ## efficacy rises mainly along GP2. True OBDC (2, 3) in both scenarios.
  # low-toxicity variant of the pair
  C1 = list(
    label = "Case study 1 (NCT03014076, GM-CSF + GP2): low toxicity",
    pT = mk(c(0.01, 0.03, 0.05,
              0.05, 0.07, 0.10), nrow = 2),
    pE = mk(c(0.10, 0.15, 0.30,
              0.12, 0.25, 0.60), nrow = 2)
  ),
  # intermediate variant: toxicity rises faster along GP2, higher efficacy ceiling
  C2 = list(
    label = "Case study 2 (NCT03014076, GM-CSF + GP2): intermediate",
    pT = mk(c(0.01, 0.05, 0.15,
              0.10, 0.15, 0.20), nrow = 2),
    pE = mk(c(0.05, 0.15, 0.35,
              0.15, 0.25, 0.75), nrow = 2)
  ),
  ## Case 2 (NCT06380309): fruquintinib 3/5 mg x IDOV-SAFE 1e9/3e9/1e10 PFU, 2x3.
  ## Viral levels (1e9/3e9/1e10/3e10 PFU) are the registered escalation; 3e10 is
  ## excluded as the toxicity ceiling. Fruquintinib 3/5 mg are label/routine
  ## clinical doses - the registry names the combination but no levels. Toripalimab
  ## is fixed at 240 mg, so it is not a dimension. Trial still recruiting with no
  ## public dose-toxicity/efficacy data: pT/pE are mechanism-based, not reported.
  ## IV viral load gives flu-like symptoms and transaminase rise, VEGFR inhibition
  ## gives hypertension and hand-foot syndrome, so toxicity climbs on both axes;
  ## efficacy climbs with viral load (intratumoral infection, immune activation)
  ## and with VEGFR inhibition (microenvironment remodelling).
  ## True OBDC (1, 3) in both scenarios: lowest TKI + highest viral load, i.e. not
  ## a grid corner - it keeps efficacy while backing off the most toxic TKI level.
  # C3: the top-top cell crosses the 0.397 toxicity reference, which pushes the
  # optimum back to (1, 3); utility advantage ~0.09
  C3 = list(
    label = "Case study 3 (NCT06380309, fruquintinib + IDOV-SAFE): 2x3 grid",
    pT = mk(c(0.05, 0.08, 0.13,
              0.15, 0.22, 0.45), nrow = 2),
    pE = mk(c(0.15, 0.22, 0.38,
              0.22, 0.32, 0.55), nrow = 2)
  ),
  # C4: same OBDC (1, 3) by a different route - no cell crosses 0.397, but
  # fruquintinib efficacy plateaus at 3 mg (5 mg adds toxicity only, pE rows match
  # at the top viral level), so the utility advantage shrinks to ~0.04 vs C3's
  # ~0.09 and the design has a harder call
  C4 = list(
    label = "Case study 4 (NCT06380309, fruquintinib + IDOV-SAFE): efficacy plateau",
    pT = mk(c(0.05, 0.08, 0.12,
              0.09, 0.15, 0.24), nrow = 2),
    pE = mk(c(0.16, 0.26, 0.38,
              0.22, 0.35, 0.38), nrow = 2)
  ),
  ## Case 3 (NCT07587827): ODI-2001 1/2/4 mg x ipilimumab 2.5/5/10 mg, 3x3
  ## (already square, so not transposed; drug A = vaccine, drug B = ipilimumab).
  ## Ipilimumab levels and the 4 mg DNA vaccine dose are from the registry; 1 and
  ## 2 mg vaccine are hypothetical lower levels probed for this platform. MVA is
  ## fixed at 1e7 pfu, and the registry defines DLT/MTD on the ipilimumab axis
  ## only. Trial not yet started: pT/pE are mechanism-based, not reported.
  ## Toxicity is driven by ipilimumab (anti-CTLA4 irAEs rise steeply with dose,
  ## crossing the 0.397 reference at high ipi above vaccine 2 mg) and rises only
  ## mildly along vaccine; efficacy rises on both axes because CTLA4 blockade
  ## amplifies the neoantigen-specific T-cell response.
  ## True OBDC (3, 1) in both scenarios: max vaccine, min ipilimumab. Each ipi step
  ## buys ~+0.03 efficacy but costs 0.33 x dp (plus 1.09 x p overshoot penalty at
  ## the high levels), so net utility falls. (3, 1) is unreachable by 1-D ipi
  ## escalation, which would stop at 10 mg - the opposite of the true optimum.
  # C5: monotone efficacy on both axes; utility advantage ~0.11
  C5 = list(
    label = "Case study 5 (NCT07587827, ODI-2001 + ipilimumab): 3x3 grid",
    pT = mk(c(0.06, 0.12, 0.22,
              0.09, 0.20, 0.45,
              0.13, 0.42, 0.60), nrow = 3),
    pE = mk(c(0.12, 0.15, 0.18,
              0.20, 0.24, 0.27,
              0.33, 0.38, 0.41), nrow = 3)
  ),
  # C6: same OBDC (3, 1), but vaccine 2 -> 4 mg adds almost no efficacy (0.27 ->
  # 0.32 at the lowest ipi level; neoantigen delivery near plateau) while toxicity
  # keeps climbing, so the utility advantage falls below 0.04 vs C5's ~0.11 - the
  # pair differs in surface steepness, not in OBDC location
  C6 = list(
    label = "Case study 6 (NCT07587827, ODI-2001 + ipilimumab): flat vaccine gain",
    pT = mk(c(0.07, 0.16, 0.30,
              0.11, 0.26, 0.50,
              0.15, 0.44, 0.62), nrow = 3),
    pE = mk(c(0.14, 0.16, 0.18,
              0.27, 0.29, 0.31,
              0.32, 0.34, 0.36), nrow = 3)
  )
)

# mirrors build.pTpE() in scenarios.R
build.case.study.pTpE <- function(name) {
  key <- as.character(name)
  if (!key %in% names(case.study.scenarios)) {
    stop(sprintf("Unknown case study scenario: %s", key))
  }
  case.study.scenarios[[key]]
}

# opt-in merge: main.R still defaults to S1-S12 from scenarios.R
register.case.study.scenarios <- function(envir = globalenv()) {
  scenarios <- get("scenarios", envir = envir)
  scenarios[names(case.study.scenarios)] <- case.study.scenarios
  assign("scenarios", scenarios, envir = envir)
  assign(
    "nrow_map",
    vapply(scenarios, function(s) nrow(s$pT), integer(1)),
    envir = envir
  )
  assign(
    "ncol_map",
    vapply(scenarios, function(s) ncol(s$pT), integer(1)),
    envir = envir
  )
  invisible(names(scenarios))
}

## True OBDC per case-study scenario, via get.true.utility() -> utility.fn(): every
## combination within OBDC_UTILITY_TOLERANCE of the max utility counts, the same
## rule run.scenario() applies, so no scenario may pre-specify its OBDC.
## name = NULL -> named list for all; a single name -> one 2-column index matrix
## (column 1 = drug A level, column 2 = drug B level).
case.study.obdc <- function(name = NULL, verbose = TRUE) {
  keys <- if (is.null(name)) {
    names(case.study.scenarios)
  } else {
    as.character(name)
  }
  unknown <- setdiff(keys, names(case.study.scenarios))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown case study scenario: %s", paste(unknown, collapse = ", ")))
  }
  obdcs <- lapply(keys, function(key) {
    scenario <- case.study.scenarios[[key]]
    get.true.obdc(scenario$pT, scenario$pE)
  })
  names(obdcs) <- keys
  if (verbose) {
    for (key in keys) {
      found <- obdcs[[key]]
      cat(sprintf(
        "%s (%s): true OBDC = %s\n",
        key, case.study.scenarios[[key]]$label,
        paste(sprintf("(%d, %d)", found[, 1], found[, 2]), collapse = ", ")
      ))
    }
  }
  if (length(keys) == 1L) obdcs[[1]] else obdcs
}
