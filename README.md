# Comb-STEIN

Reference implementation of **Comb-STEIN**, an interval design for phase I/II
dose-finding trials that combine a therapeutic cancer vaccine with a concomitant
anticancer agent. The design extends the one-dimensional STEIN framework to a
two-dimensional dose-combination grid and identifies the **optimal biological
dose combination (OBDC)** under a utility-based toxicity–efficacy trade-off.

The code reproduces every simulation result in the manuscript
*Optimizing Dose Combinations of Therapeutic Cancer Vaccines: A Novel Adaptive
Interval Design for Early-Phase Clinical Trials*.

## What the design does

Three allocation strategies share the same monitoring, escalation and final
selection rules; they differ only in where the exploration starts.

| Strategy | Starting point |
|---|---|
| `Comb-STEIN-L` | lowest dose combination of the grid |
| `Comb-STEIN-H` | highest dose combination of the grid |
| `Runin-STEIN` | a four-cohort run-in stage, then the same optimisation stage |

All three use

* Beta-binomial toxicity and efficacy monitoring with STEIN interval rules,
* local (one-step) movement on the two-dimensional grid,
* a utility function `U = q - w1*p - w2*p*I(p > phi)` for the final selection,
  with `w1 = 0.33`, `w2 = 1.09` and `phi = calc.phi(0.3, 0.5) ≈ 0.3971`.

A combination is an OBDC when its expected utility is within `0.03` of the
largest value in the grid, so a scenario can contain more than one OBDC.

A time-to-event (TITE) extension lets decisions use outcomes that are not yet
fully observed through follow-up-time-weighted effective sample sizes, with an
optional enrolment suspension rule.

## Repository layout

```
Comb-STEIN-Code/
  README.md                this file
  config.R                 design parameters, defaults, decision boundaries
  main.R                   command-line entry point
  scenarios.R              12 simulation scenarios (S1-S12), 3x4 / 3x3 / 2x3 grids
  scenarios_case_study.R   6 scenarios (C1-C6) built from three real trials
  isotonic_estimators.R    isotonic regression of the toxicity / efficacy surfaces
  dose_rules.R             escalation, stopping and final selection rules
  simulate_stein.R         trial simulator (complete outcomes)
  simulate_runin.R         run-in stage simulator
  time_to_event.R          TITE extension
  run_scenarios.R          batch runner / parallel driver
  plot/                    figure scripts (300 dpi PNG)
output/
  run_<timestamp>/         one directory per invocation: CSV summaries + logs
  figures/                 rendered figures
```

## Requirements

* R >= 4.2
* `parallel` (base), plus the packages used by the plotting scripts
  (`ggplot2`, `patchwork`)

## Running the simulations

```bash
# all 12 scenarios, 1000 trials each (the manuscript setting)
Rscript Comb-STEIN-Code/main.R

# selected scenarios with a smaller number of trials
Rscript Comb-STEIN-Code/main.R S1 S2 S3 50

# time-to-event (delayed outcome) version
Rscript Comb-STEIN-Code/main.R S1 50 --tite --accrual=1.5 --tox-window=30 --eff-window=60 --suspension=0.5

# control parallel workers, silence debug output
Rscript Comb-STEIN-Code/main.R S1 50 --workers=4 --quiet
```

Every invocation creates a new timestamped directory under `output/` containing
the per-scenario CSV summaries and a run log, so runs never overwrite one another.

The manuscript used these settings: assessment windows `AW = (1, 2)` months,
accrual 1.5 patients per month, suspension proportion `r = 0.5`, cohort size 3,
1000 simulated trials per scenario and method, random seed 2026.

## Case studies

`Comb-STEIN-Code/scenarios_case_study.R` builds six scenarios from three registered trials:

| Scenarios | Agents | Grid | Registration |
|---|---|---|---|
| C1-C2 | GM-CSF + GP2 | 2x3 | NCT03014076 |
| C3-C4 | fruquintinib + IDOV-SAFE | 2x3 | NCT06380309 |
| C5-C6 | ODI-2001 + ipilimumab | 3x3 | NCT07587827 |

Dose levels come from the trial registrations. The cell-wise toxicity and
efficacy probabilities are **scenarios constructed from the mechanisms of
action** of the agents, not results reported by the trials.

## Figures

```bash
Rscript Comb-STEIN-Code/plot/plot_delay_comparison.R    # main text Figures 3 and 4
Rscript Comb-STEIN-Code/plot/plot_sensitivity.R         # supplementary Figures S1-S3
Rscript Comb-STEIN-Code/plot/plot_trial_trajectory.R    # supplementary Figure S4
```

## Reproducibility

The random seed is fixed in `Comb-STEIN-Code/config.R` (`DEFAULT_SEED <- 2026`) and every
run records the seed and the full configuration in its output directory.

## License

MIT — see [LICENSE](LICENSE).
