# Comb-STEIN

Reference implementation of **Comb-STEIN**, an interval design for phase I/II
dose-finding trials that combine a therapeutic cancer vaccine with a concomitant
anticancer agent. The design extends the one-dimensional STEIN framework to a
two-dimensional dose-combination grid and identifies the **optimal biological
dose combination (OBDC)** under a utility-based toxicity–efficacy trade-off.

The code reproduces every simulation result in the manuscript
*Optimizing Dose Combinations of Therapeutic Cancer Vaccines: A Novel Adaptive
Interval Design for Early-Phase Clinical Trials*.

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

## License

MIT — see [LICENSE](LICENSE).
