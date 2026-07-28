# Numerical Outputs

All files in this directory are machine readable outputs from the Julia
analysis. `manifest.json` records package versions, source revisions, random
seeds, sample counts, and major numerical settings.

## Directory structure

- `aggregate/` contains case summaries, validation diagnostics, recovery
  results, critical region geometry, confidence bound summaries, and
  redispatch tail evaluations.
- `raw/experiment_steps.csv` records every step of every sequential bus
  selection trial.
- `raw/certificate_observation_steps.csv` records the trials used for the
  reported emissions certificate.

One `observation` is a noisy marginal emissions estimate at a selected bus. A
paired central difference requires two dispatch simulations. Columns ending in
`emissions_variance` refer to the proxy
`v = marginal_emissions' * demand_error_proxy * marginal_emissions`, measured
in squared tonnes of CO2 per hour. Threshold columns are measured in tonnes of
CO2 per hour.

## Primary outputs

- `case_summary.csv` records system size, congestion basis rank, exact recovery
  requirements, nominal emissions, active set margin, and validation status.
- `derivative_validation.csv` compares the automatic dispatch sensitivity with
  the independently assembled fixed active set KKT derivative and finite
  differences.
- `marginal_emissions_recovery.csv` records exact field recovery from selected
  buses.
- `experiment_efficiency.csv` compares variance informed, leverage, and uniform
  bus selection.
- `critical_region_exit.csv` and `critical_region_exit_curve.csv` report the
  analytical exit bound and Monte Carlo exit frequency.
- `emissions_tail_curve.csv` compares local and region exit adjusted thresholds
  with redispatch outcomes.
- `region_geometry.csv` reports certificate availability and the effect of
  correlated demand errors across systems.

Files with names ending in `_plot.csv`, per-policy experiment files, per-basis
files, and per-correlation tail files are compact subsets created by
`scripts/prepare_plot_data.jl`. They contain no additional simulation results.
