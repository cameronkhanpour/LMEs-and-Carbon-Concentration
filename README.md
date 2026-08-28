# Congestion Structure and Exceedance Bounds for Locational Marginal Emissions

This repository contains the Julia research code and numerical results for the
paper Congestion Structure and Exceedance Bounds for Locational Marginal
Emissions.

The project studies two questions. First, how many dispatch simulations are
needed to recover all locational marginal emissions (LMEs) when model
derivatives are unavailable? Second, how can the resulting local LME vector be
used under uncertain demand when a change in the active dispatch constraints
may invalidate it?

Within a fixed active set of a strictly convex DC optimal power flow, the LME
vector lies in the span of the uniform vector and the power transfer
distribution factor rows of the binding lines. The dimension `r` of this
congestion basis is at most one more than the number of binding lines. Thus,
`r` suitably selected nodal perturbations can recover an `n` bus LME vector.
The code also combines this representation with a demand forecast error model
and the local critical region to bound operating emissions exceedance.

Across the ten retained benchmark systems, which range from 14 to 1,354 buses,
the congestion rank ranges from 2 to 15. On the 300 bus system, 12 selected
central differences use 24 dispatch simulations to recover 300 LME values,
instead of the 600 simulations required when every bus is perturbed. These are
results at the operating points recorded in this repository, not general
guarantees for systems of the same size.

## Research status and scope

This is a reproducible research artifact, not an operational market tool. The
implementation uses a transparent DC optimal power flow as the dispatch
evaluator. LME recovery uses operating emissions evaluated at perturbed demand
points. Direct solver sensitivities are computed only to validate the recovered
values.

The current analysis assumes:

- a single period, strictly convex DC optimal power flow;
- a regular nominal active set;
- selected demand perturbations that remain in the nominal critical region;
- known binding line PTDF rows for construction of the congestion basis; and
- a declared subgaussian or Gaussian demand forecast error model.

The emission factors represent direct operating CO2 emissions. They do not
include life cycle greenhouse gas emissions. Extensions to AC power flow,
unit commitment, and multiple time periods are outside the current scope.

## Installation

[Install Julia](https://julialang.org/downloads/) 1.10 or later, then clone the
repository and instantiate its pinned environment:

```bash
git clone https://github.com/cameronkhanpour/carbonation.git
cd carbonation
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

The first run downloads solver artifacts and precompiles the environment, so it
will take longer than subsequent runs. The saved results were generated with
Julia 1.12.4. Package versions and the pinned PowerDiff revision are recorded in
[`Manifest.toml`](Manifest.toml).

## Quick start

Run a short end to end analysis of the 14 bus system:

```bash
julia --project=. scripts/smoke.jl
```

The script searches a small, declared set of operating points and reports the
binding line count, congestion rank, derivative validation error, finite
difference error, congestion span residual, and critical region exit
probability. A successful run finds one binding line and congestion rank two;
small floating point differences between platforms are expected.

Run the test suite with:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

## Reproducing the numerical results

The complete pipeline is run in the following order:

```bash
julia --project=. scripts/find_operating_points.jl
julia --project=. scripts/run_experiments.jl
julia --project=. scripts/ridge_sweep.jl
julia --project=. scripts/region_geometry.jl
julia --project=. scripts/prepare_plot_data.jl
```

The full pipeline performs repeated optimization and Monte Carlo experiments
across 14 systems, so it is substantially more expensive than the smoke test.
It rewrites the generated CSV files under `results/`.

- `find_operating_points.jl` searches a fixed grid of load and line limit
  scales for regular operating points with binding line constraints.
- `run_experiments.jl` runs the recovery, validation, active set, sequential
  design, and redispatch experiments.
- `ridge_sweep.jl` evaluates sensitivity to the confidence sequence ridge
  parameter.
- `region_geometry.jl` evaluates bound availability and correlated demand
  errors across the retained systems.
- `prepare_plot_data.jl` creates compact CSV files used by the figures.

The saved configuration, package versions, random seeds, sample counts, and
case choices are recorded in [`results/manifest.json`](results/manifest.json).
Compact summaries are in [`results/aggregate/`](results/aggregate/), while
trial level outputs are in [`results/raw/`](results/raw/).

## Data and units

The repository includes 14 PGLib-OPF MATPOWER cases and direct operating
emission factors derived from PGLib-CO2. Network sources, versions, case
licenses, the inspected emission factor revision, and the inclusion criteria
are documented in [`data/PROVENANCE.md`](data/PROVENANCE.md).

Power system calculations use per unit quantities internally. Reported demand
perturbations and forecast errors are in MW. Total operating emissions are in
metric tonnes of CO2 per hour, and LMEs are in metric tonnes of CO2 per MWh for
the one hour dispatch interval.

The implementation builds on
[PowerIO.jl](https://github.com/eigenergy/PowerIO.jl) for MATPOWER parsing and
[PowerDiff.jl](https://github.com/grid-opt-alg-lab/PowerDiff.jl) for DC dispatch
and independent sensitivity validation.

## Repository layout

- [`src/OperatingEmissionsCertificates.jl`](src/OperatingEmissionsCertificates.jl)
  contains the Julia module and exported analysis routines.
- [`scripts/`](scripts/) contains the smoke test and reproducible experiment
  entry points.
- [`test/runtests.jl`](test/runtests.jl) contains algebraic and integration
  tests.
- [`data/`](data/) contains the benchmark networks and their provenance.
- [`results/`](results/) contains the saved numerical outputs and experiment
  manifest.

## Citation

If you use this code or its numerical outputs, please cite the accompanying
paper. Complete citation metadata will be added when a public preprint is
available.
