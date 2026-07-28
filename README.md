# Operating Emissions Certificates

This repository contains a Julia implementation for learning locational
marginal emissions from selected demand perturbations and for bounding the
upper tail of operating emissions under demand forecast error.

The implementation separates three quantities that are easy to conflate:

1. A paired perturbation observation is one central finite difference at a
   selected bus. It normally requires two dispatch simulations.
2. The emissions variance proxy is
   `v = marginal_emissions' * demand_error_proxy * marginal_emissions`.
   It controls the local one sided emissions threshold.
3. The critical region exit bound accounts for demand errors that change the
   active dispatch constraints. A certificate is available only while this
   term leaves room in the requested failure probability.

Within a fixed active set, the marginal emissions vector lies in the span of a
uniform vector and the PTDF rows of the binding lines. The dimension of this
congestion basis determines the number of selected perturbations required for
exact recovery. A variance informed bus selection rule can instead bound the
emissions variance before the complete vector is identified.

The numerical pipeline uses a transparent strictly convex DC optimal power
flow as a controlled evaluator. PowerIO parses the MATPOWER cases, and
PowerDiff supplies the dispatch solves and an independent automatic
sensitivity used only for validation. The recovery procedure itself estimates
marginal emissions from dispatch evaluations at perturbed demands.

All internal power quantities are per unit. Reported demand errors are in MW,
total operating emissions are in metric tonnes of CO2 per hour, and marginal
emissions are in metric tonnes of CO2 per MWh. The generator factors represent
direct operating emissions rather than life cycle greenhouse gas emissions.

## Installation

Install Julia 1.10 or later, then instantiate the pinned environment:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Validation

Run the automated tests:

```powershell
julia --project=. -e 'using Pkg; Pkg.test()'
```

For a short end to end check on the 14 bus system:

```powershell
julia --project=. scripts/smoke.jl
```

## Reproduce the numerical outputs

Run the scripts in this order:

```powershell
julia --project=. scripts/find_operating_points.jl
julia --project=. scripts/run_experiments.jl
julia --project=. scripts/ridge_sweep.jl
julia --project=. scripts/region_geometry.jl
julia --project=. scripts/prepare_plot_data.jl
```

`find_operating_points.jl` searches a fixed grid of load and line limit scales.
`run_experiments.jl` performs the main recovery, confidence bound, active set,
and redispatch evaluations. `ridge_sweep.jl` evaluates the confidence sequence
regularization. `region_geometry.jl` evaluates certificate availability and
correlated forecast errors across systems. `prepare_plot_data.jl` creates
compact machine readable subsets for visualization.

The generated environment and experiment settings are recorded in
`results/manifest.json`. See `results/README.md` for output definitions and
`data/PROVENANCE.md` for network and emission factor provenance.

## Numerical scope

The configured analysis screens 14 PGLib-OPF cases ranging from 14 to 1,354
buses. The primary recovery case is the 300 bus system. The redispatch tail
evaluation uses the 57 bus system, and the critical region exit curve uses the
118 bus system. Case selection, random seeds, uncertainty scales, sample
counts, and regularization parameters are declared in the scripts and recorded
in the result manifest.

## Repository map

- `src/OperatingEmissionsCertificates.jl` contains the implementation.
- `test/runtests.jl` contains algebraic and integration tests.
- `scripts/` contains reproducible analysis entry points.
- `data/meshed/` contains the benchmark networks used by the scripts.
- `results/` contains machine readable numerical outputs and their manifest.
