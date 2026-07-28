using CSV
using DataFrames
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const AGGREGATE_DIR = joinpath(ROOT, "results", "aggregate")

short_case_label(name) = replace(
    name,
    "pglib_opf_case" => "",
    "_ieee_rts" => "",
    "_ieee_dtc" => "",
    "_ieee" => "",
    "_pegase" => "",
    "_pserc" => "",
    "_sdet" => "",
    "_goc" => "",
    "_snem" => "",
)

# Write one compact trajectory per bus selection policy.
experiment_efficiency = CSV.read(
    joinpath(AGGREGATE_DIR, "experiment_efficiency.csv"),
    DataFrame,
)
for policy in unique(experiment_efficiency.policy)
    CSV.write(
        joinpath(AGGREGATE_DIR, "experiment_" * policy * ".csv"),
        experiment_efficiency[experiment_efficiency.policy .== policy, :],
    )
end

# Separate the physical and generic graph basis diagnostics.
basis_ablation = CSV.read(
    joinpath(AGGREGATE_DIR, "basis_ablation.csv"),
    DataFrame,
)
for basis_name in unique(basis_ablation.basis)
    slug = replace(lowercase(basis_name), " " => "_")
    CSV.write(
        joinpath(AGGREGATE_DIR, "basis_" * slug * ".csv"),
        basis_ablation[basis_ablation.basis .== basis_name, :],
    )
end

# Add short labels while preserving the complete case names in the source
# summary.
case_summary = CSV.read(
    joinpath(AGGREGATE_DIR, "case_summary.csv"),
    DataFrame,
)
included_cases = case_summary[case_summary.included, :]
included_cases.case_short = short_case_label.(included_cases.case_name)
CSV.write(
    joinpath(AGGREGATE_DIR, "scaling_plot.csv"),
    included_cases,
)

# Logarithmic plots cannot display an empirical probability of exactly zero.
# Half an event is used only as a display floor; the original observations
# remain available in critical_region_exit_curve.csv.
exit_curve = CSV.read(
    joinpath(AGGREGATE_DIR, "critical_region_exit_curve.csv"),
    DataFrame,
)
exit_curve.plot_empirical_probability = max.(
    exit_curve.empirical_exit_probability,
    0.5 ./ exit_curve.samples,
)
exit_curve.plot_probability_bound = max.(
    exit_curve.exit_probability_bound,
    0.5 ./ exit_curve.samples,
)
CSV.write(
    joinpath(AGGREGATE_DIR, "critical_region_exit_curve_plot.csv"),
    exit_curve,
)

poisson_source = CSV.read(
    joinpath(AGGREGATE_DIR, "poisson_source.csv"),
    DataFrame,
)
CSV.write(
    joinpath(AGGREGATE_DIR, "poisson_binding_endpoints.csv"),
    poisson_source[poisson_source.binding_endpoint, :],
)
CSV.write(
    joinpath(AGGREGATE_DIR, "poisson_other_buses.csv"),
    poisson_source[.!poisson_source.binding_endpoint, :],
)

correlated_tail = CSV.read(
    joinpath(AGGREGATE_DIR, "correlated_tail_curve.csv"),
    DataFrame,
)
for correlation in unique(correlated_tail.demand_correlation)
    slug = replace(@sprintf("%.1f", correlation), "." => "p")
    CSV.write(
        joinpath(AGGREGATE_DIR, "correlated_tail_" * slug * ".csv"),
        correlated_tail[correlated_tail.demand_correlation .== correlation, :],
    )
end
