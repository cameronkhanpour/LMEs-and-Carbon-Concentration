using OperatingEmissionsCertificates
using CSV
using DataFrames
using JSON3
using LinearAlgebra
using Statistics

# Two questions this script answers, both of which need only the critical
# region and the uncertainty proxy rather than repeated dispatch solves.
#
# 1. Is the usable forecast error range a small system phenomenon? For every
#    retained system we locate the largest uncertainty scale at which the exit
#    bound still leaves room in the failure budget.
# 2. Does the correlation finding generalize? For the nearest facet of every
#    system we record the aggregate response ratio that predicts, exactly, how
#    its standardized distance moves when the demand errors share a common
#    component.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const AGGREGATE_DIR = joinpath(ROOT, "results", "aggregate")
const TOTAL_FAILURE_PROBABILITY = 0.05
const LEARNING_FAILURE_PROBABILITY = 0.005
const REMAINING_BUDGET =
    TOTAL_FAILURE_PROBABILITY - LEARNING_FAILURE_PROBABILITY
const REPORTED_CORRELATION = 0.5

short_label(name) = replace(
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

operating_points = CSV.read(
    joinpath(AGGREGATE_DIR, "operating_points.csv"),
    DataFrame,
)
case_summary = CSV.read(joinpath(AGGREGATE_DIR, "case_summary.csv"), DataFrame)
retained = Set(case_summary.case_name[case_summary.included])

rows = NamedTuple[]
for point in eachrow(operating_points)
    point.case_name in retained || continue
    result = analyze_case(
        joinpath(ROOT, "data", "meshed", point.case_name * ".m");
        load_scale=point.load_scale,
        line_limit_scale=point.line_limit_scale,
        cost_ridge=1e-4,
        uncertainty_fraction=0.01,
        exit_probability_samples=0,
        finite_difference_count=0,
        seed=2026,
    )

    threshold_scale = certificate_availability_scale(
        result.region.F,
        result.region.offsets,
        result.Sigma;
        budget=REMAINING_BUDGET,
    )

    # How much of the union bound comes from the single nearest facet. A value
    # near one means the bound is a single tail probability and is nearly
    # tight, while a small value means many facets sit at comparable distance
    # and the union double counts them.
    exit_bound, contributions, _ = gaussian_critical_region_exit_bound(
        result.region.F,
        result.region.offsets,
        result.Sigma,
    )
    total_contribution = sum(contributions)
    facet_concentration = total_contribution > 0 ?
        maximum(contributions) / total_contribution : NaN

    ratios = facet_aggregate_ratio(result.region.F, result.standard_deviation)
    closest = argmin(result.kappas)
    closest_label = result.region.labels[closest]
    aggregate_ratio = ratios[closest]
    # The exact scaling law for the standardized distance of a single facet.
    predicted_kappa = result.kappas[closest] /
        sqrt(1 - REPORTED_CORRELATION +
             REPORTED_CORRELATION * aggregate_ratio)

    correlated_sigma, _ = build_correlated_uncertainty_proxy(
        result.solved.demand,
        result.solved.base_mva;
        fraction=0.01,
        correlation=REPORTED_CORRELATION,
    )
    _, _, correlated_kappas = gaussian_critical_region_exit_bound(
        result.region.F,
        result.region.offsets,
        correlated_sigma,
    )
    correlated_threshold_scale = certificate_availability_scale(
        result.region.F,
        result.region.offsets,
        correlated_sigma;
        budget=REMAINING_BUDGET,
    )

    push!(rows, (;
        case_name=result.case_name,
        case_short=short_label(result.case_name),
        buses=result.solved.network.n,
        facets=length(result.region.offsets),
        closest_facet=closest_label,
        facet_kind=startswith(closest_label, "primal_gen") ||
                   startswith(closest_label, "dual_gen") ?
                   "generator" : "line",
        kappa_independent=result.kappas[closest],
        exit_bound,
        facet_concentration,
        aggregate_ratio,
        kappa_correlated=correlated_kappas[closest],
        kappa_correlated_predicted=predicted_kappa,
        availability_percent_independent=100 * 0.01 * threshold_scale,
        availability_percent_correlated=100 * 0.01 * correlated_threshold_scale,
    ))
    GC.gc()
end

geometry = DataFrame(rows)
sort!(geometry, :buses)
CSV.write(joinpath(AGGREGATE_DIR, "region_geometry.csv"), geometry)
display(geometry)

println()
println("aggregate ratio above one predicts a loss of margin under correlation")
for row in eachrow(geometry)
    direction = row.kappa_correlated < row.kappa_independent ? "worse" : "better"
    println(
        "  ", rpad(row.case_short, 6),
        " facet=", rpad(row.facet_kind, 10),
        " a=", round(row.aggregate_ratio; digits=2),
        " kappa ", round(row.kappa_independent; digits=3),
        " -> ", round(row.kappa_correlated; digits=3),
        " (", direction, ")",
    )
end

manifest_path = joinpath(ROOT, "results", "manifest.json")
if isfile(manifest_path)
    parsed = JSON3.read(read(manifest_path, String))
    manifest = Dict{String,Any}(String(k) => v for (k, v) in pairs(parsed))
    manifest["region_geometry_budget"] = REMAINING_BUDGET
    manifest["region_geometry_correlation"] = REPORTED_CORRELATION
    open(manifest_path, "w") do io
        JSON3.pretty(io, manifest)
    end
end
