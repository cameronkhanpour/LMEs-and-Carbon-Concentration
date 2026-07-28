using OperatingEmissionsCertificates
using CSV
using DataFrames
using JSON3
using LinearAlgebra
using Statistics

# The confidence radius grows like R*sqrt(r*log(1/lambda)) through the log
# determinant term and like sqrt(lambda)*theta_bar through the declared norm
# bound, so the reported threshold has an interior optimum in lambda. This
# script locates it for the certificate case rather than assuming it.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const AGGREGATE_DIR = joinpath(ROOT, "results", "aggregate")
const LEARNING_FAILURE_PROBABILITY = 0.005
const CERTIFICATE_OBSERVATIONS = 20
const RIDGE_SWEEP = (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1.0)

operating_points = CSV.read(
    joinpath(AGGREGATE_DIR, "operating_points.csv"),
    DataFrame,
)
point = only(
    operating_points[
        operating_points.case_name .== "pglib_opf_case57_ieee",
        :,
    ],
)
tail_case = analyze_case(
    joinpath(ROOT, "data", "meshed", point.case_name * ".m");
    load_scale=point.load_scale,
    line_limit_scale=point.line_limit_scale,
    cost_ridge=1e-4,
    uncertainty_fraction=0.01,
    exit_probability_samples=0,
    finite_difference_count=0,
    seed=2026,
)

rows = NamedTuple[]
for ridge in RIDGE_SWEEP
    steps = run_sequential_perturbation_experiment(
        tail_case.features.Phi,
        tail_case.theta,
        tail_case.Sigma;
        policies=("variance_informed",),
        seeds=0:99,
        observation_budget=CERTIFICATE_OBSERVATIONS,
        observation_noise_std=0.02,
        ridge_parameter=ridge,
        estimation_failure_probability=LEARNING_FAILURE_PROBABILITY,
    )
    final = steps[
        steps.observation_index .== CERTIFICATE_OBSERVATIONS,
        :,
    ]
    push!(rows, (;
        case_name=tail_case.case_name,
        ridge,
        feature_rank=tail_case.features.feature_rank,
        parameter_norm=norm(tail_case.theta),
        declared_parameter_bound=5sqrt(tail_case.solved.network.n),
        median_variance_bound_ratio=median(
            final.emissions_variance_upper_bound ./
            final.true_emissions_variance,
        ),
        median_threshold_ratio=median(
            sqrt.(
                final.emissions_variance_upper_bound ./
                final.true_emissions_variance,
            ),
        ),
    ))
end
ridge_optimum = DataFrame(rows)
CSV.write(joinpath(AGGREGATE_DIR, "ridge_optimum.csv"), ridge_optimum)
display(ridge_optimum)

# Record this grid alongside the pipeline manifest so the reported optimum is
# reproducible from the recorded provenance alone.
manifest_path = joinpath(ROOT, "results", "manifest.json")
if isfile(manifest_path)
    parsed = JSON3.read(read(manifest_path, String))
    manifest = Dict{String,Any}(String(k) => v for (k, v) in pairs(parsed))
    manifest["ridge_optimum_sweep"] = join(RIDGE_SWEEP, ",")
    manifest["ridge_optimum_case"] = tail_case.case_name
    open(manifest_path, "w") do io
        JSON3.pretty(io, manifest)
    end
end

println("\nmarginal emissions on ", tail_case.case_name)
println("  buses          ", tail_case.solved.network.n)
println("  total demand   ", sum(tail_case.solved.demand) * tail_case.solved.base_mva, " MW")
println("  nominal emissions ", tail_case.nominal_emissions, " tonnes CO2/h")
println("  lme min/median/max ",
    minimum(tail_case.ell), " / ", median(tail_case.ell), " / ",
    maximum(tail_case.ell), " tonnes CO2/MWh")
println("  theta norm     ", norm(tail_case.theta))
println(
    "  emissions change std dev   ",
    sqrt(tail_case.emissions_variance),
    " tonnes CO2/h",
)
