using OperatingEmissionsCertificates
using CSV
using DataFrames
using Dates
using JSON3
using LinearAlgebra
using PowerDiff
using PowerIO
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const AGGREGATE_DIR = joinpath(ROOT, "results", "aggregate")
const RAW_DIR = joinpath(ROOT, "results", "raw")
const LEARNING_FAILURE_PROBABILITY = 0.005
const TOTAL_FAILURE_PROBABILITY = 0.05
const CERTIFICATE_OBSERVATIONS = 20
const RIDGE_PARAMETER = 1e-8
const RIDGE_SWEEP = (1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1.0)
const DEMAND_CORRELATIONS = (0.3, 0.5, 0.8)
const RETAINED_CASES = Set([
    "pglib_opf_case57_ieee",
    "pglib_opf_case118_ieee",
    "pglib_opf_case300_ieee",
])
const MONTE_CARLO_CASES = Set([
    "pglib_opf_case14_ieee",
    "pglib_opf_case30_ieee",
    "pglib_opf_case57_ieee",
    "pglib_opf_case73_ieee_rts",
    "pglib_opf_case118_ieee",
    "pglib_opf_case162_ieee_dtc",
    "pglib_opf_case179_goc",
    "pglib_opf_case197_snem",
    "pglib_opf_case240_pserc",
    "pglib_opf_case300_ieee",
])

mkpath(AGGREGATE_DIR)
mkpath(RAW_DIR)

operating_points = CSV.read(
    joinpath(AGGREGATE_DIR, "operating_points.csv"),
    DataFrame,
)

analyses = Dict{String,Any}()
summary_rows = NamedTuple[]
derivative_rows = NamedTuple[]
critical_region_exit_rows = NamedTuple[]

for row in eachrow(operating_points)
    case_path = joinpath(ROOT, "data", "meshed", row.case_name * ".m")
    monte_carlo_samples =
        row.case_name in MONTE_CARLO_CASES ? 2_000 : 0
    finite_difference_count =
        row.case_name in MONTE_CARLO_CASES ? 6 : 0
    start_time = time()
    result = analyze_case(
        case_path;
        load_scale=row.load_scale,
        line_limit_scale=row.line_limit_scale,
        cost_ridge=1e-4,
        uncertainty_fraction=0.01,
        exit_probability_samples=monte_carlo_samples,
        finite_difference_count,
        seed=2026,
    )
    runtime_seconds = time() - start_time
    row.case_name in RETAINED_CASES && (analyses[row.case_name] = result)

    active_generator_constraints =
        length(result.active.gen_upper) + length(result.active.gen_lower)
    Phi = result.features.Phi
    rank_factorization = qr(transpose(Phi), ColumnNorm())
    recovery_buses =
        sort!(rank_factorization.p[1:result.features.feature_rank])
    recovered_theta = Phi[recovery_buses, :] \ result.ell[recovery_buses]
    noiseless_recovery_error =
        norm(Phi * recovered_theta - result.ell) /
        max(norm(result.ell), eps(Float64))
    included =
        result.nominal_emissions > 1e-8 &&
        result.derivative_error <= 1e-8 &&
        result.span_residual <= 1e-8 &&
        result.poisson.residual <= 1e-8
    exclusion_reason = if result.nominal_emissions <= 1e-8
        "zero operating emissions under direct factors"
    elseif result.derivative_error > 1e-8
        "sensitivity validation failed"
    elseif result.span_residual > 1e-8
        "congestion span validation failed"
    elseif result.poisson.residual > 1e-8
        "graph Poisson validation failed"
    else
        ""
    end

    push!(summary_rows, (;
        case_name=row.case_name,
        included,
        exclusion_reason,
        buses=result.solved.network.n,
        branches=result.solved.network.m,
        generators=result.solved.network.k,
        line_limit_scale=row.line_limit_scale,
        active_lines=length(result.active.line_indices),
        active_generator_constraints,
        feature_rank=result.features.feature_rank,
        selected_buses_for_exact_recovery=result.features.feature_rank,
        buses_in_full_coordinate_audit=result.solved.network.n,
        bus_selection_reduction=
            result.solved.network.n / result.features.feature_rank,
        noiseless_recovery_error,
        nominal_emissions_tph=result.nominal_emissions,
        emissions_change_std_tph=sqrt(result.emissions_variance),
        active_set_margin=result.closest_kappa,
        exit_probability_bound=result.exit_probability_bound,
        empirical_exit_probability=result.empirical_exit_probability,
        kkt_condition=result.fixed.condition_number,
        runtime_seconds,
    ))
    push!(derivative_rows, (;
        case_name=row.case_name,
        automatic_vs_active_set_kkt=result.derivative_error,
        finite_difference_error=result.finite_difference_error,
        congestion_span_residual=result.span_residual,
        poisson_residual=result.poisson.residual,
        poisson_leakage=result.poisson.leakage,
        active_feasibility_residual=result.fixed.local_error,
    ))
    push!(critical_region_exit_rows, (;
        case_name=row.case_name,
        facets=length(result.region.offsets),
        closest_facet=result.closest_label,
        active_set_margin=result.closest_kappa,
        exit_probability_bound=result.exit_probability_bound,
        empirical_exit_probability=result.empirical_exit_probability,
        monte_carlo_samples,
    ))
    row.case_name in RETAINED_CASES || GC.gc()
end

case_summary = DataFrame(summary_rows)
derivative_validation = DataFrame(derivative_rows)
critical_region_exit = DataFrame(critical_region_exit_rows)
CSV.write(joinpath(AGGREGATE_DIR, "case_summary.csv"), case_summary)
CSV.write(
    joinpath(AGGREGATE_DIR, "derivative_validation.csv"),
    derivative_validation,
)
CSV.write(
    joinpath(AGGREGATE_DIR, "critical_region_exit.csv"),
    critical_region_exit,
)

audit_case = analyses["pglib_opf_case300_ieee"]
congestion_sweep_rows = NamedTuple[]
for line_limit_scale in (1.3, 1.2, 1.1, 1.0, 0.9)
    swept_case = analyze_case(
        joinpath(ROOT, "data", "meshed", "pglib_opf_case300_ieee.m");
        load_scale=1.0,
        line_limit_scale,
        cost_ridge=1e-4,
        uncertainty_fraction=0.01,
        exit_probability_samples=0,
        finite_difference_count=0,
        seed=2026,
    )
    push!(congestion_sweep_rows, (;
        case_name=swept_case.case_name,
        line_limit_scale,
        active_lines=length(swept_case.active.line_indices),
        feature_rank=swept_case.features.feature_rank,
        buses=swept_case.solved.network.n,
        reduction=swept_case.solved.network.n /
                  swept_case.features.feature_rank,
    ))
end
congestion_sweep = DataFrame(congestion_sweep_rows)
CSV.write(
    joinpath(AGGREGATE_DIR, "congestion_sweep.csv"),
    congestion_sweep,
)

experiment_steps = run_sequential_perturbation_experiment(
    audit_case.features.Phi,
    audit_case.theta,
    audit_case.Sigma;
    policies=("variance_informed", "leverage", "uniform"),
    seeds=0:99,
    observation_budget=80,
    observation_noise_std=0.02,
    ridge_parameter=RIDGE_PARAMETER,
    estimation_failure_probability=LEARNING_FAILURE_PROBABILITY,
)
experiment_steps.case_name = fill(audit_case.case_name, nrow(experiment_steps))
experiment_steps.perturbation_observations =
    experiment_steps.observation_index
experiment_steps.cached_base_dispatch_runs =
    experiment_steps.observation_index .+ 1
experiment_steps.paired_dispatch_runs =
    2 .* experiment_steps.observation_index
experiment_steps.threshold_inflation_ratio =
    sqrt.(experiment_steps.variance_bound_ratio)
CSV.write(joinpath(RAW_DIR, "experiment_steps.csv"), experiment_steps)

experiment_efficiency = combine(
    groupby(experiment_steps, [:policy, :observation_index]),
    :variance_bound_ratio => median => :median_variance_bound_ratio,
    :variance_bound_ratio => (x -> quantile(x, 0.10)) =>
        :q10_variance_bound_ratio,
    :variance_bound_ratio => (x -> quantile(x, 0.90)) =>
        :q90_variance_bound_ratio,
    :threshold_inflation_ratio => median =>
        :median_threshold_inflation_ratio,
    :threshold_inflation_ratio => (x -> quantile(x, 0.10)) =>
        :q10_threshold_inflation_ratio,
    :threshold_inflation_ratio => (x -> quantile(x, 0.90)) =>
        :q90_threshold_inflation_ratio,
    :marginal_emissions_relative_error => median =>
        :median_marginal_emissions_error,
    :upper_bound_covers_truth => mean => :empirical_bound_coverage,
)
sort!(experiment_efficiency, [:policy, :observation_index])
CSV.write(
    joinpath(AGGREGATE_DIR, "experiment_efficiency.csv"),
    experiment_efficiency,
)

recovery = recover_marginal_emissions(
    audit_case.solved,
    audit_case.features,
    audit_case.ell;
    step_mw=1.0,
)
recovery_summary = DataFrame([(
    case_name=audit_case.case_name,
    buses=audit_case.solved.network.n,
    feature_rank=audit_case.features.feature_rank,
    paired_perturbation_observations=
        recovery.paired_perturbation_observations,
    selected_dispatch_runs=recovery.selected_dispatch_runs,
    full_dispatch_runs=recovery.full_dispatch_runs,
    relative_field_error=recovery.relative_error,
)])
CSV.write(
    joinpath(AGGREGATE_DIR, "marginal_emissions_recovery.csv"),
    recovery_summary,
)
CSV.write(
    joinpath(AGGREGATE_DIR, "selected_bus_observations.csv"),
    DataFrame(
        bus=recovery.selected_buses,
        observed_marginal_emissions=
            recovery.observations[recovery.selected_buses],
        exact_marginal_emissions=
            audit_case.ell[recovery.selected_buses],
    ),
)

step_study = perturbation_step_study(
    audit_case.solved,
    audit_case.features,
    audit_case.ell,
    audit_case.region;
    step_sizes_mw=(0.1, 1.0, 10.0, 50.0, 100.0, 250.0, 500.0),
)
step_study.case_name = fill(audit_case.case_name, nrow(step_study))
CSV.write(joinpath(AGGREGATE_DIR, "perturbation_step.csv"), step_study)

binding_study = binding_set_study(
    audit_case.solved,
    audit_case.ptdf,
    audit_case.active,
    audit_case.ell;
    extra_line_counts=(2, 5),
)
binding_study.case_name = fill(audit_case.case_name, nrow(binding_study))
CSV.write(joinpath(AGGREGATE_DIR, "binding_set_study.csv"), binding_study)

basis_results = basis_ablation(
    audit_case.solved.network,
    audit_case.ell,
    audit_case.Sigma,
    audit_case.features.Phi;
    max_dimension=50,
)
basis_results.case_name = fill(audit_case.case_name, nrow(basis_results))
CSV.write(joinpath(AGGREGATE_DIR, "basis_ablation.csv"), basis_results)

endpoint_mask = falses(audit_case.solved.network.n)
endpoint_mask[audit_case.poisson.endpoints] .= true
poisson_source = DataFrame(
    bus=1:audit_case.solved.network.n,
    source=audit_case.poisson.source,
    absolute_source=abs.(audit_case.poisson.source),
    binding_endpoint=endpoint_mask,
)
CSV.write(joinpath(AGGREGATE_DIR, "poisson_source.csv"), poisson_source)

exit_curve_case = analyses["pglib_opf_case118_ieee"]
exit_curve = critical_region_exit_curve(
    exit_curve_case.region.F,
    exit_curve_case.region.offsets,
    exit_curve_case.Sigma;
    scale_factors=range(0.35, 2.5; length=18),
    samples=50_000,
    seed=2026,
)
exit_curve.case_name = fill(exit_curve_case.case_name, nrow(exit_curve))
CSV.write(
    joinpath(AGGREGATE_DIR, "critical_region_exit_curve.csv"),
    exit_curve,
)

tail_case = analyses["pglib_opf_case57_ieee"]
const TAIL_SCALES = [0.20, 0.30, 0.40, 0.50, 0.70, 1.00]

"""
Attach the threshold that an analyst obtains from selected perturbations alone.

The estimation step spends `LEARNING_FAILURE_PROBABILITY` and the region exit
term spends `exit_bound`, so the emissions tail receives whatever remains of
the total failure probability.
"""
function attach_selected_thresholds!(path, certificate_steps)
    medians, lower, upper = Float64[], Float64[], Float64[]
    for row in eachrow(path)
        remaining =
            TOTAL_FAILURE_PROBABILITY -
            LEARNING_FAILURE_PROBABILITY -
            row.exit_bound
        if remaining > 0
            thresholds = sqrt.(
                2 .* row.uncertainty_scale^2 .*
                certificate_steps.emissions_variance_upper_bound .*
                log(1 / remaining),
            )
            push!(medians, median(thresholds))
            push!(lower, quantile(thresholds, 0.10))
            push!(upper, quantile(thresholds, 0.90))
        else
            push!(medians, NaN)
            push!(lower, NaN)
            push!(upper, NaN)
        end
    end
    path.estimated_threshold_median = medians
    path.estimated_threshold_q10 = lower
    path.estimated_threshold_q90 = upper
    path.learning_failure_probability =
        fill(LEARNING_FAILURE_PROBABILITY, nrow(path))
    path.certificate_observations =
        fill(CERTIFICATE_OBSERVATIONS, nrow(path))
    path.certificate_dispatch_runs =
        fill(2CERTIFICATE_OBSERVATIONS, nrow(path))
    path.load_forecast_std_percent = 100 .* 0.01 .* path.uncertainty_scale
    path.case_name = fill(tail_case.case_name, nrow(path))
    return path
end

certificate_observation_steps =
    run_sequential_perturbation_experiment(
    tail_case.features.Phi,
    tail_case.theta,
    tail_case.Sigma;
    policies=("variance_informed",),
    seeds=0:99,
    observation_budget=CERTIFICATE_OBSERVATIONS,
    observation_noise_std=0.02,
    ridge_parameter=RIDGE_PARAMETER,
    estimation_failure_probability=LEARNING_FAILURE_PROBABILITY,
)
certificate_observation_steps.case_name =
    fill(tail_case.case_name, nrow(certificate_observation_steps))
CSV.write(
    joinpath(RAW_DIR, "certificate_observation_steps.csv"),
    certificate_observation_steps,
)
certificate_steps = certificate_observation_steps[
    certificate_observation_steps.observation_index .==
    CERTIFICATE_OBSERVATIONS,
    :,
]
tail_path = evaluate_redispatch_tail(
    tail_case.solved,
    tail_case.ell,
    tail_case.Sigma,
    tail_case.region.F,
    tail_case.region.offsets;
    scale_factors=TAIL_SCALES,
    samples=5_000,
    failure_probability=TOTAL_FAILURE_PROBABILITY,
    seed=2026,
    solve_all_samples=false,
    inside_validation_samples=100,
    exit_bound_method=:gaussian,
)
attach_selected_thresholds!(tail_path, certificate_steps)
tail_path.demand_correlation = zeros(nrow(tail_path))
CSV.write(joinpath(AGGREGATE_DIR, "emissions_tail_curve.csv"), tail_path)

# A common demand error component leaves every bus standard deviation
# unchanged but can either reduce or increase the active set margin. The sign
# and magnitude are determined by each facet's aggregate response ratio, not
# by a generic assumption that correlated errors are easier or harder.
correlated_paths = DataFrame[]
for correlation in DEMAND_CORRELATIONS
    correlated_sigma, _ = build_correlated_uncertainty_proxy(
        tail_case.solved.demand,
        tail_case.solved.base_mva;
        fraction=0.01,
        correlation,
    )
    correlated_steps = run_sequential_perturbation_experiment(
        tail_case.features.Phi,
        tail_case.theta,
        correlated_sigma;
        policies=("variance_informed",),
        seeds=0:99,
        observation_budget=CERTIFICATE_OBSERVATIONS,
        observation_noise_std=0.02,
        ridge_parameter=RIDGE_PARAMETER,
        estimation_failure_probability=LEARNING_FAILURE_PROBABILITY,
    )
    correlated_path = evaluate_redispatch_tail(
        tail_case.solved,
        tail_case.ell,
        correlated_sigma,
        tail_case.region.F,
        tail_case.region.offsets;
        scale_factors=TAIL_SCALES,
        samples=5_000,
        failure_probability=TOTAL_FAILURE_PROBABILITY,
        seed=2026,
        solve_all_samples=false,
        inside_validation_samples=100,
        exit_bound_method=:gaussian,
    )
    attach_selected_thresholds!(
        correlated_path,
        correlated_steps[
            correlated_steps.observation_index .==
            CERTIFICATE_OBSERVATIONS,
            :,
        ],
    )
    correlated_path.demand_correlation =
        fill(correlation, nrow(correlated_path))
    push!(correlated_paths, correlated_path)
end
correlated_tail = vcat(tail_path, correlated_paths...)
CSV.write(
    joinpath(AGGREGATE_DIR, "correlated_tail_curve.csv"),
    correlated_tail,
)

# The confidence radius grows with log(det V_t / det(lambda I)), so a small
# ridge inflates the reported threshold for reasons unrelated to the data.
ridge_rows = NamedTuple[]
for ridge in RIDGE_SWEEP
    ridge_steps = run_sequential_perturbation_experiment(
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
    final_steps = ridge_steps[
        ridge_steps.observation_index .== CERTIFICATE_OBSERVATIONS,
        :,
    ]
    for row in eachrow(tail_path)
        remaining =
            TOTAL_FAILURE_PROBABILITY -
            LEARNING_FAILURE_PROBABILITY -
            row.exit_bound
        thresholds = remaining > 0 ?
            sqrt.(
                2 .* row.uncertainty_scale^2 .*
                final_steps.emissions_variance_upper_bound .*
                log(1 / remaining),
            ) : [NaN]
        push!(ridge_rows, (;
            case_name=tail_case.case_name,
            ridge,
            uncertainty_scale=row.uncertainty_scale,
            load_forecast_std_percent=row.load_forecast_std_percent,
            certificate_available=remaining > 0,
            median_variance_bound_ratio=median(
                final_steps.emissions_variance_upper_bound ./
                final_steps.true_emissions_variance,
            ),
            median_threshold=median(thresholds),
            empirical_quantile=row.empirical_quantile,
            known_variance_threshold=row.region_exit_adjusted_threshold,
        ))
    end
end
ridge_sensitivity = DataFrame(ridge_rows)
CSV.write(
    joinpath(AGGREGATE_DIR, "ridge_sensitivity.csv"),
    ridge_sensitivity,
)

lme_values = DataFrame(
    bus=collect(1:audit_case.solved.network.n),
    bus_id=audit_case.solved.network.id_map.bus_ids,
    marginal_emissions=audit_case.ell,
    projected_marginal_emissions=
        audit_case.features.Phi * audit_case.theta,
    demand_std_mw=audit_case.standard_deviation,
)
CSV.write(joinpath(AGGREGATE_DIR, "lme_values.csv"), lme_values)

manifest = Dict(
    "generated_at" => string(now()),
    "julia_version" => string(VERSION),
    "powerio_version" => string(Base.pkgversion(PowerIO)),
    "powerdiff_version" => string(Base.pkgversion(PowerDiff)),
    "powerdiff_revision" => "bebfce7e66e9afe2fbd3083fb690c0e8308ce163",
    "pglib_co2_factor_source_revision" =>
        "f5954d677ff8ec7744b3b80a6292a2d044429e6f",
    "experiment_seeds" => "0:99",
    "observation_noise_std_t_per_mwh" => 0.02,
    "sequential_experiment_case" => audit_case.case_name,
    "sequential_observation_budget" => 80,
    "congestion_sweep_case" => audit_case.case_name,
    "congestion_sweep_line_limit_scales" => "1.3,1.2,1.1,1.0,0.9",
    "perturbation_step_sizes_mw" => "0.1,1,10,50,100,250,500",
    "binding_set_extra_line_counts" => "2,5",
    "ridge_parameter" => RIDGE_PARAMETER,
    "ridge_sweep" => join(RIDGE_SWEEP, ","),
    "demand_correlations" => join(DEMAND_CORRELATIONS, ","),
    "learning_failure_probability" => LEARNING_FAILURE_PROBABILITY,
    "certificate_observations" => CERTIFICATE_OBSERVATIONS,
    "certificate_dispatch_runs" => 2CERTIFICATE_OBSERVATIONS,
    "certificate_case" => tail_case.case_name,
    "critical_region_exit_curve_monte_carlo_samples" => 50_000,
    "emissions_tail_case" => tail_case.case_name,
    "emissions_tail_samples_per_scale" => 5_000,
    "emissions_tail_inside_validation_solves_per_scale" => 100,
    "emissions_tail_failure_probability" => TOTAL_FAILURE_PROBABILITY,
    "emissions_tail_exit_bound" => "Gaussian facet union bound",
    "scaling_cases" => nrow(case_summary),
    "cost_ridge_per_unit" => 1e-4,
    "uncertainty_std_fraction" => 0.01,
)
open(joinpath(ROOT, "results", "manifest.json"), "w") do io
    JSON3.pretty(io, manifest)
end

display(case_summary)
display(derivative_validation)
display(critical_region_exit)
display(recovery_summary)
display(step_study)
display(binding_study)
selected_steps = in.(
    experiment_efficiency.observation_index,
    Ref([5, 10, 20, 40, 80]),
)
display(experiment_efficiency[selected_steps, :])
