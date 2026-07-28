using OperatingEmissionsCertificates
using LinearAlgebra
using PowerIO
using Random
using Test

@testset "PGLib-CO2 factors" begin
    @test PGLIB_CO2_FACTORS["COW"] > PGLIB_CO2_FACTORS["NG"] > 0
    @test PGLIB_CO2_FACTORS["SYNC"] == 0
end

@testset "Quadratic emissions variance upper bound" begin
    rng = MersenneTwister(17)
    for _ in 1:100
        r = 5
        V_root = randn(rng, r, r)
        M_root = randn(rng, r, r)
        V = V_root' * V_root + I
        M = M_root' * M_root
        theta_hat = randn(rng, r)
        beta = 0.7
        direction = randn(rng, r)
        direction ./= sqrt(dot(direction, V * direction))
        theta = theta_hat + beta * rand(rng) * direction
        upper = quadratic_variance_upper_bound(theta_hat, V, M, beta)
        @test dot(theta, M * theta) <= upper + 1e-10
    end
end

@testset "PowerIO and PowerDiff integration" begin
    case_path = joinpath(
        @__DIR__,
        "..",
        "data",
        "meshed",
        "pglib_opf_case14_ieee.m",
    )
    parsed = PowerIO.parse_file(case_path)
    @test PowerIO.n_buses(parsed) == 14

    result = analyze_case(
        case_path;
        line_limit_scale=0.6,
        exit_probability_samples=1_000,
        finite_difference_count=5,
    )
    @test length(result.active.line_indices) == 1
    @test result.features.feature_rank == 2
    @test result.derivative_error <= 1e-10
    @test result.span_residual <= 1e-10
    @test result.poisson.residual <= 1e-9
    @test result.poisson.leakage <= 1e-9
    endpoint_features = result.features.Phi[result.poisson.endpoints, :]
    @test rank(endpoint_features) == result.features.feature_rank
    endpoint_theta = endpoint_features \ result.ell[result.poisson.endpoints]
    @test norm(result.features.Phi * endpoint_theta - result.ell) <= 1e-9
    @test result.finite_difference_error <= 1e-4
    @test result.empirical_exit_probability <=
          result.exit_probability_bound + 1e-3

    ablation = basis_ablation(
        result.solved.network,
        result.ell,
        result.Sigma,
        result.features.Phi;
        max_dimension=5,
    )
    active_rows = ablation[ablation.basis .== "active PTDF", :]
    @test active_rows.sigma_error[end] <= 1e-10

    recovery = recover_marginal_emissions(
        result.solved,
        result.features,
        result.ell,
    )
    @test length(recovery.selected_buses) == result.features.feature_rank
    @test recovery.relative_error <= 1e-4

    curve = critical_region_exit_curve(
        result.region.F,
        result.region.offsets,
        result.Sigma;
        scale_factors=[0.8, 1.0],
        samples=1_000,
    )
    @test all(
        curve.empirical_exit_probability .<=
        curve.exit_probability_bound .+ 0.01,
    )
    gaussian_exit_bound, _, _ = gaussian_critical_region_exit_bound(
        result.region.F,
        result.region.offsets,
        result.Sigma,
    )
    @test gaussian_exit_bound <= result.exit_probability_bound

    tail_curve = evaluate_redispatch_tail(
        result.solved,
        result.ell,
        result.Sigma,
        result.region.F,
        result.region.offsets;
        scale_factors=[1.0],
        samples=50,
        failure_probability=0.10,
        exit_bound_method=:gaussian,
    )
    @test size(tail_curve, 1) == 1
    @test tail_curve.dispatch_solves[1] == 50
    @test tail_curve.certificate_available[1]
    @test tail_curve.maximum_inside_linearization_error[1] <= 1e-4
    @test tail_curve.adjusted_violation_rate[1] <= 0.10

    experiment_log = run_sequential_perturbation_experiment(
        result.features.Phi,
        result.theta,
        result.Sigma;
        policies=("variance_informed", "uniform"),
        seeds=0:2,
        observation_budget=5,
    )
    @test size(experiment_log, 1) == 30
    @test all(
        experiment_log.emissions_variance_upper_bound .>=
        experiment_log.true_emissions_variance,
    )

    steps = perturbation_step_study(
        result.solved,
        result.features,
        result.ell,
        result.region;
        step_sizes_mw=(1.0, 500.0),
    )
    # A small step stays inside the region and recovers the field, while a step
    # large enough to leave it loses accuracy or the dispatch outright.
    @test steps.inside_region[1]
    @test steps.relative_field_error[1] <= 1e-4
    @test !steps.inside_region[2]
    @test !steps.solved_perturbations[2] ||
          steps.relative_field_error[2] > steps.relative_field_error[1]

    declared = binding_set_study(
        result.solved,
        result.ptdf,
        result.active,
        result.ell;
        extra_line_counts=(2,),
    )
    exact = only(declared[declared.declared_set .== "exact", :])
    @test exact.recovery_error <= 1e-9
    superset = only(declared[startswith.(declared.declared_set, "add"), :])
    # Naming lines that are not binding enlarges the basis but stays exact.
    @test superset.feature_rank > exact.feature_rank
    @test superset.recovery_error <= 1e-9
end

@testset "Correlated uncertainty proxy" begin
    demand_pu = [0.5, 1.0, 1.5]
    base_mva = 100.0
    independent, sd = build_uncertainty_proxy(demand_pu, base_mva)
    correlated, correlated_sd = build_correlated_uncertainty_proxy(
        demand_pu,
        base_mva;
        correlation=0.5,
    )
    @test sd ≈ correlated_sd
    # Correlation changes the off diagonal entries only, so every marginal
    # variance is preserved and the proxy stays positive semidefinite.
    @test diag(correlated) ≈ diag(independent)
    @test correlated[1, 2] ≈ 0.5 * sd[1] * sd[2]
    @test minimum(eigvals(Symmetric(Matrix(correlated)))) >= -1e-12
end
