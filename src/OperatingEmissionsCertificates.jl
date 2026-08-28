module OperatingEmissionsCertificates

using CSV
using DataFrames
using JSON3
using LinearAlgebra
using PowerDiff
using PowerIO
using Printf
using Random
using SparseArrays
using SpecialFunctions
using Statistics

export PGLIB_CO2_FACTORS
export analyze_case
export basis_ablation
export certificate_availability_scale
export binding_set_study
export build_correlated_uncertainty_proxy
export build_uncertainty_proxy
export facet_aggregate_ratio
export operating_emissions
export recover_marginal_emissions
export emission_factors
export evaluate_redispatch_tail
export find_operating_point
export graph_poisson_diagnostics
export perturbation_step_study
export quadratic_variance_upper_bound
export run_sequential_perturbation_experiment
export critical_region_exit_bound
export gaussian_critical_region_exit_bound
export critical_region_exit_curve

"""
Direct operating emission factors in metric tonnes of CO₂ per MWh.

The values reproduce the lookup table distributed with PGLib-CO2. They describe
the benchmark's operating emissions, not life cycle greenhouse gas emissions.
"""
const PGLIB_CO2_FACTORS = Dict(
    "ANT" => 0.9095,
    "COW" => 0.8204,
    "PEL" => 0.7001,
    "NG" => 0.5173,
    "CCGT" => 0.3621,
    "ICE" => 0.6030,
    "THERMAL" => 0.6874,
    "NUC" => 0.0,
    "NUCLEAR" => 0.0,
    "RE" => 0.0,
    "HYD" => 0.0,
    "HYDRO" => 0.0,
    "SOLAR" => 0.0,
    "WIND" => 0.0,
    "SYNC" => 0.0,
    "SYNC_COND" => 0.0,
    "N/A" => 0.0,
)

const DEFAULT_PRIMAL_TOL = 2e-5
const DEFAULT_DUAL_TOL = 1e-6

"""
Read generator fuel labels from comments in a MATPOWER generator table.

PGLib-OPF stores labels such as `NG`, `COW`, and `SYNC` after each generator
row. Returning labels by source row preserves alignment after PowerIO filters
out of service records.
"""
function _matpower_generator_fuels(path::AbstractString)
    fuels = String[]
    in_table = false
    for line in eachline(path)
        if occursin(r"^\s*mpc\.gen\s*=\s*\[", line)
            in_table = true
            continue
        end
        in_table || continue
        occursin(r"^\s*\];", line) && break
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "%")) && continue
        occursin(';', line) || continue
        label_match = match(r";\s*%\s*([A-Za-z0-9_/-]+)\s*$", line)
        push!(fuels, isnothing(label_match) ? "N/A" : uppercase(label_match.captures[1]))
    end
    isempty(fuels) && throw(ArgumentError("no MATPOWER generator rows found in $path"))
    return fuels
end

"""
Return an emission factor and fuel label for each generator retained by
PowerDiff. Source row identifiers from PowerIO prevent indexing errors when a
case contains disabled generators.
"""
function emission_factors(path::AbstractString, network::PowerDiff.DCNetwork)
    source_fuels = _matpower_generator_fuels(path)
    factors = zeros(network.k)
    fuels = Vector{String}(undef, network.k)
    for (local_index, source_index) in enumerate(network.id_map.gen_ids)
        source_index <= length(source_fuels) ||
            throw(ArgumentError("generator source row $source_index is absent from $path"))
        fuel = source_fuels[source_index]
        haskey(PGLIB_CO2_FACTORS, fuel) ||
            throw(ArgumentError("no PGLib-CO2 factor is configured for fuel label $fuel"))
        fuels[local_index] = fuel
        factors[local_index] = PGLIB_CO2_FACTORS[fuel]
    end
    return factors, fuels
end

"""Total hourly operating emissions in metric tonnes of CO2."""
operating_emissions(pg_pu, factors, base_mva) =
    base_mva * dot(factors, pg_pu)

function _copy_network_with_settings(
    network::PowerDiff.DCNetwork;
    line_limit_scale::Float64,
    cost_ridge::Float64,
    flow_regularization::Float64,
)
    cost_ridge > 0 || throw(ArgumentError("cost_ridge must be positive"))
    line_limit_scale > 0 || throw(ArgumentError("line_limit_scale must be positive"))
    return PowerDiff.DCNetwork(
        network.n,
        network.m,
        network.k,
        network.A,
        network.G_inc,
        network.b;
        sw=copy(network.sw),
        fmax=line_limit_scale .* network.fmax,
        gmax=copy(network.gmax),
        gmin=copy(network.gmin),
        angmax=copy(network.angmax),
        angmin=copy(network.angmin),
        cq=network.cq .+ cost_ridge,
        cl=copy(network.cl),
        c_shed=copy(network.c_shed),
        demand=copy(network.demand),
        pg_init=copy(network.pg_init),
        ref_bus=network.ref_bus,
        tau=flow_regularization,
    )
end

function _solve_case(
    path::AbstractString;
    load_scale::Float64=1.0,
    line_limit_scale::Float64=1.0,
    cost_ridge::Float64=1e-4,
    flow_regularization::Float64=0.0,
)
    parsed = PowerIO.parse_file(path)
    base = PowerDiff.DCNetwork(parsed; tau=flow_regularization)
    network = _copy_network_with_settings(
        base;
        line_limit_scale=line_limit_scale,
        cost_ridge=cost_ridge,
        flow_regularization=flow_regularization,
    )
    demand = load_scale .* network.demand
    problem = PowerDiff.DCOPFProblem(network, demand)
    solution = PowerDiff.solve!(problem)
    base_mva = Float64(PowerIO.base_mva(parsed))
    factors, fuels = emission_factors(path, network)
    return (; parsed, network, demand, problem, solution, base_mva, factors, fuels)
end

function _active_constraints(
    network::PowerDiff.DCNetwork,
    solution::PowerDiff.DCOPFSolution;
    primal_tol::Float64=DEFAULT_PRIMAL_TOL,
    dual_tol::Float64=DEFAULT_DUAL_TOL,
)
    m, k = network.m, network.k
    line_upper_slack = network.fmax .- solution.f
    line_lower_slack = network.fmax .+ solution.f
    gen_upper_slack = network.gmax .- solution.pg
    gen_lower_slack = solution.pg .- network.gmin

    line_upper = findall((line_upper_slack .<= primal_tol) .&
                         (solution.lam_ub .>= dual_tol))
    line_lower = findall((line_lower_slack .<= primal_tol) .&
                         (solution.lam_lb .>= dual_tol))
    gen_upper = findall((gen_upper_slack .<= primal_tol) .&
                        (solution.rho_ub .>= dual_tol))
    gen_lower = findall((gen_lower_slack .<= primal_tol) .&
                        (solution.rho_lb .>= dual_tol))

    all_slacks = vcat(
        line_upper_slack,
        line_lower_slack,
        gen_upper_slack,
        gen_lower_slack,
    )
    indices = vcat(
        line_upper,
        m .+ line_lower,
        2m .+ gen_upper,
        2m + k .+ gen_lower,
    )
    line_indices = sort!(unique!(vcat(line_upper, line_lower)))
    return (;
        indices,
        line_indices,
        line_upper,
        line_lower,
        gen_upper,
        gen_lower,
        all_slacks,
    )
end

function _ptdf(network::PowerDiff.DCNetwork)
    zero_state = PowerDiff.DCPowerFlowState(
        network,
        zeros(network.n),
        zeros(network.n),
    )
    return PowerDiff.ptdf_matrix(zero_state)
end

function _inequality_matrices(network::PowerDiff.DCNetwork, ptdf)
    generator_map = Matrix(network.G_inc)
    line_generator = ptdf * generator_map
    k = network.k
    generator_identity = Matrix{Float64}(I, k, k)
    A = vcat(
        line_generator,
        -line_generator,
        generator_identity,
        -generator_identity,
    )
    b = vcat(
        network.fmax,
        network.fmax,
        network.gmax,
        -network.gmin,
    )
    E = vcat(
        ptdf,
        -ptdf,
        zeros(2k, network.n),
    )
    labels = vcat(
        ["line_upper_$j" for j in 1:network.m],
        ["line_lower_$j" for j in 1:network.m],
        ["gen_upper_$j" for j in 1:k],
        ["gen_lower_$j" for j in 1:k],
    )
    return (; A, b, E, labels)
end

function _independent_rows(K::AbstractMatrix; tolerance::Float64=1e-10)
    factorization = qr(transpose(K), ColumnNorm())
    diagonal = abs.(diag(factorization.R))
    isempty(diagonal) && return Int[]
    rank_K = count(>=(tolerance * maximum(diagonal)), diagonal)
    return sort!(factorization.p[1:rank_K])
end

function _fixed_active_model(
    network::PowerDiff.DCNetwork,
    demand,
    solution,
    ptdf,
    active;
    rank_tolerance::Float64=1e-10,
)
    inequalities = _inequality_matrices(network, ptdf)
    balance_K = ones(1, network.k)
    balance_D = ones(1, network.n)
    K_all = vcat(balance_K, inequalities.A[active.indices, :])
    D_all = vcat(balance_D, inequalities.E[active.indices, :])

    independent = _independent_rows(K_all; tolerance=rank_tolerance)
    1 in independent || throw(ArgumentError("balance row was removed as dependent"))
    K = K_all[independent, :]
    D = D_all[independent, :]
    active_positions = [i - 1 for i in independent if i > 1]
    selected_active = active.indices[active_positions]

    q_diagonal = 2 .* network.cq
    minimum(q_diagonal) > 0 ||
        throw(ArgumentError("the dispatch Hessian is not positive definite"))
    Q_inverse_Kt = Diagonal(1 ./ q_diagonal) * transpose(K)
    W = Symmetric(K * Q_inverse_Kt)
    X = W \ D
    G = Q_inverse_Kt * X
    N = -X

    gradient = q_diagonal .* solution.pg .+ network.cl
    multipliers = -(W \ (K * (gradient ./ q_diagonal)))
    inequality_multipliers = multipliers[2:end]

    local_error = norm(
        K * solution.pg - vcat(sum(demand),
                               inequalities.b[selected_active] .+
                               inequalities.E[selected_active, :] * demand),
        Inf,
    )

    return (;
        G,
        N,
        K,
        D,
        W,
        multipliers,
        inequality_multipliers,
        selected_active,
        inequalities,
        local_error,
        condition_number=cond(Matrix(W)),
    )
end

function _congestion_features(ptdf, active_lines; rank_tolerance::Float64=1e-10)
    n = size(ptdf, 2)
    raw = isempty(active_lines) ?
        ones(n, 1) :
        hcat(ones(n), transpose(ptdf[active_lines, :]))
    decomposition = svd(raw)
    threshold = rank_tolerance * first(decomposition.S)
    feature_rank = count(>=(threshold), decomposition.S)
    Phi = decomposition.U[:, 1:feature_rank]
    return (; Phi, raw, feature_rank, singular_values=decomposition.S)
end

"""
Verify the graph Poisson form of the marginal emissions field.

For a connected DC network, the weighted Laplacian applied to the centered
field must lie in the span of incidence vectors for binding lines. The
resulting source is therefore supported only at endpoints of those lines.
"""
function graph_poisson_diagnostics(network, ell, active_lines)
    weights = -network.b .* network.sw
    laplacian = sparse(network.A' * Diagonal(weights) * network.A)
    centered = ell .- mean(ell)
    source = Vector(laplacian * centered)
    # The absolute normalization avoids amplifying roundoff when the field is
    # nearly constant and the true Poisson source is zero.
    scale = max(1.0, norm(source))

    if isempty(active_lines)
        fitted_source = zeros(length(source))
        endpoints = Int[]
    else
        dictionary = Matrix(transpose(network.A[active_lines, :]))
        fitted_source = dictionary * (dictionary \ source)
        endpoint_mask = vec(
            sum(abs.(network.A[active_lines, :]); dims=1),
        ) .> 0
        endpoints = findall(endpoint_mask)
    end

    outside = setdiff(eachindex(source), endpoints)
    residual = norm(source - fitted_source) / scale
    leakage = isempty(outside) ? 0.0 : norm(source[outside]) / scale
    return (;
        laplacian,
        source,
        fitted_source,
        endpoints,
        residual,
        leakage,
    )
end

"""
Compare the physical active-PTDF representation with a generic graph Fourier
basis. Errors are reported in both Euclidean and uncertainty-weighted norms.
"""
function basis_ablation(
    network,
    ell,
    Sigma,
    active_features;
    max_dimension::Int=40,
)
    n = length(ell)
    largest_dimension = min(max_dimension, n)
    graph_laplacian = Symmetric(Matrix(network.A' * network.A))
    graph_eigenvectors = eigen(graph_laplacian).vectors[:, 1:largest_dimension]
    risk_scale = max(dot(ell, Sigma * ell), eps(Float64))
    l2_scale = max(norm(ell), eps(Float64))
    rows = NamedTuple[]

    function add_row!(basis_name, dimension, basis)
        estimate = basis * (transpose(basis) * ell)
        error = ell - estimate
        push!(rows, (;
            basis=basis_name,
            dimension,
            l2_error=norm(error) / l2_scale,
            sigma_error=sqrt(max(0.0, dot(error, Sigma * error)) / risk_scale),
        ))
    end

    active_rank = size(active_features, 2)
    for dimension in 1:active_rank
        add_row!(
            "active PTDF",
            dimension,
            active_features[:, 1:dimension],
        )
    end
    for dimension in 1:largest_dimension
        add_row!(
            "graph Fourier",
            dimension,
            graph_eigenvectors[:, 1:dimension],
        )
    end
    return DataFrame(rows)
end

function _critical_region(
    network,
    demand,
    solution,
    active,
    fixed_model,
    base_mva;
    offset_tolerance::Float64=1e-9,
    normal_tolerance::Float64=1e-12,
)
    inequalities = fixed_model.inequalities
    inactive = setdiff(1:size(inequalities.A, 1), fixed_model.selected_active)
    primal_F_pu = inequalities.A[inactive, :] * fixed_model.G -
                  inequalities.E[inactive, :]
    primal_s = active.all_slacks[inactive]
    primal_labels = ["primal_" * inequalities.labels[i] for i in inactive]

    dual_F_pu = -fixed_model.N[2:end, :]
    dual_s = fixed_model.inequality_multipliers
    dual_labels = ["dual_" * inequalities.labels[i]
                   for i in fixed_model.selected_active]

    F_pu = vcat(primal_F_pu, dual_F_pu)
    offsets = vcat(primal_s, dual_s)
    labels = vcat(primal_labels, dual_labels)
    keep = [
        offsets[i] > offset_tolerance &&
        norm(view(F_pu, i, :)) > normal_tolerance
        for i in eachindex(offsets)
    ]
    # Demand uncertainty is represented in MW in the statistical layer.
    F_mw = F_pu[keep, :] ./ base_mva
    return (; F=F_mw, offsets=offsets[keep], labels=labels[keep])
end

"""
Construct a diagonal sub-Gaussian proxy in MW². Each bus receives a standard
deviation equal to `fraction` of nominal demand, subject to a small floor.
"""
function build_uncertainty_proxy(
    demand_pu,
    base_mva;
    fraction::Float64=0.01,
    floor_mw::Float64=0.25,
)
    demand_mw = base_mva .* demand_pu
    standard_deviation = fraction .* max.(abs.(demand_mw), floor_mw)
    return Diagonal(standard_deviation .^ 2), standard_deviation
end

"""
Construct a sub-Gaussian proxy in MW² whose nodal demand forecast errors share
a common system wide component. Every bus keeps the marginal standard deviation
produced by `build_uncertainty_proxy`, so the only change relative to that proxy
is the correlation `rho` between any two buses.
"""
function build_correlated_uncertainty_proxy(
    demand_pu,
    base_mva;
    fraction::Float64=0.01,
    floor_mw::Float64=0.25,
    correlation::Float64=0.5,
)
    0 <= correlation < 1 ||
        throw(ArgumentError("correlation must lie in [0, 1)"))
    _, standard_deviation = build_uncertainty_proxy(
        demand_pu,
        base_mva;
        fraction,
        floor_mw,
    )
    common = correlation .* (standard_deviation * transpose(standard_deviation))
    independent = (1 - correlation) .* Diagonal(standard_deviation .^ 2)
    return Symmetric(common + independent), standard_deviation
end

"""Upper quantile of the standard normal distribution."""
normal_quantile(p) = sqrt(2) * erfinv(2p - 1)

"""
Report how strongly each facet of the critical region responds to a common
demand component.

Under the proxy `rho * s * s' + (1 - rho) * Diagonal(s.^2)` with marginal
standard deviations `s`, facet `j` satisfies
`F_j' * Sigma_rho * F_j = (1 - rho + rho * a_j) * F_j' * Sigma_0 * F_j`, where
the aggregate response ratio `a_j = (s' F_j)^2 / sum_i (s_i F_ji)^2` lies in
`[0, n]`. The standardized distance therefore scales exactly as
`kappa_j(rho) = kappa_j(0) / sqrt(1 - rho + rho * a_j)`. Facets with `a_j > 1`
lose margin under correlation, while facets with `a_j < 1` gain margin. The
constraint type alone does not determine which case applies.
"""
function facet_aggregate_ratio(F, standard_deviation)
    ratios = zeros(size(F, 1))
    for j in axes(F, 1)
        weighted = standard_deviation .* view(F, j, :)
        denominator = sum(abs2, weighted)
        ratios[j] = denominator > 0 ? sum(weighted)^2 / denominator : 0.0
    end
    return ratios
end

"""
Find the largest uncertainty scale at which the region exit bound still leaves
room in the failure budget.

The exit bound increases with the scale, so a bisection on
`gaussian_critical_region_exit_bound(F, offsets, scale^2 * Sigma)` locates
the crossing.
Returning `NaN` means the bound already exceeds `budget` at `minimum_scale`.
"""
function certificate_availability_scale(
    F,
    offsets,
    Sigma;
    budget::Float64=0.045,
    minimum_scale::Float64=1e-3,
    maximum_scale::Float64=4.0,
    tolerance::Float64=1e-4,
)
    0 < budget < 1 || throw(ArgumentError("budget must lie in (0, 1)"))
    0 < minimum_scale < maximum_scale ||
        throw(ArgumentError("scale bracket must be positive and ordered"))
    exit_bound(scale) =
        first(gaussian_critical_region_exit_bound(
            F,
            offsets,
            scale^2 .* Sigma,
        ))

    exit_bound(minimum_scale) >= budget && return NaN
    exit_bound(maximum_scale) < budget && return maximum_scale

    low, high = minimum_scale, maximum_scale
    while high - low > tolerance
        middle = 0.5 * (low + high)
        if exit_bound(middle) < budget
            low = middle
        else
            high = middle
        end
    end
    return low
end

"""
Compute the facet union bound for escape from a critical region.
"""
function critical_region_exit_bound(F, offsets, Sigma)
    contributions = zeros(length(offsets))
    kappas = fill(Inf, length(offsets))
    for i in eachindex(offsets)
        variance = dot(view(F, i, :), Sigma * view(F, i, :))
        if variance > eps(Float64)
            kappas[i] = offsets[i] / sqrt(variance)
            contributions[i] = exp(-0.5 * kappas[i]^2)
        end
    end
    return min(1.0, sum(contributions)), contributions, kappas
end

"""
Compute the tighter facet union bound available for Gaussian demand errors.
"""
function gaussian_critical_region_exit_bound(F, offsets, Sigma)
    contributions = zeros(length(offsets))
    kappas = fill(Inf, length(offsets))
    for i in eachindex(offsets)
        variance = dot(view(F, i, :), Sigma * view(F, i, :))
        if variance > eps(Float64)
            kappas[i] = offsets[i] / sqrt(variance)
            contributions[i] =
                0.5 * erfc(kappas[i] / sqrt(2.0))
        end
    end
    return min(1.0, sum(contributions)), contributions, kappas
end

function _empirical_escape_probability(F, offsets, Sigma, samples, rng)
    isempty(offsets) && return 0.0
    samples <= 0 && return NaN
    covariance_root = if Sigma isa Diagonal
        Diagonal(sqrt.(max.(diag(Sigma), 0.0)))
    else
        cholesky(Symmetric(Matrix(Sigma) + 1e-14I)).L
    end
    escapes = 0
    chunk_size = min(256, samples)
    for first_sample in 1:chunk_size:samples
        count = min(chunk_size, samples - first_sample + 1)
        xi = covariance_root * randn(rng, size(Sigma, 1), count)
        violations = F * xi .> offsets
        escapes += sum(vec(any(violations; dims=1)))
    end
    return escapes / samples
end

function _escape_scores(F, offsets, Sigma, samples, rng)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    isempty(offsets) && return fill(-Inf, samples)
    covariance_root = if Sigma isa Diagonal
        Diagonal(sqrt.(max.(diag(Sigma), 0.0)))
    else
        cholesky(Symmetric(Matrix(Sigma) + 1e-14I)).L
    end
    scores = zeros(samples)
    chunk_size = min(256, samples)
    for first_sample in 1:chunk_size:samples
        count = min(chunk_size, samples - first_sample + 1)
        indices = first_sample:(first_sample + count - 1)
        xi = covariance_root * randn(rng, size(Sigma, 1), count)
        normalized = (F * xi) ./ offsets
        scores[indices] .= vec(maximum(normalized; dims=1))
    end
    return scores
end

function _wilson_interval(successes, samples; z=1.959963984540054)
    proportion = successes / samples
    denominator = 1 + z^2 / samples
    center = (proportion + z^2 / (2samples)) / denominator
    radius = z / denominator * sqrt(
        proportion * (1 - proportion) / samples + z^2 / (4samples^2),
    )
    return max(0.0, center - radius), min(1.0, center + radius)
end

"""
Evaluate the critical-region escape bound over a declared uncertainty scale
path. One shared set of Gaussian draws is reused at every scale.
"""
function critical_region_exit_curve(
    F,
    offsets,
    Sigma;
    scale_factors=range(0.4, 2.0; length=13),
    samples::Int=50_000,
    seed::Int=2026,
)
    rng = MersenneTwister(seed)
    scores = _escape_scores(F, offsets, Sigma, samples, rng)
    rows = NamedTuple[]
    for scale in scale_factors
        scale > 0 || throw(ArgumentError("scale factors must be positive"))
        scaled_proxy = scale^2 .* Sigma
        bound, _, kappas = critical_region_exit_bound(
            F,
            offsets,
            scaled_proxy,
        )
        successes = sum(scale .* scores .> 1)
        empirical = successes / samples
        lower, upper = _wilson_interval(successes, samples)
        push!(rows, (;
            uncertainty_scale=scale,
            closest_kappa=minimum(kappas),
            exit_probability_bound=bound,
            empirical_exit_probability=empirical,
            empirical_lower=lower,
            empirical_upper=upper,
            samples,
        ))
    end
    return DataFrame(rows)
end

"""
Evaluate the complete operating-emissions bound against repeated DC dispatch.

For each uncertainty scale, this routine applies shared Gaussian demand errors,
solves the dispatch problem for every realization, and compares the resulting
emissions changes with two thresholds. `fixed_region_threshold` uses only the
nominal marginal emissions vector. `region_exit_adjusted_threshold` allocates
part of the requested failure probability to leaving the nominal critical
region.
It is reported as `NaN` when the critical-region exit bound already exceeds
the requested failure probability. If `solve_all_samples` is false, dispatch
is solved for every realization outside the critical region and for a declared
number inside it. The exact affine continuation is used for the remaining
inside-region realizations.
"""
function evaluate_redispatch_tail(
    solved,
    ell,
    Sigma,
    F,
    offsets;
    scale_factors=(0.35, 0.45, 0.55, 0.65, 0.8, 1.0),
    samples::Int=2_000,
    failure_probability::Float64=0.05,
    seed::Int=2026,
    solve_all_samples::Bool=true,
    inside_validation_samples::Int=0,
    exit_bound_method::Symbol=:subgaussian,
)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    inside_validation_samples >= 0 ||
        throw(ArgumentError("inside_validation_samples must be nonnegative"))
    0 < failure_probability < 1 ||
        throw(ArgumentError("failure_probability must lie in (0, 1)"))
    exit_bound_method in (:subgaussian, :gaussian) ||
        throw(ArgumentError(
            "exit_bound_method must be :subgaussian or :gaussian",
        ))
    covariance_root = if Sigma isa Diagonal
        Diagonal(sqrt.(max.(diag(Sigma), 0.0)))
    else
        cholesky(Symmetric(Matrix(Sigma) + 1e-14I)).L
    end
    rng = MersenneTwister(seed)
    standardized_errors =
        covariance_root * randn(rng, size(Sigma, 1), samples)
    nominal_emissions = operating_emissions(
        solved.solution.pg,
        solved.factors,
        solved.base_mva,
    )
    base_emissions_variance = dot(ell, Sigma * ell)
    rows = NamedTuple[]

    for scale in scale_factors
        scale > 0 || throw(ArgumentError("scale factors must be positive"))
        demand_errors = scale .* standardized_errors
        linear_changes = vec(transpose(ell) * demand_errors)
        actual_changes = copy(linear_changes)
        escaped = if isempty(offsets)
            falses(samples)
        else
            vec(any(F * demand_errors .> offsets; dims=1))
        end

        solved_samples = solve_all_samples ? trues(samples) : copy(escaped)
        if !solve_all_samples && inside_validation_samples > 0
            inside_indices = findall(.!escaped)
            validation_count =
                min(inside_validation_samples, length(inside_indices))
            solved_samples[inside_indices[1:validation_count]] .= true
        end
        for sample in findall(solved_samples)
            demand = solved.demand .+
                     view(demand_errors, :, sample) ./ solved.base_mva
            solution = PowerDiff.solve!(
                PowerDiff.DCOPFProblem(solved.network, demand),
            )
            actual_changes[sample] = operating_emissions(
                solution.pg,
                solved.factors,
                solved.base_mva,
            ) - nominal_emissions
        end

        scaled_proxy = scale^2 .* Sigma
        exit_bound, _, kappas = exit_bound_method == :gaussian ?
            gaussian_critical_region_exit_bound(F, offsets, scaled_proxy) :
            critical_region_exit_bound(F, offsets, scaled_proxy)
        emissions_variance = scale^2 * base_emissions_variance
        fixed_region_threshold = sqrt(
            2emissions_variance * log(1 / failure_probability),
        )
        certificate_available = exit_bound < failure_probability
        region_exit_adjusted_threshold = certificate_available ?
            sqrt(
                2emissions_variance *
                log(1 / (failure_probability - exit_bound)),
            ) :
            NaN
        # The Chernoff threshold applies to any sub-Gaussian forecast error.
        # For Gaussian errors, the exact normal quantile isolates the margin
        # introduced by the weaker distributional assumption.
        gaussian_fixed_region_threshold =
            sqrt(emissions_variance) *
            normal_quantile(1 - failure_probability)
        gaussian_region_exit_adjusted_threshold = certificate_available ?
            sqrt(emissions_variance) *
            normal_quantile(1 - (failure_probability - exit_bound)) :
            NaN
        gaussian_violations = certificate_available ?
            sum(
                actual_changes .>=
                gaussian_region_exit_adjusted_threshold,
            ) : 0
        fixed_violations = sum(actual_changes .>= fixed_region_threshold)
        fixed_lower, fixed_upper = _wilson_interval(
            fixed_violations,
            samples,
        )
        adjusted_violations = certificate_available ?
            sum(actual_changes .>= region_exit_adjusted_threshold) : 0
        adjusted_lower, adjusted_upper = certificate_available ?
            _wilson_interval(adjusted_violations, samples) : (NaN, NaN)
        inside = .!escaped
        validated_inside = inside .& solved_samples
        inside_errors = abs.(
            actual_changes[validated_inside] .-
            linear_changes[validated_inside],
        )
        outside_errors =
            abs.(actual_changes[escaped] .- linear_changes[escaped])
        fixed_violation_exit_fraction = fixed_violations == 0 ? 0.0 :
            sum(escaped .& (actual_changes .>= fixed_region_threshold)) /
            fixed_violations

        push!(rows, (;
            uncertainty_scale=scale,
            exit_bound_method=String(exit_bound_method),
            closest_kappa=isempty(kappas) ? Inf : minimum(kappas),
            empirical_exit=mean(escaped),
            exit_bound,
            failure_probability,
            emissions_change_standard_deviation=sqrt(emissions_variance),
            empirical_quantile=quantile(
                actual_changes,
                1 - failure_probability,
            ),
            linear_quantile=quantile(
                linear_changes,
                1 - failure_probability,
            ),
            fixed_region_threshold,
            region_exit_adjusted_threshold,
            gaussian_fixed_region_threshold,
            gaussian_region_exit_adjusted_threshold,
            gaussian_violation_rate=
                certificate_available ? gaussian_violations / samples : NaN,
            certificate_available,
            fixed_violation_rate=fixed_violations / samples,
            fixed_violation_lower=fixed_lower,
            fixed_violation_upper=fixed_upper,
            adjusted_violation_rate=
                certificate_available ? adjusted_violations / samples : NaN,
            adjusted_violation_lower=adjusted_lower,
            adjusted_violation_upper=adjusted_upper,
            fixed_violation_exit_fraction,
            maximum_inside_linearization_error=
                isempty(inside_errors) ? NaN : maximum(inside_errors),
            median_outside_linearization_error=
                isempty(outside_errors) ? NaN : median(outside_errors),
            samples,
            dispatch_solves=sum(solved_samples),
            inside_validation_solves=sum(validated_inside),
        ))
    end
    return DataFrame(rows)
end

"""
Closed form upper confidence bound for `theta' * M * theta` over an ellipsoid
`norm(theta - theta_hat, V) <= beta`.
"""
function quadratic_variance_upper_bound(theta_hat, V, M, beta)
    beta >= 0 || throw(ArgumentError("beta must be nonnegative"))
    factor = cholesky(Symmetric(V))
    L = factor.L
    estimate = dot(theta_hat, M * theta_hat)
    cross = 2beta * norm(L \ (M * theta_hat))
    weighted = (L \ M) / transpose(L)
    max_eigenvalue = max(0.0, eigmax(Symmetric(weighted)))
    return max(0.0, estimate + cross + beta^2 * max_eigenvalue)
end

function _confidence_radius(
    V,
    ridge_parameter,
    observation_noise_bound,
    parameter_norm_bound,
    estimation_failure_probability,
)
    eigenvalues = eigvals(Symmetric(V))
    logdet_ratio =
        sum(log, eigenvalues) -
        length(eigenvalues) * log(ridge_parameter)
    return observation_noise_bound * sqrt(
        max(
            0.0,
            logdet_ratio + 2log(1 / estimation_failure_probability),
        ),
    ) + sqrt(ridge_parameter) * parameter_norm_bound
end

function _select_perturbation_bus(Phi, V, M, policy, rng)
    n = size(Phi, 1)
    V_inverse = inv(Symmetric(V))
    if policy == "uniform"
        return rand(rng, 1:n)
    elseif policy == "leverage"
        scores = [dot(view(Phi, i, :), V_inverse * view(Phi, i, :))
                  for i in 1:n]
        return argmax(scores)
    elseif policy == "variance_informed"
        scores = zeros(n)
        for i in 1:n
            phi = view(Phi, i, :)
            direction = V_inverse * phi
            scores[i] = dot(direction, M * direction) /
                        (1 + dot(phi, direction))
        end
        return argmax(scores)
    end
    throw(ArgumentError("unknown bus selection policy $policy"))
end

function _single_perturbation_run(
    Phi,
    theta,
    Sigma;
    policy::String,
    seed::Int,
    observation_budget::Int,
    observation_noise_std::Float64,
    ridge_parameter::Float64,
    estimation_failure_probability::Float64,
    parameter_norm_bound::Float64,
)
    rng = MersenneTwister(seed)
    feature_rank = size(Phi, 2)
    M = Symmetric(transpose(Phi) * Sigma * Phi)
    true_emissions_variance = dot(theta, M * theta)
    V = ridge_parameter .* Matrix{Float64}(I, feature_rank, feature_rank)
    response_sum = zeros(feature_rank)
    rows = NamedTuple[]

    for observation_index in 1:observation_budget
        bus = _select_perturbation_bus(Phi, V, M, policy, rng)
        phi = Vector(view(Phi, bus, :))
        observation =
            dot(phi, theta) + observation_noise_std * randn(rng)
        V .+= phi * transpose(phi)
        response_sum .+= phi .* observation
        theta_hat = V \ response_sum
        beta = _confidence_radius(
            V,
            ridge_parameter,
            observation_noise_std,
            parameter_norm_bound,
            estimation_failure_probability,
        )
        estimated_emissions_variance = dot(theta_hat, M * theta_hat)
        emissions_variance_upper_bound =
            quadratic_variance_upper_bound(theta_hat, V, M, beta)
        ell_error = norm(Phi * (theta_hat - theta)) /
                    max(norm(Phi * theta), eps(Float64))
        push!(rows, (;
            policy,
            seed,
            observation_index,
            selected_bus=bus,
            observation,
            true_emissions_variance,
            estimated_emissions_variance,
            emissions_variance_upper_bound,
            variance_bound_ratio=
                emissions_variance_upper_bound /
                max(true_emissions_variance, eps(Float64)),
            marginal_emissions_relative_error=ell_error,
            upper_bound_covers_truth=
                emissions_variance_upper_bound + 1e-12 >=
                true_emissions_variance,
        ))
    end
    return rows
end

"""
Run repeated noisy paired perturbation experiments.

Each observation represents one central finite difference at a selected bus.
The estimation procedure sees only the noisy marginal emissions observation,
not the coefficient vector `theta`.
"""
function run_sequential_perturbation_experiment(
    Phi,
    theta,
    Sigma;
    policies=("variance_informed", "leverage", "uniform"),
    seeds=0:99,
    observation_budget::Int=60,
    observation_noise_std::Float64=0.02,
    ridge_parameter::Float64=1e-8,
    estimation_failure_probability::Float64=0.05,
    parameter_norm_bound::Float64=5sqrt(size(Phi, 1)),
)
    observation_budget > 0 ||
        throw(ArgumentError("observation_budget must be positive"))
    observation_noise_std >= 0 ||
        throw(ArgumentError("observation_noise_std must be nonnegative"))
    ridge_parameter > 0 ||
        throw(ArgumentError("ridge_parameter must be positive"))
    0 < estimation_failure_probability < 1 ||
        throw(ArgumentError(
            "estimation_failure_probability must lie in (0, 1)",
        ))
    parameter_norm_bound >= 0 ||
        throw(ArgumentError("parameter_norm_bound must be nonnegative"))
    all_rows = NamedTuple[]
    for policy in policies, seed in seeds
        append!(
            all_rows,
            _single_perturbation_run(
                Phi,
                theta,
                Sigma;
                policy,
                seed,
                observation_budget,
                observation_noise_std,
                ridge_parameter,
                estimation_failure_probability,
                parameter_norm_bound,
            ),
        )
    end
    return DataFrame(all_rows)
end

function _finite_difference_columns(
    solved,
    columns;
    step_mw::Float64=1.0,
)
    n = solved.network.n
    derivative = fill(NaN, n)
    step_pu = step_mw / solved.base_mva
    for i in columns
        upper_demand = copy(solved.demand)
        lower_demand = copy(solved.demand)
        upper_demand[i] += step_pu
        lower_demand[i] -= step_pu
        upper_solution = PowerDiff.solve!(
            PowerDiff.DCOPFProblem(solved.network, upper_demand),
        )
        lower_solution = PowerDiff.solve!(
            PowerDiff.DCOPFProblem(solved.network, lower_demand),
        )
        upper_emissions = operating_emissions(
            upper_solution.pg,
            solved.factors,
            solved.base_mva,
        )
        lower_emissions = operating_emissions(
            lower_solution.pg,
            solved.factors,
            solved.base_mva,
        )
        derivative[i] =
            (upper_emissions - lower_emissions) / (2step_mw)
    end
    return derivative
end

"""
Recover the complete marginal emissions field from a rank-revealing set of
central finite difference perturbation experiments.

Each selected bus uses two dispatch evaluations, at `d0 ± step_mw`. A full
coordinate finite difference audit would require the same pair at every bus.
"""
function recover_marginal_emissions(
    solved,
    features,
    ell;
    step_mw::Float64=1.0,
)
    Phi = features.Phi
    rank = size(Phi, 2)
    factorization = qr(transpose(Phi), ColumnNorm())
    selected_buses = sort!(factorization.p[1:rank])
    observations = _finite_difference_columns(
        solved,
        selected_buses;
        step_mw,
    )
    selected_features = Phi[selected_buses, :]
    theta_hat = selected_features \ observations[selected_buses]
    ell_hat = Phi * theta_hat
    relative_error = norm(ell_hat - ell) / max(norm(ell), eps(Float64))
    return (;
        selected_buses,
        observations,
        theta_hat,
        ell_hat,
        relative_error,
        paired_perturbation_observations=rank,
        selected_dispatch_runs=2rank,
        full_dispatch_runs=2length(ell),
    )
end

"""
Measure how recovery degrades as the paired perturbation step grows.

Recovery assumes that both perturbed demands stay in the nominal critical
region, where the dispatch response is affine. This routine repeats the
selected central differences at a range of step sizes and reports both the
recovered field error and whether every selected perturbation is provably
inside the region. A perturbation of `step_mw` at bus `i` stays inside when
every facet slack remains positive for both signs of that perturbation.
"""
function perturbation_step_study(
    solved,
    features,
    ell,
    region;
    step_sizes_mw=(0.1, 1.0, 10.0, 50.0, 100.0, 250.0, 500.0),
)
    Phi = features.Phi
    feature_rank = size(Phi, 2)
    factorization = qr(transpose(Phi), ColumnNorm())
    selected_buses = sort!(factorization.p[1:feature_rank])
    rows = NamedTuple[]

    for step_mw in step_sizes_mw
        step = Float64(step_mw)
        step > 0 || throw(ArgumentError("step sizes must be positive"))
        # A step large enough to leave the region can also leave the feasible
        # set, in which case the dispatch model returns no answer at all.
        solved_perturbations = true
        relative_field_error = NaN
        try
            observations = _finite_difference_columns(
                solved,
                selected_buses;
                step_mw=step,
            )
            theta_hat = Phi[selected_buses, :] \ observations[selected_buses]
            relative_field_error = norm(Phi * theta_hat - ell) /
                                   max(norm(ell), eps(Float64))
        catch error
            solved_perturbations = false
            @debug "perturbation step failed" step error
        end

        inside = true
        worst_facet_use = 0.0
        if !isempty(region.offsets)
            for bus in selected_buses
                use = maximum(
                    step .* abs.(view(region.F, :, bus)) ./ region.offsets,
                )
                worst_facet_use = max(worst_facet_use, use)
            end
            inside = worst_facet_use < 1
        end

        push!(rows, (;
            step_mw=step,
            inside_region=inside,
            worst_facet_use,
            solved_perturbations,
            relative_field_error,
            dispatch_runs=2feature_rank,
        ))
    end
    return DataFrame(rows)
end

"""
Measure the effect of a misdeclared binding line set on the congestion basis.

The analyst receives the binding line identities from the model owner and may
receive a set that omits a binding line or that names lines which are not
binding at the nominal operating point. Each declared set produces its own
basis, its own rank, and therefore its own number of dispatch runs.
"""
function binding_set_study(
    solved,
    ptdf,
    active,
    ell;
    extra_line_counts=(2, 5),
)
    binding = collect(active.line_indices)
    isempty(binding) &&
        throw(ArgumentError("the operating point has no binding lines"))
    slack = min.(
        solved.network.fmax .- solved.solution.f,
        solved.network.fmax .+ solved.solution.f,
    )
    nonbinding = [j for j in sortperm(slack) if !(j in binding)]
    rows = NamedTuple[]

    function add_row!(label, lines)
        declared = sort!(unique(collect(lines)))
        features = _congestion_features(ptdf, declared)
        Phi = features.Phi
        span_residual = norm(ell - Phi * (transpose(Phi) * ell)) /
                        max(norm(ell), eps(Float64))
        factorization = qr(transpose(Phi), ColumnNorm())
        selected = sort!(factorization.p[1:features.feature_rank])
        theta_hat = Phi[selected, :] \ ell[selected]
        recovery_error = norm(Phi * theta_hat - ell) /
                         max(norm(ell), eps(Float64))
        push!(rows, (;
            declared_set=label,
            declared_lines=length(declared),
            feature_rank=features.feature_rank,
            dispatch_runs=2features.feature_rank,
            span_residual,
            recovery_error,
        ))
    end

    add_row!("exact", binding)
    for line in binding
        add_row!("omit one binding line", setdiff(binding, [line]))
    end
    for count in extra_line_counts
        count <= length(nonbinding) || continue
        add_row!(
            "add $count slack lines",
            vcat(binding, nonbinding[1:count]),
        )
    end
    return DataFrame(rows)
end

"""
Analyze one MATPOWER case at a documented operating point.
"""
function analyze_case(
    path::AbstractString;
    load_scale::Float64=1.0,
    line_limit_scale::Float64=1.0,
    cost_ridge::Float64=1e-4,
    uncertainty_fraction::Float64=0.01,
    exit_probability_samples::Int=20_000,
    finite_difference_count::Int=12,
    seed::Int=2026,
)
    solved = _solve_case(
        path;
        load_scale,
        line_limit_scale,
        cost_ridge,
        flow_regularization=0.0,
    )
    maximum(solved.solution.psh) <= 1e-7 ||
        throw(ArgumentError("nominal point uses load shedding"))

    ptdf = _ptdf(solved.network)
    active = _active_constraints(solved.network, solved.solution)
    fixed = _fixed_active_model(
        solved.network,
        solved.demand,
        solved.solution,
        ptdf,
        active,
    )

    dispatch_sensitivity = Matrix(PowerDiff.calc_sensitivity(
        solved.problem,
        :pg,
        :d,
    ))
    ell = transpose(dispatch_sensitivity) * solved.factors
    ell_fixed = transpose(fixed.G) * solved.factors
    derivative_error = norm(dispatch_sensitivity - fixed.G) /
                       max(norm(dispatch_sensitivity), eps(Float64))

    features = _congestion_features(ptdf, active.line_indices)
    theta = transpose(features.Phi) * ell
    span_residual = norm(ell - features.Phi * theta) /
                    max(norm(ell), eps(Float64))

    Sigma, standard_deviation = build_uncertainty_proxy(
        solved.demand,
        solved.base_mva;
        fraction=uncertainty_fraction,
    )
    emissions_variance = dot(ell, Sigma * ell)
    region = _critical_region(
        solved.network,
        solved.demand,
        solved.solution,
        active,
        fixed,
        solved.base_mva,
    )
    exit_probability_bound, facet_contributions, kappas =
        critical_region_exit_bound(
            region.F,
            region.offsets,
            Sigma,
        )
    rng = MersenneTwister(seed)
    empirical_exit_probability = _empirical_escape_probability(
        region.F,
        region.offsets,
        Sigma,
        exit_probability_samples,
        rng,
    )

    n_fd = min(finite_difference_count, solved.network.n)
    if n_fd > 0
        fd_columns = unique!(
            round.(Int, range(1, solved.network.n; length=n_fd)),
        )
        ell_fd = _finite_difference_columns(solved, fd_columns)
        finite_difference_error =
            norm(ell[fd_columns] - ell_fd[fd_columns]) /
            max(norm(ell[fd_columns]), eps(Float64))
    else
        fd_columns = Int[]
        ell_fd = fill(NaN, solved.network.n)
        finite_difference_error = NaN
    end

    poisson = graph_poisson_diagnostics(
        solved.network,
        ell,
        active.line_indices,
    )

    nominal_emissions = operating_emissions(
        solved.solution.pg,
        solved.factors,
        solved.base_mva,
    )
    closest_index = isempty(kappas) ? 0 : argmin(kappas)
    closest_label = closest_index == 0 ? "none" : region.labels[closest_index]
    closest_kappa = closest_index == 0 ? Inf : kappas[closest_index]

    return (;
        case_name=splitext(basename(path))[1],
        path=abspath(path),
        solved,
        ptdf,
        active,
        fixed,
        dispatch_sensitivity,
        ell,
        ell_fixed,
        ell_fd,
        fd_columns,
        features,
        poisson,
        theta,
        Sigma,
        standard_deviation,
        region,
        emissions_variance,
        exit_probability_bound,
        facet_contributions,
        kappas,
        empirical_exit_probability,
        nominal_emissions,
        closest_label,
        closest_kappa,
        derivative_error,
        finite_difference_error,
        span_residual,
    )
end

"""
Search a fixed, declared grid of line limit scales and return the first feasible
operating point with at least `minimum_active_lines` binding lines.
"""
function find_operating_point(
    path::AbstractString;
    line_limit_scales=(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
    load_scales=(1.0, 1.1, 1.2),
    minimum_active_lines::Int=1,
    maximum_active_lines::Int=12,
    cost_ridge::Float64=1e-4,
)
    attempts = NamedTuple[]
    for load_scale in load_scales, line_limit_scale in line_limit_scales
        try
            solved = _solve_case(
                path;
                load_scale,
                line_limit_scale,
                cost_ridge,
                flow_regularization=0.0,
            )
            active = _active_constraints(solved.network, solved.solution)
            shedding = maximum(solved.solution.psh)
            push!(attempts, (;
                load_scale,
                line_limit_scale,
                active_lines=length(active.line_indices),
                shedding,
                feasible=shedding <= 1e-7,
            ))
            if shedding <= 1e-7 &&
               minimum_active_lines <= length(active.line_indices) <= maximum_active_lines
                return (;
                    load_scale,
                    line_limit_scale,
                    attempts=DataFrame(attempts),
                )
            end
        catch error
            push!(attempts, (;
                load_scale,
                line_limit_scale,
                active_lines=-1,
                shedding=Inf,
                feasible=false,
            ))
            @debug "operating point failed" path load_scale line_limit_scale error
        end
    end
    throw(ArgumentError(
        "no operating point met the active line target; attempts=$(DataFrame(attempts))",
    ))
end

end
