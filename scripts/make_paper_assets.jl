using CSV
using DataFrames
using Printf
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const AGGREGATE_DIR = joinpath(ROOT, "results", "aggregate")
const TABLE_DIR = joinpath(ROOT, "results", "tables")
mkpath(TABLE_DIR)

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

experiments = CSV.read(
    joinpath(AGGREGATE_DIR, "experiment_efficiency.csv"),
    DataFrame,
)
for policy in unique(experiments.policy)
    CSV.write(
        joinpath(AGGREGATE_DIR, "experiment_" * policy * ".csv"),
        experiments[experiments.policy .== policy, :],
    )
end

basis = CSV.read(joinpath(AGGREGATE_DIR, "basis_ablation.csv"), DataFrame)
for basis_name in unique(basis.basis)
    slug = replace(lowercase(basis_name), " " => "_")
    CSV.write(
        joinpath(AGGREGATE_DIR, "basis_" * slug * ".csv"),
        basis[basis.basis .== basis_name, :],
    )
end

cases = CSV.read(joinpath(AGGREGATE_DIR, "case_summary.csv"), DataFrame)
derivatives = CSV.read(
    joinpath(AGGREGATE_DIR, "derivative_validation.csv"),
    DataFrame,
)
included_cases = cases[cases.included, :]
included_cases.case_short = short_label.(included_cases.case_name)
CSV.write(joinpath(AGGREGATE_DIR, "scaling_plot.csv"), included_cases)

switching_path = CSV.read(
    joinpath(AGGREGATE_DIR, "switching_curve.csv"),
    DataFrame,
)
switching_path.plot_empirical = max.(
    switching_path.empirical_switch,
    0.5 ./ switching_path.samples,
)
switching_path.plot_bound = max.(
    switching_path.switch_bound,
    0.5 ./ switching_path.samples,
)
CSV.write(
    joinpath(AGGREGATE_DIR, "switching_curve_plot.csv"),
    switching_path,
)

sources = CSV.read(joinpath(AGGREGATE_DIR, "poisson_source.csv"), DataFrame)
CSV.write(
    joinpath(AGGREGATE_DIR, "poisson_binding_endpoints.csv"),
    sources[sources.binding_endpoint, :],
)
CSV.write(
    joinpath(AGGREGATE_DIR, "poisson_other_buses.csv"),
    sources[.!sources.binding_endpoint, :],
)

selected_labels = Set(["14", "118", "300", "588", "1354"])
selected = included_cases[in.(included_cases.case_short, Ref(selected_labels)), :]
sort!(selected, :buses)
open(joinpath(TABLE_DIR, "case_summary.tex"), "w") do io
    println(io, "\\begin{tabular}{lrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "System & buses \$n\$ & binding lines \$q\$ & rank \$r\$ & \$n/r\$ \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(selected)
        @printf(
            io,
            "%s & %d & %d & %d & %.1f \\\\\n",
            row.case_short,
            row.buses,
            row.active_lines,
            row.feature_rank,
            row.experiment_reduction,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

open(joinpath(TABLE_DIR, "full_case_summary.tex"), "w") do io
    println(io, "\\begin{tabular}{lrrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "System & \$n\$ & \$q\$ & \$r\$ & \$n/r\$ & " *
        "\$E_{\\mathrm{span}}\$ & \$E_{\\mathrm{P}}\$ \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(included_cases)
        validation = only(
            derivatives[derivatives.case_name .== row.case_name, :],
        )
        @printf(
            io,
            "%s & %d & %d & %d & %.1f & %.1e & %.1e \\\\\n",
            row.case_short,
            row.buses,
            row.active_lines,
            row.feature_rank,
            row.experiment_reduction,
            validation.congestion_span_residual,
            validation.poisson_residual,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

excluded = cases[.!cases.included, :]
open(joinpath(TABLE_DIR, "excluded_cases.tex"), "w") do io
    println(io, "\\begin{tabular}{ll}")
    println(io, "\\toprule")
    println(io, "System & exclusion reason \\\\")
    println(io, "\\midrule")
    for row in eachrow(excluded)
        @printf(
            io,
            "%s & %s \\\\\n",
            short_label(row.case_name),
            row.exclusion_reason,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

congestion_sweep = CSV.read(
    joinpath(AGGREGATE_DIR, "congestion_sweep.csv"),
    DataFrame,
)
open(joinpath(TABLE_DIR, "congestion_sweep.tex"), "w") do io
    println(io, "\\begin{tabular}{rrrr}")
    println(io, "\\toprule")
    println(io, "line limit multiplier & binding lines \$q\$ & rank \$r\$ & \$n/r\$ \\\\")
    println(io, "\\midrule")
    for row in eachrow(congestion_sweep)
        @printf(
            io,
            "%.1f & %d & %d & %.1f \\\\\n",
            row.line_limit_scale,
            row.active_lines,
            row.feature_rank,
            row.reduction,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

step_study = CSV.read(
    joinpath(AGGREGATE_DIR, "perturbation_step.csv"),
    DataFrame,
)
open(joinpath(TABLE_DIR, "perturbation_step.tex"), "w") do io
    println(io, "\\begin{tabular}{rlrr}")
    println(io, "\\toprule")
    println(
        io,
        "step \$h\$ (MW) & inside region & facet use & recovery error \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(step_study)
        error_text = row.solved_perturbations ?
            @sprintf("%.1e", row.relative_field_error) : "infeasible"
        @printf(
            io,
            "%.1f & %s & %.2f & %s \\\\\n",
            row.step_mw,
            row.inside_region ? "yes" : "no",
            row.worst_facet_use,
            error_text,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

binding_study = CSV.read(
    joinpath(AGGREGATE_DIR, "binding_set_study.csv"),
    DataFrame,
)
binding_summary = combine(
    groupby(binding_study, :declared_set),
    nrow => :variants,
    :feature_rank => minimum => :min_rank,
    :feature_rank => maximum => :max_rank,
    :recovery_error => median => :median_recovery_error,
    :recovery_error => maximum => :max_recovery_error,
)
open(joinpath(TABLE_DIR, "binding_set_study.tex"), "w") do io
    println(io, "\\begin{tabular}{lrlrr}")
    println(io, "\\toprule")
    println(
        io,
        "declared binding set & variants & rank \$r\$ & " *
        "median error & worst error \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(binding_summary)
        rank_text = row.min_rank == row.max_rank ?
            string(row.min_rank) :
            string(row.min_rank) * " to " * string(row.max_rank)
        @printf(
            io,
            "%s & %d & %s & %.1e & %.1e \\\\\n",
            row.declared_set,
            row.variants,
            rank_text,
            row.median_recovery_error,
            row.max_recovery_error,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

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
open(joinpath(TABLE_DIR, "correlated_tail.tex"), "w") do io
    println(io, "\\begin{tabular}{rrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "\$\\rho\$ & bus error std. (\\% load) & \$\\kappa_{\\min}\$ & " *
        "exit freq. & \$\\varepsilon_{\\rm cr}\$ & threshold \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(correlated_tail)
        threshold = isfinite(row.selected_query_threshold_median) ?
            @sprintf("%.2f", row.selected_query_threshold_median) : "--"
        @printf(
            io,
            "%.1f & %.2f & %.2f & %.3f & %.3f & %s \\\\\n",
            row.demand_correlation,
            row.load_forecast_std_percent,
            row.closest_kappa,
            row.empirical_exit,
            row.exit_bound,
            threshold,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

ridge_optimum = CSV.read(
    joinpath(AGGREGATE_DIR, "ridge_optimum.csv"),
    DataFrame,
)
open(joinpath(TABLE_DIR, "ridge_optimum.tex"), "w") do io
    println(io, "\\begin{tabular}{rrr}")
    println(io, "\\toprule")
    println(
        io,
        "\$\\lambda\$ & \$\\overline v_t/v\$ & threshold inflation \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(ridge_optimum)
        @printf(
            io,
            "\$10^{%d}\$ & %.3f & %.3f \\\\\n",
            round(Int, log10(row.ridge)),
            row.median_risk_ratio,
            row.median_threshold_ratio,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

switch_tightness = CSV.read(
    joinpath(AGGREGATE_DIR, "switching.csv"),
    DataFrame,
)
geometry = CSV.read(joinpath(AGGREGATE_DIR, "region_geometry.csv"), DataFrame)
tightness = innerjoin(
    geometry,
    select(switch_tightness, :case_name, :switch_bound, :empirical_switch),
    on=:case_name,
)
sort!(tightness, :buses)
open(joinpath(TABLE_DIR, "exit_bound_tightness.tex"), "w") do io
    println(io, "\\begin{tabular}{lrrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "System & facets & \$\\kappa_{\\min}\$ & sub-G. & Gauss. & " *
        "exit freq. & conc. \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(tightness)
        empirical = ismissing(row.empirical_switch) ||
                    isnan(row.empirical_switch) ?
            "--" : @sprintf("%.3f", row.empirical_switch)
        @printf(
            io,
            "%s & %d & %.2f & %.3f & %.3f & %s & %.2f \\\\\n",
            row.case_short,
            row.facets,
            row.kappa_independent,
            row.switch_bound,
            row.exit_bound,
            empirical,
            row.facet_concentration,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

open(joinpath(TABLE_DIR, "region_geometry.tex"), "w") do io
    println(io, "\\begin{tabular}{lrrlrr}")
    println(io, "\\toprule")
    println(
        io,
        "System & \$\\kappa_{\\min}\$ & \$a_\\star\$ & nearest facet & " *
        "avail.\\ \\(\\rho{=}0\\) & avail.\\ \\(\\rho{=}0.5\\) \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(geometry)
        available(x) = x >= 3.99 ? "\$>\$4.0" : @sprintf("%.2f", x)
        @printf(
            io,
            "%s & %.2f & %.1f & %s & %s & %s \\\\\n",
            row.case_short,
            row.kappa_independent,
            row.aggregate_ratio,
            row.facet_kind,
            available(row.availability_percent_independent),
            available(row.availability_percent_correlated),
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

tail_results = CSV.read(
    joinpath(AGGREGATE_DIR, "emissions_tail_curve.csv"),
    DataFrame,
)
open(joinpath(TABLE_DIR, "emissions_tail_results.tex"), "w") do io
    println(io, "\\begin{tabular}{rrrrrrrr}")
    println(io, "\\toprule")
    println(
        io,
        "bus error std. (\\% load) & \$\\kappa_{\\min}\$ & exit freq. & " *
        "\$\\varepsilon_{\\rm cr}\$ & empirical \$q_{.95}\$ & " *
        "known \$v\$ \$z_{.95}\$ & selected perturb. \$z_{.95}\$ & solves \\\\",
    )
    println(io, "\\midrule")
    for row in eachrow(tail_results)
        selected_query_threshold =
            isfinite(row.selected_query_threshold_median) ?
            @sprintf("%.2f", row.selected_query_threshold_median) : "--"
        @printf(
            io,
            "%.2f & %.2f & %.3f & %.3f & %.2f & %.2f & %s & %d \\\\\n",
            row.load_forecast_std_percent,
            row.closest_kappa,
            row.empirical_exit,
            row.exit_bound,
            row.empirical_quantile,
            row.fixed_region_threshold,
            selected_query_threshold,
            row.dispatch_solves,
        )
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end
