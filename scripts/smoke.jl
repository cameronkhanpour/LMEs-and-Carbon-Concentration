using OperatingEmissionsCertificates

case_path = joinpath(@__DIR__, "..", "data", "meshed", "pglib_opf_case14_ieee.m")
operating_point = find_operating_point(
    case_path;
    line_limit_scales=(1.0, 0.8, 0.6, 0.5, 0.4, 0.3),
    load_scales=(1.0, 1.1, 1.2),
)
display(operating_point.attempts)
analysis = analyze_case(
    case_path;
    load_scale=operating_point.load_scale,
    line_limit_scale=operating_point.line_limit_scale,
    exit_probability_samples=1_000,
    finite_difference_count=5,
)
println((
    case_name=analysis.case_name,
    active_lines=length(analysis.active.line_indices),
    feature_rank=analysis.features.feature_rank,
    derivative_error=analysis.derivative_error,
    finite_difference_error=analysis.finite_difference_error,
    span_residual=analysis.span_residual,
    exit_probability_bound=analysis.exit_probability_bound,
    empirical_exit_probability=analysis.empirical_exit_probability,
))
