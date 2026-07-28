using OperatingEmissionsCertificates
using CSV
using DataFrames

root = normpath(joinpath(@__DIR__, ".."))
case_names = [
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case30_ieee.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case162_ieee_dtc.m",
    "pglib_opf_case179_goc.m",
    "pglib_opf_case197_snem.m",
    "pglib_opf_case240_pserc.m",
    "pglib_opf_case300_ieee.m",
    "pglib_opf_case500_goc.m",
    "pglib_opf_case588_sdet.m",
    "pglib_opf_case793_goc.m",
    "pglib_opf_case1354_pegase.m",
]

rows = NamedTuple[]
for case_name in case_names
    path = joinpath(root, "data", "meshed", case_name)
    operating_point = find_operating_point(
        path;
        line_limit_scales=(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
        load_scales=(1.0, 1.1, 1.2),
        minimum_active_lines=1,
        maximum_active_lines=40,
    )
    push!(rows, (;
        case_name=splitext(case_name)[1],
        load_scale=operating_point.load_scale,
        line_limit_scale=operating_point.line_limit_scale,
        attempts=nrow(operating_point.attempts),
    ))
    display(operating_point.attempts)
end

output = DataFrame(rows)
mkpath(joinpath(root, "results", "aggregate"))
CSV.write(joinpath(root, "results", "aggregate", "operating_points.csv"), output)
display(output)
