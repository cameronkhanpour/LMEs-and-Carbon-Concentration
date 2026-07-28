# Data Provenance

## PGLib-OPF network cases

The numerical pipeline reads the following MATPOWER cases:

- `meshed/pglib_opf_case14_ieee.m`
- `meshed/pglib_opf_case30_ieee.m`
- `meshed/pglib_opf_case57_ieee.m`
- `meshed/pglib_opf_case73_ieee_rts.m`
- `meshed/pglib_opf_case118_ieee.m`
- `meshed/pglib_opf_case162_ieee_dtc.m`
- `meshed/pglib_opf_case179_goc.m`
- `meshed/pglib_opf_case197_snem.m`
- `meshed/pglib_opf_case240_pserc.m`
- `meshed/pglib_opf_case300_ieee.m`
- `meshed/pglib_opf_case500_goc.m`
- `meshed/pglib_opf_case588_sdet.m`
- `meshed/pglib_opf_case793_goc.m`
- `meshed/pglib_opf_case1354_pegase.m`

Their headers identify them as PGLib-OPF version 23.07 benchmark cases from
<https://github.com/power-grid-lib/pglib-opf>. The case headers identify the
underlying systems and state a Creative Commons Attribution 4.0 license. The
analysis does not rewrite the case files. PowerIO parses the in service network
records and normalizes power quantities to per unit.

The prespecified inclusion checks retain cases 14, 30, 57, 118, 162, 179, 240,
300, 588, and 1354. Cases 73, 197, 500, and 793 are retained in the screening
output but excluded from subsequent comparisons because their selected
operating points have zero direct operating emissions under the factors below.

## Operating emission factors

Generator comments in the PGLib files provide fuel labels such as `COW`, `NG`,
`PEL`, and `SYNC`. `src/OperatingEmissionsCertificates.jl` maps those labels to
the direct CO2 intensity lookup distributed in:

- PGLib-CO2: <https://github.com/jacobcho0103/PGLib-CO2>
- inspected revision: `f5954d677ff8ec7744b3b80a6292a2d044429e6f`
- inspection date: 2026-07-23

The factors are metric tonnes of CO2 per MWh. They represent operating
emissions only. The implementation assigns zero to synchronous condensers,
nuclear, hydro, wind, and solar resources, consistent with that direct
emissions boundary. These values are not life cycle greenhouse gas
intensities.
