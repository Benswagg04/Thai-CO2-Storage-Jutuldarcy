# Step 3b: tables for the verification report, from output/03b/*.jls.
# Run from the project folder:  julia --project=. scripts/03b_report.jl
# Writes results/03b_tables.md (figures: scripts/make_figures.jl).

include(joinpath(@__DIR__, "..", "src", "case.jl"))
using Serialization, Printf

dir = joinpath(@__DIR__, "..", "output", "03b")
load(n) = (f = joinpath(dir, "$n.jls"); isfile(f) ? deserialize(f) : nothing)
at(r, y) = findfirst(==(y), r.years)
io = IOBuffer()
md(s) = (println(io, s); println(s))
f1(x) = @sprintf("%.1f", x); f2(x) = @sprintf("%.2f", x); f3(x) = @sprintf("%.3f", x)

# 1. Boundary transport and conservation residual.
md("## Boundary transport and CO2 conservation (tonnes)\n")
md("| Case | Year | Injected | In reservoir | In well | Out via top (free / dissolved) | Out via sides | Residual | Max abs residual | Ministeps (cut) | Newton its | Wall (min) |")
md("|---|---|---|---|---|---|---|---|---|---|---|---|")
for n in ("flux_pe5.5_leaky", "flux_pe5.5_leaky_tight", "flux_pe1.2_mid", "flux_pe1.2_mid_tight")
    r = load(n); isnothing(r) && continue
    for y in (25, 100)
        i = at(r, y)
        md("| $n | $y | $(f2(r.injected[i])) | $(f2(r.in_res[i])) | $(f3(r.well[i])) | $(f3(r.top_free[i])) / $(f3(r.top_diss[i])) | " *
           "$(f3(r.side_free[i] + r.side_diss[i])) | $(f3(r.residual[i])) | $(f3(r.max_abs_residual)) | " *
           "$(r.stats.ministeps_ok) ($(r.stats.ministeps_failed)) | $(r.stats.newton) | $(f1(r.stats.wall_min)) |")
    end
end

# 2. Fault invasion vs crossing the seal.
md("\n## Fault invasion vs crossing the seal (tonnes)\n")
md("| Case | Year | Free CO2 in fault | Dissolved CO2 in fault | Free in overburden | Dissolved in overburden | Cumulative net crossing of the seal | Net out of storage formation |")
md("|---|---|---|---|---|---|---|---|")
for n in ("flux_pe5.5_leaky", "flux_pe1.2_mid")
    r = load(n); isnothing(r) && continue
    for y in (25, 100)
        i = at(r, y)
        md("| $n | $y | $(f1(r.fault_free[i])) | $(f1(r.fault_diss[i])) | $(f1(r.beyond_free[i])) | $(f1(r.beyond_diss[i])) | $(f1(r.crossed_seal[i])) | $(f1(r.out_of_storage[i])) |")
    end
end

# 3. Injection rate at fixed total mass.
md("\n## Injection rate at fixed total mass (37,500 t), fault Pe 1.4 bar, k 0.1/1 mD\n")
md("| Rate (t/yr) | Injected by simulator (t) | Peak entry drive (bar) | Peak Δp (bar) | Out of storage: end of injection / +25 yr / yr 100 (t) | Crossed seal: end / +25 / 100 (t) |")
md("|---|---|---|---|---|---|")
rates = [(750, 50), (1500, 25), (2500, 15)]
for (rate, yrs_inj) in rates
    r = load("rate$(rate)_pe1.4_mid"); isnothing(r) && continue
    a, b, c = at(r, yrs_inj), at(r, yrs_inj + 25), at(r, 100)
    d = filter(isfinite, r.drive)
    md("| $rate × $yrs_inj yr | $(f1(r.injected[end])) | $(f2(maximum(d))) | $(f2(maximum(r.dp_max))) | " *
       "$(f1(r.out_of_storage[a])) / $(f1(r.out_of_storage[b])) / $(f1(r.out_of_storage[c])) | " *
       "$(f1(r.crossed_seal[a])) / $(f1(r.crossed_seal[b])) / $(f1(r.crossed_seal[c])) |")
end

# 4. Time step and grid.
md("\n## Time-step and grid sensitivity (fault k 0.1/1 mD)\n")
md("| Case | Cells | Out of storage, yr 25 / 100 (t) | Crossed seal, yr 25 / 100 (t) | Free CO2 in fault, yr 100 (t) | Max abs residual (t) | Ministeps (cut) | Wall (min) |")
md("|---|---|---|---|---|---|---|---|")
grid_cases = ["dt90_pe1.2_mid", "dt30_pe1.2_mid", "grid10x2_pe1.2_mid", "grid10x1_pe1.2_mid",
              "grid5x2_pe1.2_mid", "grid5x1_pe1.2_mid", "grid10x2_pe0.8_mid", "grid5x1_pe0.8_mid"]
for n in grid_cases
    r = load(n); isnothing(r) && continue
    sec = build_section(r.spec)
    md("| $n | $(number_of_cells(sec.mesh)) | $(f1(r.out_of_storage[at(r, 25)])) / $(f1(r.out_of_storage[end])) | " *
       "$(f1(r.crossed_seal[at(r, 25)])) / $(f1(r.crossed_seal[end])) | $(f1(r.fault_free[end])) | " *
       "$(f3(r.max_abs_residual)) | $(r.stats.ministeps_ok) ($(r.stats.ministeps_failed)) | $(f1(r.stats.wall_min)) |")
end
mkpath(joinpath(@__DIR__, "..", "results"))
write(joinpath(@__DIR__, "..", "results", "03b_tables.md"), String(take!(io)))

println("
Wrote results/03b_tables.md (figures: julia --project=. scripts/make_figures.jl 03b)")
