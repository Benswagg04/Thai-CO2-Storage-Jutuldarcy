# Step 2 figures from output/02_*.jls (run scripts/02_run_cases.jl first).
# Run from the project folder:  julia --project=. scripts/02_plot.jl

include(joinpath(@__DIR__, "..", "src", "case.jl"))
using CairoMakie, Serialization

outdir = joinpath(@__DIR__, "..", "output")
figdir = joinpath(@__DIR__, "..", "figures")
mkpath(figdir)
load(name) = deserialize(joinpath(outdir, "$name.jls"))

cases = [
    ("02_nopc_clay", "No Pc, clay-rich fault"),
    ("02_pc_clay", "Pc, clay-rich fault (Pe 5.5 bar)"),
    ("02_pc_sandy", "Pc, sand-rich fault (Pe 0.55 bar)"),
]
res = Dict(n => load(n) for (n, _) in cases)
fine = load("02_pc_sandy_fine")
years(r) = r.t ./ yr

section_axis(pos, title) = Axis3(pos, title = title,
    zreversed = true, azimuth = 1.5π, elevation = 0.0, aspect = (1.0, 0.05, 0.25),
    xlabel = "x (m)", zlabel = "Depth (m)", protrusions = (60, 10, 10, 40),
    ylabelvisible = false, yticklabelsvisible = false, yticksvisible = false,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)

# 1. CO2 saturation at year 100 (1.2 km core only).
fig = Figure(size = (1100, 900))
for (row, (n, label)) in enumerate(cases)
    r = res[n]
    sec = build_section(r.spec)
    ax = section_axis(fig[row, 1], "$label: year 100")
    plot_cell_data!(ax, sec.mesh, r.states[end][:Saturations][2, :], cells = findall(sec.core),
        colormap = :inferno, colorrange = (0, 0.8), shading = NoShading)
    plot_mesh_edges!(ax, sec.mesh, cells = findall(sec.facies .== FAULT_FACIES), color = :gray70)
end
Colorbar(fig[1:3, 2], colormap = :inferno, colorrange = (0, 0.8), label = "S_CO₂")
save(joinpath(figdir, "02_saturation_year100.png"), fig)

# 2. Where the CO2 outside the storage formation is at year 100, free vs dissolved.
cats = [(:fault, "Fault zone"), (:caprock, "Caprock"), (:overburden, "Overburden")]
fig = Figure(size = (1100, 420))
ax = Axis(fig[1, 1], title = "CO₂ outside the storage formation at year 100 (% of injected)",
    xticks = (1:4, [[l for (_, l) in cats]..., "Left model\n(through top)"]), ylabel = "% of injected")
w = 0.25
colors = Makie.wong_colors()
for (j, (n, label)) in enumerate(cases)
    r = res[n]
    inj = r.inventory[end][:injected] / 1000
    xs = (1:4) .+ (j - 2) * w
    free = [100 * r.acc.free[k][end] / inj for (k, _) in cats]
    diss = [100 * r.acc.diss[k][end] / inj for (k, _) in cats]
    push!(free, 100 * r.exited[end] / inj); push!(diss, 0.0)
    barplot!(ax, xs, free .+ diss, width = w, color = (colors[j], 0.35), strokecolor = colors[j], strokewidth = 1.5)
    barplot!(ax, xs, free, width = w, color = colors[j], label = label)
end
axislegend(ax, position = :lt, labelsize = 10)
Label(fig[2, 1], "Solid: free (CO₂-rich phase). Light: dissolved in brine. CO₂ that left the model is reported as one total.", fontsize = 11)
save(joinpath(figdir, "02_accounting_year100.png"), fig)

# 3. Pressure at the fault: entry drive vs entry pressure, and injection overpressure.
fig = Figure(size = (1100, 420))
ax1 = Axis(fig[1, 1], title = "CO₂ entry drive at the fault", xlabel = "Years",
    ylabel = "p_CO₂(sand) − p_brine(fault) (bar)")
for (j, (n, label)) in enumerate(cases[2:3])
    r = res[n]
    sec = build_section(r.spec)
    pe = entry_pressures(sec, r.cs)[FAULT_FACIES]
    lines!(ax1, years(r), r.diag.drive, color = colors[j + 1], linewidth = 2, label = label)
    hlines!(ax1, [pe], color = colors[j + 1], linestyle = :dash)
end
vlines!(ax1, [25], color = :gray, linestyle = :dot)
axislegend(ax1, position = :rc, labelsize = 10)
text!(ax1, 60, 5.0, text = "dashed: entry pressure", fontsize = 11)
ax2 = Axis(fig[1, 2], title = "Pressure rise from injection", xlabel = "Years", ylabel = "Δp (bar)")
for (j, (n, label)) in enumerate(cases)
    lines!(ax2, years(res[n]), res[n].diag.dp_max, color = colors[j], linewidth = 2, label = "$label: max anywhere")
end
lines!(ax2, years(res["02_pc_clay"]), res["02_pc_clay"].diag.overpressure_fault, color = :black,
    linestyle = :dash, label = "Pc, clay-rich: at the fault")
axislegend(ax2, position = :rt, labelsize = 9)
save(joinpath(figdir, "02_fault_pressure.png"), fig)

# 4. Grid refinement: sand-rich fault on the 10 m × 2 m and 5 m × 1 m grids.
fig = Figure(size = (1100, 420))
outside(r) = begin
    inj = [i[:injected] / 1000 for i in r.inventory]
    stor = r.acc.total[:fw_storage] .+ r.acc.total[:hw_storage]
    100 .* (inj .- stor) ./ max.(inj, 1e-9)
end
ax1 = Axis(fig[1, 1], title = "Sand-rich fault: CO₂ outside the storage formation",
    xlabel = "Years", ylabel = "% of injected")
coarse = res["02_pc_sandy"]
lines!(ax1, years(coarse), outside(coarse), linewidth = 2, label = "10 m × 2 m grid")
lines!(ax1, years(fine), outside(fine), linewidth = 2, linestyle = :dash, label = "5 m × 1 m grid")
axislegend(ax1, position = :lt)
ax2 = Axis(fig[1, 2], title = "Sand-rich fault: entry drive", xlabel = "Years", ylabel = "bar")
lines!(ax2, years(coarse), coarse.diag.drive, linewidth = 2, label = "10 m × 2 m grid")
lines!(ax2, years(fine), fine.diag.drive, linewidth = 2, linestyle = :dash, label = "5 m × 1 m grid")
hlines!(ax2, [0.548], color = :gray, linestyle = :dot)
axislegend(ax2, position = :rb)
save(joinpath(figdir, "02_grid_refinement.png"), fig)
println("Saved figures/02_saturation_year100.png, 02_accounting_year100.png, 02_fault_pressure.png, 02_grid_refinement.png")
