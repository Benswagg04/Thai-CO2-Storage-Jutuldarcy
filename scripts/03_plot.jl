# Step 3 figures from output/03_sweep/*.jls (run scripts/03_sweep.jl first).
# Run from the project folder:  julia --project=. scripts/03_plot.jl

include(joinpath(@__DIR__, "..", "src", "case.jl"))
using CairoMakie, Serialization

outdir = joinpath(@__DIR__, "..", "output", "03_sweep")
figdir = joinpath(@__DIR__, "..", "figures")
pe_values = [0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0, 5.5]
k_levels = [(0.01, 0.1), (0.1, 1.0), (1.0, 10.0)]
runname(pe, k) = "pe$(pe)_ka$(k[1])_kl$(k[2])"
klabel(k) = "k = $(k[1]) / $(k[2]) mD (across / along)"

R = Dict{Tuple{Float64, Tuple{Float64, Float64}}, Any}()
for k in k_levels, pe in pe_values
    f = joinpath(outdir, "$(runname(pe, k)).jls")
    isfile(f) && (R[(pe, k)] = deserialize(f))
end
println("Loaded $(length(R)) of $(length(pe_values) * length(k_levels)) runs")

outside_pct(r) = begin
    inj = r.inventory[end][:injected] / 1000
    100 .* ([i[:injected] / 1000 for i in r.inventory] .-
            (r.acc.total[:fw_storage] .+ r.acc.total[:hw_storage])) ./ inj
end
colors = Makie.wong_colors()

# Entry drive with a sealing fault (highest Pe): what the fault has to withstand.
# Buoyancy alone after injection (year 100) and the peak at the end of injection.
sealed = [R[(5.5, k)] for k in k_levels if haskey(R, (5.5, k))]
drive_post = isempty(sealed) ? NaN : maximum(r.diag.drive[end] for r in sealed)
drive_peak = isempty(sealed) ? NaN : maximum(maximum(filter(isfinite, r.diag.drive)) for r in sealed)

# 1. Key figure: leakage vs fault entry pressure.
fig = Figure(size = (1150, 460))
for (col, (step, title)) in enumerate([(25, "End of injection (year 25)"), (100, "Year 100")])
    ax = Axis(fig[1, col], title = title, xscale = log10,
        xticks = [0.4, 0.6, 1.0, 1.7, 3.0, 5.5], xlabel = "Fault capillary entry pressure (bar)",
        ylabel = "CO₂ outside the storage formation (% of injected)")
    if isfinite(drive_post)
        vspan!(ax, drive_post, drive_peak, color = (:gray, 0.2))
        vlines!(ax, [drive_post, drive_peak], color = :gray, linestyle = :dash)
    end
    for (j, k) in enumerate(k_levels)
        pes = [pe for pe in pe_values if haskey(R, (pe, k))]
        isempty(pes) && continue
        y = [outside_pct(R[(pe, k)])[step] for pe in pes]
        scatterlines!(ax, pes, y, color = colors[j], label = klabel(k), linewidth = 2)
    end
    col == 1 && axislegend(ax, position = :rt, labelsize = 10)
end
Label(fig[2, 1:2], "Grey band: entry drive at the fault with a sealing fault, from buoyancy alone after injection " *
    "($(round(drive_post, digits = 2)) bar) to the end-of-injection peak ($(round(drive_peak, digits = 2)) bar).", fontsize = 11)
save(joinpath(figdir, "03_leakage_vs_entry_pressure.png"), fig)

# 2. Time series for the middle permeability level.
kmid = k_levels[2]
fig = Figure(size = (1150, 440))
ax1 = Axis(fig[1, 1], title = "CO₂ outside the storage formation, $(klabel(kmid))",
    xlabel = "Years", ylabel = "% of injected")
ax2 = Axis(fig[1, 2], title = "Entry drive at the fault (dashed: entry pressure)",
    xlabel = "Years", ylabel = "bar")
pes = [pe for pe in pe_values if haskey(R, (pe, kmid)) && pe <= 2.0]
cmap = cgrad(:viridis, max(length(pes), 2), categorical = true)
for (i, pe) in enumerate(pes)
    r = R[(pe, kmid)]
    lines!(ax1, r.t ./ yr, outside_pct(r), color = cmap[i], linewidth = 2, label = "Pe = $pe bar")
    lines!(ax2, r.t ./ yr, r.diag.drive, color = cmap[i], linewidth = 2)
    hlines!(ax2, [pe], color = cmap[i], linestyle = :dash, linewidth = 1)
end
vlines!(ax1, [25], color = :gray, linestyle = :dot); vlines!(ax2, [25], color = :gray, linestyle = :dot)
Legend(fig[1, 3], ax1, labelsize = 10)
save(joinpath(figdir, "03_time_series.png"), fig)

# 3. Final CO2 saturation for three entry pressures at the middle permeability.
section_axis(pos, title) = Axis3(pos, title = title,
    zreversed = true, azimuth = 1.5π, elevation = 0.0, aspect = (1.0, 0.05, 0.25),
    xlabel = "x (m)", zlabel = "Depth (m)", protrusions = (60, 10, 10, 40),
    ylabelvisible = false, yticklabelsvisible = false, yticksvisible = false,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)
picks = [pe for pe in (0.6, 1.2, 2.0) if haskey(R, (pe, kmid))]
if !isempty(picks)
    fig = Figure(size = (1100, 300 * length(picks)))
    for (row, pe) in enumerate(picks)
        r = R[(pe, kmid)]
        sec = build_section(r.spec)
        ax = section_axis(fig[row, 1], "Fault Pe = $pe bar, $(klabel(kmid)): year 100")
        plot_cell_data!(ax, sec.mesh, r.states[end][:Saturations][2, :], cells = findall(sec.core),
            colormap = :inferno, colorrange = (0, 0.8), shading = NoShading)
        plot_mesh_edges!(ax, sec.mesh, cells = findall(sec.facies .== FAULT_FACIES), color = :gray70)
    end
    Colorbar(fig[1:length(picks), 2], colormap = :inferno, colorrange = (0, 0.8), label = "S_CO₂")
    save(joinpath(figdir, "03_saturation_examples.png"), fig)
end
println("Saved figures/03_leakage_vs_entry_pressure.png, 03_time_series.png, 03_saturation_examples.png")
