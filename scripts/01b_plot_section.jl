# Step 1b: figure of the faulted cross-section (facies and permeability).
# Run from the project folder:  julia --project=. scripts/01b_plot_section.jl

include(joinpath(@__DIR__, "..", "src", "section.jl"))
using CairoMakie

sec = build_section()
mesh = sec.mesh
mkpath(joinpath(@__DIR__, "..", "figures"))

labels = [[F.name for F in FACIES]..., "Fault zone"]
colors = [:khaki, :gray40, :goldenrod, :gray20, :firebrick]

# A 2D section viewed from the side: Axis3 looking along y, with lighting off so
# colours are not darkened and the y axis hidden.
section_axis(pos, title) = Axis3(pos, title = title,
    zreversed = true, azimuth = 1.5π, elevation = 0.0, aspect = (1.0, 0.05, 0.25),
    xlabel = "x (m)", zlabel = "Depth (m)", protrusions = (60, 10, 10, 40),
    ylabelvisible = false, yticklabelsvisible = false, yticksvisible = false,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)

fig = Figure(size = (1100, 650))
ax1 = section_axis(fig[1, 1], "Facies (hanging wall downthrown 30 m)")
core = findall(sec.core)  # leave out the wide padding columns
plot_cell_data!(ax1, mesh, sec.facies, cells = core, colormap = cgrad(colors, categorical = true), colorrange = (0.5, 5.5), shading = NoShading)
Legend(fig[1, 2], [PolyElement(color = c) for c in colors], labels)

ax2 = section_axis(fig[2, 1], "Horizontal permeability, log10(mD)")
plt = plot_cell_data!(ax2, mesh, log10.(sec.perm[1, :] ./ mD), cells = core, colormap = :viridis, shading = NoShading)
Colorbar(fig[2, 2], plt)

out = joinpath(@__DIR__, "..", "figures", "01_section.png")
save(out, fig)
println("Saved ", out)
