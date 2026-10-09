# Step 1: build the faulted cross-section and check it before running any flow.
# Run from the project folder:  julia --project=. scripts/01_build_section.jl

include(joinpath(@__DIR__, "..", "src", "section.jl"))

spec = section_spec()
sec = build_section(spec)
mesh = sec.mesh

println("Cells: ", number_of_cells(mesh), "  (", sec.nxt, " × 1 × ", sec.nz, ", including ",
    spec.pad_cells, " padding columns per side)")
println("Core: ", spec.Lx, " m wide, cells ", spec.dx, " m × ", spec.Ly, " m × ", spec.dz, " m")
println("Throw: ", spec.throw, " m,  dip: ", spec.dip_deg, "°")
println()

print_section(sec, well_x = 345.0)
println()

println("Footwall storage sand juxtaposed against (thickness, m):")
for (name, h) in sort(collect(juxtaposition_table(sec)), by = last, rev = true)
    println("  ", rpad(name, 20), h)
end
println()

# Pore volume by facies in the core: how much room each unit has for CO2.
domain = reservoir_domain(mesh, permeability = sec.perm, porosity = sec.poro)
pv = domain[:volumes] .* sec.poro
names = [[F.name for F in FACIES]..., "Fault zone"]
println("Pore volume by facies in the 1.2 km core (m³):")
for f in 1:FAULT_FACIES
    println("  ", rpad(names[f], 20), round(sum(pv[(sec.facies .== f) .& sec.core]), sigdigits = 3))
end

# Depth check: the tilt should make the right side shallower.
z = domain[:cell_centroids][3, sec.core]
println()
println("Core: shallowest cell centre ", round(minimum(z), digits = 1), " m, deepest ", round(maximum(z), digits = 1), " m")
