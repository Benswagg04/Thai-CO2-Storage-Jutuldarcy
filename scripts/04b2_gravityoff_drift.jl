# Quantify the water redistribution in the gravity-off no-well test (stage B2a):
# total water inventory, volume of water that moved across the initial
# contact, and the shift of the gas–water contact.
# Run from the project folder:  julia --project=. scripts/04b2_gravityoff_drift.jl

include(joinpath(@__DIR__, "04b2_gas_2d.jl"))

m = build_model(gravity = false)
s0 = initial_state(m, gravity = false)
off = setup_reservoir_forces(m.model, control = Dict(:Well => DisabledControl()))
_, s_init, _ = simulate_reservoir(s0, m.model, [1.0], parameters = m.parameters, forces = off, info_level = -1)
ws, states, t = run_schedule(m, s0, [(10, off)])

pv = m.domain[:volumes] .* m.domain[:porosity]
sw0 = s0[:Reservoir][:ImmiscibleSaturation]
sw1 = states[end][:Saturations][1, :]
z = m.domain[:cell_centroids][3, :]
iw = 3   # water component index (CO2, C1, Water)
W0 = sum(s_init[1][:TotalMasses][iw, :]); W1 = sum(states[end][:TotalMasses][iw, :])
above = z .< GWC
moved_up = sum(max.(sw1[above] .- sw0[above], 0) .* pv[above])
@printf("Total water mass: %.3f t → %.3f t (change %.2e of inventory)\n", W0 / 1000, W1 / 1000, (W1 - W0) / W0)
@printf("Water volume gained by cells above the initial contact: %.3f m³ (%.2e of water pore volume)\n",
    moved_up, moved_up / sum(sw0 .* pv))
@printf("Max |ΔSw| %.2e; cells with |ΔSw| > 1e-4: %d of %d\n", maximum(abs.(sw1 .- sw0)), count(>(1e-4), abs.(sw1 .- sw0)), length(sw0))
# Contact position: per column, depth where Sw crosses 0.5 (linear interpolation).
function contact_depths(sw)
    out = Float64[]
    for I in 1:NX
        cells = [c for c in eachindex(sw) if cell_ijk(m.mesh, c)[1] == I]
        sort!(cells, by = c -> z[c])
        k = findfirst(c -> sw[c] >= 0.5, cells)
        (isnothing(k) || k == 1) && continue
        a, b = cells[k-1], cells[k]
        push!(out, z[a] + (z[b] - z[a]) * (0.5 - sw[a]) / (sw[b] - sw[a]))
    end
    return out
end
d0, d1 = contact_depths(sw0), contact_depths(sw1)
@printf("Gas–water contact (Sw = 0.5 crossing): max shift over %d columns %.3f m\n", length(d0), maximum(abs.(d1 .- d0)))
