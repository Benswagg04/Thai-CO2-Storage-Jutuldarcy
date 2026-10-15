# Step 4C: verify gas–water capillary pressure in the three-phase
# compositional model (immiscible brine + CO2–CH4 Peng–Robinson), JutulDarcy 0.3.10.
#
# In this system the reference pressure is the hydrocarbon-LIQUID phase (index
# 2), which never exists here. JutulDarcy computes phase pressures as
#   p_water = p + Pc_ow(Sw),   p_gas = p + Pc_og(Sg).
# For a gas–water system we set Pc_ow ≡ 0 and Pc_og(Sg) = Pc_gw(Sw = 1 − Sg),
# so p_gas − p_water = Pc_gw (Brooks–Corey, entry pressure Pe at Sw = 1).
#
# Test 1 (slope): a 60 m column, sharp gas-over-brine contact, relaxes for
#   100 years. At capillary–gravity equilibrium Pc_gw(Sw) must decrease linearly
#   with depth at slope (ρ_brine − ρ_gas)·g, and gas must be absent where Pc < Pe.
# Test 2 (barrier): a 4 m layer with Pe = 1 bar above a gas column of 10 m
#   (buoyancy ≈ 0.9 bar: must hold) or 20 m (≈ 1.8 bar: must break through).
# Run from the project folder:  julia --project=. scripts/04c_pc_check.jl

include(joinpath(@__DIR__, "..", "src", "gas.jl"))
using Printf

const mD = si_unit(:darcy) / 1000
const GS = gas_spec()
const T_RES = convert_to_si(GS.T_C, :Celsius)
const LAMBDA = 2.0

"Brooks–Corey gas–water capillary pressure as a function of Sw (Pa)."
pc_gw(sw, pe_bar; swr = GS.swc) = pe_bar * bar * max((sw - swr) / (1 - swr), 1e-3)^(-1 / LAMBDA)

"Pc tables for (Pc_ow, Pc_og) per region, as JutulDarcy expects for three phases."
function pc_tables(pe_bar_per_region)
    sg = vcat(collect(range(0.0, 0.7, 71)), collect(range(0.7, 1.0, 151))[2:end])
    zero_tab = Jutul.get_1d_interpolator([0.0, 1.0], [0.0, 0.0])
    pcow = Tuple(zero_tab for _ in pe_bar_per_region)
    pcog = Tuple(Jutul.get_1d_interpolator(sg, pc_gw.(1 .- sg, pe)) for pe in pe_bar_per_region)
    return (pcow, pcog)
end

"Vertical column, 1 m cells; region 1 = sand, region 2 = barrier (cells in `barrier`)."
function column_model(; nz = 60, pe_sand = 0.2, pe_barrier = 1.0, barrier = Int[])
    mesh = reservoir_mesh((1, 1, nz), (10.0, 10.0, Float64(nz)), origin = [0.0, 0.0, 1500.0])
    domain = reservoir_domain(mesh, permeability = 100mD, porosity = 0.2, temperature = T_RES)
    sys = gas_system(GS)
    model = setup_reservoir_model(domain, sys)
    set_brine_properties!(model, GS)
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0, 2.0], [GS.swc, 0.0, 0.0]))
    regions = ones(Int, nz); regions[barrier] .= 2
    pc = JutulDarcy.SimpleCapillaryPressure(pc_tables([pe_sand, pe_barrier]), regions = regions)
    set_secondary_variables!(reservoir_model(model), CapillaryPressure = pc)
    push!(reservoir_model(model).output_variables, :TotalMasses, :Saturations, :PhaseMassDensities)
    return (; mesh, domain, model, regions, pe = [pe_sand, pe_barrier])
end

"Sharp initial state: gas (Sw = swc) in `gas_cells`, brine elsewhere; brine-hydrostatic pressure."
function sharp_state(m, gas_cells)
    nz = number_of_cells(m.mesh)
    z = m.domain[:cell_centroids][3, :]
    p = 152bar .+ 990.0 * 9.81 .* (z .- 1520.0)
    sw = ones(nz); sw[gas_cells] .= GS.swc
    return setup_reservoir_state(m.model, Pressure = p, OverallMoleFractions = [0.0, 1.0], ImmiscibleSaturation = sw)
end

"""
Gas column in capillary–gravity equilibrium: the free-water level (Pc = 0) at
the base of `gas_cells`; at height h above it, Pc = Δρ g h and Sw follows the
Brooks–Corey inverse, Sw = Swr + (1 − Swr)(Pe/(Δρ g h))^λ (Sw = 1 where Δρ g h < Pe).
Starting gas at Sw = Swr instead would mean Pc at its maximum (the table cap,
≈ 32 Pe), i.e. a large artificial capillary overpressure.
"""
function equilibrium_column_state(m, gas_cells; drho = 900.0)
    nz = number_of_cells(m.mesh)
    z = m.domain[:cell_centroids][3, :]
    p = 152bar .+ 990.0 * 9.81 .* (z .- 1520.0)
    z_fwl = maximum(z[gas_cells]) + 0.5               # base of the deepest gas cell
    sw = ones(nz)
    for c in gas_cells
        pc = drho * 9.81 * (z_fwl - z[c])
        pe = m.pe[m.regions[c]] * bar
        sw[c] = pc <= pe ? 1.0 : GS.swc + (1 - GS.swc) * (pe / pc)^LAMBDA
    end
    return setup_reservoir_state(m.model, Pressure = p, OverallMoleFractions = [0.0, 1.0], ImmiscibleSaturation = sw)
end

relax(m, s0, years) = simulate_reservoir(s0, m.model, fill(1.0yr, years), info_level = -1,
    max_timestep = 30day, relaxation = true, timestep_max_increase = 2.0)

println("== Test 1: capillary–gravity equilibrium slope (Pe = 0.2 bar, 60 m column, top 25 m gas)")
let
    m = column_model()
    s0 = sharp_state(m, 1:25)
    ws, states, t = relax(m, s0, 100)
    s = states[end]
    z = m.domain[:cell_centroids][3, :]
    sw = s[:Saturations][1, :]
    pc = [pc_gw(sw[i], m.pe[m.regions[i]]) for i in eachindex(sw)]
    rho_w, rho_g = s[:PhaseMassDensities][1, :], s[:PhaseMassDensities][3, :]
    # Fit Pc vs depth over the transition zone (gas present and Sw not at its lower cap).
    tz = findall(i -> 0.25 < sw[i] < 0.99, eachindex(sw))
    A = hcat(ones(length(tz)), z[tz]); coef = A \ pc[tz]
    drho = sum(rho_w[tz] .- rho_g[tz]) / length(tz)
    @printf("  transition-zone cells: %d (depths %.1f–%.1f m)\n", length(tz), minimum(z[tz]), maximum(z[tz]))
    @printf("  fitted dPc/dz = %.1f Pa/m; expected −(ρw − ρg)g = %.1f Pa/m (Δρ = %.1f kg/m³) → ratio %.4f\n",
        coef[2], -drho * 9.81, drho, coef[2] / (-drho * 9.81))
    gas_cells = findall(<(0.999), sw)
    deepest = gas_cells[argmax(z[gas_cells])]
    @printf("  deepest gas-bearing cell: depth %.1f m, Sw = %.4f, Pc = %.3f bar (entry pressure %.2f bar)\n",
        z[deepest], sw[deepest], pc[deepest] / bar, m.pe[1])
    dsw_last = maximum(abs.(states[end][:Saturations][1, :] .- states[end-10][:Saturations][1, :]))
    @printf("  max |ΔSw| over the last 10 years: %.2e (equilibrium reached?)\n", dsw_last)
end

for (label, gas_top, expect) in (("10 m column", 31, "HOLD"), ("20 m column", 21, "BREAK THROUGH"))
    println("== Test 2: barrier (cells 27–30, Pe = 1.0 bar) above a $label (expected: $expect)")
    m = column_model(barrier = 27:30)
    s0 = equilibrium_column_state(m, gas_top:(gas_top + (label == "10 m column" ? 9 : 19)))
    z = m.domain[:cell_centroids][3, :]
    drho_est = 990.0 - 90.0
    H = label == "10 m column" ? 10.0 : 20.0
    @printf("  buoyancy of the column ≈ (ρw − ρg) g H = %.2f bar vs barrier Pe 1.00 bar\n", drho_est * 9.81 * H / bar)
    ws, states, t = relax(m, s0, 50)
    s = states[end]
    above = findall(<(27), 1:length(z))
    gas_above = sum(s[:TotalMasses][2, above]) / 1000
    gas_total = sum(s[:TotalMasses][2, :]) / 1000
    @printf("  CH4 above the barrier after 50 yr: %.3f t of %.3f t (%.2f%%)\n", gas_above, gas_total, 100 * gas_above / gas_total)
end
