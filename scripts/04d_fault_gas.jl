# Stage B2c: faulted depleted gas field.
#   equilibrium : capillary–gravity-consistent methane trap; 10 years without
#                 wells, with methane-based (scale 1) and CO2-scaled (scale 0.6)
#                 capillary curves. Must stay at rest; if the CO2-scaled curves
#                 cannot hold the initial methane column, that is flagged.
# Run from the project folder:  julia --project=. scripts/04d_fault_gas.jl equilibrium

include(joinpath(@__DIR__, "..", "src", "gasfault.jl"))
using Serialization

const GS = gas_spec()
const FS = fault_gas_spec()

function equilibrium_test(scale; years = 10)
    sec = build_section()
    m = fault_gas_model(sec, GS, FS; pe_scale = scale)
    s0 = fault_gas_initial_state(sec, m, GS, FS)
    off = setup_reservoir_forces(m.model, control = Dict(:Well => DisabledControl()))
    _, s_init, _ = simulate_reservoir(s0, m.model, [1.0], parameters = m.parameters, forces = off, info_level = -1)
    t0 = time()
    res = simulate_reservoir(s0, m.model, fill(1.0yr, years); parameters = m.parameters, forces = off,
        info_level = -1, SOLVER_KW...)
    ws, states, t = res
    wall = (time() - t0) / 60
    mini = [mm for r in res.result.reports for mm in r[:ministeps]]

    z = m.domain[:cell_centroids][3, :]
    trap = trap_cells(sec, m.domain, FS)
    beyond = fault_gas_beyond_seal(sec)
    fault = findall(sec.facies .== FAULT_FACIES)
    r0 = s0[:Reservoir]
    sw0 = r0[:ImmiscibleSaturation]; sw1 = states[end][:Saturations][1, :]
    M0 = s_init[1][:TotalMasses]; M1 = states[end][:TotalMasses]
    ch4(M, cells) = sum(M[2, cells]) / 1000
    # Column height and buoyancy at the fault: deepest trap cell next to the fault.
    x_sand = maximum(sec.xc[c] for c in trap)
    col = [c for c in trap if sec.xc[c] == x_sand]
    h_fault = FS.z_fwl - minimum(z[col]) + 1.0   # cell centre to top of cell (dz/2 = 1 m)
    rg = gas_density(gas_eos(GS), r0[:Pressure][col[1]], T_GAS, [0.0, 1.0])
    buoy = (GS.brine_density - rg) * 9.81 * h_fault / bar

    println("== Pc scale $scale: entry pressures (bar) ", join(round.(m.pe, sigdigits = 3), ", "),
        " [overburden, caprock, storage sand, underburden, fault]")
    @printf("   trap: %d cells; gas column at the fault ≈ %.1f m → buoyancy ≈ %.2f bar vs fault Pe %.2f bar (ratio %.2f)\n",
        length(trap), h_fault, buoy, m.pe[FAULT_FACIES], m.pe[FAULT_FACIES] / buoy)
    @printf("   initial CH4: total %.1f t, in trap %.1f t; initial p %.1f–%.1f bar\n",
        ch4(M0, 1:m.nc), ch4(M0, trap), minimum(r0[:Pressure]) / bar, maximum(r0[:Pressure]) / bar)
    @printf("   after %d yr without wells (%.1f min, %d ok / %d cut ministeps):\n", years, wall,
        count(mm -> mm[:success], mini), count(mm -> !mm[:success], mini))
    @printf("     max |Δp| %.2e bar, max |ΔSw| %.2e, CH4 total change %.2e t\n",
        maximum(abs.(states[end][:Pressure] .- r0[:Pressure])) / bar, maximum(abs.(sw1 .- sw0)),
        ch4(M1, 1:m.nc) - ch4(M0, 1:m.nc))
    @printf("     CH4 in fault: %.4f → %.4f t; beyond the seal: %.4f → %.4f t\n",
        ch4(M0, fault), ch4(M1, fault), ch4(M0, beyond), ch4(M1, beyond))
    return (; scale, pe = m.pe, buoy, h_fault, wall, dp = maximum(abs.(states[end][:Pressure] .- r0[:Pressure])) / bar,
        dsw = maximum(abs.(sw1 .- sw0)), ch4_fault = (ch4(M0, fault), ch4(M1, fault)), ch4_beyond = (ch4(M0, beyond), ch4(M1, beyond)))
end

group = isempty(ARGS) ? "equilibrium" : ARGS[1]
if group == "equilibrium"
    out = [equilibrium_test(s) for s in (1.0, FS.co2_ratio)]
    mkpath(joinpath(@__DIR__, "..", "output", "04d")); serialize(joinpath(@__DIR__, "..", "output", "04d", "equilibrium.jls"), out)
end

"""
Numerical equilibration: relax the hand-built capillary–gravity state for
`relax_years` without wells, then run a second no-well period from the relaxed
state, which must be at rest. Reports where the first-period changes occurred.
"""
function relaxed_state(scale; relax_years = 10, check_years = 10)
    sec = build_section()
    m = fault_gas_model(sec, GS, FS; pe_scale = scale)
    s0 = fault_gas_initial_state(sec, m, GS, FS)
    off = setup_reservoir_forces(m.model, control = Dict(:Well => DisabledControl()))
    run(s, n) = simulate_reservoir(s, m.model, fill(1.0yr, n); parameters = m.parameters, forces = off, info_level = -1, SOLVER_KW...)
    res1 = run(s0, relax_years)
    st1 = res1.result.states[end]                          # full state (reservoir + well)
    z = m.domain[:cell_centroids][3, :]
    dp = abs.(res1.states[end][:Pressure] .- s0[:Reservoir][:Pressure])
    dsw = abs.(res1.states[end][:Saturations][1, :] .- s0[:Reservoir][:ImmiscibleSaturation])
    names = [[F.name for F in FACIES]..., "Fault zone"]
    cp, cs = argmax(dp), argmax(dsw)
    @printf("== scale %.1f relaxation: max |Δp| %.3f bar at %s (depth %.1f m, x %.0f m); max |ΔSw| %.4f at %s (depth %.1f m, x %.0f m)\n",
        scale, dp[cp] / bar, names[sec.facies[cp]], z[cp], sec.xc[cp], dsw[cs], names[sec.facies[cs]], z[cs], sec.xc[cs])
    @printf("   cells with |ΔSw| > 1e-3: %d (all in trap? %s)\n", count(>(1e-3), dsw),
        all(c -> c in trap_cells(sec, m.domain, FS), findall(>(1e-3), dsw)))
    res2 = simulate_reservoir(st1, m.model, fill(1.0yr, check_years); parameters = m.parameters, forces = off, info_level = -1, SOLVER_KW...)
    dp2 = maximum(abs.(res2.states[end][:Pressure] .- res1.states[end][:Pressure])) / bar
    dsw2 = maximum(abs.(res2.states[end][:Saturations][1, :] .- res1.states[end][:Saturations][1, :]))
    @printf("   second %d-yr no-well period from the relaxed state: max |Δp| %.2e bar, max |ΔSw| %.2e\n", check_years, dp2, dsw2)
    return (; sec, m, state = st1, dp2, dsw2)
end

if group == "relax"
    for s in (1.0, FS.co2_ratio); relaxed_state(s); end
end
