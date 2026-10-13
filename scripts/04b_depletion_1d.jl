# Step 4B, stage 1: one-dimensional depleted gas reservoir.
#   10 yr gas production at a fixed mass rate (with a bottom-hole pressure
#   floor) → 1 yr shut-in → 10 yr CO2 injection → 10 yr shut-in.
# Checks: component mass conservation (exact, per ministep), p/Z material
# balance during depletion, phase labelling, Peng–Robinson CO2 density, and
# grid sensitivity of the CO2–CH4 mixing zone.
# Run from the project folder:  julia --project=. scripts/04b_depletion_1d.jl

include(joinpath(@__DIR__, "..", "src", "gas.jl"))
using CairoMakie, Printf, Serialization
import JutulDarcy.CO2Properties: co2_brine_property_tables

const mD = si_unit(:darcy) / 1000
const GS = gas_spec()
const T_RES = convert_to_si(GS.T_C, :Celsius)
const EOS = gas_eos(GS)

function run_1d(nx; L = 1000.0, W = 100.0, H = 20.0, k = 100.0, phi = 0.2,
        prod_t_per_yr = 2000.0, bhp_min_bar = 30.0, inj_t_per_yr = 8000.0,
        years = (deplete = 10, shut1 = 1, inject = 10, shut2 = 10), max_dt = 30day)
    mesh = reservoir_mesh((nx, 1, 1), (L, W, H))
    domain = reservoir_domain(mesh, permeability = k * mD, porosity = phi, temperature = T_RES)
    inj = setup_well(domain, [(1, 1, 1)], name = :Injector, simple_well = true)
    prod = setup_well(domain, [(nx, 1, 1)], name = :Producer, simple_well = true)
    sys = gas_system(GS)
    model = setup_reservoir_model(domain, sys, wells = [inj, prod])
    set_brine_properties!(model, GS)
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0, 2.0], [GS.swc, 0.0, 0.0])
    replace_variables!(model, RelativePermeabilities = kr)
    push!(reservoir_model(model).output_variables, :TotalMasses, :Saturations, :PhaseMassDensities)
    state0 = setup_reservoir_state(model, Pressure = GS.p_init_bar * bar,
        OverallMoleFractions = [0.0, 1.0], ImmiscibleSaturation = GS.swc)
    parameters = setup_parameters(model)

    off = DisabledControl()
    f_dep = setup_reservoir_forces(model, control = Dict(:Injector => off, :Producer => gas_producer_control(prod_t_per_yr, EOS)),
        limits = Dict(:Producer => (bhp = bhp_min_bar * bar,)))
    f_shut = setup_reservoir_forces(model, control = Dict(:Injector => off, :Producer => off))
    f_inj = setup_reservoir_forces(model, control = Dict(:Producer => off,
        :Injector => co2_injector_control(inj_t_per_yr, T_RES)))

    # Initial component masses: one 1-second step with all wells shut (no flux).
    _, s_init, _ = simulate_reservoir(state0, model, [1.0], parameters = parameters, forces = f_shut, info_level = -1)
    M0 = s_init[1][:TotalMasses]

    phases = [(years.deplete, f_dep), (years.shut1, f_shut), (years.inject, f_inj), (years.shut2, f_shut)]
    dt = vcat([fill(1.0yr, n) for (n, _) in phases]...)
    forces = vcat([fill(f, n) for (n, f) in phases]...)
    res = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces,
        info_level = -1, output_substates = true, max_timestep = max_dt)
    ws, states, t = res
    _, dts, report_ix = Jutul.expand_to_ministeps(res.result)
    return (; model, domain, states, t, dts, report_ix, ws, nx, phi, M0, years)
end

"Component mass balance: cumulative well injection/production integrated over ministeps."
function mass_balance(r)
    names = reservoir_model(r.model).system.components   # ["CO2", "C1", "Water"]
    out = Dict{String, Any}()
    for (j, n) in enumerate(names)
        qin = zeros(length(r.dts)); qout = zeros(length(r.dts))
        for (wk, w) in pairs(r.ws.wells)
            key = Symbol("$(n)_mass_rate")
            haskey(w, key) || continue
            q = w[key]
            @assert length(q) == length(r.dts)
            qin .+= max.(q, 0.0) .* r.dts
            qout .+= max.(-q, 0.0) .* r.dts
        end
        M = [sum(s[:TotalMasses][j, :]) for s in r.states]
        M0 = sum(r.M0[j, :])
        inj = cumsum(qin); prod = cumsum(qout)
        out[n] = (; M0, M, inj, prod, resid = M .- M0 .- (inj .- prod))
    end
    return names, out
end

"p/Z of the pore-volume-averaged pressure vs the closed-tank material balance line."
function pz_check(r, mb)
    pv = r.domain[:volumes] .* r.phi
    mw = EOS.mixture.properties[2].mw
    G = mb["C1"].M0 / mw
    pZ_i = GS.p_init_bar * bar / gas_Z(EOS, GS.p_init_bar * bar, T_RES, [0.0, 1.0])
    rows = []
    for y in 1:r.years.deplete
        i = findlast(j -> r.report_ix[j] == y, eachindex(r.dts))
        pbar = sum(r.states[i][:Pressure] .* pv) / sum(pv)
        pz_sim = pbar / gas_Z(EOS, pbar, T_RES, [0.0, 1.0])
        pz_mb = pZ_i * (1 - (mb["C1"].prod[i] / mw) / G)
        push!(rows, (y, pbar / bar, mb["C1"].prod[i] / 1000, pz_sim, pz_mb))
    end
    return rows
end

function main()
    println("Peng–Robinson CO2 density vs JutulDarcy's CO2 table (NIST-based), 95 °C:")
    tab = co2_brine_property_tables()
    for pb in (40.0, 80.0, 120.0, 150.0)
        rpr = gas_density(EOS, pb * bar, T_RES, [1.0, 0.0])
        rref = tab[:density](pb * bar, T_RES)[2]
        @printf("  %5.0f bar: PR %6.1f, table %6.1f kg/m³ (%+.1f%%); CH4 (PR) %5.1f kg/m³\n",
            pb, rpr, rref, 100 * (rpr - rref) / rref, gas_density(EOS, pb * bar, T_RES, [0.0, 1.0]))
    end

    results = Dict{Int, Any}()
    for nx in (100, 400)
        t0 = time()
        r = run_1d(nx)
        names, mb = mass_balance(r)
        results[nx] = (r, mb)
        @printf("\n== nx = %d: %.1f min, %d ministeps\n", nx, (time() - t0) / 60, length(r.dts))
        for n in names
            m = mb[n]
            @printf("  %-5s initial %9.1f t, injected %8.1f t, produced %8.1f t, final %9.1f t, max |residual| %.2e t (%.1e of initial+injected)\n",
                n, m.M0 / 1000, m.inj[end] / 1000, m.prod[end] / 1000, m.M[end] / 1000,
                maximum(abs.(m.resid)) / 1000, maximum(abs.(m.resid)) / (m.M0 + m.inj[end]))
        end
        S = hcat([s[:Saturations] for s in r.states]...)
        @printf("  max hydrocarbon-liquid saturation over the run: %.2e (0 = single gas phase)\n", maximum(S[2, :]))
    end

    r, mb = results[100]
    println("\np/Z material balance during depletion (nx = 100, pore-volume-averaged pressure):")
    println("  year  p_avg (bar)  CH4 produced (t)   p/Z simulated   p/Z material balance   diff %")
    for (y, pb, gp, pzs, pzm) in pz_check(r, mb)
        @printf("  %4d   %8.2f      %8.1f         %9.4e       %9.4e           %+.3f\n", y, pb, gp, pzs, pzm, 100 * (pzs - pzm) / pzm)
    end

    # Mixing zone at the end of injection: front position (z_CO2 = 0.5) and the
    # 10–90% width, by linear interpolation along x.
    println("\nCO2–CH4 mixing zone at the end of injection (year 21):")
    for nx in (100, 400)
        r, _ = results[nx]
        i = findlast(k -> r.report_ix[k] == r.years.deplete + r.years.shut1 + r.years.inject, eachindex(r.dts))
        x = r.domain[:cell_centroids][1, :]
        z = r.states[i][:OverallMoleFractions][1, :]
        xat(v) = (j = findfirst(<(v), z); x[j-1] + (x[j] - x[j-1]) * (z[j-1] - v) / (z[j-1] - z[j]))
        @printf("  nx = %3d (dx = %4.1f m): front at %.1f m, 90%%→10%% width %.1f m (%.1f cells)\n",
            nx, x[2] - x[1], xat(0.5), xat(0.1) - xat(0.9), (xat(0.1) - xat(0.9)) / (x[2] - x[1]))
        pv = r.domain[:volumes] .* r.phi
        @printf("           average pressure at year 21: %.2f bar\n", sum(r.states[i][:Pressure] .* pv) / sum(pv) / bar)
    end

    # Light results for the figure (drawn by scripts/make_figures.jl, so it can be
    # restyled without rerunning): pressure history, p/Z check, composition profiles.
    fig_data = Dict{Int, Any}()
    for nx in (100, 400)
        r, mb = results[nx]
        pv = r.domain[:volumes] .* r.phi
        pavg = [sum(s[:Pressure] .* pv) / sum(pv) / bar for s in r.states]
        x = r.domain[:cell_centroids][1, :]
        profiles = Dict(y => r.states[findlast(k -> r.report_ix[k] == y, eachindex(r.dts))][:OverallMoleFractions][1, :]
            for y in (r.years.deplete + r.years.shut1 + r.years.inject, sum(r.years)))
        fig_data[nx] = (; t_years = r.t ./ yr, pavg, x, profiles, years = r.years,
            pz = nx == 100 ? pz_check(r, mb) : nothing)
    end
    mkpath(joinpath(@__DIR__, "..", "output", "04b"))
    serialize(joinpath(@__DIR__, "..", "output", "04b", "depletion_1d_figdata.jls"), fig_data)
    println("\nSaved output/04b/depletion_1d_figdata.jls (figure: julia --project=. scripts/make_figures.jl 04b_)")
    return results
end

results = (@isdefined SKIP_MAIN) ? nothing : main()
