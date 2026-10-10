# Step 4B, stage 2: small 2D tilted gas reservoir with a gas–water contact.
#   a) "stability": equilibrium initial state, 10 years without wells. Nothing
#      should move, with gravity on (hydrostatic initial pressure) and with
#      gravity off (uniform initial pressure).
#   b) "gravity": production then CO2 injection at the crest, then 30 years
#      shut-in, run twice with identical inputs except gravity on/off.
# Run from the project folder:
#   julia --project=. scripts/04b2_gas_2d.jl stability
#   julia --project=. scripts/04b2_gas_2d.jl gravity

include(joinpath(@__DIR__, "..", "src", "gas.jl"))
using CairoMakie, Printf, Serialization

const mD = si_unit(:darcy) / 1000
const GS = gas_spec()
const T_RES = convert_to_si(GS.T_C, :Celsius)
const EOS = gas_eos(GS)

# Geometry: 600 m × 100 m × 60 m slab, top at 1500 m on the left, rising 3° to
# the right (crest on the right). Gas–water contact at 1520 m.
const NX, NZ = 60, 30
const GWC = 1520.0
const P_GWC_BAR = 152.0   # brine pressure at the contact (≈ hydrostatic from sea level)

function build_model(; gravity = true)
    mesh = reservoir_mesh((NX, 1, NZ), (600.0, 100.0, 60.0), origin = [0.0, 0.0, 1500.0])
    slope = tand(3.0)
    for (i, pt) in enumerate(mesh.node_points)
        mesh.node_points[i] = pt - [0.0, 0.0, slope * pt[1]]
    end
    nc = number_of_cells(mesh)
    perm = zeros(3, nc); perm[1, :] .= perm[2, :] .= 100mD; perm[3, :] .= 10mD
    domain = reservoir_domain(mesh, permeability = perm, porosity = 0.2, temperature = T_RES)
    crest = [(NX, 1, k) for k in 1:4]                 # top 8 m at the crest
    well = setup_well(domain, crest, name = :Well, simple_well = true)
    sys = gas_system(GS)
    model = setup_reservoir_model(domain, sys, wells = [well])
    set_brine_properties!(model, GS)
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0, 2.0], [GS.swc, 0.0, 0.0])
    replace_variables!(model, RelativePermeabilities = kr)
    push!(reservoir_model(model).output_variables, :TotalMasses, :Saturations, :PhaseMassDensities)
    parameters = setup_parameters(model)
    if !gravity
        parameters[:Reservoir][:TwoPointGravityDifference] .= 0.0
    end
    return (; mesh, domain, model, parameters, nc)
end

"""
Initial state. Gravity on: hydrostatic equilibration with the gas–water contact
(no hydrocarbon-liquid zone: both contacts at the GWC). Gravity off: same
saturations and composition, uniform pressure (the equilibrium without gravity).
"""
function initial_state(m; gravity = true)
    eql = EquilibriumRegion(m.model, P_GWC_BAR * bar, GWC, woc = GWC, goc = GWC,
        composition = [0.0, 1.0], temperature = T_RES)
    s0 = setup_reservoir_state(m.model, eql)
    # JutulDarcy's equilibration takes connate water from tabulated
    # (ReservoirRelativePermeabilities) curves only; with Brooks–Corey it leaves
    # Sw = 0 in the gas zone. Set connate water there. It is immobile (kr_w = 0
    # at Sw = swc) and there is no capillary pressure, so pressure equilibrium is
    # unchanged.
    sw = s0[:Reservoir][:ImmiscibleSaturation]
    sw[sw .< 1.0] .= max.(sw[sw .< 1.0], GS.swc)
    if !gravity
        p = s0[:Reservoir][:Pressure]
        pv = m.domain[:volumes] .* m.domain[:porosity]
        p .= sum(p .* pv) / sum(pv)
    end
    return s0
end

function run_schedule(m, s0, phases; max_dt = 30day)
    dt = vcat([fill(1.0yr, n) for (n, _) in phases]...)
    forces = vcat([fill(f, n) for (n, f) in phases]...)
    return simulate_reservoir(s0, m.model, dt, parameters = m.parameters, forces = forces,
        info_level = -1, max_timestep = max_dt)
end

function stability()
    for gravity in (true, false)
        m = build_model(gravity = gravity)
        s0 = initial_state(m, gravity = gravity)
        off = setup_reservoir_forces(m.model, control = Dict(:Well => DisabledControl()))
        ws, states, t = run_schedule(m, s0, [(10, off)])
        r0 = s0[:Reservoir]
        sw0 = r0[:ImmiscibleSaturation]
        dp = maximum(abs.(states[end][:Pressure] .- r0[:Pressure])) / bar
        dsw = maximum(abs.(states[end][:Saturations][1, :] .- sw0))
        dz = maximum(abs.(states[end][:OverallMoleFractions] .- r0[:OverallMoleFractions]))
        gas = count(<(0.99), sw0)
        @printf("Gravity %-3s: %d gas cells (Sw < 0.99) of %d; initial p %.2f–%.2f bar. After 10 yr without wells: max |Δp| = %.2e bar, max |ΔSw| = %.2e, max |Δz| = %.2e\n",
            gravity ? "on" : "off", gas, m.nc, minimum(r0[:Pressure]) / bar, maximum(r0[:Pressure]) / bar, dp, dsw, dz)
        if gravity
            z = m.domain[:cell_centroids][3, :]
            sw_above = maximum(sw0[z .< GWC - 1.0]); sw_below = minimum(sw0[z .> GWC + 1.0])
            @printf("            Sw above the contact: max %.3f (connate %.2f); below: min %.3f\n", sw_above, GS.swc, sw_below)
        end
    end
end

"""
Production at the crest well, shut-in, CO2 injection through the same well,
then a long shut-in. Identical inputs except gravity.
"""
function gravity_case(gravity; prod_t_per_yr = 1500.0, inj_t_per_yr = 6000.0,
        years = (deplete = 10, shut1 = 1, inject = 10, shut2 = 30))
    m = build_model(gravity = gravity)
    s0 = initial_state(m, gravity = gravity)
    off = DisabledControl()
    f_dep = setup_reservoir_forces(m.model, control = Dict(:Well => gas_producer_control(prod_t_per_yr, EOS)),
        limits = Dict(:Well => (bhp = 30bar,)))
    f_shut = setup_reservoir_forces(m.model, control = Dict(:Well => off))
    f_inj = setup_reservoir_forces(m.model, control = Dict(:Well => co2_injector_control(inj_t_per_yr, T_RES)))
    _, s_init, _ = simulate_reservoir(s0, m.model, [1.0], parameters = m.parameters, forces = f_shut, info_level = -1)
    phases = [(years.deplete, f_dep), (years.shut1, f_shut), (years.inject, f_inj), (years.shut2, f_shut)]
    dt = vcat([fill(1.0yr, n) for (n, _) in phases]...)
    forces = vcat([fill(f, n) for (n, f) in phases]...)
    t0 = time()
    res = simulate_reservoir(s0, m.model, dt, parameters = m.parameters, forces = forces,
        info_level = -1, max_timestep = 30day, output_substates = true)
    wall = (time() - t0) / 60
    ws, states, t = res
    _, dts, report_ix = Jutul.expand_to_ministeps(res.result)
    reps = res.result.reports
    mini = [mm for r in reps for mm in r[:ministeps]]
    return (; m, s0, M0 = s_init[1][:TotalMasses], ws, states, t, dts, report_ix, years, wall,
        n_ok = count(mm -> mm[:success], mini), n_cut = count(mm -> !mm[:success], mini))
end

function report_gravity(r, label)
    names = reservoir_model(r.m.model).system.components   # CO2, C1, Water
    z = r.m.domain[:cell_centroids][3, :]
    pv = r.m.domain[:volumes] .* r.m.domain[:porosity]
    @printf("\n== %s: %.1f min, %d ministeps accepted, %d cut\n", label, r.wall, r.n_ok, r.n_cut)
    # Component mass balance over ministeps.
    for (j, n) in enumerate(names)
        net = zeros(length(r.dts))
        for (wk, w) in pairs(r.ws.wells)
            key = Symbol("$(n)_mass_rate")
            haskey(w, key) && (net .+= w[key] .* r.dts)
        end
        M = [sum(s[:TotalMasses][j, :]) for s in r.states]
        resid = M .- sum(r.M0[j, :]) .- cumsum(net)
        @printf("  %-5s initial %9.1f t, net well %+10.1f t, final %9.1f t, max |residual| %.2e t\n",
            n, sum(r.M0[j, :]) / 1000, sum(net) / 1000, M[end] / 1000, maximum(abs.(resid)) / 1000)
    end
    # Mass-weighted mean depth of each gas component, and of water above the initial contact.
    ends = [findlast(==(y), r.report_ix) for y in 1:maximum(r.report_ix)]
    zc(s, j) = sum(s[:TotalMasses][j, :] .* z) / sum(s[:TotalMasses][j, :])
    rows = []
    for y in (10, 21, 31, 41, 51)
        i = ends[y]; s = r.states[i]
        pavg = sum(s[:Pressure] .* pv) / sum(pv) / bar
        push!(rows, (y, pavg, zc(s, 1), zc(s, 2)))
        @printf("  year %2d: p_avg %6.2f bar; mean depth CO2 %7.2f m, CH4 %7.2f m (CO2 − CH4 = %+6.2f m)\n",
            y, pavg, zc(s, 1), zc(s, 2), zc(s, 1) - zc(s, 2))
    end
    return (; ends, z, rows)
end

function gravity_compare()
    out = Dict{Bool, Any}()
    for g in (true, false)
        r = gravity_case(g)
        out[g] = (r, report_gravity(r, g ? "gravity on" : "gravity off"))
    end
    # Figure: CO2 mole fraction at end of injection and at the end, gravity on/off.
    fig = Figure(size = (1200, 620))
    for (col, g) in enumerate((true, false)), (row, y) in enumerate((21, 51))
        r, rep = out[g]
        ax = Axis3(fig[row, col], title = "Gravity $(g ? "on" : "off"), year $y: CO₂ mole fraction in gas",
            zreversed = true, azimuth = 1.5π, elevation = 0.0, aspect = (1.0, 0.05, 0.25),
            ylabelvisible = false, yticklabelsvisible = false, yticksvisible = false,
            xgridvisible = false, ygridvisible = false, zgridvisible = false, xlabel = "x (m)", zlabel = "Depth (m)")
        s = r.states[rep.ends[y]]
        gas = findall(<(0.99), s[:Saturations][1, :])
        plot_cell_data!(ax, r.m.mesh, s[:OverallMoleFractions][1, :], cells = gas,
            colormap = :viridis, colorrange = (0, 1), shading = NoShading)
    end
    Colorbar(fig[1:2, 3], colormap = :viridis, colorrange = (0, 1), label = "z_CO₂ (gas cells only)")
    save(joinpath(@__DIR__, "..", "figures", "04b2_gravity.png"), fig)

    fig = Figure(size = (700, 420))
    ax = Axis(fig[1, 1], title = "Mean depth of CO₂ minus mean depth of CH₄", xlabel = "Years", ylabel = "m (positive: CO₂ deeper)")
    for (g, c) in ((true, :firebrick), (false, :steelblue))
        r, rep = out[g]
        zc(s, j) = sum(s[:TotalMasses][j, :] .* rep.z) / sum(s[:TotalMasses][j, :])
        ys = 11:maximum(r.report_ix)
        lines!(ax, collect(ys), [zc(r.states[rep.ends[y]], 1) - zc(r.states[rep.ends[y]], 2) for y in ys],
            color = c, linewidth = 2, label = "gravity $(g ? "on" : "off")")
    end
    vlines!(ax, [21], color = :gray, linestyle = :dot)
    axislegend(ax, position = :lt)
    save(joinpath(@__DIR__, "..", "figures", "04b2_segregation.png"), fig)
    println("\nSaved figures/04b2_gravity.png, figures/04b2_segregation.png")
end

group = isempty(ARGS) ? "stability" : ARGS[1]
group == "stability" && stability()
group == "gravity" && gravity_compare()
