# Step 4D: faulted depleted gas field.
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

# ---------------------------------------------------------------------------
# Main runs: one depletion history, then three 30-year branches from the same
# depleted state (a: methane curves + injection; b: CO2-scaled + injection;
# c: CO2-scaled, no injection = control for the property switch).
# ---------------------------------------------------------------------------
const OUT4D = joinpath(@__DIR__, "..", "output", "04d")

"Cells used for diagnostics: regions, and the matched-depth pair across the fault at the trap crest."
function diagnostics_cells(sec, m)
    z = m.domain[:cell_centroids][3, :]
    trap = trap_cells(sec, m.domain, FS)
    x_sand = maximum(sec.xc[c] for c in trap)
    crest = argmin([c in trap && sec.xc[c] == x_sand ? z[c] : Inf for c in eachindex(z)])  # top trap cell at the fault
    x_fault = minimum(sec.xc[c] for c in eachindex(z) if sec.facies[c] == FAULT_FACIES)
    x_fault2 = maximum(sec.xc[c] for c in eachindex(z) if sec.facies[c] == FAULT_FACIES)
    level(c) = cell_ijk(sec.mesh, c)[3]
    fault_cell = only(c for c in eachindex(z) if sec.xc[c] == x_fault && level(c) == level(crest))
    x_hw = minimum(sec.xc[c] for c in eachindex(z) if sec.xc[c] > x_fault2)
    hw_cell = only(c for c in eachindex(z) if sec.xc[c] == x_hw && level(c) == level(crest))
    fw_store = [c for c in eachindex(z) if sec.facies[c] == 3 && sec.xc[c] < sec.spec.fault_x[1]]
    return (; trap, crest, fault_cell, hw_cell, z,
        regions = (fw_storage = fw_store,
                   hw_storage = [c for c in eachindex(z) if sec.facies[c] == 3 && sec.xc[c] > sec.spec.fault_x[2]],
                   fault = findall(sec.facies .== FAULT_FACIES),
                   caprock = findall(sec.facies .== 2),
                   beyond_seal = fault_gas_beyond_seal(sec)))
end

"""
Post-ministep hook: cumulative well mass flow per component (t, injection
positive), from the facility's total mass rate × ministep length. Injection is
pure CO2 (the injected stream; the well's own contents are not, so they cannot
be used for the split). Production is split by the composition of the well's
contents, i.e. the produced mixture. Exact in time, unlike end-of-step rate ×
step length.
"""
function well_flow_hook()
    acc = zeros(3); t_now = Ref(0.0)
    log = Tuple{Float64, Vector{Float64}}[]
    function hook(done, report, sim, dt, forces, max_iter, cfg)
        if done
            q = Jutul.value(only(sim.storage[:Facility].state.TotalSurfaceMassRate))
            if q > 0
                x = [1.0, 0.0, 0.0]
            else
                x = vec(sum(Jutul.value.(sim.storage[:Well].state.TotalMasses), dims = 2)); x ./= sum(x)
            end
            acc .+= q * dt / 1000 .* x
            t_now[] += dt
            push!(log, (t_now[], copy(acc)))
        end
        return (done, report)
    end
    return hook, log
end
cumulative_at(log, t) = (i = findlast(e -> e[1] <= t * (1 + 1e-12), log); i === nothing ? zeros(3) : log[i][2])

"""
Per-report-step summary of one run segment. Conservation per component:
(reservoir + well contents) − (same at the segment start) − cumulative well flow.
"""
function summarise_segment(sec, m, d, res, t_offset, state0, flow_log)
    ws, states, t = res
    tot(M) = vec(sum(M, dims = 2)) ./ 1000
    start = tot(state0[:Reservoir][:TotalMasses]) .+ tot(state0[:Well][:TotalMasses])
    pe_fault = m.pe[FAULT_FACIES]
    pcog = m.model[:Reservoir].secondary_variables[:CapillaryPressure].pc[2]
    pc_sand(sg) = JutulDarcy.table_by_region(pcog, 3)(sg)
    rows = []
    for (i, s) in enumerate(states)
        M = s[:TotalMasses]
        mass(j, cells) = sum(M[j, cells]) / 1000
        p = s[:Pressure]; Sg = s[:Saturations][3, :]
        # Entry drive = brine-pressure difference across the sand–fault interface + Pc at the crest.
        dp_brine = (p[d.crest] - p[d.fault_cell]) / bar
        pc_crest = pc_sand(Sg[d.crest]) / bar
        drive = dp_brine + pc_crest
        pv = m.domain[:volumes][d.regions.fw_storage] .* m.domain[:porosity][d.regions.fw_storage]
        well_node = tot(res.result.states[i][:Well][:TotalMasses])      # mass held in the well itself
        flow = cumulative_at(flow_log, t[i])
        push!(rows, (; year = (t_offset + t[i]) / yr,
            p_trap = sum(p[d.regions.fw_storage] .* pv) / sum(pv) / bar,
            p_crest = p[d.crest] / bar, p_fault = p[d.fault_cell] / bar, p_hw_matched = p[d.hw_cell] / bar,
            drive, dp_brine, pc_crest, pe_fault,
            sg_max_fault = maximum(Sg[d.regions.fault]),
            co2 = map(cells -> mass(1, cells), d.regions), ch4 = map(cells -> mass(2, cells), d.regions),
            co2_total = mass(1, 1:size(M, 2)), ch4_total = mass(2, 1:size(M, 2)), water_total = mass(3, 1:size(M, 2)),
            well_node, well_flow = flow,
            residual = tot(M) .+ well_node .- start .- flow))
    end
    # Well mass flows per report step (rate at the end of each step × step length).
    dts = diff(vcat(0.0, t))
    well = Dict(n => sum(w[Symbol("$(n)_mass_rate")] .* dts for (k, w) in pairs(ws.wells) if haskey(w, Symbol("$(n)_mass_rate"))) / 1000
        for n in ("CO2", "C1", "Water"))
    mini = [mm for r in res.result.reports for mm in r[:ministeps]]
    stats = (; ok = count(mm -> mm[:success], mini), cut = count(mm -> !mm[:success], mini))
    return (; rows, well, stats)
end

function run_segment(m, state0, phases)
    dt = vcat([fill(1.0yr, n) for (n, _) in phases]...)
    forces = vcat([fill(f, n) for (n, f) in phases]...)
    t0 = time()
    hook, flow_log = well_flow_hook()
    res = simulate_reservoir(state0, m.model, dt; parameters = m.parameters, forces = forces, info_level = -1,
        post_ministep_hook = hook, SOLVER_KW...)
    return res, (time() - t0) / 60, flow_log
end

"""
Depletion once, then branches from the same depleted state.
`set = :main`: constant-rate injection (a, b) and the no-injection control (c).
`set = :bhp`: the same injection with a bottom-hole-pressure cap equal to the
original brine pressure at the top perforation, as an operator would apply to a
depleted field (a_bhp, b_bhp).
"""
function main_runs(; set = :main, out = OUT4D, fault_subcells = 1, branch_names = nothing, prod_t_per_yr = 300.0, inj_t_per_yr = 1500.0, years = (deplete = 10, shut1 = 1, inject = 10, shut2 = 20))
    mkpath(out)
    sec = build_section(section_spec(; fault_subcells))
    eos = gas_eos(GS)
    m1 = fault_gas_model(sec, GS, FS; pe_scale = 1.0)
    m2 = fault_gas_model(sec, GS, FS; pe_scale = FS.co2_ratio)
    d = diagnostics_cells(sec, m1)
    off(m) = setup_reservoir_forces(m.model, control = Dict(:Well => DisabledControl()))
    prod(m) = setup_reservoir_forces(m.model, control = Dict(:Well => gas_producer_control(prod_t_per_yr, eos)),
        limits = Dict(:Well => (bhp = 30bar,)))
    inj(m) = setup_reservoir_forces(m.model, control = Dict(:Well => co2_injector_control(inj_t_per_yr, T_GAS)))
    p_orig = fault_gas_initial_state(sec, m1, GS, FS)[:Reservoir][:Pressure]
    p_cap = p_orig[physical_representation(m1.model[:Well].domain).perforations.reservoir[1]]
    inj_cap(m) = setup_reservoir_forces(m.model, control = Dict(:Well => co2_injector_control(inj_t_per_yr, T_GAS)),
        limits = Dict(:Well => (bhp = p_cap,)))

    f_dep = joinpath(out, "depletion.jls")
    if isfile(f_dep)
        dep = deserialize(f_dep); println("loaded depletion")
    else
        relax, = run_segment(m1, fault_gas_initial_state(sec, m1, GS, FS), [(10, off(m1))])
        s_relaxed = relax.result.states[end]
        res, wall, flow_log = run_segment(m1, s_relaxed, [(years.deplete, prod(m1)), (years.shut1, off(m1))])
        dep = (; state = res.result.states[end], summary = summarise_segment(sec, m1, d, res, 0.0, s_relaxed, flow_log), wall,
            M0 = relax.states[end][:TotalMasses])
        serialize(f_dep, dep)
        @printf("depletion: %.1f min, ministeps %d ok / %d cut\n", wall, dep.summary.stats.ok, dep.summary.stats.cut)
    end
    t_dep = (years.deplete + years.shut1) * yr
    branches = set == :main ?
        [("a_methane_inject", m1, [(years.inject, inj(m1)), (years.shut2, off(m1))]),
         ("b_co2scaled_inject", m2, [(years.inject, inj(m2)), (years.shut2, off(m2))]),
         ("c_co2scaled_control", m2, [(years.inject + years.shut2, off(m2))])] :
        [("a_bhp_methane_inject", m1, [(years.inject, inj_cap(m1)), (years.shut2, off(m1))]),
         ("b_bhp_co2scaled_inject", m2, [(years.inject, inj_cap(m2)), (years.shut2, off(m2))])]
    set == :bhp && @printf("injector BHP cap: %.2f bar (original brine pressure at the top perforation)\n", p_cap / bar)
    isnothing(branch_names) || filter!(b -> b[1] in branch_names, branches)
    for (name, m, phases) in branches
        println("start $name at ", Libc.strftime("%H:%M", time())); flush(stdout)
        f = joinpath(out, "$name.jls")
        isfile(f) && (println("skip $name"); continue)
        try
            res, wall, flow_log = run_segment(m, dep.state, phases)
            sm = summarise_segment(sec, m, d, res, t_dep, dep.state, flow_log)
            maps = Dict(Int(round((t_dep + res.time[i]) / yr)) => (Sg = res.states[i][:Saturations][3, :], zco2 = res.states[i][:OverallMoleFractions][1, :])
                for i in eachindex(res.states) if round(Int, (t_dep + res.time[i]) / yr) in (years.deplete + years.shut1 + years.inject, 41))
            serialize(f, (; name, summary = sm, wall, maps))
            r = sm.rows[end]
            @printf("%s: %.1f min, ministeps %d ok / %d cut\n", name, wall, sm.stats.ok, sm.stats.cut)
            @printf("   year %.0f: p_trap %.1f bar; CO2 in fault %.2f t, beyond seal %.2f t; CH4 in fault %.2f t, beyond seal %.2f t; max Sg in fault %.3f\n",
                r.year, r.p_trap, r.co2.fault, r.co2.beyond_seal, r.ch4.fault, r.ch4.beyond_seal, r.sg_max_fault)
        catch err
            println("FAILED $name: ", sprint(showerror, err)[1:min(end, 400)])
            flush(stdout)
        end
        flush(stdout); GC.gc()
    end
end

group == "runs" && main_runs()
group == "runs_bhp" && main_runs(set = :bhp)

# Fault-resolution test: the fault zone (10 m) split into 5 columns of 2 m, so the
# brine pressure next to the sand–fault interface is resolved (see notes/04d_fault_gas.md).
group == "runs_fault5" && main_runs(out = OUT4D * "_fault5", fault_subcells = 5,
    branch_names = ["a_methane_inject", "b_co2scaled_inject"])
group == "runs_fault10" && main_runs(out = OUT4D * "_fault10", fault_subcells = 10,
    branch_names = ["a_methane_inject", "b_co2scaled_inject"])
