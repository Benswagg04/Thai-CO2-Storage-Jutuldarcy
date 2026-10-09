# Run one case end to end: equilibrium check, simulation, CO2 accounting,
# fault entry-pressure diagnostics, saved results.

include(joinpath(@__DIR__, "case.jl"))
using Serialization, Printf

"""
Pressure diagnostics at the fault for each report step:
- `drive`: max over the storage-sand cells touching the fault of
  p_CO2(sand) − p_brine(fault cell at the same level), in bar. CO2 can enter
  the fault when this exceeds the fault entry pressure. It includes buoyancy
  and injection overpressure together.
- `overpressure_fault`: brine pressure rise in the fault cell where `drive` is
  largest (injection effect alone), in bar.
- `dp_max`: largest pressure rise anywhere in the model, in bar.
"""
function fault_diagnostics(sec, states, p0, pc_tables)
    (; mesh, facies, xc, spec) = sec
    x_sand = maximum(xc[c] for c in eachindex(xc) if xc[c] < spec.fault_x[1])
    x_fault = minimum(xc[c] for c in eachindex(xc) if facies[c] == FAULT_FACIES)
    level(c) = cell_ijk(mesh, c)[3]
    fault_at = Dict(level(c) => c for c in eachindex(xc) if facies[c] == FAULT_FACIES && xc[c] == x_fault)
    pairs_ = [(c, fault_at[level(c)]) for c in eachindex(xc) if xc[c] == x_sand && facies[c] == 3]
    pc_sand(sg) = isnothing(pc_tables) ? 0.0 : pc_tables[3](sg)
    drive = Float64[]; over = Float64[]; dp_max = Float64[]
    for s in states
        p, S = s[:Pressure], s[:Saturations]
        best, ov = -Inf, 0.0
        for (cs_, cf) in pairs_
            S[2, cs_] > 0.01 || continue
            d = (p[cs_] + pc_sand(S[2, cs_]) - p[cf]) / bar
            if d > best
                best, ov = d, (p[cf] - p0[cf]) / bar
            end
        end
        push!(drive, isfinite(best) ? best : NaN)
        push!(over, ov)
        push!(dp_max, maximum(p .- p0) / bar)
    end
    return (; drive, overpressure_fault = over, dp_max)
end

"""
    run_case(name, sec, cs; output_dir, save_states = :all)

`save_states = :all` saves saturation and pressure for every report step;
`:final` saves only the last step (for sweeps, ~20 MB less per run).
"""
function run_case(name, sec, cs; output_dir = joinpath(@__DIR__, "..", "output"), save_states = :all,
        sim_kw = (;))
    case = setup_case(sec, cs)
    (; model, state0, parameters, forces, forces_shut, dt, regions, pc_tables, edge) = case
    p0 = state0[:Reservoir][:Pressure]
    println("== Case $name ==  ($(number_of_cells(sec.mesh)) cells, open top: $(cs.open_top), capillary: $(cs.capillary))")
    if cs.capillary
        names = [[F.name for F in FACIES]..., "Fault zone"]
        println("Entry pressures (bar): ", join(["$(n) $(round(p, sigdigits = 3))" for (n, p) in zip(names, entry_pressures(sec, cs))], ", "))
    end

    # 1. Equilibrium check: with no injection, nothing should move.
    _, states_eq, _ = simulate_reservoir(state0, model, fill(1.0yr, 5),
        parameters = parameters, forces = forces_shut, info_level = -1)
    @printf("Equilibrium check: max pressure change after 5 years without wells = %.4f bar\n",
        maximum(abs.(states_eq[end][:Pressure] .- p0)) / bar)

    # 2. Injection and migration.
    t_start = time()
    ws, states, t = simulate_reservoir(state0, model, dt;
        parameters = parameters, forces = forces, info_level = -1,
        merge((max_timestep = 90day,), sim_kw)...)
    @printf("Simulation wall time: %.1f min\n", (time() - t_start) / 60)

    # 3. CO2 accounting. Regions partition the model; CO2 that has left the
    #    model is injected minus what is inside it.
    inventory = co2_inventory(model, ws, states, t)
    acc = co2_by_region(states, regions)
    exited = [i[:outside_domain] / 1000 for i in inventory]
    edge_co2 = [sum(s[:TotalMasses][2, edge]) / 1000 for s in states]
    diag = fault_diagnostics(sec, states, p0, pc_tables)

    for step in (cs.inject_years, length(dt))
        inj = inventory[step][:injected] / 1000
        println("Year $step: CO2 (t)        total     free  dissolved")
        for k in keys(regions)
            @printf("  %-12s %12.0f %8.0f %10.0f\n", k, acc.total[k][step], acc.free[k][step], acc.diss[k][step])
        end
        # "left model" is defined as injected − in model, so it also absorbs any
        # mass-balance error from the solver tolerance. In fully sealed cases it
        # is that error (a few tonnes at most); check it against CO2 actually
        # reaching the boundaries before reading it as leakage.
        in_model = sum(acc.total[k][step] for k in keys(regions))
        @printf("  %-12s %12.1f   (= injected %.1f − in model %.1f)\n", "left model", exited[step], inj, in_model)
    end
    @printf("Max CO2 ever in side-boundary columns: %.3f t\n", maximum(edge_co2))
    @printf("Trapping at year %d (t): mobile %.0f, residual %.0f, dissolved %.0f\n", length(dt),
        inventory[end][:mobile] / 1000, inventory[end][:residual] / 1000, inventory[end][:dissolved] / 1000)
    @printf("Max pressure build-up during injection: %.2f bar\n", maximum(diag.dp_max[1:cs.inject_years]))
    if cs.capillary
        pe_f = entry_pressures(sec, cs)[FAULT_FACIES]
        @printf("Fault entry drive (bar): yr5 %.2f, yr25 %.2f, yr100 %.2f vs fault Pe %.2f\n",
            diag.drive[5], diag.drive[cs.inject_years], diag.drive[end], pe_f)
    end

    mkpath(output_dir)
    out = joinpath(output_dir, "$name.jls")
    kept = save_states == :all ? states : states[end:end]
    serialize(out, (; t, states = [Dict(:Saturations => s[:Saturations], :Pressure => s[:Pressure]) for s in kept],
        inventory, acc, exited, edge_co2, diag, p0, cs, spec = sec.spec))
    println("Saved $out\n")
    return (; ws, states, t, inventory, acc, exited, diag)
end
