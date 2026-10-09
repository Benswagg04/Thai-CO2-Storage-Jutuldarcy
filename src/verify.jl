# Verification run: like run_case, but integrates the CO2 flow across the model
# boundaries over every internal time step ("ministep") and computes a CO2
# conservation residual independently of it:
#
#   injected(t) = CO2 in reservoir(t) + CO2 in well(t) + CO2 out through boundaries(t) + residual(t)
#
# Boundary flux: JutulDarcy's own `compute_bc_mass_fluxes`, evaluated on the
# converged state at the end of each ministep and multiplied by the ministep
# length. This matches the backward-Euler time discretisation the solver uses.
# It is accumulated in a `post_ministep_hook`, so full intermediate states never
# have to be kept in memory (they would need several GB on the finest grid).
#
# From the JutulDarcy 0.3.10 source: the boundary flux is driven by the brine
# (reference) pressure difference only, without gravity across the half cell or
# capillary pressure ("# TODO: Capillary pressure" in multicomponent/bc.jl).

include(joinpath(@__DIR__, "run.jl"))

"Cells whose CO2 has passed the seal: overburden sand, plus fault cells above the footwall caprock top."
function beyond_seal_cells(sec)
    ztop = sec.spec.layer_bottoms[1]
    return [c for c in eachindex(sec.facies)
            if sec.facies[c] == 1 || (sec.facies[c] == FAULT_FACIES && sec.zr[c] < ztop)]
end

"""
Accumulator for boundary CO2 transport (kg), updated after every converged
ministep. Keys: top/side × free/dissolved, plus the part leaving from
beyond-seal cells and from storage-formation cells.
"""
function boundary_flux_hook(sys, bc, is_top_bc, seal, store)
    gmap = Jutul.TrivialGlobalMap()
    acc = Dict(k => 0.0 for k in (:top_free, :top_diss, :side_free, :side_diss, :from_beyond_seal, :from_storage))
    log = Tuple{Float64, Dict{Symbol, Float64}}[]   # (time, cumulative) after each ministep
    t_now = Ref(0.0)
    function hook(done, report, sim, dt, forces, max_iter, cfg)
        if done
            s = Jutul.get_output_state(sim.storage, sim.model)[:Reservoir]
            st = (; (Symbol(k) => v for (k, v) in s)...)
            for (j, b) in enumerate(bc)
                q = JutulDarcy.compute_bc_mass_fluxes(sys, b, gmap, st)
                size(q, 1) == 1 && continue      # inflow: brine only (fractional flow [1, 0])
                q_diss, q_free = q[1, 2], q[2, 2]  # liquid row: dissolved CO2; vapour row: free CO2
                if is_top_bc[j]
                    acc[:top_free] += q_free * dt; acc[:top_diss] += q_diss * dt
                else
                    acc[:side_free] += q_free * dt; acc[:side_diss] += q_diss * dt
                end
                b.cell in seal && (acc[:from_beyond_seal] += (q_free + q_diss) * dt)
                b.cell in store && (acc[:from_storage] += (q_free + q_diss) * dt)
            end
            t_now[] += dt
            push!(log, (t_now[], copy(acc)))
        end
        return (done, report)
    end
    return hook, log
end

"Cumulative value (tonnes) of `key` at time t from the hook log (last entry at or before t)."
function cumulative_at(log, key, t)
    i = findlast(e -> e[1] <= t * (1 + 1e-9), log)
    return isnothing(i) ? 0.0 : log[i][2][key] / 1000
end

function verify_case(name, sec, cs; sim_kw = (;), output_dir = joinpath(@__DIR__, "..", "output", "03b"))
    case = setup_case(sec, cs)
    (; model, state0, parameters, forces, dt, regions, pc_tables, bc, top) = case
    p0 = state0[:Reservoir][:Pressure]
    println("== Verification case $name ($(number_of_cells(sec.mesh)) cells) sim_kw = $sim_kw")

    # setup_case creates the side BCs first and the top BCs (if any) last, so
    # the position in `bc` identifies the boundary (the two top corner cells
    # carry one BC of each kind).
    n_top = cs.open_top ? length(top) : 0
    is_top_bc = [j > length(bc) - n_top for j in eachindex(bc)]
    seal = Set(beyond_seal_cells(sec))
    store = Set(vcat(regions.fw_storage, regions.hw_storage))
    hook, flog = boundary_flux_hook(reservoir_model(model).system, bc, is_top_bc, seal, store)

    t_start = time()
    res = simulate_reservoir(state0, model, dt; parameters = parameters, forces = forces,
        info_level = -1, post_ministep_hook = hook, merge((max_timestep = 90day,), sim_kw)...)
    wall = (time() - t_start) / 60
    ws, states, t = res

    # Solver statistics from the reports.
    reps = res.result.reports
    mini = [m for r in reps for m in r[:ministeps]]
    n_fail = count(m -> !m[:success], mini)
    n_newton = sum(m -> haskey(m, :steps) ? length(m[:steps]) : m[:stats].newtons, mini)
    stats = (; wall_min = wall, report_steps = length(reps), ministeps_ok = length(mini) - n_fail,
        ministeps_failed = n_fail, newton = n_newton)
    @printf("Wall time %.1f min, %d ministeps accepted, %d cut\n", wall, stats.ministeps_ok, n_fail)

    # Inventories at report steps (whole years).
    yrs = round.(Int, t ./ yr)
    inventory = co2_inventory(model, ws, states, t)
    injected = [i[:injected] / 1000 for i in inventory]
    in_res = [sum(s[:TotalMasses][2, :]) / 1000 for s in states]
    full = res.result.states
    well = [haskey(s, :Injector) && haskey(s[:Injector], :TotalMasses) ?
            sum(s[:Injector][:TotalMasses][2, :]) / 1000 : 0.0 for s in full]
    cum(k) = [cumulative_at(flog, k, ti) for ti in t]
    B = cum(:top_free) .+ cum(:top_diss) .+ cum(:side_free) .+ cum(:side_diss)
    residual = injected .- in_res .- well .- B

    acc = co2_by_region(states, regions)
    beyond = beyond_seal_cells(sec)
    M_beyond = [sum(s[:TotalMasses][2, beyond]) / 1000 for s in states]
    M_store = acc.total[:fw_storage] .+ acc.total[:hw_storage]
    # Cumulative net CO2 that has crossed into the beyond-seal region (what it
    # holds now + what has left it through the model boundaries), and out of
    # the storage formation (injected − held − exited directly from storage
    # cells; this one also absorbs the conservation residual).
    crossed_seal = M_beyond .+ cum(:from_beyond_seal)
    out_of_storage = injected .- M_store .- cum(:from_storage)
    diag = fault_diagnostics(sec, states, p0, pc_tables)

    summary = (; years = yrs, injected, in_res, well,
        top_free = cum(:top_free), top_diss = cum(:top_diss),
        side_free = cum(:side_free), side_diss = cum(:side_diss),
        residual, max_abs_residual = maximum(abs.(residual)),
        crossed_seal, out_of_storage,
        fault_free = acc.free[:fault], fault_diss = acc.diss[:fault],
        beyond_free = acc.free[:overburden], beyond_diss = acc.diss[:overburden],
        drive = diag.drive, dp_max = diag.dp_max,
        stats, sim_kw = string(sim_kw), cs, spec = sec.spec)

    for y in unique((cs.inject_years, yrs[end]))
        i = findfirst(==(y), yrs)
        @printf("Year %d: injected %.2f t = reservoir %.2f + well %.3f + out via top %.3f (free %.3f, dissolved %.3f) + sides %.3f + residual %.3f\n",
            y, injected[i], in_res[i], well[i], summary.top_free[i] + summary.top_diss[i],
            summary.top_free[i], summary.top_diss[i], summary.side_free[i] + summary.side_diss[i], residual[i])
        @printf("         crossed the seal (cumulative) %.3f t; out of storage formation %.3f t\n",
            crossed_seal[i], out_of_storage[i])
    end
    @printf("Max |residual| at report steps: %.3f t (%.2e of injected)\n",
        summary.max_abs_residual, summary.max_abs_residual / injected[end])
    @printf("Solver: %d report steps, %d ministeps (%d cut), %d Newton iterations\n",
        stats.report_steps, stats.ministeps_ok, n_fail, stats.newton)

    mkpath(output_dir)
    serialize(joinpath(output_dir, "$name.jls"), summary)
    println("Saved output/03b/$name.jls\n")
    return summary
end
