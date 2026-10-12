# Stage B2c: cost of the mass-balance tolerance in branch (a) (methane curves,
# 1,500 t/yr for `inj_years`, then shut-in for `shut_years`), from the depleted state.
# Prints one progress line per simulated year and stops after `wall_limit` minutes.
# Needs output/04d/depletion.jls (scripts/04d_fault_gas.jl runs).
# Run from the project folder:
#   julia --project=. scripts/04d_tolerance.jl tol_mb [inj_years] [shut_years] [wall_limit_min]

include(joinpath(@__DIR__, "..", "src", "gasfault.jl"))
using Serialization

GS = gas_spec(); FS = fault_gas_spec()
tol_mb = parse(Float64, ARGS[1])
inj_years = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
shut_years = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
wall_limit = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 15.0

dep = deserialize(joinpath(@__DIR__, "..", "output", "04d", "depletion.jls"))
sec = build_section()
m1 = fault_gas_model(sec, GS, FS; pe_scale = 1.0)
inj = setup_reservoir_forces(m1.model, control = Dict(:Well => co2_injector_control(1500.0, T_GAS)))
off = setup_reservoir_forces(m1.model, control = Dict(:Well => DisabledControl()))
forces = vcat(fill(inj, inj_years), fill(off, shut_years))

t_start = time(); n_ok = Ref(0); n_cut = Ref(0); t_sim = Ref(0.0); dts = Float64[]; next_year = Ref(1)
function hook(done, report, sim, dt, forces, max_iter, cfg)
    if done
        n_ok[] += 1; t_sim[] += dt; push!(dts, dt)
        if t_sim[] >= next_year[] * yr * (1 - 1e-9)
            @printf("  year %2d done: wall %5.1f min, ministeps %d ok / %d cut so far\n", next_year[], (time() - t_start) / 60, n_ok[], n_cut[])
            flush(stdout); next_year[] += 1
        end
    else
        n_cut[] += 1
    end
    if (time() - t_start) / 60 > wall_limit
        @printf("STOP at wall limit: simulated %.3f yr, %d ok / %d cut\n", t_sim[] / yr, n_ok[], n_cut[]); flush(stdout)
        error("wall limit")
    end
    return (done, report)
end
@printf("tol_mb %.0e: %d yr injection + %d yr shut-in, wall limit %.0f min\n", tol_mb, inj_years, shut_years, wall_limit); flush(stdout)
try
    simulate_reservoir(dep.state, m1.model, fill(1.0yr, length(forces)); parameters = m1.parameters, forces = forces,
        info_level = -1, post_ministep_hook = hook, merge(SOLVER_KW, (; tol_mb))...)
catch err
    occursin("wall limit", sprint(showerror, err)) || rethrow()
end
@printf("tol_mb %.0e: %.1f min wall, simulated %.3f yr, ministeps %d ok / %d cut, median dt %.2f day, min dt %.3g day\n",
    tol_mb, (time() - t_start) / 60, t_sim[] / yr, n_ok[], n_cut[], isempty(dts) ? NaN : sort(dts)[cld(length(dts), 2)] / day,
    isempty(dts) ? NaN : minimum(dts) / day)
