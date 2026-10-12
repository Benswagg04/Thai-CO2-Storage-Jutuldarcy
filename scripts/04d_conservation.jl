# Stage B2c conservation check: rerun the no-injection control (branch c) from the
# depleted state and track the mass held inside the shut-in well as well as in the
# reservoir. With closed boundaries and no well flow, reservoir + well must be constant.
# Optional argument: tol_mb (default 1e-7, the JutulDarcy default; SOLVER_KW uses 1e-10). The MB check is
# relative to the mass of the whole model (~1.27e10 kg, mostly brine).
# Run from the project folder:  julia --project=. scripts/04d_conservation.jl [tol_mb]

include(joinpath(@__DIR__, "..", "src", "gasfault.jl"))
using Serialization

GS = gas_spec(); FS = fault_gas_spec()
dep = deserialize(joinpath(@__DIR__, "..", "output", "04d", "depletion.jls"))
sec = build_section()
m2 = fault_gas_model(sec, GS, FS; pe_scale = FS.co2_ratio)
tol_mb = isempty(ARGS) ? 1e-7 : parse(Float64, ARGS[1])
off = setup_reservoir_forces(m2.model, control = Dict(:Well => DisabledControl()))
res = simulate_reservoir(dep.state, m2.model, fill(1.0yr, 30); parameters = m2.parameters, forces = off,
    info_level = -1, merge(SOLVER_KW, (; tol_mb))...)
mini = [mm for r in res.result.reports for mm in r[:ministeps]]
@printf("tol_mb = %.0e: ministeps %d ok / %d cut\n", tol_mb, count(mm -> mm[:success], mini), count(mm -> !mm[:success], mini))

tot(s) = vec(sum(s[:TotalMasses], dims = 2)) ./ 1000      # t per component (CO2, CH4, water)
r0, w0 = tot(dep.state[:Reservoir]), tot(dep.state[:Well])
@printf("start: reservoir CH4 %.3f t; well CH4 %.4f t, water %.4f t\n", r0[2], w0[2], w0[3])
println(" yr | reservoir ΔCH4 | well ΔCH4 | total ΔCH4 | reservoir Δwater | well Δwater | total Δwater (t)")
for i in (1, 2, 5, 10, 20, 30)
    r, w = tot(res.result.states[i][:Reservoir]), tot(res.result.states[i][:Well])
    @printf("%3d | %9.4f | %9.4f | %9.4f | %9.4f | %9.4f | %9.4f\n", i,
        r[2] - r0[2], w[2] - w0[2], r[2] + w[2] - r0[2] - w0[2], r[3] - r0[3], w[3] - w0[3], r[3] + w[3] - r0[3] - w0[3])
end
