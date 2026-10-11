# Step 4B, stage 2: solver-tuning experiment on the gravity-on 2D case.
# Tolerances are kept at their defaults; only time stepping and nonlinear
# update controls change. Each run records cost, where steps are cut, mass
# conservation and the physical results (a time-step sensitivity check).
# Restartable: finished runs in output/04b2_tuning/ are skipped.
# Run from the project folder:  julia --project=. scripts/04b2_tuning.jl

include(joinpath(@__DIR__, "04b2_gas_2d.jl"))

const CASES = [
    ("baseline_dt30", (;), nothing),
    ("A_dt10", (max_timestep = 10day,), nothing),
    ("B_relaxation", (relaxation = true,), nothing),
    ("C_growth2", (timestep_max_increase = 2.0,), nothing),
    ("D_dzmax0.05", (;), 0.05),
]
outdir = joinpath(@__DIR__, "..", "output", "04b2_tuning")
mkpath(outdir)

function summarise(name, r)
    z = r.m.domain[:cell_centroids][3, :]
    pv = r.m.domain[:volumes] .* r.m.domain[:porosity]
    ends = [findlast(==(y), r.report_ix) for y in 1:maximum(r.report_ix)]
    zc(s, j) = sum(s[:TotalMasses][j, :] .* z) / sum(s[:TotalMasses][j, :])
    seg = Dict(y => zc(r.states[ends[y]], 1) - zc(r.states[ends[y]], 2) for y in (21, 31, 51))
    pavg = sum(r.states[end][:Pressure] .* pv) / sum(pv) / bar
    resid = Float64[]
    for (j, n) in enumerate(reservoir_model(r.m.model).system.components)
        net = zeros(length(r.dts))
        for (wk, w) in pairs(r.ws.wells)
            key = Symbol("$(n)_mass_rate"); haskey(w, key) && (net .+= w[key] .* r.dts)
        end
        M = [sum(s[:TotalMasses][j, :]) for s in r.states]
        push!(resid, maximum(abs.(M .- sum(r.M0[j, :]) .- cumsum(net))) / max(sum(r.M0[j, :]), sum(abs.(net))))
    end
    y = r.years
    periods = (deplete = 1:y.deplete, shut1 = y.deplete .+ (1:y.shut1),
        inject = (y.deplete + y.shut1) .+ (1:y.inject), shut2 = (y.deplete + y.shut1 + y.inject) .+ (1:y.shut2))
    cuts = map(p -> sum(r.cuts_per_year[p]), periods)
    return (; name, wall = r.wall, ok = r.n_ok, cut = r.n_cut, newton = r.newton, min_dt_days = r.min_dt / day,
        cuts, seg, pavg, max_rel_resid = maximum(resid))
end

for (name, kw, dz) in CASES
    f = joinpath(outdir, "$name.jls")
    isfile(f) && (println("skip $name"); continue)
    println("== $name: sim_kw = $kw, dz_max = $(something(dz, "default 0.1"))")
    try
        r = gravity_case(true; sim_kw = kw, dz_max = dz)
        s = summarise(name, r)
        serialize(f, s)
        @printf("   %.1f min; ministeps %d ok / %d cut; Newton %d; min dt %.2f d\n", s.wall, s.ok, s.cut, s.newton, s.min_dt_days)
        @printf("   cuts: deplete %d, shut-in 1 %d, inject %d, shut-in 2 %d\n", s.cuts.deplete, s.cuts.shut1, s.cuts.inject, s.cuts.shut2)
        @printf("   CO2−CH4 mean depth: yr21 %+.2f, yr31 %+.2f, yr51 %+.2f m; p_avg yr51 %.2f bar; max rel. mass residual %.1e\n",
            s.seg[21], s.seg[31], s.seg[51], s.pavg, s.max_rel_resid)
    catch err
        println("   FAILED: ", sprint(showerror, err)[1:min(end, 300)])
    end
    GC.gc()
end
