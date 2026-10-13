# Step 3c: fault-resolution check of the Step 3 sweep on representative scenarios.
# Stage 4B-2c showed that one cell across the 10 m fault overestimates leakage in
# the depleted-field case. Here six Step 3 scenarios (all three regimes) are rerun
# with the fault split into 5 columns of 2 m; everything else (grid, fluid, schedule,
# solver settings) is as in Step 3, so the only difference is fault resolution.
# The 1-column results are read from the Step 3 sweep (output/03_sweep).
#
# Restartable: finished runs (output/03c/<name>.jls) are skipped.
# Run from the project folder:  julia --project=. scripts/03c_fault_resolution.jl

include(joinpath(@__DIR__, "..", "src", "run.jl"))

scenarios = [   # (pe_fault_bar, (k_across, k_along) mD, Step 3 regime)
    (0.6, (0.1, 1.0), "persistent leakage"),
    (1.0, (0.1, 1.0), "near threshold"),
    (1.2, (0.1, 1.0), "injection-triggered"),
    (1.4, (0.1, 1.0), "injection-triggered (weak)"),
    (2.0, (0.1, 1.0), "capillary seal"),
    (1.2, (0.01, 0.1), "injection-triggered, tight fault"),
]
nsub = 5
dir1 = joinpath(@__DIR__, "..", "output", "03_sweep")
dir5 = joinpath(@__DIR__, "..", "output", "03c"); mkpath(dir5)
runname(pe, k) = "pe$(pe)_ka$(k[1])_kl$(k[2])"

for (pe, k, _) in scenarios
    name = runname(pe, k)
    isfile(joinpath(dir5, "$name.jls")) && (println("skip $name (done)"); continue)
    println("start $name at ", Libc.strftime("%H:%M", time())); flush(stdout)
    sec = build_section(section_spec(k_across = k[1], k_along = k[2], fault_subcells = nsub))
    cs = case_spec(capillary = true, pe_fault_bar = pe)
    t0 = time()
    try
        run_case(name, sec, cs, output_dir = dir5, save_states = :final)
        @printf("done %s in %.1f min\n", name, (time() - t0) / 60)
    catch err
        println("FAILED $name: ", sprint(showerror, err)[1:min(end, 300)])
    end
    flush(stdout); GC.gc()
end

# Same metrics as results/03_sweep_summary.csv, for 1 and 5 fault columns.
outside(r) = [i[:injected] / 1000 for i in r.inventory] .- (r.acc.total[:fw_storage] .+ r.acc.total[:hw_storage])
# threshold year: first year with > 0.1% of the injected CO2 outside the storage sand (free +
# dissolved), a reporting threshold. First invasion: first year with > 0.01 t of free CO2 in the fault.
const OUTSIDE_REGIONS = (:fault, :caprock, :overburden, :underburden)
function metrics(r)
    inj = r.inventory[end][:injected] / 1000
    o = outside(r)
    thr = findfirst(>(0.001 * inj), o)
    inv = findfirst(>(0.01), r.acc.free[:fault])
    d = filter(isfinite, r.diag.drive)
    return (; out25 = 100 * o[25] / inj, out100 = 100 * o[end] / inj,
        thr = isnothing(thr) ? "none" : string(thr), inv = isnothing(inv) ? "none" : string(inv),
        free = 100 * sum(r.acc.free[k][end] for k in OUTSIDE_REGIONS) / inj,
        diss = 100 * sum(r.acc.diss[k][end] for k in OUTSIDE_REGIONS) / inj,
        above = 100 * (r.acc.total[:overburden][end] + r.exited[end]) / inj, drive_max = maximum(d))
end
rows = ["pe_fault_bar,k_across_mD,k_along_mD,step3_regime,fault_columns,outside_pct_yr25,outside_pct_yr100,threshold_0p1pct_year,first_free_co2_in_fault_year,free_outside_pct_yr100,dissolved_outside_pct_yr100,above_seal_pct_yr100,drive_max_bar"]
println("\nPe  | k across | regime                            | cols | outside % yr 25 / 100 | 0.1% yr | 1st free in fault | free / dissolved outside % | above seal % | max drive")
for (pe, k, regime) in scenarios, (ncol, dir) in ((1, dir1), (nsub, dir5))
    f = joinpath(dir, "$(runname(pe, k)).jls"); isfile(f) || continue
    m = metrics(deserialize(f))
    @printf("%.1f | %5.2f | %-33s | %2d | %7.3f / %7.3f | %5s | %5s | %7.3f / %6.3f | %7.3f | %.3f\n", pe, k[1], regime, ncol,
        m.out25, m.out100, m.thr, m.inv, m.free, m.diss, m.above, m.drive_max)
    push!(rows, join([pe, k[1], k[2], "\"$regime\"", ncol, round(m.out25, digits = 3), round(m.out100, digits = 3),
        m.thr, m.inv, round(m.free, digits = 3), round(m.diss, digits = 3), round(m.above, digits = 3), round(m.drive_max, digits = 3)], ","))
end
out = joinpath(@__DIR__, "..", "results", "03c_fault_resolution.csv")
write(out, join(rows, "\n") * "\n")
println("written ", out)
