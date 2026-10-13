# Step 3: fault-property sweep. Fault capillary entry pressure and fault
# permeability are varied independently (in Step 2 they were tied together by
# Leverett scaling, so the sand-rich vs clay-rich comparison changed both).
#
# 12 entry pressures × 3 permeability levels = 36 runs on the 10 m × 2 m grid,
# which Step 2c showed is adequate for the leakage threshold and for total
# leakage out of the storage formation (not for trapping fractions).
#
# Restartable: finished runs (output/03_sweep/<name>.jls) are skipped.
# Run from the project folder:  julia --project=. scripts/03_sweep.jl

include(joinpath(@__DIR__, "..", "src", "run.jl"))

# Entry pressures (bar), closely spaced around the expected threshold: the
# entry drive at the fault is ~1.1 bar from buoyancy alone and ~1.7 bar at the
# end of injection (Step 2c, clay-rich fault).
pe_values = [0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0, 5.5]

# Fault permeability levels, mD (across, along).
k_levels = [(0.01, 0.1), (0.1, 1.0), (1.0, 10.0)]

outdir = joinpath(@__DIR__, "..", "output", "03_sweep")
mkpath(outdir)
runname(pe, k) = "pe$(pe)_ka$(k[1])_kl$(k[2])"

failed = String[]
for k in k_levels, pe in pe_values
    name = runname(pe, k)
    isfile(joinpath(outdir, "$name.jls")) && (println("skip $name (done)"); continue)
    sec = build_section(section_spec(k_across = k[1], k_along = k[2]))
    cs = case_spec(capillary = true, pe_fault_bar = pe)
    try
        run_case(name, sec, cs, output_dir = outdir, save_states = :final)
    catch err
        println("FAILED $name: ", sprint(showerror, err)[1:min(end, 300)])
        push!(failed, name)
    end
    GC.gc()
end

# Summary table (CSV), one row per run.
outside(r) = [i[:injected] / 1000 for i in r.inventory] .-
             (r.acc.total[:fw_storage] .+ r.acc.total[:hw_storage])
# threshold_0p1pct_year: first year with > 0.1% of the injected CO2 outside the storage sand
# (free + dissolved); a reporting threshold, not the first entry into the fault.
# first_free_co2_in_fault_year: first year with > 0.01 t of free (CO2-rich phase) CO2 in the fault.
# free/dissolved_outside: CO2 in the fault, caprock, overburden and underburden (CO2 that left
# through the top is reported separately as left_top).
const OUTSIDE_REGIONS = (:fault, :caprock, :overburden, :underburden)
first_free_in_fault(r) = (i = findfirst(>(0.01), r.acc.free[:fault]); isnothing(i) ? "none" : i)
free_outside(r) = sum(r.acc.free[k][end] for k in OUTSIDE_REGIONS)
diss_outside(r) = sum(r.acc.diss[k][end] for k in OUTSIDE_REGIONS)
rows = ["pe_fault_bar,k_across_mD,k_along_mD,outside_pct_yr25,outside_pct_yr100,threshold_0p1pct_year,first_free_co2_in_fault_year,free_outside_pct_yr100,dissolved_outside_pct_yr100,above_seal_pct_yr100,left_top_pct_yr100,drive_max_bar,drive_yr25_bar,drive_yr100_bar,dp_max_bar"]
for k in k_levels, pe in pe_values
    f = joinpath(outdir, "$(runname(pe, k)).jls")
    isfile(f) || continue
    r = deserialize(f)
    inj = r.inventory[end][:injected] / 1000
    o = outside(r)
    onset = findfirst(>(0.001 * inj), o)   # 0.1% reporting threshold (free + dissolved)
    above = r.acc.total[:overburden][end] + r.exited[end]
    d = filter(isfinite, r.diag.drive)
    push!(rows, join([pe, k[1], k[2],
        round(100 * o[25] / inj, digits = 3), round(100 * o[end] / inj, digits = 3),
        isnothing(onset) ? "none" : onset, first_free_in_fault(r),
        round(100 * free_outside(r) / inj, digits = 3), round(100 * diss_outside(r) / inj, digits = 3),
        round(100 * above / inj, digits = 3), round(100 * r.exited[end] / inj, digits = 3),
        round(maximum(d), digits = 3), round(r.diag.drive[25], digits = 3), round(r.diag.drive[end], digits = 3),
        round(maximum(r.diag.dp_max), digits = 3)], ","))
end
write(joinpath(@__DIR__, "..", "results", "03_sweep_summary.csv"), join(rows, "\n") * "\n")
println("\nWrote results/03_sweep_summary.csv ($(length(rows) - 1) runs)")
isempty(failed) || println("Failed runs: ", join(failed, ", "))
