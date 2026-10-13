# Stage B2c report: tables and figures from output/04d/*.jls
# (run scripts/04d_fault_gas.jl runs, and runs_bhp, first).
# Run from the project folder:  julia --project=. scripts/04d_report.jl

include(joinpath(@__DIR__, "..", "src", "gasfault.jl"))
using CairoMakie, Serialization, Printf

const OUT4D = joinpath(@__DIR__, "..", "output", "04d")
figdir = joinpath(@__DIR__, "..", "figures"); mkpath(figdir)
resdir = joinpath(@__DIR__, "..", "results"); mkpath(resdir)

branches = [
    ("a_methane_inject", "(a) methane curves, 1,500 t/yr"),
    ("b_co2scaled_inject", "(b) CO₂-scaled curves, 1,500 t/yr"),
    ("c_co2scaled_control", "(c) CO₂-scaled curves, no injection"),
    ("a_bhp_methane_inject", "(a′) methane curves, BHP-capped"),
    ("b_bhp_co2scaled_inject", "(b′) CO₂-scaled curves, BHP-capped"),
]
dep = deserialize(joinpath(OUT4D, "depletion.jls"))
B = Dict(n => deserialize(joinpath(OUT4D, "$n.jls")) for (n, _) in branches if isfile(joinpath(OUT4D, "$n.jls")))
avail = [(n, l) for (n, l) in branches if haskey(B, n)]

# ---- Depletion table -------------------------------------------------------
println("Depletion (10 yr at a 300 t/yr standard-volume target, then 1 yr shut-in)")
println("year | crest p (bar) | HW matched p (bar) | entry drive (bar) | CH4 in place (t)")
for r in dep.summary.rows
    @printf("%4.0f | %6.1f | %6.1f | %7.2f | %7.1f\n", r.year, r.p_crest, r.p_hw_matched, r.drive, r.ch4_total)
end

# ---- Branch summary ---------------------------------------------------------
yr_at(rows, y) = rows[argmin([abs(r.year - y) for r in rows])]
open(joinpath(resdir, "04d_branches.csv"), "w") do io
    println(io, "branch,co2_injected_t,peak_crest_bar,peak_drive_bar,pe_fault_bar,first_year_drive_gt_pe,",
        "co2_fault_y21_t,co2_beyond_y21_t,co2_beyond_y41_t,ch4_beyond_y21_t,ch4_beyond_y41_t,max_sg_fault,",
        "co2_max_resid_t,ch4_max_resid_t,ministeps_ok,ministeps_cut,wall_min")
    for (n, l) in avail
        s = B[n].summary; rows = s.rows
        # Cumulative well flow (ministep-integrated). In the branches the only well flow is
        # CO2 injection; runs made before the hook fix split it by the well's contents, so
        # the CO2 flow is taken as the sum over components (no change for later runs).
        flow_co2(r) = sum(r.well_flow)
        resid(r) = r.residual .+ r.well_flow .- [flow_co2(r), 0.0, 0.0]
        inj = flow_co2(rows[end])
        first_y = something(findfirst(r -> r.drive > r.pe_fault, rows), 0)
        fy = first_y == 0 ? "never" : @sprintf("%.0f", rows[first_y].year)
        # Conservation: largest |(reservoir + well contents) − start − well flow| over all report steps.
        co2_res = maximum(r -> abs(resid(r)[1]), rows)
        ch4_res = maximum(r -> abs(resid(r)[2]), rows)
        r21, r41 = yr_at(rows, 21), yr_at(rows, 41)
        @printf(io, "%s,%.1f,%.2f,%.3f,%.2f,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%d,%d,%.1f\n",
            n, inj, maximum(r -> r.p_crest, rows), maximum(r -> r.drive, rows), rows[1].pe_fault, fy,
            r21.co2.fault, r21.co2.beyond_seal, r41.co2.beyond_seal, r21.ch4.beyond_seal, r41.ch4.beyond_seal,
            maximum(r -> r.sg_max_fault, rows), co2_res, ch4_res, s.stats.ok, s.stats.cut, B[n].wall)
    end
end
print(read(joinpath(resdir, "04d_branches.csv"), String))
@printf("depletion: max |residual| CO2 %.4f, CH4 %.4f, water %.4f t; CH4 produced %.1f t\n",
    maximum(r -> abs(r.residual[1]), dep.summary.rows), maximum(r -> abs(r.residual[2]), dep.summary.rows),
    maximum(r -> abs(r.residual[3]), dep.summary.rows), -dep.summary.rows[end].well_flow[2])

# ---- What drives gas into the fault: brine-pressure difference vs capillary pressure ----
println("\nEntry drive = Δp_brine (crest − fault cell) + Pc(crest)   [bar]")
for (n, l) in avail
    println(l)
    for y in (11, 12, 13, 15, 21, 25, 31, 41)
        r = yr_at(B[n].summary.rows, y)
        abs(r.year - y) < 0.5 || continue
        @printf("  yr %2.0f: crest %6.1f, fault %6.1f, HW %6.1f | Δp_brine %6.2f + Pc %5.2f = %6.2f (Pe %.1f)\n",
            r.year, r.p_crest, r.p_fault, r.p_hw_matched, r.dp_brine, r.pc_crest, r.drive, r.pe_fault)
    end
end
for y in (0, 11)
    r = y == 0 ? dep.summary.rows[1] : dep.summary.rows[end]
    @printf("depletion yr %2.0f: Δp_brine %6.2f + Pc %5.2f = %6.2f\n", r.year, r.dp_brine, r.pc_crest, r.drive)
end

# ---- Tolerance sensitivity: same runs with JutulDarcy's default tol_mb = 1e-7 ----
old = joinpath(@__DIR__, "..", "output", "04d_tolmb1e-7")
if isdir(old)
    println("\nTolerance sensitivity (CO2 / CH4 beyond the seal at year 41, t): tol_mb 1e-9 vs 1e-7")
    for (n, l) in avail
        f = joinpath(old, "$n.jls"); isfile(f) || continue
        r_new, r_old = B[n].summary.rows[end], deserialize(f).summary.rows[end]
        @printf("  %-40s CO2 %8.1f vs %8.1f | CH4 %7.2f vs %7.2f\n", l, r_new.co2.beyond_seal, r_old.co2.beyond_seal,
            r_new.ch4.beyond_seal, r_old.ch4.beyond_seal)
    end
end

# ---- Fault resolution: 1, 5 and 10 columns across the 10 m fault zone ----
res_dirs = [(1, OUT4D), (5, OUT4D * "_fault5"), (10, OUT4D * "_fault10")]
res_rows = []
println("\nFault resolution (columns across the 10 m fault; first fault-cell centre at 5, 1, 0.5 m from the sand)")
println("branch | columns | peak drive | Δp_brine yr 15 | CO2 beyond seal yr 21 / 41 | CH4 beyond seal yr 41 | max Sg in fault")
for (n, l) in branches[1:2], (ncol, dir) in res_dirs
    f = joinpath(dir, "$n.jls"); isfile(f) || continue
    rows = deserialize(f).summary.rows
    r15, r21, r41 = yr_at(rows, 15), yr_at(rows, 21), yr_at(rows, 41)
    push!(res_rows, (; n, ncol, rows))
    @printf("%-20s | %2d | %6.2f | %6.2f | %8.1f / %8.1f | %7.2f | %.3f\n", n, ncol, maximum(r -> r.drive, rows),
        r15.dp_brine, r21.co2.beyond_seal, r41.co2.beyond_seal, r41.ch4.beyond_seal, maximum(r -> r.sg_max_fault, rows))
end

# ---- Year-by-year results for every run (plain CSV; output/*.jls are not version-controlled) ----
open(joinpath(resdir, "04d_timeseries.csv"), "w") do io
    regs = (:fw_storage, :hw_storage, :fault, :caprock, :beyond_seal)
    println(io, "run,fault_columns,year,p_crest_bar,p_fault_bar,p_hw_matched_bar,dp_brine_bar,pc_crest_bar,drive_bar,pe_fault_bar,sg_max_fault,",
        join(["co2_$(r)_t" for r in regs], ","), ",", join(["ch4_$(r)_t" for r in regs], ","),
        ",co2_total_t,ch4_total_t,well_flow_total_t,co2_residual_t,ch4_residual_t")
    function dump(name, ncol, rows; branch)
        for r in rows
            fl = branch ? sum(r.well_flow) : r.well_flow[2]              # branches: CO2 injected; depletion: CH4 (negative = produced)
            res = branch ? r.residual .+ r.well_flow .- [sum(r.well_flow), 0.0, 0.0] : r.residual
            @printf(io, "%s,%d,%.0f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%s,%s,%.3f,%.3f,%.3f,%.4f,%.4f\n",
                name, ncol, r.year, r.p_crest, r.p_fault, r.p_hw_matched, r.dp_brine, r.pc_crest, r.drive, r.pe_fault, r.sg_max_fault,
                join([@sprintf("%.3f", getproperty(r.co2, g)) for g in regs], ","), join([@sprintf("%.3f", getproperty(r.ch4, g)) for g in regs], ","),
                r.co2_total, r.ch4_total, fl, res[1], res[2])
        end
    end
    dump("depletion", 1, dep.summary.rows; branch = false)
    for (n, _) in branches, (ncol, dir) in res_dirs
        f = joinpath(dir, "$n.jls"); isfile(f) && dump(n, ncol, deserialize(f).summary.rows; branch = true)
    end
end
println("written results/04d_timeseries.csv")

# ---- Crossing of the sand–fault interface ----
# Gas leaves the footwall storage sand almost only through the fault face (the caprock takes ≤ 1.2 t,
# printed below), so CO2 outside the storage sand ≈ cumulative net crossing of that interface.
# First invasion = first report year with > 0.01 t of CO2 in the fault (yearly resolution).
println("\nCrossing of the sand–fault interface (CO2, t) and first invasion year")
println("run | columns | first year CO2 in fault | crossed by yr 15 / 21 / 41 | max CO2 in caprock")
for (n, _) in branches[1:2], (ncol, dir) in res_dirs
    f = joinpath(dir, "$n.jls"); isfile(f) || continue
    rows = deserialize(f).summary.rows
    crossed(r) = r.co2_total - r.co2.fw_storage - r.co2.hw_storage
    i = findfirst(r -> r.co2.fault > 0.01, rows)
    @printf("%-20s | %2d | %s | %7.1f / %7.1f / %7.1f | %.3f\n", n, ncol, isnothing(i) ? "none" : @sprintf("%.0f", rows[i].year),
        crossed(yr_at(rows, 15)), crossed(yr_at(rows, 21)), crossed(yr_at(rows, 41)), maximum(r -> r.co2.caprock, rows))
end

# ---- Figure 1: time series ------------------------------------------------------
colors = Makie.wong_colors()
fig = Figure(size = (1200, 800))
axp = Axis(fig[1, 1], title = "Pressure at the trap crest (solid) and opposite it across the fault (dashed)", xlabel = "year", ylabel = "bar")
axd = Axis(fig[1, 2], title = "Entry drive p_gas(crest) − p_brine(fault)", xlabel = "year", ylabel = "bar")
axc = Axis(fig[2, 1], title = "CO₂ beyond the seal", xlabel = "year", ylabel = "t")
axm = Axis(fig[2, 2], title = "CH₄ beyond the seal", xlabel = "year", ylabel = "t")
yd = [r.year for r in dep.summary.rows]
lines!(axp, yd, [r.p_crest for r in dep.summary.rows], color = :black, label = "depletion")
lines!(axp, yd, [r.p_hw_matched for r in dep.summary.rows], color = :black, linestyle = :dash)
lines!(axd, yd, [r.drive for r in dep.summary.rows], color = :black)
for (k, (n, l)) in enumerate(avail)
    rows = B[n].summary.rows; y = [r.year for r in rows]
    c = colors[k]
    lines!(axp, y, [r.p_crest for r in rows], color = c, label = l)
    lines!(axp, y, [r.p_hw_matched for r in rows], color = c, linestyle = :dash)
    lines!(axd, y, [r.drive for r in rows], color = c)
    lines!(axc, y, [r.co2.beyond_seal for r in rows], color = c)
    lines!(axm, y, [r.ch4.beyond_seal for r in rows], color = c)
end
p0 = dep.summary.rows[1]
hlines!(axd, [2.5], color = :gray40, linestyle = :dot); text!(axd, 0.5, 2.5, text = "fault Pe, methane curves (2.5)", fontsize = 10, align = (:left, :bottom))
hlines!(axd, [1.5], color = :gray60, linestyle = :dot); text!(axd, 0.5, 1.5, text = "fault Pe, CO₂-scaled (1.5)", fontsize = 10, align = (:left, :top))
vspan!(axd, 11, 21, color = (:gray, 0.1)); vspan!(axc, 11, 21, color = (:gray, 0.1)); vspan!(axm, 11, 21, color = (:gray, 0.1))
Legend(fig[3, 1:2], axp, orientation = :horizontal, nbanks = 2, labelsize = 11, framevisible = false)
Label(fig[4, 1:2], "Years 0–10 production, 10–11 shut-in, 11–21 CO₂ injection (shaded), 21–41 shut-in. Dissolution is not modelled.", fontsize = 11)
save(joinpath(figdir, "04d_timeseries.png"), fig)

# ---- Figure 2: gas saturation near the fault ----------------------------------------
sec = build_section()
near = findall(c -> sec.core[c] && 450 <= sec.xc[c] <= 950 && sec.zr[c] <= 90, eachindex(sec.facies))
section_axis(pos, title) = Axis3(pos, title = title, titlesize = 12, titlegap = 30,
    zreversed = true, azimuth = 1.5π, elevation = 0.0, aspect = (1.0, 0.05, 0.3),
    xlabel = "x (m)", zlabel = "Depth (m)", protrusions = (60, 10, 10, 30),
    ylabelvisible = false, yticklabelsvisible = false, yticksvisible = false,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)
mapped = [(n, l) for (n, l) in avail if !isempty(B[n].maps)]
fig = Figure(size = (1200, 260 * length(mapped) + 60))
for (row, (n, l)) in enumerate(mapped), (col, y) in enumerate(sort(collect(keys(B[n].maps))))
    ax = section_axis(fig[row, col], "$l: year $y")
    plot_cell_data!(ax, sec.mesh, B[n].maps[y].Sg, cells = near, colormap = :inferno, colorrange = (0, 0.8), shading = NoShading)
    plot_mesh_edges!(ax, sec.mesh, cells = findall(sec.facies .== FAULT_FACIES), color = :gray70)
end
Colorbar(fig[1:length(mapped), 3], colormap = :inferno, colorrange = (0, 0.8), label = "gas saturation")
save(joinpath(figdir, "04d_gas_maps.png"), fig)
println("figures written")

# ---- Figure 3: fault resolution --------------------------------------------------
if any(r -> r.ncol > 1, res_rows)
    fig = Figure(size = (1200, 420))
    axd = Axis(fig[1, 1], title = "Entry drive vs fault entry pressure", xlabel = "year", ylabel = "bar")
    axc = Axis(fig[1, 2], title = "CO₂ beyond the seal", xlabel = "year", ylabel = "t")
    styles = Dict(1 => :solid, 5 => :dash, 10 => :dot)
    for r in res_rows
        c = r.n == "a_methane_inject" ? colors[1] : colors[2]
        y = [x.year for x in r.rows]
        lab = (r.n == "a_methane_inject" ? "(a) methane curves" : "(b) CO₂-scaled") * ", $(r.ncol) column" * (r.ncol > 1 ? "s" : "")
        lines!(axd, y, [x.drive for x in r.rows], color = c, linestyle = styles[r.ncol], label = lab)
        lines!(axc, y, [x.co2.beyond_seal for x in r.rows], color = c, linestyle = styles[r.ncol])
    end
    hlines!(axd, [2.5], color = colors[1], linewidth = 0.8); hlines!(axd, [1.5], color = colors[2], linewidth = 0.8)
    vspan!(axd, 11, 21, color = (:gray, 0.1)); vspan!(axc, 11, 21, color = (:gray, 0.1))
    Legend(fig[1, 3], axd, labelsize = 10, framevisible = false)
    Label(fig[2, 1:3], "Thin horizontal lines: fault entry pressure (2.5 bar methane curves, 1.5 bar CO₂-scaled). Shaded: CO₂ injection.", fontsize = 11)
    save(joinpath(figdir, "04d_fault_resolution.png"), fig)
    println("fault-resolution figure written")
end
