# All figures in figures/, from saved results (no simulations are run here,
# except Step 4A's flash tables, which are read from results/).
# Run from the project folder:
#   julia --project=. scripts/make_figures.jl            # every figure
#   julia --project=. scripts/make_figures.jl 01 04d     # figures whose name starts with these prefixes

include(joinpath(@__DIR__, "..", "src", "run.jl"))        # section, case, run helpers
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
using Serialization, Printf, DelimitedFiles
set_theme!(PUB_THEME)

const ROOT = joinpath(@__DIR__, "..")
const FIG = joinpath(ROOT, "figures")
const OUT = joinpath(ROOT, "output")
loadout(parts...) = deserialize(joinpath(OUT, parts...))
yrs(r) = r.t ./ yr

# --------------------------------------------------------------------------- Step 1
function fig_01_section()
    sec = build_section()
    polys = section_polygons(sec)
    xl, dl = (0.0, 1200.0), (1455.0, 1625.0)
    cells = window_cells(sec, xl, dl)
    fig = pub_figure(height = 360)
    d(x, z) = sec.spec.top_depth + z - sec.slope * x       # depth of a point z below the top surface
    # (a) facies
    ax = section_axis(fig[1, 1]; xlims = xl, dlims = dl, ve = 3, title = "Facies; hanging wall (right) downthrown 30 m")
    section_facies!(ax, sec, polys; cells)
    fault_outline!(ax, sec, polys)
    names = [[F.name for F in FACIES]..., "Fault zone (10 m)"]
    Legend(fig[1, 2], [PolyElement(color = c, strokecolor = (:black, 0.3), strokewidth = 0.4) for c in FACIES_COLORS], names,
        tellheight = false, valign = :center)
    # Injector of Steps 2–3.
    wx = 345.0
    lines!(ax, [wx, wx], [d(wx, -6.0), d(wx, 84.0)], color = :black, linewidth = 0.8)
    lines!(ax, [wx, wx], [d(wx, 64.0), d(wx, 84.0)], color = :black, linewidth = 2.2)
    text!(ax, wx + 8, d(wx, 74.0), text = "injector, Steps 2–3", fontsize = 7, align = (:left, :center))
    # The sand–sand window: footwall storage sand facing hanging-wall overburden sand.
    fx = sec.spec.fault_x[2]
    wtop, wbot = d(fx, 44.0), d(fx, 54.0)
    lines!(ax, [fx + 3, fx + 3], [wtop, wbot], color = OI.vermillion, linewidth = 2)
    lines!(ax, [fx + 5, 860.0], [(wtop + wbot) / 2, d(860.0, 8.0)], color = OI.vermillion, linewidth = 0.6)
    text!(ax, 864.0, d(860.0, 8.0), text = "10 m sand–sand window", fontsize = 7, color = OI.vermillion, align = (:left, :center))
    panel_label!(fig[1, 1], "a")
    # (b) horizontal permeability
    ax = section_axis(fig[2, 1]; xlims = xl, dlims = dl, ve = 3, title = "Horizontal permeability")
    logk = log10.(sec.perm[1, :] ./ mD)
    section_field!(ax, polys, logk; cells, colormap = CMAP_LOGK, colorrange = (-3, 2.5))
    fault_outline!(ax, sec, polys)
    Colorbar(fig[2, 2], colormap = CMAP_LOGK, colorrange = (-3, 2.5), label = "log₁₀ k (mD)", height = Relative(0.7),
        tellheight = false, halign = :left)
    panel_label!(fig[2, 1], "b")
    save_fig(joinpath(FIG, "01_section.png"), fig)
end

# ------------------------------------------------------- shared: saturation maps
"""
Stacked cross-sections of a gas-saturation-like field over faded geology.
`panels` = [(sec, values, title)]; only cells above `threshold` are coloured.
"""
function field_panels(panels, path; xl, dl, ve, colormap = CMAP_SG, colorrange = (0, 0.8),
        label = "CO₂ saturation", threshold = 0.005, rowh = 120, ncols = 1, width = PT_FULL, letters = 'a':'z', caption = "")
    nrows = cld(length(panels), ncols)
    fig = pub_figure(; width, height = nrows * rowh + 30)
    for (i, (sec, vals, title)) in enumerate(panels)
        row, col = fldmod1(i, ncols)
        polys = section_polygons(sec)
        cells = window_cells(sec, xl, dl)
        ax = section_axis(fig[row, col]; xlims = xl, dlims = dl, ve, title, ve_in_title = ncols == 1,
            xlabelvisible = row == nrows, xticklabelsvisible = row == nrows,
            ylabelvisible = col == 1)
        section_facies!(ax, sec, polys; cells, fade = 0.55)
        gas = [c for c in cells if vals[c] > threshold]
        section_field!(ax, polys, vals; cells = gas, colormap, colorrange)
        fault_outline!(ax, sec, polys)
        panel_label!(fig[row, col], string(letters[i]))
    end
    Colorbar(fig[1:nrows, ncols + 1], colormap = colormap, colorrange = colorrange, label = label,
        height = Relative(min(1.0, 0.25 * nrows + 0.3)), tellheight = false)
    isempty(caption) || Label(fig[nrows + 1, 1:ncols], caption, fontsize = 7, color = :grey30, tellwidth = false)
    rowgap!(fig.layout, 6)
    save_fig(path, fig)
end

"Facies legend entries, for figures that show faded geology behind a field."
faded_facies_note(fig_pos) = Label(fig_pos, "Background: facies (faded), as in Fig. 01; fault zone outlined.", fontsize = 7, color = :grey35, halign = :left, tellwidth = false)

# --------------------------------------------------------------------------- Step 2
const STEP2_CASES = [("02_nopc_clay", "No capillary pressure, clay-rich fault", OI.grey),
                     ("02_pc_clay", "Capillary pressure, clay-rich fault (Pe 5.5 bar)", OI.blue),
                     ("02_pc_sandy", "Capillary pressure, sand-rich fault (Pe 0.55 bar)", OI.vermillion)]

function fig_02_saturation()
    panels = [(build_section(r.spec), r.states[end][:Saturations][2, :], "$label, year 100")
              for (n, label, _) in STEP2_CASES for r in (loadout("$n.jls"),)]
    field_panels(panels, joinpath(FIG, "02_saturation_year100.png"); xl = (250.0, 1200.0), dl = (1458.0, 1585.0), ve = 2.5, rowh = 125)
end

function fig_02_accounting()
    cats = [(:fault, "Fault zone"), (:caprock, "Caprock"), (:overburden, "Overburden")]
    fig = pub_figure(height = 190)
    ax = Axis(fig[1, 1], xticks = (1:4, [[l for (_, l) in cats]..., "Left the model\n(through the top)"]),
        ylabel = "CO₂ (% of injected)", title = "CO₂ outside the storage formation at year 100")
    w = 0.26
    for (j, (n, label, c)) in enumerate(STEP2_CASES)
        r = loadout("$n.jls")
        inj = r.inventory[end][:injected] / 1000
        xs = (1:4) .+ (j - 2) * (w + 0.02)
        free = [100 * r.acc.free[k][end] / inj for (k, _) in cats]; push!(free, 100 * r.exited[end] / inj)
        diss = [100 * r.acc.diss[k][end] / inj for (k, _) in cats]; push!(diss, 0.0)
        barplot!(ax, xs, free .+ diss, width = w, color = (c, 0.35), strokewidth = 0)
        barplot!(ax, xs, free, width = w, color = c, strokewidth = 0)
    end
    ylims!(ax, 0, nothing)
    elems = [[PolyElement(color = c) for (_, _, c) in STEP2_CASES]; [PolyElement(color = :grey40), PolyElement(color = (:grey40, 0.35))]]
    labs = [[l for (_, l, _) in STEP2_CASES]; ["free (CO₂-rich phase)", "dissolved in brine"]]
    labs[1:3] = ["No Pc, clay-rich fault", "Pc, clay-rich fault (Pe 5.5 bar)", "Pc, sand-rich fault (Pe 0.55 bar)"]
    Legend(fig[1, 2][1, 1], elems, labs, framevisible = false, valign = :bottom)
    Label(fig[1, 2][2, 1], "Pc, clay-rich fault: 0 in
every category (no bars)", fontsize = 7, color = OI.blue, halign = :left, valign = :top, justification = :left)
    save_fig(joinpath(FIG, "02_accounting_year100.png"), fig)
end

function fig_02_fault_pressure()
    res = Dict(n => loadout("$n.jls") for (n, _, _) in STEP2_CASES)
    fig = pub_figure(height = 250)
    ax1 = Axis(fig[1, 1], title = "Entry drive at the fault", xlabel = "Time (years)",
        ylabel = "p_CO₂(sand) − p_brine(fault) (bar)")
    inj_span!(ax1, 0, 25)
    for (n, label, c) in STEP2_CASES[2:3]
        r = res[n]; pe = entry_pressures(build_section(r.spec), r.cs)[FAULT_FACIES]
        lines!(ax1, yrs(r), r.diag.drive, color = c)
        hlines!(ax1, [pe], color = c, linestyle = :dash, linewidth = 0.9)
        # Label below its line, on the right, where no curve passes.
        text!(ax1, 98, pe, text = "Pe = $(round(pe, digits = 2)) bar", color = c, fontsize = 7, align = (:right, :top), offset = (0, -1))
    end
    ylims!(ax1, 0, 6.2)
    panel_label!(fig[1, 1], "a")
    ax2 = Axis(fig[1, 2], title = "Pressure rise from injection", xlabel = "Time (years)", ylabel = "Δp (bar)")
    inj_span!(ax2, 0, 25)
    for (n, _, c) in STEP2_CASES
        lines!(ax2, yrs(res[n]), res[n].diag.dp_max, color = c)
    end
    lines!(ax2, yrs(res["02_pc_clay"]), res["02_pc_clay"].diag.overpressure_fault, color = OI.blue, linestyle = :dot)
    panel_label!(fig[1, 2], "b")
    Legend(fig[2, 1:2], [LineElement(color = c) for (_, _, c) in STEP2_CASES],
        ["No Pc, clay-rich fault", "Pc, clay-rich fault (Pe 5.5 bar)", "Pc, sand-rich fault (Pe 0.55 bar)"],
        orientation = :horizontal, tellwidth = false)
    Label(fig[3, 1:2], "a: dashed, fault entry pressure.   b: solid, maximum Δp anywhere; dotted, Δp in the clay-rich fault (Pc case).   Shaded: injection.",
        fontsize = 7, color = :grey30, tellwidth = false)
    rowgap!(fig.layout, 2)
    save_fig(joinpath(FIG, "02_fault_pressure.png"), fig)
end

outside_pct_series(r) = (inj = [i[:injected] / 1000 for i in r.inventory];
    100 .* (inj .- (r.acc.total[:fw_storage] .+ r.acc.total[:hw_storage])) ./ max.(inj, 1e-9))

function fig_02_grid()
    coarse, fine = loadout("02_pc_sandy.jls"), loadout("02_pc_sandy_fine.jls")
    pe = entry_pressures(build_section(coarse.spec), coarse.cs)[FAULT_FACIES]
    fig = pub_figure(height = 200)
    ax1 = Axis(fig[1, 1], title = "Sand-rich fault: CO₂ outside the storage formation", xlabel = "Time (years)", ylabel = "% of injected")
    ax2 = Axis(fig[1, 2], title = "Sand-rich fault: entry drive", xlabel = "Time (years)", ylabel = "p_CO₂(sand) − p_brine(fault) (bar)")
    for ax in (ax1, ax2); inj_span!(ax, 0, 25); end
    for (r, ls, lab) in ((coarse, :solid, "10 m × 2 m cells"), (fine, :dash, "5 m × 1 m cells"))
        lines!(ax1, yrs(r), outside_pct_series(r), color = OI.vermillion, linestyle = ls, label = lab)
        lines!(ax2, yrs(r), r.diag.drive, color = OI.vermillion, linestyle = ls)
    end
    hlines!(ax2, [pe], color = :grey40, linestyle = :dot, linewidth = 0.9)
    text!(ax2, 99, pe, text = "fault Pe = $(round(pe, digits = 2)) bar", fontsize = 7, color = :grey30, align = (:right, :top))
    axislegend(ax1, position = :rb)
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    save_fig(joinpath(FIG, "02_grid_refinement.png"), fig)
end

# --------------------------------------------------------------------------- Step 3
const PE_VALUES = [0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0, 5.5]
const K_LEVELS = [((0.01, 0.1), OI.blue), ((0.1, 1.0), OI.orange), ((1.0, 10.0), OI.green)]
runname3(pe, k) = "pe$(pe)_ka$(k[1])_kl$(k[2])"
klabel(k) = "k = $(k[1]) / $(k[2]) mD (across / along)"
load3(pe, k) = (f = joinpath(OUT, "03_sweep", "$(runname3(pe, k)).jls"); isfile(f) ? deserialize(f) : nothing)

function fig_03_leakage()
    R = Dict((pe, k) => load3(pe, k) for pe in PE_VALUES for (k, _) in K_LEVELS)
    sealed = [R[(5.5, k)] for (k, _) in K_LEVELS]
    post = maximum(r.diag.drive[end] for r in sealed)
    peak = maximum(maximum(filter(isfinite, r.diag.drive)) for r in sealed)
    fig = pub_figure(height = 215)
    for (col, (step, title)) in enumerate([(25, "End of injection (year 25)"), (100, "Year 100")])
        ax = Axis(fig[1, col], title = title, xscale = log10, xticks = [0.4, 0.6, 1.0, 1.7, 3.0, 5.5],
            xlabel = "Fault capillary entry pressure (bar)", ylabel = col == 1 ? "CO₂ outside storage (% of injected)" : "")
        vspan!(ax, post, peak, color = (:grey, 0.15))
        col == 1 && text!(ax, sqrt(post * peak), 0.97, text = "entry drive\nwith a sealing\nfault", space = :relative,
            fontsize = 7, color = :grey30, align = (:center, :top))
        for (k, c) in K_LEVELS
            pes = [pe for pe in PE_VALUES if !isnothing(R[(pe, k)])]
            scatterlines!(ax, pes, [outside_pct_series(R[(pe, k)])[step] for pe in pes], color = c, markersize = 4)
        end
        ylims!(ax, -0.5, nothing)
        panel_label!(fig[1, col], string('a' + col - 1))
    end
    Legend(fig[2, 1:2], [LineElement(color = c) for (_, c) in K_LEVELS], ["$(k[1]) / $(k[2]) mD" for (k, _) in K_LEVELS],
        "Fault permeability (across / along):", orientation = :horizontal, tellwidth = false,
        titleposition = :left, titlefont = :regular, titlesize = 8)
    save_fig(joinpath(FIG, "03_leakage_vs_entry_pressure.png"), fig)
end

function fig_03_time_series()
    kmid = K_LEVELS[2][1]
    pes = [0.4, 0.8, 1.0, 1.2, 1.4, 2.0]
    cm = cgrad(:viridis, length(pes), categorical = true)[1:end]
    cols = [get(cgrad(:viridis), 0.85 * (i - 1) / (length(pes) - 1)) for i in eachindex(pes)]
    fig = pub_figure(height = 215)
    ax1 = Axis(fig[1, 1], title = "CO₂ outside the storage formation", xlabel = "Time (years)", ylabel = "% of injected")
    ax2 = Axis(fig[1, 2], title = "Entry drive at the fault", xlabel = "Time (years)", ylabel = "p_CO₂(sand) − p_brine(fault) (bar)")
    for ax in (ax1, ax2); inj_span!(ax, 0, 25); end
    ends = Float64[]
    for (pe, c) in zip(pes, cols)
        r = load3(pe, kmid)
        y = outside_pct_series(r)
        lines!(ax1, yrs(r), y, color = c)
        push!(ends, y[end])
        lines!(ax2, yrs(r), r.diag.drive, color = c)
        hlines!(ax2, [pe], color = c, linestyle = :dash, linewidth = 0.7)
    end
    # End-of-line labels, spread so they do not overlap (gap ≈ one text line).
    ytop = maximum(ends) * 1.12
    for (pe, c, yl) in zip(pes, cols, spread_labels(ends, 0.075 * ytop))
        text!(ax1, 102, yl, text = "$pe", fontsize = 7, color = c, align = (:left, :center))
    end
    text!(ax1, 102, ytop, text = "Pe (bar)", fontsize = 7, align = (:left, :center))
    xlims!(ax1, 0, 118); ylims!(ax1, -0.2 * ytop, ytop * 1.05)
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    Label(fig[2, 1:2], "Fault $(klabel(kmid)). Line colour: fault entry pressure Pe (labelled in a; dashed lines in b). Shaded: injection.",
        fontsize = 7, color = :grey30, tellwidth = false)
    save_fig(joinpath(FIG, "03_time_series.png"), fig)
end

function fig_03_saturation()
    kmid = K_LEVELS[2][1]
    panels = [(build_section(r.spec), r.states[end][:Saturations][2, :], "Fault Pe = $pe bar, year 100")
              for pe in (0.6, 1.2, 2.0) for r in (load3(pe, kmid),)]
    field_panels(panels, joinpath(FIG, "03_saturation_examples.png"); xl = (250.0, 1200.0), dl = (1458.0, 1585.0), ve = 2.5, rowh = 125)
end

# -------------------------------------------------------------------------- Step 3b
load3b(n) = (f = joinpath(OUT, "03b", "$n.jls"); isfile(f) ? deserialize(f) : nothing)

function fig_03b_grid()
    cases = [("grid10x2_pe1.2_mid", "10 × 2 m (base)", OI.black, :solid), ("grid10x1_pe1.2_mid", "10 × 1 m", OI.blue, :solid),
             ("grid5x2_pe1.2_mid", "5 × 2 m", OI.orange, :solid), ("grid5x1_pe1.2_mid", "5 × 1 m", OI.green, :solid),
             ("dt30_pe1.2_mid", "10 × 2 m, max Δt 30 days", OI.purple, :dash)]
    fig = pub_figure(height = 215)
    ax1 = Axis(fig[1, 1], title = "CO₂ out of the storage formation", xlabel = "Time (years)", ylabel = "CO₂ (t)")
    ax2 = Axis(fig[1, 2], title = "Cumulative net crossing of the seal", xlabel = "Time (years)", ylabel = "CO₂ (t)")
    for ax in (ax1, ax2); inj_span!(ax, 0, 25); end
    for (n, lab, c, ls) in cases
        r = load3b(n); isnothing(r) && continue
        lines!(ax1, r.years, r.out_of_storage, color = c, linestyle = ls, label = lab)
        lines!(ax2, r.years, r.crossed_seal, color = c, linestyle = ls)
    end
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    Legend(fig[2, 1:2], ax1, "Cell size (fault Pe 1.2 bar, k 0.1 / 1 mD)", orientation = :horizontal, nbanks = 2, tellwidth = false, titleposition = :top, titlefont = :regular, titlesize = 8)
    save_fig(joinpath(FIG, "03b_grid_timestep.png"), fig)
end

function fig_03b_rate()
    rates = [(750, 50, OI.blue), (1500, 25, OI.orange), (2500, 15, OI.green)]
    fig = pub_figure(height = 235)
    ax1 = Axis(fig[1, 1], title = "Same total mass (37,500 t), different rate", xlabel = "Time (years)", ylabel = "CO₂ out of storage (t)")
    ax2 = Axis(fig[1, 2], title = "Entry drive at the fault", xlabel = "Time (years)", ylabel = "p_CO₂(sand) − p_brine(fault) (bar)")
    for (rate, yi, c) in rates
        r = load3b("rate$(rate)_pe1.4_mid"); isnothing(r) && continue
        lines!(ax1, r.years, r.out_of_storage, color = c)
        lines!(ax2, r.years, r.drive, color = c)
        for ax in (ax1, ax2); vlines!(ax, [yi], color = c, linestyle = :dot, linewidth = 0.8); end
    end
    hlines!(ax2, [1.4], color = :grey40, linestyle = :dash, linewidth = 0.9)
    text!(ax2, 99, 1.4, text = "fault Pe = 1.4 bar", fontsize = 7, color = :grey30, align = (:right, :bottom))
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    Legend(fig[2, 1:2], [[LineElement(color = c) for (_, _, c) in rates]; LineElement(color = :grey40, linestyle = :dot)],
        [["$rate t/yr for $yi years" for (rate, yi, _) in rates]; "end of injection"], orientation = :horizontal, tellwidth = false)
    save_fig(joinpath(FIG, "03b_rate.png"), fig)
end

# --------------------------------------------------------------------------- Step 4A
function readcsv(name)
    d, h = readdlm(joinpath(ROOT, "results", name), ',', header = true)
    return (; d, h = vec(String.(h)))
end
col(t, name) = Float64.(t.d[:, findfirst(==(name), t.h)])
const SAL_COLORS = [get(cgrad(:deep), v) for v in (0.35, 0.6, 0.88)]   # 0, 30,000, 55,000 mg/L: light to dark
fmt_sal(s) = s == 0 ? "0 mg/L" : string(round(Int, s ÷ 1000), ",000 mg/L")

function fig_04a_co2_brine()
    t = readcsv("04a_co2_brine.csv")
    sal, p = col(t, "salinity_mgL"), col(t, "p_bar")
    fig = pub_figure(height = 215)
    ax1 = Axis(fig[1, 1], title = "CO₂ solubility in brine, 95 °C", xlabel = "Pressure (bar)", ylabel = "CO₂ (mol per kg water)")
    ax2 = Axis(fig[1, 2], title = "Brine (liquid) density, 95 °C", xlabel = "Pressure (bar)", ylabel = "Density (kg/m³)")
    for (s, c) in zip(sort(unique(sal)), SAL_COLORS)
        i = findall(==(s), sal)
        lines!(ax1, p[i], col(t, "mCO2_table")[i], color = c)
        scatterlines!(ax1, p[i], col(t, "mCO2_SW")[i], color = c, linestyle = :dash, markersize = 4)
        text!(ax1, p[i][end] + 4, col(t, "mCO2_table")[i][end], text = fmt_sal(s), color = c, fontsize = 7, align = (:left, :center))
        lines!(ax2, p[i], col(t, "rho_brine_table")[i], color = c)
        scatterlines!(ax2, p[i], col(t, "rho_liq_SW")[i], color = c, linestyle = :dash, markersize = 4)
    end
    xlims!(ax1, 40, 400)
    text!(ax2, 175, 840, text = "Søreide–Whitson (Peng–Robinson,\nno volume shift): 15–19% low", fontsize = 7, color = :grey30, align = (:center, :bottom))
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    Legend(fig[2, 1:2], [LineElement(color = :grey30), [LineElement(color = :grey30, linestyle = :dash), MarkerElement(color = :grey30, marker = :circle, markersize = 4)]],
        ["JutulDarcy CO₂–brine correlations (used in Steps 2–3)", "Søreide–Whitson equation of state"], orientation = :horizontal, tellwidth = false)
    Label(fig[3, 1:2], "Colour: NaCl-equivalent salinity (labelled in a).", fontsize = 7, color = :grey30, tellwidth = false)
    rowgap!(fig.layout, 2)
    save_fig(joinpath(FIG, "04a_co2_brine.png"), fig)
end

function fig_04a_ch4_water()
    t = readcsv("04a_ch4_water.csv")
    sal, p = col(t, "salinity_mgL"), col(t, "p_bar")
    fig = pub_figure(width = 320, height = 190)
    ax = Axis(fig[1, 1], title = "CH₄ solubility, 95 °C (Søreide–Whitson)", xlabel = "Pressure (bar)", ylabel = "CH₄ (mol per kg water)")
    for (s, c) in zip(sort(unique(sal)), SAL_COLORS)
        i = findall(==(s), sal)
        lines!(ax, p[i], col(t, "mCH4_published_Tr2")[i], color = c)
        lines!(ax, p[i], col(t, "mCH4_implemented_Tr")[i], color = c, linestyle = :dash)
        text!(ax, p[i][end] + 4, col(t, "mCH4_published_Tr2")[i][end], text = fmt_sal(s), color = c, fontsize = 7, align = (:left, :center))
    end
    xlims!(ax, 40, 395); axislegend(ax, [LineElement(color = :grey30), LineElement(color = :grey30, linestyle = :dash)], ["published form (T_r²)", "MultiComponentFlash 1.1.19 (T_r)"], position = (1.0, 0.24))
    text!(ax, 110, 0.016, text = "≈ 10× lower", fontsize = 7, color = :grey30, align = (:left, :bottom))
    save_fig(joinpath(FIG, "04a_ch4_water.png"), fig)
end

# --------------------------------------------------------------------------- Step 4B
function fig_04b_depletion()
    f = joinpath(OUT, "04b", "depletion_1d_figdata.jls")
    isfile(f) || (println("skip 04b: run scripts/04b_depletion_1d.jl first"); return)
    D = deserialize(f)
    y = D[100].years
    t1, t2, t3 = y.deplete, y.deplete + y.shut1, y.deplete + y.shut1 + y.inject
    styles = Dict(100 => (OI.blue, :solid, "100 cells (10 m)"), 400 => (OI.vermillion, :dash, "400 cells (2.5 m)"))
    fig = pub_figure(height = 200)
    ax1 = Axis(fig[1, 1], title = "Average pressure", xlabel = "Time (years)", ylabel = "Pressure (bar)")
    vspan!(ax1, t2, t3, color = (:grey, 0.12))
    for nx in (100, 400)
        c, ls, _ = styles[nx]; lines!(ax1, D[nx].t_years, D[nx].pavg, color = c, linestyle = ls)
    end
    text!(ax1, t1 / 2, 0.04, text = "CH₄\nproduction", space = :relative, fontsize = 7, color = :grey30, align = (:center, :bottom))
    text!(ax1, (t2 + t3) / 2 / sum(y), 0.04, text = "CO₂\ninjection", space = :relative, fontsize = 7, color = :grey30, align = (:center, :bottom))
    # (b) p/Z material balance
    ax2 = Axis(fig[1, 2], title = "p/Z check during production", xlabel = "CH₄ produced (kt)", ylabel = "p/Z (bar)")
    rows = D[100].pz
    gp = [r[3] / 1000 for r in rows]; pzs = [r[4] for r in rows] ./ bar; pzm = [r[5] for r in rows] ./ bar
    lines!(ax2, gp, pzm, color = :grey40, label = "material balance")
    scatter!(ax2, gp, pzs, color = OI.blue, markersize = 5, label = "simulated (100 cells)")
    text!(ax2, 0.04, 0.05, text = "max difference 0.12%\n(= connate-water expansion)", space = :relative, fontsize = 7, color = :grey30, align = (:left, :bottom))
    # (c) composition profile at the end of injection
    ax3 = Axis(fig[1, 3], title = "Gas composition, year $t3", xlabel = "Distance from injector (m)", ylabel = "CO₂ mole fraction in gas")
    for nx in (100, 400)
        c, ls, _ = styles[nx]; lines!(ax3, D[nx].x, D[nx].profiles[t3], color = c, linestyle = ls)
    end
    for (i, l) in enumerate(("a", "b", "c")); panel_label!(fig[1, i], l); end
    Legend(fig[2, 1:3], [[LineElement(color = styles[n][1], linestyle = styles[n][2]) for n in (100, 400)];
            [LineElement(color = :grey40), MarkerElement(color = OI.blue, marker = :circle, markersize = 5)]],
        [[styles[n][3] for n in (100, 400)]; ["p/Z material balance (b)", "p/Z simulated, 100 cells (b)"]],
        orientation = :horizontal, tellwidth = false)
    save_fig(joinpath(FIG, "04b_depletion_1d.png"), fig)
end

# -------------------------------------------------------------------------- Step 4B2
function fig_04b2()
    f = joinpath(OUT, "04b2", "gravity_figdata.jls")
    isfile(f) || (println("skip 04b2: run scripts/04b2_gas_2d.jl gravity first"); return)
    G = deserialize(f)
    polys = tilted_cell_polygons(collect(range(0, G.L, G.NX + 1)), G.H / G.NZ, G.NZ, G.top, G.slope)
    xl, dl = (0.0, G.L), (G.top - G.slope * G.L - 3, G.top + G.H + 3)
    fig = pub_figure(height = 215)
    for (col, g) in enumerate((true, false)), (row, y) in enumerate((21, 51))
        m = G.data[g].maps[y]
        ax = section_axis(fig[row, col]; xlims = xl, dlims = dl, ve = 3, title = "Gravity $(g ? "on" : "off"), year $y",
            xlabelvisible = row == 2, xticklabelsvisible = row == 2, ylabelvisible = col == 1)
        water = findall(>=(0.99), m.sw); gas = findall(<(0.99), m.sw)
        poly!(ax, polys[water], color = RGBf(0.93, 0.93, 0.93), strokecolor = RGBf(0.93, 0.93, 0.93), strokewidth = 0.25)
        section_field!(ax, polys, m.zco2; cells = gas, colormap = CMAP_ZCO2, colorrange = (0, 1))
        panel_label!(fig[row, col], string("abcd"[2 * (col - 1) + row]))
    end
    Colorbar(fig[1:2, 3], colormap = CMAP_ZCO2, colorrange = (0, 1), label = "CO₂ mole fraction in gas", height = Relative(0.8))
    Label(fig[3, 1:2], "Grey: water zone. CO₂ injected at the crest (right) after methane production; years 11–21 injection, then shut-in.",
        fontsize = 7, color = :grey30, tellwidth = false)
    rowgap!(fig.layout, 4)
    save_fig(joinpath(FIG, "04b2_gravity.png"), fig)

    fig = pub_figure(width = 300, height = 200)
    ax = Axis(fig[1, 1], title = "Separation of CO₂ and CH₄", xlabel = "Time (years)", ylabel = "Mean depth CO₂ − CH₄ (m)")
    vspan!(ax, 11, 21, color = (:grey, 0.12)); hlines!(ax, [0], color = :grey60, linewidth = 0.6)
    for (g, c) in ((true, OI.blue), (false, OI.vermillion))
        d = G.data[g]; lines!(ax, d.years, d.separation, color = c)
        text!(ax, d.years[end], d.separation[end], text = "gravity $(g ? "on" : "off")", color = c, fontsize = 7, align = (:right, g ? :bottom : :top), offset = (0, g ? 2 : -2))
    end
    text!(ax, 0.98, 0.56, text = "positive: CO₂ deeper", space = :relative, fontsize = 7, color = :grey30, align = (:right, :bottom))
    ylims!(ax, -21, 18)
    save_fig(joinpath(FIG, "04b2_segregation.png"), fig)
end

# --------------------------------------------------------------------------- Stage 4B-2c (04d)
const BRANCHES = [("a_methane_inject", "(a) methane curves", OI.blue), ("b_co2scaled_inject", "(b) CO₂-scaled curves", OI.vermillion),
                  ("c_co2scaled_control", "(c) control, no injection", OI.green), ("a_bhp_methane_inject", "(a′) methane, capped", OI.sky),
                  ("b_bhp_co2scaled_inject", "(b′) CO₂-scaled, capped", OI.orange)]
load4d(n; dir = "04d") = (f = joinpath(OUT, dir, "$n.jls"); isfile(f) ? deserialize(f) : nothing)
rowv(rows, f) = [f(r) for r in rows]

function label_ends!(ax, xs_end, ys_end, labels, colors; gap, dx = 0.6)
    for (y, l, c) in zip(spread_labels(ys_end, gap), labels, colors)
        text!(ax, xs_end + dx, y, text = l, color = c, fontsize = 7, align = (:left, :center))
    end
end

function fig_04d_timeseries()
    dep = load4d("depletion").summary.rows
    B = [(n, l, c, load4d(n).summary.rows) for (n, l, c) in BRANCHES]
    yd = rowv(dep, r -> r.year)
    fig = pub_figure(height = 330)
    ax1 = Axis(fig[1, 1], title = "Pressure at the trap crest", xlabel = "Time (years)", ylabel = "Pressure (bar)")
    ax2 = Axis(fig[1, 2], title = "Entry drive at the fault", xlabel = "Time (years)", ylabel = "p_gas − p_brine (bar)")
    ax3 = Axis(fig[2, 1], title = "CO₂ beyond the seal", xlabel = "Time (years)", ylabel = "CO₂ (t)")
    ax4 = Axis(fig[2, 2], title = "CH₄ beyond the seal", xlabel = "Time (years)", ylabel = "CH₄ (t)")
    for ax in (ax1, ax2, ax3, ax4); vspan!(ax, 11, 21, color = (:grey, 0.12)); end
    lines!(ax1, yd, rowv(dep, r -> r.p_crest), color = :black); lines!(ax1, yd, rowv(dep, r -> r.p_hw_matched), color = :black, linestyle = :dash)
    lines!(ax2, yd, rowv(dep, r -> r.drive), color = :black)
    for (n, l, c, rows) in B
        y = rowv(rows, r -> r.year)
        if n in ("a_methane_inject", "c_co2scaled_control", "a_bhp_methane_inject")    # (b), (b′) overlap (a), (a′) in pressure
            lines!(ax1, y, rowv(rows, r -> r.p_crest), color = c); lines!(ax1, y, rowv(rows, r -> r.p_hw_matched), color = c, linestyle = :dash)
        end
        lines!(ax2, y, rowv(rows, r -> r.drive), color = c)
        lines!(ax3, y, rowv(rows, r -> r.co2.beyond_seal), color = c)
        lines!(ax4, y, rowv(rows, r -> r.ch4.beyond_seal), color = c)
    end
    hlines!(ax2, [2.5], color = OI.blue, linestyle = :dot, linewidth = 0.9); hlines!(ax2, [1.5], color = OI.vermillion, linestyle = :dot, linewidth = 0.9)
    text!(ax2, 0.5, 2.5, text = "Pe 2.5", fontsize = 7, color = OI.blue, align = (:left, :bottom))
    text!(ax2, 0.5, 1.5, text = "Pe 1.5", fontsize = 7, color = OI.vermillion, align = (:left, :top), offset = (0, -1))
    hlines!(ax1, [153.1], color = :grey50, linewidth = 0.6)
    text!(ax1, 1, 153.1, text = "original crest pressure", fontsize = 7, color = :grey30, align = (:left, :bottom))
    for (ax, f, gap) in ((ax3, r -> r.co2.beyond_seal, 190), (ax4, r -> r.ch4.beyond_seal, 20))
        ends = [f(rows[end]) for (_, _, _, rows) in B]
        label_ends!(ax, 41, ends, [split(l)[1] for (_, l, _, _) in B], [c for (_, _, c, _) in B]; gap)
        xlims!(ax, 0, 46)
    end
    for (i, ax) in enumerate((ax1, ax2, ax3, ax4)); panel_label!(fig[fldmod1(i, 2)...], string("abcd"[i])); end
    Legend(fig[3, 1:2], [[LineElement(color = :black)]; [LineElement(color = c) for (_, _, c) in BRANCHES]],
        [["depletion (years 0–11)"]; [l for (_, l, _) in BRANCHES]], orientation = :horizontal, nbanks = 2, tellwidth = false)
    Label(fig[4, 1:2], "a: solid, trap crest; dashed, the cell opposite it across the fault ((b), (b′) overlap (a), (a′) and are not drawn).
b: dotted, fault entry pressure (blue: methane curves; vermillion: CO₂-scaled).  Shaded: CO₂ injection.  Dissolution is not modelled.",
        fontsize = 7, color = :grey30, tellwidth = false)
    rowgap!(fig.layout, 4)
    save_fig(joinpath(FIG, "04d_timeseries.png"), fig)
end

function fig_04d_maps()
    sec = build_section()
    panels = []
    for (n, l, _) in BRANCHES, y in (21, 41)
        b = load4d(n)
        push!(panels, (sec, b.maps[y].Sg, "$l, yr $y"))
    end
    field_panels(panels, joinpath(FIG, "04d_gas_maps.png"); xl = (450.0, 950.0), dl = (1468.0, 1545.0), ve = 2.5,
        ncols = 2, rowh = 100, label = "Gas saturation", letters = 'a':'j', width = 500,
        caption = "Left: end of injection (year 21). Right: year 41. Vertical exaggeration ×2.5. Background: facies (faded); fault zone outlined.")
end

function fig_04d_resolution()
    fig = pub_figure(height = 225)
    ax1 = Axis(fig[1, 1], title = "Entry drive at the fault", xlabel = "Time (years)", ylabel = "p_gas(crest) − p_brine(fault) (bar)")
    ax2 = Axis(fig[1, 2], title = "CO₂ beyond the seal", xlabel = "Time (years)", ylabel = "CO₂ (t)")
    for ax in (ax1, ax2); vspan!(ax, 11, 21, color = (:grey, 0.12)); end
    styles = [(1, "04d", :solid), (5, "04d_fault5", :dash), (10, "04d_fault10", :dot)]
    for (n, l, c) in BRANCHES[1:2], (ncol, dir, ls) in styles
        b = load4d(n; dir); isnothing(b) && continue
        rows = b.summary.rows; y = rowv(rows, r -> r.year)
        lines!(ax1, y, rowv(rows, r -> r.drive), color = c, linestyle = ls)
        lines!(ax2, y, rowv(rows, r -> r.co2.beyond_seal), color = c, linestyle = ls)
    end
    hlines!(ax1, [2.5], color = OI.blue, linewidth = 0.6); hlines!(ax1, [1.5], color = OI.vermillion, linewidth = 0.6)
    text!(ax1, 41, 2.5, text = "Pe 2.5", fontsize = 7, color = OI.blue, align = (:right, :bottom))
    text!(ax1, 41, 1.5, text = "Pe 1.5", fontsize = 7, color = OI.vermillion, align = (:right, :bottom))
    panel_label!(fig[1, 1], "a"); panel_label!(fig[1, 2], "b")
    Legend(fig[2, 1:2], [[LineElement(color = c) for (_, _, c) in BRANCHES[1:2]]; [LineElement(color = :grey30, linestyle = ls) for (_, _, ls) in styles]],
        [[l for (_, l, _) in BRANCHES[1:2]]; ["$n column$(n > 1 ? "s" : "") across the fault" for (n, _, _) in styles]],
        orientation = :horizontal, nbanks = 2, tellwidth = false)
    save_fig(joinpath(FIG, "04d_fault_resolution.png"), fig)
end

# --------------------------------------------------------------------------- driver
const FIGURES = [
    "04a_co2_brine" => fig_04a_co2_brine, "04a_ch4_water" => fig_04a_ch4_water,
    "04b_depletion" => fig_04b_depletion, "04b2" => fig_04b2,
    "04d_timeseries" => fig_04d_timeseries, "04d_maps" => fig_04d_maps, "04d_resolution" => fig_04d_resolution,
    "01_section" => fig_01_section,
    "02_saturation" => fig_02_saturation, "02_accounting" => fig_02_accounting,
    "02_fault_pressure" => fig_02_fault_pressure, "02_grid" => fig_02_grid,
    "03_leakage" => fig_03_leakage, "03_time_series" => fig_03_time_series, "03_saturation" => fig_03_saturation,
    "03b_grid" => fig_03b_grid, "03b_rate" => fig_03b_rate,
]
sel = isempty(ARGS) ? FIGURES : filter(p -> any(a -> startswith(p.first, a), ARGS), FIGURES)
for (name, f) in sel
    f()
end
