# Shared figure style and helpers (used by scripts/make_figures.jl).
#
# - One publication theme for every figure: no top/right spines, outward ticks,
#   consistent font sizes; figures sized in points for a page and saved at high
#   resolution.
# - Fixed categorical palette (Okabe–Ito, designed to remain distinguishable for
#   common colour-vision deficiencies). Each scenario keeps its colour in every
#   figure; variants of a scenario use the same hue with a different line style.
# - Perceptually uniform cmocean colormaps (Thyng et al. 2016) for fields.
# - Cross-sections drawn in 2D with each cell as its true (tilted) quadrilateral,
#   depth increasing downward, and the vertical exaggeration stated on the axis.

using CairoMakie

const PT_FULL = 470    # full text width, points (~6.5 in)
const PX_PER_PT = 4    # saved resolution: ~1,900 px across a full-width figure

const OI = (blue = colorant"#0072B2", vermillion = colorant"#D55E00", green = colorant"#009E73",
    sky = colorant"#56B4E9", orange = colorant"#E69F00", purple = colorant"#CC79A7",
    yellow = colorant"#F0E442", black = colorant"#000000", grey = colorant"#7F7F7F")

const CMAP_SG = :dense      # gas saturation: near-white (no gas) to dark
const CMAP_ZCO2 = :deep     # CO2 fraction in the gas
const CMAP_LOGK = Reverse(:turbid)   # log10 permeability: low k dark (shale), high k light (sand)

# Facies colours for geological context (light, neutral; the fault darker).
const FACIES_COLORS = [colorant"#EFE3C2", colorant"#B8B3AA", colorant"#DCC27E", colorant"#8F8A82", colorant"#5B4A3F"]

const PUB_THEME = Theme(
    fontsize = 9, figure_padding = (4, 8, 4, 4),
    Axis = (xgridvisible = false, ygridvisible = false, topspinevisible = false, rightspinevisible = false,
        spinewidth = 0.7, xtickwidth = 0.7, ytickwidth = 0.7, xticksize = 3, yticksize = 3,
        xtickalign = 0, ytickalign = 0, titlesize = 9, titlefont = :regular, titlealign = :left,
        xlabelsize = 9, ylabelsize = 9, xticklabelsize = 8, yticklabelsize = 8, xlabelpadding = 2, ylabelpadding = 2),
    Legend = (framevisible = false, labelsize = 8, patchsize = (16, 8), rowgap = 0, colgap = 8, padding = (2, 2, 2, 2)),
    Colorbar = (labelsize = 9, ticklabelsize = 8, spinewidth = 0.5, tickwidth = 0.7, ticksize = 3, width = 7),
    Lines = (linewidth = 1.4,), Scatter = (markersize = 5,),
)

"Figure with the publication theme applied (width in points)."
pub_figure(; width = PT_FULL, height) = Figure(size = (width, height))

save_fig(path, fig) = (save(path, fig; px_per_unit = PX_PER_PT); println("wrote ", path))

"Bold panel letter at the top-left corner of a layout cell."
panel_label!(pos, s) = Label(pos[1, 1, TopLeft()], s, fontsize = 10, font = :bold, padding = (0, 6, 2, 0), halign = :right)

"Shade the CO2-injection period."
inj_span!(ax, t0, t1) = vspan!(ax, t0, t1, color = (:grey, 0.12))

# ---------------------------------------------------------------- cross-sections

"""
Quadrilateral of every cell of a tilted Cartesian section (one cell in y):
x edges `xe`, layer thickness `dz`, top depth `top`, slope (depth decreases by
`slope` per metre of x). Cell order matches Jutul (x fastest, then depth).
"""
function tilted_cell_polygons(xe, dz, nz, top, slope)
    d(x, zrel) = top + zrel - slope * x
    nx = length(xe) - 1
    polys = Vector{Vector{Point2f}}(undef, nx * nz)
    for k in 1:nz, i in 1:nx
        zt, zb = (k - 1) * dz, k * dz
        polys[i + nx * (k - 1)] = Point2f[(xe[i], d(xe[i], zt)), (xe[i+1], d(xe[i+1], zt)),
            (xe[i+1], d(xe[i+1], zb)), (xe[i], d(xe[i], zb))]
    end
    return polys
end
section_polygons(sec) = tilted_cell_polygons(sec.x_edges, sec.spec.dz, sec.nz, sec.spec.top_depth, sec.slope)

"""
Axis for a cross-section with depth increasing downward and a fixed vertical
exaggeration `ve` over the window `xlims` × `dlims` (depths, top first).
"""
function section_axis(pos; xlims, dlims, ve = 5, title = "", xlabel = "x (m)", ylabel = "Depth (m)", ve_in_title = true, kw...)
    w, h = xlims[2] - xlims[1], dlims[2] - dlims[1]
    # The vertical exaggeration goes in the title, so the depth label stays short.
    ttl = (ve == 1 || !ve_in_title) ? title : (isempty(title) ? "vert. exag. ×$ve" : "$title  (vert. exag. ×$ve)")
    ax = Axis(pos; title = ttl, xlabel, ylabel,
        yreversed = true, aspect = AxisAspect(w / (h * ve)), limits = (xlims, dlims), kw...)
    return ax
end

"""
Vertical positions for end-of-line labels: the original `ys`, pushed apart so
neighbours are at least `gap` apart (same units as `ys`), preserving order.
"""
function spread_labels(ys, gap)
    p = sortperm(ys); y = float.(ys[p])
    for _ in 1:50, i in 2:length(y)
        d = y[i] - y[i-1]
        d < gap && (s = (gap - d) / 2; y[i-1] -= s; y[i] += s)
    end
    out = similar(y); out[p] = y
    return out
end

"Colours of `values` on `colormap` over `colorrange` (clamped)."
function map_colors(values, colormap, colorrange)
    cm = Makie.to_colormap(colormap)          # handles Reverse(...) as well as names
    n = length(cm)
    idx(v) = 1 + round(Int, (n - 1) * clamp((v - colorrange[1]) / (colorrange[2] - colorrange[1]), 0, 1))
    return [cm[idx(v)] for v in values]
end

# Cells are stroked in their own fill colour so that anti-aliasing does not leave
# faint seams between neighbouring cells.
"Draw cells (filled polygons) coloured by `values`, optionally only `cells`."
function section_field!(ax, polys, values; cells = eachindex(values), colormap, colorrange)
    cols = map_colors(values[cells], colormap, colorrange)
    poly!(ax, polys[cells], color = cols, strokecolor = cols, strokewidth = 0.25)
end

"Facies as geological context, optionally faded (`fade` 0–1 towards white)."
function section_facies!(ax, sec, polys; cells = eachindex(sec.facies), fade = 0.0)
    fadec(c) = RGBf(c.r + (1 - c.r) * fade, c.g + (1 - c.g) * fade, c.b + (1 - c.b) * fade)
    cols = [fadec(RGBf(FACIES_COLORS[f])) for f in sec.facies[cells]]
    poly!(ax, polys[cells], color = cols, strokecolor = cols, strokewidth = 0.25)
end

"Outline of the fault zone (left and right edges), as thin lines."
function fault_outline!(ax, sec, polys; color = (:black, 0.6), linewidth = 0.5)
    fc = findall(sec.facies .== FAULT_FACIES)
    xs = (minimum(p[1][1] for p in polys[fc]), maximum(p[2][1] for p in polys[fc]))
    for x in xs
        cs = [c for c in fc if any(q -> abs(q[1] - x) < 1e-6, polys[c])]
        d = sort(unique(round.([q[2] for c in cs for q in polys[c] if abs(q[1] - x) < 1e-6], digits = 3)))
        lines!(ax, [x, x], [first(d), last(d)], color = color, linewidth = linewidth)
    end
end

"Cells whose centre lies in the window (used to limit drawing to the visible area)."
window_cells(sec, xlims, dlims) = findall(c -> xlims[1] - 20 <= sec.xc[c] <= xlims[2] + 20, eachindex(sec.xc))
