# Faulted 2D cross-section: a tilted fault block with one vertical fault.
#
# Geometry choices (see notes/01_geometry.md for the reasoning):
# - 2D vertical slice (nx × 1 × nz). Fast enough to run many fault scenarios on a laptop.
# - Beds rise toward the fault (tilted fault-block trap), so buoyant CO2 migrates
#   up-dip *to* the fault instead of away from it.
# - Throw is represented by shifting the facies on the hanging-wall side. The
#   grid itself stays continuous; what changes across the fault is which rock
#   sits opposite which (juxtaposition).
# - The fault is a zone of its own cells with an anisotropic permeability:
#   k_across (x) controls cross-fault flow, k_along (z) controls up-fault flow.
# - Everything is defined in metres, so the grid can be refined without
#   changing the geology (e.g. the fault stays 10 m wide).
# - Optional padding columns of growing width on both sides move the lateral
#   boundaries kilometres away from the area of interest.

using Jutul, JutulDarcy

const mD = si_unit(:darcy) / 1000

# Facies 1-4 are the stratigraphy, 5 is the fault zone.
const FACIES = (
    (name = "Overburden sand",   symbol = 'o', k = 200.0, kv_kh = 0.1, phi = 0.22),
    (name = "Caprock shale",     symbol = '-', k = 1e-3,  kv_kh = 1.0, phi = 0.05),
    (name = "Storage sand",      symbol = 's', k = 300.0, kv_kh = 0.1, phi = 0.25),
    (name = "Underburden shale", symbol = '=', k = 1e-3,  kv_kh = 1.0, phi = 0.05),
)
const FAULT_FACIES = 5

"""
    section_spec(; kwarg...)

All geometric inputs in one place, in metres. `layer_bottoms` gives the depth
below the top surface of the base of each footwall layer (top to bottom).
`fault_x` is the x-interval occupied by the fault zone.
"""
section_spec(; dx = 10.0, dz = 2.0, Lx = 1200.0, Ly = 100.0, Lz = 120.0,
    top_depth = 1500.0, dip_deg = 2.0,
    fault_x = (690.0, 700.0), throw = 30.0,
    layer_bottoms = (24.0, 44.0, 84.0, 120.0),
    k_across = 0.01, k_along = 0.1, phi_fault = 0.10,
    pad_cells = 7, pad_first = 20.0,
) = (; dx, dz, Lx, Ly, Lz, top_depth, dip_deg, fault_x, throw, layer_bottoms,
       k_across, k_along, phi_fault, pad_cells, pad_first)

function footwall_facies(zr, layer_bottoms)
    for (f, zb) in enumerate(layer_bottoms)
        zr < zb && return f
    end
    return length(layer_bottoms)
end

"Facies at x (m) and depth zr below the top surface (m). The hanging wall (right side) is downthrown by `throw`."
function facies_at(x, zr, spec)
    spec.fault_x[1] < x < spec.fault_x[2] && return FAULT_FACIES
    x < spec.fault_x[1] && return footwall_facies(zr, spec.layer_bottoms)
    # Hanging wall: stratigraphy shifted down; above it is overburden sand.
    zs = zr - spec.throw
    zs <= 0 && return 1
    return footwall_facies(zs, spec.layer_bottoms)
end

function build_section(spec = section_spec())
    (; dx, dz, Lx, Ly, Lz, top_depth, dip_deg, pad_cells, pad_first) = spec
    nx = round(Int, Lx / dx)
    nz = round(Int, Lz / dz)
    # Padding widths double each column (20, 40, ..., 1280 m by default, 2.54 km
    # per side), independent of the core cell size.
    pads = [pad_first * 2.0^(k - 1) for k in 1:pad_cells]
    dxs = vcat(reverse(pads), fill(dx, nx), pads)
    x0 = -sum(pads; init = 0.0)
    nxt = length(dxs)
    mesh = reservoir_mesh((nxt, 1, nz), (dxs, [Ly], fill(dz, nz)), origin = [x0, 0.0, top_depth])

    # Tilt: z is depth (positive down), so subtracting makes beds shallower to the right.
    slope = tand(dip_deg)
    pts = mesh.node_points
    for (i, pt) in enumerate(pts)
        pts[i] = pt - [0.0, 0.0, slope * pt[1]]
    end

    x_edges = x0 .+ vcat(0.0, cumsum(dxs))
    nc = number_of_cells(mesh)
    xc = zeros(nc)    # cell-centre x (unaffected by the tilt)
    zr = zeros(nc)    # depth of cell centre below the top surface
    facies = zeros(Int, nc)
    perm = zeros(3, nc)
    poro = zeros(nc)
    for c in 1:nc
        I, _, K = cell_ijk(mesh, c)
        xc[c] = 0.5 * (x_edges[I] + x_edges[I + 1])
        zr[c] = (K - 0.5) * dz
        f = facies_at(xc[c], zr[c], spec)
        facies[c] = f
        if f == FAULT_FACIES
            perm[1, c] = spec.k_across * mD   # across the fault (x)
            perm[2, c] = spec.k_along * mD
            perm[3, c] = spec.k_along * mD    # along the fault (z)
            poro[c] = spec.phi_fault
        else
            F = FACIES[f]
            perm[1, c] = perm[2, c] = F.k * mD
            perm[3, c] = F.k * F.kv_kh * mD
            poro[c] = F.phi
        end
    end
    core = (xc .> 0.0) .& (xc .< Lx)
    return (; mesh, facies, perm, poro, spec, xc, zr, core, nxt, nz)
end

"Column index (including padding) of the cell column that contains x."
column_at(sec, x) = cell_ijk(sec.mesh, argmin(abs.(sec.xc .- x)))[1]

"Print the core facies as a character map (every `stride`-th column)."
function print_section(sec; stride = 2, well_x = nothing)
    (; mesh, facies, spec, nz) = sec
    sym(f) = f == FAULT_FACIES ? '|' : FACIES[f].symbol
    pad = spec.pad_cells
    nx = round(Int, spec.Lx / spec.dx)
    grid = fill(' ', nz, sec.nxt)
    for c in eachindex(facies)
        I, _, K = cell_ijk(mesh, c)
        grid[K, I] = sym(facies[c])
    end
    fault_cols = unique(cell_ijk(mesh, c)[1] for c in findall(facies .== FAULT_FACIES))
    cols = sort(unique(vcat(collect(pad+1:stride:pad+nx), fault_cols)))
    wcol = isnothing(well_x) ? -1 : column_at(sec, well_x)
    for K in 1:nz
        row = [(I == wcol && grid[K, I] == 's') ? 'W' : grid[K, I] for I in cols]
        println(lpad(K, 3), " ", String(row))
    end
    println("Legend: o overburden sand, - caprock, s storage sand, = underburden, | fault, W well")
end

"Thickness (m) of footwall storage sand that faces each hanging-wall facies across the fault."
function juxtaposition_table(sec)
    (; spec, nz) = sec
    out = Dict{String, Float64}()
    x_hw = spec.fault_x[2] + 1e-6
    for K in 1:nz
        zr = (K - 0.5) * spec.dz
        footwall_facies(zr, spec.layer_bottoms) == 3 || continue
        name = FACIES[facies_at(x_hw, zr, spec)].name
        out[name] = get(out, name, 0.0) + spec.dz
    end
    return out
end
