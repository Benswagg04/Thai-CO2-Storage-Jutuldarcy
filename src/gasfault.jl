# Step 4D: methane-filled footwall trap against the Step 3 fault, for depletion
# and CO2 re-injection. Fluid model: immiscible brine + CO2–CH4 Peng–Robinson
# (src/gas.jl); dissolution of gas in brine is NOT modelled.
#
# Geometry: the exact Step 1–3 cross-section (src/section.jl), 8,040 cells,
# outer boundaries closed (the padding columns act as a finite aquifer).
# Capillary pressure (verified in notes/04c_pc_check.md): one Brooks–Corey
# gas–water curve per facies, entered as Pc_ow ≡ 0 and Pc_og(Sg) = Pc_gw(1 − Sg),
# because the reference pressure is the (absent) hydrocarbon-liquid phase.
# See notes/04d_fault_gas.md.

include(joinpath(@__DIR__, "section.jl"))
include(joinpath(@__DIR__, "gas.jl"))
using Printf

const T_GAS = convert_to_si(95.0, :Celsius)
const LAMBDA_BC = 2.0

fault_gas_spec(; z_fwl = 1534.0, p_datum_bar = 151.0, datum_depth = 1500.0,
    pe_sand = 0.05, pe_fault_ch4 = 2.5, co2_ratio = 0.6,
    well_x = 665.0, perf_depths = (44.0, 50.0)) =
    (; z_fwl, p_datum_bar, datum_depth, pe_sand, pe_fault_ch4, co2_ratio, well_x, perf_depths)

"""
Entry pressures (bar) per facies: Leverett-scaled from the storage sand as in
Steps 2–3 for facies 1–4; the fault value is set directly. `scale` multiplies
all of them (1 for methane–brine; `co2_ratio` for the CO2-scaled case, since
interfacial tension and wettability are fluid properties).
"""
function fault_gas_entry_pressures(sec, fs; scale = 1.0)
    ref = FACIES[3]
    J(k, phi) = fs.pe_sand * sqrt((ref.k / ref.phi) / (k / phi))
    pe = [J(F.k, F.phi) for F in FACIES]
    push!(pe, fs.pe_fault_ch4)
    return scale .* pe
end

pc_gw(sw, pe_bar; swr) = pe_bar * bar * max((sw - swr) / (1 - swr), 1e-3)^(-1 / LAMBDA_BC)

"Brooks–Corey inverse: Sw at capillary pressure pc (Pa); 1 below the entry pressure."
sw_from_pc(pc, pe_bar; swr) = pc <= pe_bar * bar ? 1.0 : swr + (1 - swr) * (pe_bar * bar / pc)^LAMBDA_BC

function fault_gas_pc(sec, gs, pe)
    sg = vcat(collect(range(0.0, 0.7, 71)), collect(range(0.7, 1.0, 151))[2:end])
    zero_tab = Jutul.get_1d_interpolator([0.0, 1.0], [0.0, 0.0])
    pcow = Tuple(zero_tab for _ in pe)
    pcog = Tuple(Jutul.get_1d_interpolator(sg, pc_gw.(1 .- sg, p; swr = gs.swc)) for p in pe)
    return JutulDarcy.SimpleCapillaryPressure((pcow, pcog), regions = sec.facies)
end

"Model with the Step 3 geometry and one well near the crest of the footwall trap."
function fault_gas_model(sec, gs, fs; pe_scale = 1.0)
    domain = reservoir_domain(sec.mesh, permeability = sec.perm, porosity = sec.poro, temperature = T_GAS)
    nc = number_of_cells(sec.mesh)
    wcol = column_at(sec, fs.well_x)
    perf = [cell_ijk(sec.mesh, c)[3] for c in 1:nc
            if cell_ijk(sec.mesh, c)[1] == wcol && fs.perf_depths[1] < sec.zr[c] < fs.perf_depths[2]]
    well = setup_vertical_well(domain, wcol, 1, heel = minimum(perf), toe = maximum(perf), name = :Well, simple_well = true)
    sys = gas_system(gs)
    model = setup_reservoir_model(domain, sys, wells = [well])
    set_brine_properties!(model, gs)
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0, 2.0], [gs.swc, 0.0, 0.0]))
    pe = fault_gas_entry_pressures(sec, fs; scale = pe_scale)
    set_secondary_variables!(reservoir_model(model), CapillaryPressure = fault_gas_pc(sec, gs, pe))
    push!(reservoir_model(model).output_variables, :TotalMasses, :Saturations, :PhaseMassDensities)
    parameters = setup_parameters(model)
    return (; domain, model, parameters, pe, nc)
end

"""
Footwall-trap cells: storage sand on the footwall side, above the free-water level.
"""
trap_cells(sec, domain, fs) = [c for c in eachindex(sec.facies)
    if sec.facies[c] == 3 && sec.xc[c] < sec.spec.fault_x[1] && domain[:cell_centroids][3, c] < fs.z_fwl]

"""
Capillary–gravity-consistent initial state.
- Brine pressure hydrostatic everywhere: dp/dz = ρw(p) g, integrated from the
  datum with the compressible brine density used in the model.
- In the footwall trap, gas (pure CH4) with Pc = (ρw − ρg) g (z_FWL − z) and Sw
  from that cell's Brooks–Corey curve; ρg from Peng–Robinson at the local pressure.
- All other cells water-filled (Sw = 1), including the hanging-wall sand above the
  free-water level opposite the fault: the fault is what keeps gas out of it.
The reference (liquid) pressure equals the brine pressure because Pc_ow ≡ 0.
"""
function fault_gas_initial_state(sec, m, gs, fs)
    z = m.domain[:cell_centroids][3, :]
    rho_w(p) = gs.brine_density * (1 + gs.brine_compressibility * (p - gs.p_init_bar * bar))
    zs = range(minimum(z) - 1, maximum(z) + 1, length = 4001)
    i0 = argmin(abs.(zs .- fs.datum_depth))
    pz = zeros(length(zs)); pz[i0] = fs.p_datum_bar * bar
    for i in i0+1:length(zs); pz[i] = pz[i-1] + rho_w(pz[i-1]) * 9.81 * (zs[i] - zs[i-1]); end
    for i in i0-1:-1:1;       pz[i] = pz[i+1] - rho_w(pz[i+1]) * 9.81 * (zs[i+1] - zs[i]); end
    interp(zq) = (j = clamp(searchsortedlast(zs, zq), 1, length(zs) - 1); pz[j] + (pz[j+1] - pz[j]) * (zq - zs[j]) / (zs[j+1] - zs[j]))
    p = interp.(z)
    sw = ones(m.nc)
    eos = gas_eos(gs)
    for c in trap_cells(sec, m.domain, fs)
        rg = gas_density(eos, p[c], T_GAS, [0.0, 1.0])
        pc = (rho_w(p[c]) - rg) * 9.81 * (fs.z_fwl - z[c])
        sw[c] = sw_from_pc(pc, m.pe[sec.facies[c]]; swr = gs.swc)
    end
    return setup_reservoir_state(m.model, Pressure = p, OverallMoleFractions = [0.0, 1.0], ImmiscibleSaturation = sw)
end

"Cells whose gas has passed the seal: overburden sand, plus fault cells above the footwall caprock top."
fault_gas_beyond_seal(sec) = [c for c in eachindex(sec.facies)
    if sec.facies[c] == 1 || (sec.facies[c] == FAULT_FACIES && sec.zr[c] < sec.spec.layer_bottoms[1])]

# tol_mb: JutulDarcy's mass-balance check is relative to the mass of the whole
# model (~1.27e10 kg, mostly brine), so the default 1e-7 allows ~1.3 t of error per
# time step in any component, including the ~1,000 t of CH4. 1e-10 brings the
# no-injection control from a 1.95 t CH4 loss to 7 kg at no extra cost
# (scripts/04d_conservation.jl), but stalls CO2 injection once gas enters the
# fault (median step 0.14 day; 9.6 of 30 yr in 40 min). 1e-9 completes branch (a)
# in 39 min vs 12 min at 1e-7 (scripts/04d_tolerance.jl, output/04d/tolerance*.log).
const SOLVER_KW = (relaxation = true, timestep_max_increase = 2.0, max_timestep = 30day, tol_mb = 1e-9)
