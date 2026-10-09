# Simulation case on top of the faulted section: fluid, initial state, boundaries,
# injector and schedule. See notes/02_baseline.md and notes/02c_verification.md
# for the reasoning.

include(joinpath(@__DIR__, "section.jl"))
import JutulDarcy: KilloughHysteresis, ReservoirRelativePermeabilities, add_relperm_parameters!, brooks_corey_relperm

const bar, day, yr = si_units(:bar, :day, :year)

"""
    nacl_mole_fraction(salinity_mg_per_L)

Convert NaCl-equivalent salinity to the salt mole fraction JutulDarcy expects
(moles NaCl / moles of NaCl + water). Treats mg/L ≈ mg/kg, which is within a
few percent at these salinities.
"""
function nacl_mole_fraction(S_mgL)
    w = S_mgL / 1e6
    n_salt = w / 58.44277e-3
    n_h2o = (1 - w) / 18.01528e-3
    return n_salt / (n_salt + n_h2o)
end

"""
Cell groups used to track where the CO2 ends up. Together they cover every
cell exactly once, so their sum is the CO2 in the model.
"""
function region_cells(sec)
    (; facies, xc, spec) = sec
    fw(c) = xc[c] < spec.fault_x[1]
    cells = 1:length(facies)
    return (
        fw_storage  = [c for c in cells if facies[c] == 3 && fw(c)],
        hw_storage  = [c for c in cells if facies[c] == 3 && !fw(c)],
        fault       = [c for c in cells if facies[c] == FAULT_FACIES],
        caprock     = [c for c in cells if facies[c] == 2],
        overburden  = [c for c in cells if facies[c] == 1],
        underburden = [c for c in cells if facies[c] == 4],
    )
end

case_spec(; temperature_C = 95.0, salinity_mgL = 30_000.0,
    p_datum_bar = 151.0, datum_depth = 1500.0,
    well_x = 345.0, perf_depths = (64.0, 84.0),
    rate_t_per_yr = 1500.0, inject_years = 25, post_years = 75,
    swr = 0.2, sgc = 0.05, sgr_max = 0.25,
    capillary = false, pe_ref_bar = 0.05, bc_lambda = 2.0, pc_se_min = 1e-3,
    pe_fault_bar = nothing, open_top = true,
) = (; temperature_C, salinity_mgL, p_datum_bar, datum_depth, well_x, perf_depths,
       rate_t_per_yr, inject_years, post_years, swr, sgc, sgr_max,
       capillary, pe_ref_bar, bc_lambda, pc_se_min, pe_fault_bar, open_top)

"""
    entry_pressures(sec, cs)

Capillary entry pressure (bar) for each facies (1-4) and the fault zone (5),
scaled from the storage sand with Leverett J: Pe ∝ sqrt(φ/k). For the fault,
the across-fault permeability is used, since CO2 has to enter the fault from
the side, unless `pe_fault_bar` sets the fault value directly.
"""
function entry_pressures(sec, cs)
    ref = FACIES[3]
    J(k, phi) = cs.pe_ref_bar * sqrt((ref.k / ref.phi) / (k / phi))
    pe = [J(F.k, F.phi) for F in FACIES]
    push!(pe, isnothing(cs.pe_fault_bar) ? J(sec.spec.k_across, sec.spec.phi_fault) : cs.pe_fault_bar)
    return pe
end

"""
Brooks-Corey drainage capillary pressure as a function of CO2 saturation:
Pc = Pe * Se^(-1/λ), Se = (1 - Sg - Swr)/(1 - Swr). Se is floored at
`pc_se_min` to keep the table finite as brine approaches its residual
saturation; with λ = 2 and 1e-3 the cap is ≈ 32 Pe. The table is denser at high
CO2 saturation, where the curve is steep.
"""
function brooks_corey_pc_table(pe_bar, cs)
    sg = vcat(collect(range(0.0, 0.7, 71)), collect(range(0.7, 1.0, 151))[2:end])
    se = max.((1 .- sg .- cs.swr) ./ (1 - cs.swr), cs.pc_se_min)
    return Jutul.get_1d_interpolator(sg, pe_bar .* se .^ (-1 / cs.bc_lambda) .* bar)
end

function setup_case(sec, cs = case_spec())
    T = convert_to_si(cs.temperature_C, :Celsius)
    domain = reservoir_domain(sec.mesh, permeability = sec.perm, porosity = sec.poro, temperature = T)
    nc = number_of_cells(sec.mesh)

    # Vertical injector perforated over a depth interval below the top surface.
    wcol = column_at(sec, cs.well_x)
    perf = [cell_ijk(sec.mesh, c)[3] for c in 1:nc
            if cell_ijk(sec.mesh, c)[1] == wcol && cs.perf_depths[1] < sec.zr[c] < cs.perf_depths[2]]
    injector = setup_vertical_well(domain, wcol, 1,
        heel = minimum(perf), toe = maximum(perf), name = :Injector, simple_well = true)

    # Two-component (H2O, CO2) model with K-values from the CO2-brine
    # correlations (default source: Saló-Salgado et al. 2024). CO2 can dissolve
    # in brine and water can evaporate into the CO2-rich phase.
    model = setup_reservoir_model(domain, :co2brine,
        wells = injector,
        salt_names = ["NaCl"],
        salt_mole_fractions = [nacl_mole_fraction(cs.salinity_mgL)],
        co2_physics = :kvalue,
    )

    # Relative permeability: Brooks-Corey, with hysteresis on the CO2 phase so
    # that CO2 left behind when brine re-enters is trapped (residual trapping).
    sl = range(0, 1, 50)
    sg = range(0, 1, 50)
    kr_brine = PhaseRelativePermeability(sl, brooks_corey_relperm.(sl, n = 3, residual = cs.swr), label = :og)
    kr_drain = PhaseRelativePermeability(sg, brooks_corey_relperm.(sg, n = 2, residual = cs.sgc, residual_total = cs.swr + cs.sgc), label = :g)
    kr_imb = PhaseRelativePermeability(sg, brooks_corey_relperm.(sg, n = 2, residual = cs.sgr_max, residual_total = cs.swr + cs.sgr_max), label = :g)
    relperm = ReservoirRelativePermeabilities(g = (kr_drain, kr_imb), og = kr_brine, hysteresis_g = KilloughHysteresis())
    replace_variables!(model, RelativePermeabilities = relperm)
    add_relperm_parameters!(model)

    # Optional capillary pressure, one curve per facies (region = facies index).
    # JutulDarcy uses brine as the reference phase, so p_CO2 = p_brine + Pc(S_CO2).
    pc_tables = nothing
    if cs.capillary
        pc_tables = Tuple(brooks_corey_pc_table(pe, cs) for pe in entry_pressures(sec, cs))
        pc = JutulDarcy.SimpleCapillaryPressure((pc_tables,), regions = sec.facies)
        set_secondary_variables!(reservoir_model(model), CapillaryPressure = pc)
    end

    # Hydrostatic initial state from a datum pressure, all brine. In this model
    # brine is the "liquid" phase and CO2 the "vapor" phase, so the contact
    # between them (goc) is placed above the model top: everything starts as brine.
    eql = EquilibriumRegion(model, cs.p_datum_bar * bar, cs.datum_depth,
        goc = 1000.0, composition = [1.0, 0.0], temperature = T)
    state0 = setup_reservoir_state(model, eql)
    parameters = setup_parameters(model)

    # Boundaries, all held at the initial (hydrostatic) pressure, brine inflow:
    # - left and right edges: the formations continue laterally (x-direction);
    # - top (optional): connection to shallower aquifers (z-direction). CO2 that
    #   leaves through it is counted as having leaked above the model.
    p0 = state0[:Reservoir][:Pressure]
    edge = [c for c in 1:nc if cell_ijk(sec.mesh, c)[1] in (1, sec.nxt)]
    bc = flow_boundary_condition(edge, domain, p0[edge], T, fractional_flow = [1.0, 0.0], dir = :x)
    top = [c for c in 1:nc if cell_ijk(sec.mesh, c)[3] == 1]
    if cs.open_top
        bc = vcat(bc, flow_boundary_condition(top, domain, p0[top], T, fractional_flow = [1.0, 0.0], dir = :z))
    end

    # Injection specified as a mass rate. TotalRateTarget is a surface volume
    # rate, so divide by the CO2 reference (surface) density.
    rho_co2_surface = JutulDarcy.reference_densities(:co2brine)[2]
    q = cs.rate_t_per_yr * 1000 / yr / rho_co2_surface
    ctrl = InjectorControl(TotalRateTarget(q), [0.0, 1.0], density = rho_co2_surface, temperature = T)
    forces_inj = setup_reservoir_forces(model, control = Dict(:Injector => ctrl), bc = bc)
    forces_shut = setup_reservoir_forces(model, bc = bc)

    dt = fill(1.0yr, cs.inject_years + cs.post_years)
    forces = vcat(fill(forces_inj, cs.inject_years), fill(forces_shut, cs.post_years))

    outvars = reservoir_model(model).output_variables
    push!(outvars, :TotalMasses)
    push!(outvars, :PhaseMobilities)   # needed to recompute boundary fluxes afterwards
    unique!(outvars)
    return (; model, state0, parameters, forces, forces_shut, dt, domain,
        regions = region_cells(sec), cs, pc_tables, edge, top, bc)
end

"""
CO2 mass (tonnes) in each region for every report step, split into free
(CO2-rich phase) and dissolved (in brine). The split uses each cell's phase
masses: S·ρ·(mass fraction of CO2) for the vapour and liquid phases.
"""
function co2_by_region(states, regions)
    total = Dict{Symbol, Vector{Float64}}()
    free = Dict{Symbol, Vector{Float64}}()
    diss = Dict{Symbol, Vector{Float64}}()
    for (name, cells) in pairs(regions)
        total[name] = zeros(length(states))
        free[name] = zeros(length(states))
        diss[name] = zeros(length(states))
        for (i, s) in enumerate(states)
            S, rho, X, Y, M = s[:Saturations], s[:PhaseMassDensities], s[:LiquidMassFractions], s[:VaporMassFractions], s[:TotalMasses]
            for c in cells
                m = M[2, c]
                a = S[2, c] * rho[2, c] * Y[2, c]   # CO2 in the CO2-rich phase
                b = S[1, c] * rho[1, c] * X[2, c]   # CO2 dissolved in brine
                f = a + b > 0 ? a / (a + b) : 0.0
                total[name][i] += m / 1000
                free[name][i] += f * m / 1000
                diss[name][i] += (1 - f) * m / 1000
            end
        end
    end
    return (; total, free, diss)
end
