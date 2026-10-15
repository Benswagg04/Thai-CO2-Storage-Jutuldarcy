# Depleted gas field fluid model (Step 4B; no dissolution, see notes/04b_depleted_gas_1d.md).
#
# Three phases: brine as a separate, immiscible aqueous phase, plus a CO2–CH4
# mixture described by the Peng–Robinson equation of state (liquid and vapour
# hydrocarbon phases; at 95 °C both CO2 and CH4 are supercritical, so the
# mixture stays in one gas-like phase).
#
# What this formulation includes: CO2–CH4 mixing and density/viscosity of the
# mixture. What it excludes: dissolution of CO2 and CH4 in brine, and
# evaporation of water into the gas. (The Søreide–Whitson route, Step 4A, is
# what would add these.) See notes/04b_depleted_gas_1d.md.

using Jutul, JutulDarcy, MultiComponentFlash

const bar, day, yr = si_units(:bar, :day, :year)

gas_spec(; T_C = 95.0, p_init_bar = 150.0, swc = 0.2, kij_co2_c1 = 0.10,
    brine_density = 990.0, brine_viscosity = 3.5e-4, brine_compressibility = 4.5e-10) =
    (; T_C, p_init_bar, swc, kij_co2_c1, brine_density, brine_viscosity, brine_compressibility)

"Peng–Robinson EOS for CO2 and CH4 (properties from MultiComponentFlash's CoolProp table)."
function gas_eos(gs)
    k = gs.kij_co2_c1
    mix = MultiComponentMixture([MolecularProperty("CarbonDioxide"), MolecularProperty("Methane")],
        names = ["CO2", "C1"], A_ij = [0.0 k; k 0.0])
    return GenericCubicEOS(mix, PengRobinson())
end

"Reference density used to convert mass rates to the surface-volume rates well controls expect."
const RHO_REF_GAS = 1.0

function gas_system(gs)
    phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
    return MultiPhaseCompositionalSystemLV(gas_eos(gs), phases,
        reference_densities = (gs.brine_density, 700.0, RHO_REF_GAS))
end

"""
Replace JutulDarcy's default water properties (reference density 1 kg/m³,
viscosity 1 mPa·s) with brine values: density `brine_density` at the initial
pressure, compressibility and viscosity as given.
"""
function set_brine_properties!(model, gs)
    wpvt = JutulDarcy.ConstMuBTable(gs.p_init_bar * bar, 1.0, gs.brine_compressibility, gs.brine_viscosity, 0.0)
    for (k, m) in pairs(model.models)
        if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
            set_secondary_variables!(m,
                PhaseViscosities = JutulDarcy.ThreePhaseLBCViscositiesLV(wpvt),
                PhaseMassDensities = JutulDarcy.ThreePhaseCompositionalDensitiesLV(wpvt))
        end
    end
    return model
end

# Well rates. JutulDarcy 0.3.10 cannot use TotalMassRateTarget for producers
# (its limit checker rejects the automatic "mrat" limit), so rates are given as
# surface volumes, which are interpreted differently for the two well types:
# - injector: surface volume × the `density` given to InjectorControl = mass;
# - producer (compositional): surface volume at standard conditions (1 atm,
#   15 °C) from a flash of the produced fluid, i.e. m³ of gas at standard
#   conditions. Found in Step 4B: a "2,000" target produced 1,361 t/yr of CH4,
#   a ratio of 1.47 m³/kg, the standard volume of 1 kg of methane.
const P_STD, T_STD = 101325.0, 288.15

"Injector control for pure CO2 at a mass rate (mixture order: CO2, C1, Water)."
co2_injector_control(rate_t_per_yr, T) =
    InjectorControl(TotalRateTarget(rate_t_per_yr * 1000 / yr / RHO_REF_GAS), [1.0, 0.0, 0.0],
        density = RHO_REF_GAS, temperature = T)

"""
Producer control giving `rate_t_per_yr` of gas with composition z (default pure
CH4), converted to standard-condition volume with the Peng–Robinson density.
"""
function gas_producer_control(rate_t_per_yr, eos; z = [0.0, 1.0])
    rho_std = gas_density(eos, P_STD, T_STD, z)
    return ProducerControl(TotalRateTarget(-rate_t_per_yr * 1000 / yr / rho_std))
end

"Gas compressibility factor Z of a CO2–CH4 mixture (vapour root of Peng–Robinson)."
function gas_Z(eos, p, T, z)
    cond = (p = p, T = T, z = z, phase = :vapor)
    return MultiComponentFlash.mixture_compressibility_factor(eos, cond)
end

"Gas mass density (kg/m³) of a CO2–CH4 mixture from Peng–Robinson."
function gas_density(eos, p, T, z)
    M = sum(z[i] * eos.mixture.properties[i].mw for i in eachindex(z))
    return p * M / (gas_Z(eos, p, T, z) * 8.314462618 * T)
end
