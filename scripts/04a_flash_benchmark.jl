# Step 4A: standalone phase-equilibrium benchmark of the Søreide–Whitson (SW)
# equation of state in MultiComponentFlash.jl (no reservoir simulation).
#
# Compares, at reservoir conditions:
#   1. CO2 solubility in NaCl brine: SW vs the CO2-brine correlations used by
#      JutulDarcy's :co2brine model (co2_source = :salo24) in Steps 2–3.
#   2. Water content of the CO2-rich phase: SW vs the same correlations.
#   3. Brine (liquid) density: SW (Peng-Robinson based, no volume shift) vs the
#      brine density correlation.
#   4. CH4 solubility with the SW hydrocarbon–water interaction coefficient as
#      implemented (third term ∝ Tr) and as published (third term ∝ Tr²).
#
# Run from the project folder:  julia --project=. scripts/04a_flash_benchmark.jl
# Writes results/04a_*.csv and figures/04a_*.png.

using MultiComponentFlash, JutulDarcy, Printf
import JutulDarcy.CO2Properties: co2_brine_property_tables

const R = 8.314462618
const M_W = 18.01528e-3
bar = 1e5

molality_nacl(S_mgL) = (S_mgL / 1e6 / 58.44277e-3) / (1 - S_mgL / 1e6)   # mol NaCl per kg water
nacl_mole_fraction(S_mgL) = (w = S_mgL / 1e6; (w / 58.44277e-3) / (w / 58.44277e-3 + (1 - w) / M_W))

# Vapour-phase (non-aqueous) BIPs. SW (1992) values for water–gas pairs, a
# typical Peng-Robinson value for CO2–CH4. Order: Water, CarbonDioxide, Methane.
const BIP_VAPOUR = [0.0 0.1896 0.4850; 0.1896 0.0 0.10; 0.4850 0.10 0.0]

function sw_eos(names, molality)
    idx = Dict("Water" => 1, "CarbonDioxide" => 2, "Methane" => 3)
    ix = [idx[n] for n in names]
    mix = MultiComponentMixture(names, A_ij = BIP_VAPOUR[ix, ix])
    return GenericCubicEOS(mix, SoreideWhitson(mix, molality = molality))
end

"Flash at (p, T, z); returns liquid and vapour mole fractions, or nothing if single phase."
function flash_xy(eos, p, T, z)
    K = initial_guess_K(eos, (p = p, T = T, z = z))
    V = flash_2ph(eos, (p = p, T = T, z = z), K)
    isnan(V) && return nothing
    x = liquid_mole_fraction.(z, K, V)
    y = vapor_mole_fraction.(z, K, V)
    return (; x, y, V)
end

"""
Mass density (kg/m³) of a phase with mole fractions x from the cubic EOS.
SW uses phase-specific coefficients, so both sets are built explicitly: the
default argument of `mixture_compressibility_factor` builds a single set and
fails for SW (MultiComponentFlash 1.1.19).
"""
function eos_density(eos, p, T, x, phase)
    cond = (p = p, T = T, z = x, phase = phase)
    fc(ph) = MultiComponentFlash.force_coefficients(eos, (p = p, T = T, z = x, phase = ph))
    forces = (liquid = fc(:liquid), vapor = fc(:vapor))
    Z = MultiComponentFlash.mixture_compressibility_factor(eos, cond, forces)
    M = sum(x[i] * eos.mixture.properties[i].mw for i in eachindex(x))
    return p * M / (Z * R * T)
end

# Reference: JutulDarcy's CO2-brine tables (Spycher et al. 2003 / Duan & Sun
# activity for salt, as implemented for Saló-Salgado et al. 2024). For two
# components, the K-values give the phase compositions directly.
function table_binary(tab, p, T)
    Kw, Kc = tab[:K].K(p, T)
    xc = (1 - Kw) / (Kc - Kw)
    rho_b, _ = tab[:density](p, T)
    return (; xc, yw = Kw * (1 - xc), rho_brine = rho_b)
end

mkpath(joinpath(@__DIR__, "..", "results")); mkpath(joinpath(@__DIR__, "..", "figures"))
T = 368.15                                  # 95 °C, as in Steps 2–3
pressures = collect(50.0:25.0:300.0)        # bar
salinities = [0.0, 30_000.0, 55_000.0]      # mg/L NaCl-equivalent

# 1–3. CO2–brine binary.
rows = ["salinity_mgL,molality,p_bar,mCO2_SW,mCO2_table,rel_diff_pct,yH2O_SW,yH2O_table,rho_liq_SW,rho_brine_table,rho_rel_diff_pct"]
co2 = Dict{Float64, Any}()
for S in salinities
    m = molality_nacl(S)
    eos = sw_eos(["Water", "CarbonDioxide"], m)
    tab = S == 0 ? co2_brine_property_tables() :
          co2_brine_property_tables(salt_names = ["NaCl"], salt_mole_fractions = [nacl_mole_fraction(S)])
    res = []
    for pb in pressures
        p = pb * bar
        f = flash_xy(eos, p, T, [0.9, 0.1])
        r = table_binary(tab, p, T)
        mSW = f.x[2] / (f.x[1] * M_W)
        mTab = r.xc / ((1 - r.xc) * M_W)
        rhoSW = eos_density(eos, p, T, f.x, :liquid)
        push!(res, (; pb, mSW, mTab, ywSW = f.y[1], ywTab = r.yw, rhoSW, rhoTab = r.rho_brine))
        push!(rows, join([S, round(m, digits = 4), pb, round(mSW, digits = 4), round(mTab, digits = 4),
            round(100 * (mSW - mTab) / mTab, digits = 2), round(f.y[1], sigdigits = 4), round(r.yw, sigdigits = 4),
            round(rhoSW, digits = 1), round(r.rho_brine, digits = 1), round(100 * (rhoSW - r.rho_brine) / r.rho_brine, digits = 2)], ","))
    end
    co2[S] = res
end
write(joinpath(@__DIR__, "..", "results", "04a_co2_brine.csv"), join(rows, "\n") * "\n")
println("CO2–brine at 95 °C (SW vs JutulDarcy CO2-brine tables):")
println(" salinity   p(bar)  mCO2 SW  mCO2 tab  diff%   yH2O SW   yH2O tab   rho SW  rho tab")
for S in salinities, r in co2[S]
    r.pb in (50.0, 100.0, 150.0, 200.0, 300.0) || continue
    @printf(" %7.0f  %6.0f   %7.3f  %7.3f  %6.1f   %.5f   %.5f   %6.1f  %6.1f\n",
        S, r.pb, r.mSW, r.mTab, 100 * (r.mSW - r.mTab) / r.mTab, r.ywSW, r.ywTab, r.rhoSW, r.rhoTab)
end

# 4. CH4–water: as implemented, then with the published Tr² term.
function ch4_curve(m)
    eos = sw_eos(["Water", "Methane"], m)
    out = Float64[]
    for pb in pressures
        f = flash_xy(eos, pb * bar, T, [0.9, 0.1])
        push!(out, isnothing(f) ? NaN : f.x[2] / (f.x[1] * M_W))
    end
    return out
end
sw_type(m) = SoreideWhitson(["Water", "Methane"], molality = m)
k_aq_ch4(m) = MultiComponentFlash.soreide_whitson_hc_aqueous_bic(sw_type(m),
    MolecularProperty("Methane").ω, T / MolecularProperty("Methane").T_c)

ch4_impl = Dict(S => ch4_curve(molality_nacl(S)) for S in salinities)
k_impl = Dict(S => k_aq_ch4(molality_nacl(S)) for S in salinities)

# Published form: k = A0(1+α0 c) + A1 Tr (1+α1 c) + A2 Tr² (1+α2 c).
@eval MultiComponentFlash function soreide_whitson_hc_aqueous_bic(sw::SoreideWhitson, acf_i, T_ri)
    c_sw = sw.molality
    a_0, a_1, a_2 = sw.A
    a_mw_0, a_mw_1, a_mw_2 = sw.A_mw
    α_0, α_1, α_2 = sw.alphas
    A_0 = a_0 + a_mw_0*sign(acf_i)*abs(acf_i)^(-0.1)
    A_1 = a_1 + a_mw_1*acf_i
    A_2 = a_2 + a_mw_2*acf_i
    return A_0*(1 + α_0*c_sw) + A_1*T_ri*(1 + α_1*c_sw) + A_2*T_ri^2*(1 + α_2*c_sw)
end
ch4_pub = Dict(S => Base.invokelatest(ch4_curve, molality_nacl(S)) for S in salinities)
k_pub = Dict(S => Base.invokelatest(k_aq_ch4, molality_nacl(S)) for S in salinities)

rows = ["salinity_mgL,p_bar,mCH4_implemented_Tr,mCH4_published_Tr2,k_aq_implemented,k_aq_published"]
println("\nCH4–water at 95 °C: aqueous BIP and solubility (mol/kg water)")
for S in salinities
    @printf(" salinity %6.0f mg/L: k_aq implemented (Tr) = %+.4f, published (Tr²) = %+.4f\n", S, k_impl[S], k_pub[S])
    for (i, pb) in enumerate(pressures)
        push!(rows, join([S, pb, round(ch4_impl[S][i], sigdigits = 4), round(ch4_pub[S][i], sigdigits = 4),
            round(k_impl[S], digits = 4), round(k_pub[S], digits = 4)], ","))
        pb in (100.0, 200.0, 300.0) && @printf("   p = %3.0f bar: implemented %.4f, published %.4f\n", pb, ch4_impl[S][i], ch4_pub[S][i])
    end
end
write(joinpath(@__DIR__, "..", "results", "04a_ch4_water.csv"), join(rows, "\n") * "\n")

println("
Wrote results/04a_co2_brine.csv, results/04a_ch4_water.csv (figures: julia --project=. scripts/make_figures.jl 04a)")
