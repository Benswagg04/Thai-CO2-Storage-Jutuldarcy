# Step 4A: Søreide–Whitson phase-equilibrium benchmark (no flow simulation)

Script: `scripts/04a_flash_benchmark.jl`. Results: `results/04a_co2_brine.csv`, `results/04a_ch4_water.csv`. Figures: `figures/04a_co2_brine.png`, `figures/04a_ch4_water.png`.

## Why this step
The depleted-gas-field model needs CO2–CH4 mixing *and* gas dissolution in brine in one formulation. In JutulDarcy 0.3.10 the candidate is the Søreide–Whitson (SW) equation of state in MultiComponentFlash.jl (1.1.19, now a direct project dependency, same version as before). SW is Peng–Robinson with:
- a salinity-dependent attraction term for water;
- **phase-specific** binary interaction coefficients (BIPs): one set for the aqueous phase (salinity- and temperature-dependent), one for the vapour phase.

MultiComponentFlash's own test suite contains **no tests for SW**, so this is the first check of it in this project. The reference for CO2 is the CO2-brine correlation set used by JutulDarcy's `:co2brine` model in Steps 2–3 (Spycher et al. 2003 fugacities with Duan & Sun-type salt activity, as implemented for Saló-Salgado et al. 2024). That set is a model, not data; experimental comparison is still to come.

## Set-up
- 95 °C (368.15 K); 50–300 bar; NaCl 0, 30,000 and 55,000 mg/L (0, 0.53 and 1.00 mol/kg).
- Components named "Water", "CarbonDioxide", "Methane" (SW recognises water and CO2 by name).
- Vapour-phase BIPs (supplied by us; the SW constructor does not set them): H2O–CO2 0.1896, H2O–CH4 0.4850 (SW 1992 values, as quoted in later papers), CO2–CH4 0.10 (a typical Peng–Robinson value, *assumed*).
- Solubility in mol per kg water, from the flashed liquid composition.

## Results

### 1. CO2 solubility (mol/kg water)
| p (bar) | 0 mg/L: SW / tables | 30,000: SW / tables | 55,000: SW / tables |
|---|---|---|---|
| 50 | 0.422 / 0.478 (−12%) | 0.415 / 0.424 (−2%) | 0.380 / 0.382 (−1%) |
| 100 | 0.743 / 0.807 (−8%) | 0.731 / 0.713 (+3%) | 0.665 / 0.642 (+4%) |
| 200 | 1.107 / 1.135 (−2%) | 1.088 / 0.999 (+9%) | 0.984 / 0.897 (+10%) |
| 300 | 1.306 / 1.302 (0%) | 1.283 / 1.144 (+12%) | 1.157 / 1.026 (+13%) |

- In pure water the two models agree within 0–12%, best at high pressure.
- **Salting-out is much weaker in SW.** At 200 bar, going from 0 to 30,000 mg/L lowers solubility by 1.7% in SW and by 12.0% in the tables (0 → 55,000 mg/L: 11.1% vs 21.0%).
  - The reason is in the CO2 aqueous BIP: at T_r = 1.21, its salinity terms nearly cancel (k ≈ −0.032 for both 0 and 0.53 mol/kg).
- Later work re-fitted the SW CO2–brine BIPs to improve this (Chabab et al. 2019, 2020, as cited in Ahmed et al. 2025). MultiComponentFlash implements the original 1992 form.
- Which model is closer to measurements is an empirical question that our own CO2–NaCl solubility work is suited to answer.

### 2. Water content of the CO2-rich phase
SW is 3–8% above the tables at all salinities (e.g. y_H2O = 0.0181 vs 0.0173 at 200 bar, fresh water). Good agreement.

### 3. Liquid (brine) density
SW gives 809–824 kg/m³ against 964–1,012 kg/m³ from the brine correlation: **16–19% too low**, and almost insensitive to salinity, because salt is not a component in SW. This is the known weakness of Peng–Robinson for liquid water.
- **Consequence:** SW densities must not be used for brine in flow simulation. That needs a volume shift, or (as in the `:co2brine` model) brine-density and viscosity correlations kept separate from the phase-equilibrium model.

### 4. CH4 solubility: a likely implementation error
SW's aqueous BIP for hydrocarbon–water pairs is, per Søreide & Whitson (1992) eq. 12, as reproduced in generalised form by Ahmed, …, Møyner et al. (2025, arXiv:2506.02582, eq. 18):

  k_aq = A0(1 + α0·c) + A1·T_r(1 + α1·c) + A2·**T_r²**(1 + α2·c)

MultiComponentFlash 1.1.19 (`src/eos/soreide_whitson.jl`, `soreide_whitson_hc_aqueous_bic`) computes the third term as `A_2*T_ri*(1 + α_2*c_sw)`, i.e. with **T_r instead of T_r²**. Its code comment says it implements "eq 12 from paper". The coefficient values themselves match.

| 95 °C, 200 bar | k_aq (CH4–H2O) | CH4 solubility (mol/kg) |
|---|---|---|
| Fresh water, as implemented (T_r) | +0.211 | 0.0106 |
| Fresh water, published form (T_r²) | −0.095 | 0.108 |
| 30,000 mg/L, as implemented | +0.232 | 0.0088 |
| 30,000 mg/L, published form | −0.076 | 0.091 |

- The implemented form gives methane solubility about **10× lower**.
- Sanity check: commonly cited methane solubility in fresh water near 100 °C and 200 bar is of order 0.1 mol/kg (mole fraction ≈ 0.002). That matches the published form. **This must be confirmed against a real dataset** (e.g. Duan & Mao 2006 model tables, or Culberson & McKetta 1951 measurements) and against SW 1992 eq. 12 itself, before anything is reported upstream.
- For a depleted *gas* field, methane dissolution and exsolution matter. The implemented form would effectively remove methane from the brine.

### 5. A second, minor API issue
`mixture_compressibility_factor(eos, cond)` fails for SW: its default argument builds one set of coefficients, while SW expects separate liquid and vapour sets. The flash itself is not affected. Workaround in the script: pass both sets explicitly.

## What this means for Step 4
- SW is a usable route for CO2–CH4–brine phase *equilibrium*, but not as it stands:
  1. the hydrocarbon–water BIP needs the T_r² term (patched in our script, for testing only);
  2. the CO2–brine salting-out is weak (original 1992 BIPs);
  3. liquid density is unusable without correction.
- Before any flow simulation (Stage A2): confirm points 1–2 against experimental data, and decide how to supply brine density and viscosity.
- Whether JutulDarcy's flow coupling handles SW's phase-specific coefficients correctly is still untested.

## Open items needing the researcher
- Confirm SW 1992 eq. 12 from the original paper (Fluid Phase Equilibria 77, 217–240).
- Experimental CO2–NaCl solubility at 60–120 °C (our compiled database is the natural source) and CH4–water/brine data.
- ~~Decide whether to report the T_r issue to the MultiComponentFlash maintainers.~~ Resolved: the term is already corrected upstream (MultiComponentFlash 2.0.1); no issue to file. See `notes/04a_upstream_issue_draft.md`.
