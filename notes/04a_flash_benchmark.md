# Step 4A: Søreide–Whitson phase-equilibrium benchmark (no flow simulation)

Script: `scripts/04a_flash_benchmark.jl`. Results: `results/04a_co2_brine.csv`, `results/04a_ch4_water.csv`. Figures: `figures/04a_co2_brine.png`, `figures/04a_ch4_water.png`.

## Why
The depleted-field model needs CO2–CH4 mixing and gas dissolution in brine in one formulation. In JutulDarcy 0.3.10 the candidate is the Søreide–Whitson (SW) equation of state in MultiComponentFlash.jl (1.1.19, which I added as a direct dependency at the same version). SW is Peng–Robinson with a salinity-dependent attraction term for water and with separate binary interaction coefficients (BIPs) for the aqueous phase (depending on salinity and temperature) and the vapour phase.

MultiComponentFlash's own test suite has no tests for SW, so before relying on it I compared it with something I trust. For CO2 the reference is the correlation set behind JutulDarcy's `:co2brine` model in Steps 2–3 (Spycher et al. 2003 fugacities with Duan & Sun-type salt activity, as implemented for Saló-Salgado et al. 2024). That is a model rather than data; comparison with measurements is still to come.

## Set-up
- 95 °C (368.15 K), 50–300 bar, NaCl 0, 30,000 and 55,000 mg/L (0, 0.53 and 1.00 mol/kg).
- Components named "Water", "CarbonDioxide" and "Methane" (SW recognises water and CO2 by name).
- The SW constructor does not set the vapour-phase BIPs, so I supplied them: H2O–CO2 0.1896 and H2O–CH4 0.4850 (the SW 1992 values, as quoted in later papers), and CO2–CH4 0.10 (a typical Peng–Robinson value, assumed).
- Solubility in mol per kg water, from the flashed liquid composition.

## Results

### CO2 solubility (mol/kg water)

| p (bar) | 0 mg/L: SW / tables | 30,000: SW / tables | 55,000: SW / tables |
|---|---|---|---|
| 50 | 0.422 / 0.478 (−12%) | 0.415 / 0.424 (−2%) | 0.380 / 0.382 (−1%) |
| 100 | 0.743 / 0.807 (−8%) | 0.731 / 0.713 (+3%) | 0.665 / 0.642 (+4%) |
| 200 | 1.107 / 1.135 (−2%) | 1.088 / 0.999 (+9%) | 0.984 / 0.897 (+10%) |
| 300 | 1.306 / 1.302 (0%) | 1.283 / 1.144 (+12%) | 1.157 / 1.026 (+13%) |

In pure water the two agree to within 0–12%, best at high pressure. The salting-out effect is much weaker in SW: at 200 bar, going from 0 to 30,000 mg/L lowers solubility by 1.7% in SW and by 12.0% in the tables (0 → 55,000 mg/L: 11.1% against 21.0%). The reason is in the CO2 aqueous BIP: at T_r = 1.21 its salinity terms nearly cancel (k ≈ −0.032 at both 0 and 0.53 mol/kg). Later work re-fitted the SW CO2–brine BIPs to improve this (Chabab et al. 2019, 2020, as cited in Ahmed et al. 2025); MultiComponentFlash uses the original 1992 form. Which model is closer to measurements is an empirical question, and CO2–NaCl solubility data are what I worked on before this project, so this is a natural comparison to make.

### Water content of the CO2-rich phase
SW is 3–8% above the tables at all salinities (for example y_H2O = 0.0181 against 0.0173 at 200 bar in fresh water), which is good agreement.

### Brine density
SW gives 809–824 kg/m³ against 964–1,012 kg/m³ from the brine correlation, 16–19% too low, and almost independent of salinity because salt is not a component in SW. This is the known weakness of Peng–Robinson for liquid water. SW densities therefore cannot be used for brine in a flow simulation without a volume shift, or (as in the `:co2brine` model) brine density and viscosity taken from separate correlations.

### CH4 solubility: a discrepancy in the installed version
SW's aqueous BIP for hydrocarbon–water pairs is, from Søreide & Whitson (1992) eq. 12, as written in generalised form by Ahmed, …, Møyner et al. (2025, arXiv:2506.02582, eq. 18):

  k_aq = A0(1 + α0·c) + A1·T_r(1 + α1·c) + A2·T_r²(1 + α2·c)

MultiComponentFlash 1.1.19 (`soreide_whitson_hc_aqueous_bic` in `src/eos/soreide_whitson.jl`) computes the third term as `A_2*T_ri*(1 + α_2*c_sw)`, with T_r instead of T_r², although its comment says it implements "eq 12 from paper". The coefficient values match.

| 95 °C, 200 bar | k_aq (CH4–H2O) | CH4 solubility (mol/kg) |
|---|---|---|
| Fresh water, as implemented (T_r) | +0.211 | 0.0106 |
| Fresh water, published form (T_r²) | −0.095 | 0.108 |
| 30,000 mg/L, as implemented | +0.232 | 0.0088 |
| 30,000 mg/L, published form | −0.076 | 0.091 |

The implemented form gives about 10 times less dissolved methane. Commonly cited methane solubility in fresh water near 100 °C and 200 bar is of order 0.1 mol/kg (mole fraction ≈ 0.002), which matches the published form, but this still needs checking against a real dataset (for example Duan & Mao 2006, or the Culberson & McKetta 1951 measurements) and against eq. 12 in the original paper. For a depleted gas field this matters, because the implemented form would effectively remove methane from the brine.

Update (11 October 2026): before reporting this, I checked the later releases. The term is still `T_ri` in 1.1.20 and is `T_ri^2` in 2.0.1, changed in pull request #31 ("Modernize docs + EOS fixes + Clapeyron comparisons", merged 13 September 2026). So the problem is already fixed upstream, and there is nothing to report. What this step found is a discrepancy between the older release this project uses and the published equation, and its size; not a new bug. MultiComponentFlash 2.x needs JutulDarcy 0.4 and Jutul 0.5, so using the fix here means upgrading the whole environment. The lesson for me: check the latest release and changelog before assuming a defect.

### A minor API issue
`mixture_compressibility_factor(eos, cond)` fails for SW, because its default argument builds one set of coefficients while SW expects separate liquid and vapour sets. The flash itself is not affected; the script passes both sets explicitly. I have not checked whether this is still the case in 2.0.1.

## What this means for Step 4
SW is a possible route to CO2–CH4–brine phase equilibrium, but not as it stands in this environment: the hydrocarbon–water BIP needs the T_r² term (patched in my script, for testing only), CO2 salting-out is weak with the original BIPs, and the liquid density is unusable without correction. Whether JutulDarcy's flow coupling handles SW's phase-specific coefficients correctly is also untested. For the flow models in Step 4 I therefore used a simpler formulation without dissolution (Step 4B onwards), and left the SW route for later.

## Open items
- Confirm eq. 12 in the original paper (Fluid Phase Equilibria 77, 217–240).
- Compare with experimental CO2–NaCl solubility at 60–120 °C and with CH4–water/brine data.
