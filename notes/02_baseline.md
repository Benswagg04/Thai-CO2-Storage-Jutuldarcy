# Step 2: Baseline CO2 injection (no capillary pressure)

## Inputs and why

| Input | Value | Reasoning |
|---|---|---|
| Fluid model | `:co2brine`, `co2_physics = :kvalue` | Two components (H2O, CO2), two phases. CO2 dissolves in brine and water evaporates into CO2, via K-values from the Saló-Salgado et al. (2024) correlations (JutulDarcy default `co2_source = :salo24`). |
| Temperature | 95 °C, isothermal | Low end of the Pattani Basin gradient (≈27 °C + 1.5 km × 45 °C/km). The correlations are validated to about 100 °C; JutulDarcy's tables extend to 200 °C but are built without range checks. |
| Salinity | 30,000 mg/L NaCl-equivalent (**hypothetical baseline**) | Hypothetical moderate salinity, chosen for illustration; not a documented basin value. Converted to an NaCl mole fraction of about 0.0094. |
| Initial pressure | Hydrostatic from 151 bar at 1,500 m | Normal hydrostatic gradient (~0.1 bar/m). Equilibrated with `EquilibriumRegion`; checked by running 5 years with no wells. |
| Boundaries | Left and right edges held at initial pressure | Formations continue beyond the slice; prevents unrealistic pressure build-up in a closed 2D box. |
| Injector | Column 35 (x ≈ 345 m), bottom 20 m of storage sand | Injecting low lets CO2 contact more brine as it rises (more dissolution and residual trapping). |
| Rate | 1,500 t/yr for 25 years | Per 100 m of along-strike width (the slice thickness). Specified as mass: surface volume rate × CO2 reference density. |
| Post-injection | 75 years | Fault-related leakage is mostly a post-injection, buoyancy-driven process. |
| Rel. perm. | Brooks-Corey; brine n = 3, Swr = 0.2; CO2 n = 2, Sgc = 0.05 (drainage), Sgr,max = 0.25 (imbibition); Killough hysteresis | Hysteresis is what produces residual trapping when brine re-enters the plume tail. |

## Checks
- Equilibrium: maximum pressure change after 5 years without wells should be close to zero.
- Mass balance: injected = in model + left model.

## Results (fault k_across = 0.01 mD, k_along = 0.1 mD; no capillary pressure)

Checks passed: the equilibrium drift was 0.0000 bar after 5 years, and the mass balance closes (37,500 t injected = 37,498 t in model + 2 t out). The run took 34 s.

| CO2 (t) | Year 25 | Year 100 |
|---|---|---|
| Footwall storage sand | 36,735 | 33,258 |
| Fault zone | 282 | 927 |
| Caprock | 446 | 1,331 |
| Overburden sand | 35 | 1,982 |
| Hanging-wall storage sand | 0 | 0 |

- Trapping at year 100: 46% mobile, 31% residual, 23% dissolved.
- Simulator densities: CO2 366 kg/m³, brine 996 kg/m³. So Δρ ≈ 630 kg/m³, and 1 bar of entry pressure holds about 16 m of CO2 (not 20 m as first estimated).
- Maximum pressure build-up: 1.2 bar. It stays small because of the open boundaries and the 300 mD sand.

Interpretation:
- CO2 rises to the caprock within ~5 years and reaches the fault by year 5.
- **Up-fault migration:** CO2 enters the fault column at the storage-sand level and climbs it. It reaches the top of the model by year 100 and spreads into the overburden on both sides of the fault (FW and HW).
- Leakage *accelerates* after injection stops, because buoyancy keeps driving it.
- **Caprock seepage:** about 1,300 t, mostly away from the fault. It enters across the whole plume footprint because a permeability-only seal has no entry threshold. This is expected to disappear with capillary pressure (Step 2b).
- About 11% of the injected CO2 is outside the storage sand by year 100. This number depends on k_along, which is a guessed value, hence the sweep in Step 3.
