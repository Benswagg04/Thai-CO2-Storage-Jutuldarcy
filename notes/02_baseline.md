# Step 2: Baseline CO2 injection (no capillary pressure)

## Inputs

| Input | Value | Reason |
|---|---|---|
| Fluid model | `:co2brine`, `co2_physics = :kvalue` | Two components (H2O, CO2) and two phases. CO2 dissolves in brine and water evaporates into the CO2, using K-values from the Saló-Salgado et al. (2024) correlations (JutulDarcy's default `co2_source = :salo24`). |
| Temperature | 95 °C, isothermal | Low end of the Pattani Basin gradient (≈27 °C + 1.5 km × 45 °C/km). The correlations are validated to about 100 °C; JutulDarcy's tables go to 200 °C but are built without range checks. |
| Salinity | 30,000 mg/L NaCl-equivalent (hypothetical) | A moderate salinity chosen for illustration, not a documented basin value. NaCl mole fraction about 0.0094. |
| Initial pressure | Hydrostatic from 151 bar at 1,500 m | Normal gradient (~0.1 bar/m), set with `EquilibriumRegion` and checked by running 5 years with no wells. |
| Boundaries | Left and right edges held at initial pressure | The formations continue beyond the slice; this avoids an unrealistic pressure build-up in a closed 2D box. |
| Injector | Column 35 (x ≈ 345 m), bottom 20 m of the storage sand | Injecting low lets the CO2 contact more brine on its way up (more dissolution and residual trapping). |
| Rate | 1,500 t/yr for 25 years | Per 100 m of along-strike width. Given as mass: surface volume rate × CO2 reference density. |
| Post-injection | 75 years | Leakage along faults is mostly driven by buoyancy after injection stops. |
| Relative permeability | Brooks–Corey; brine n = 3, Swr = 0.2; CO2 n = 2, Sgc = 0.05 (drainage), Sgr,max = 0.25 (imbibition); Killough hysteresis | Hysteresis is what traps CO2 residually when brine moves back into the tail of the plume. |

## Results (fault k_across = 0.01 mD, k_along = 0.1 mD)

The model was in equilibrium (pressure drift 0.0000 bar over 5 years without wells), and the CO2 adds up: 37,500 t injected, 37,498 t in the model and 2 t out. The run took 34 s.

| CO2 (t) | Year 25 | Year 100 |
|---|---|---|
| Footwall storage sand | 36,735 | 33,258 |
| Fault zone | 282 | 927 |
| Caprock | 446 | 1,331 |
| Overburden sand | 35 | 1,982 |
| Hanging-wall storage sand | 0 | 0 |

At year 100, 46% of the CO2 is mobile, 31% residually trapped and 23% dissolved. The simulator gives densities of 366 kg/m³ for CO2 and 996 kg/m³ for brine, so Δρ ≈ 630 kg/m³ and 1 bar of entry pressure holds about 16 m of CO2, not the 20 m I estimated in Step 1. The maximum pressure build-up is only 1.2 bar, thanks to the open boundaries and the 300 mD sand.

CO2 reaches the caprock and the fault within about 5 years. It then enters the fault column at the storage-sand level, climbs it, reaches the top of the model by year 100 and spreads into the overburden on both sides. The leakage speeds up after injection stops, because buoyancy keeps driving it.

About 1,300 t also seeps into the caprock, mostly away from the fault. With permeability as the only seal there is no entry threshold, so CO2 enters across the whole footprint of the plume; I expect this to disappear once capillary pressure is added (Step 2b).

In total about 11% of the injected CO2 is outside the storage sand by year 100. This depends on k_along, which is a guess, so Step 3 sweeps it.
