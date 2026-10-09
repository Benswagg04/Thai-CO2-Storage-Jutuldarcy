# Step 2b: Capillary pressure

## Set-up
- One Brooks-Corey drainage curve per facies: Pc = Pe·Se^(−1/λ), λ = 2, Se = (1 − S_CO2 − Swr)/(1 − Swr), capped at Se = 0.01.
- Entry pressure scaled from the storage sand (Pe = 0.05 bar at 300 mD, φ = 0.25) with Leverett J, Pe ∝ √(φ/k).
- Resulting entry pressures:

| Unit | Pe (bar) | CO2 column it can hold (Δρ ≈ 630 kg/m³) |
|---|---|---|
| Storage sand | 0.05 | – |
| Overburden sand | 0.057 | – |
| Caprock / underburden shale | 12.2 | ~200 m |
| Clay-rich fault (k_across 0.01 mD) | 5.48 | ~89 m |
| Sand-rich fault (k_across 1 mD) | 0.55 | ~9 m |

- JutulDarcy convention, confirmed in the source: brine is the reference phase, and **p_CO2 = p_brine + Pc(S_CO2)**. The table therefore starts at Pe when S_CO2 = 0.
- The fault uses the *across-fault* permeability for its entry pressure, because CO2 has to enter it from the side. A single scalar Pc per cell cannot represent anisotropic pore-throat sizes; this is a limitation.

## Checks
1. Equilibrium: 0.0000 bar drift in all cases.
2. Regression: the refactored code reproduces the Step 2 baseline to 0.1 t.
3. Positive control: the sand-rich fault leaks. So capillary pressure is not simply blocking all flow.
4. Threshold arithmetic: CO2 column against the fault vs. fault entry pressure.

| Case | Column at fault | Buoyancy (Δρ g h) | Fault Pe | Outcome |
|---|---|---|---|---|
| Clay-rich fault | 18 m | 1.11 bar | 5.48 bar | Holds |
| Sand-rich fault, yr 25 | 16 m | 0.99 bar | 0.55 bar | Leaks |
| Sand-rich fault, yr 100 | 14 m | 0.87 bar | 0.55 bar | Still leaking; column draining toward ~9 m |

## Results at year 100 (% of 37,500 t injected outside the footwall storage sand)
| Case | Outside storage sand |
|---|---|
| No Pc, clay-rich fault (Step 2) | 11% (fault chimney + caprock seepage) |
| Pc, clay-rich fault | ~0% |
| Pc, sand-rich fault (k 1/10 mD) | 46% (overburden 9,730 t; 7,120 t left the model) |

Interpretation:
- Capillary pressure removes the caprock seepage completely, as the threshold argument predicts (1 bar of buoyancy against 12 bar of entry pressure).
- With capillary pressure, fault seal becomes a **threshold** process. The fault holds completely until the CO2 column exceeds Pe/(Δρ g), then it leaks. Leakage in the sand-rich case starts at about year 8, once the column has grown tall enough.
- After injection, the leaking fault bleeds the column down toward the height it can hold (~9 m), so the leakage rate slows.
- The sand-rich fault also has 100× higher permeability, so the sand-rich vs clay-rich comparison changes two things at once. Step 3 separates them.

## Limitations noticed
- The model top is a no-flow boundary. CO2 that reaches the overburden pools against it, which is artificial. In the sand-rich case, CO2 also migrates up-dip along the top of the hanging wall and exits through the right boundary ("left the model").
- Leverett scaling and λ = 2 are assumptions; no measured capillary data are used.
