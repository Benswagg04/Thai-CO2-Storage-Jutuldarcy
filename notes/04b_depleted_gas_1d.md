# Step 4B, stage 1: one-dimensional depleted gas reservoir (Path B)

Scripts: `src/gas.jl`, `scripts/04b_depletion_1d.jl`. Figure: `figures/04b_depletion_1d.png`. Runs on the existing environment (JutulDarcy 0.3.10).

## Fluid model and what it leaves out
- Three phases: brine as a separate **immiscible** phase, plus CO2 and CH4 described by Peng–Robinson (PR) with k_CO2–CH4 = 0.10 (a typical value, assumed).
- At 95 °C both gases are supercritical, so the mixture stays in one gas-like phase. Checked: the hydrocarbon-liquid saturation is exactly 0 throughout.
- **Excluded:** CO2 and CH4 dissolving in brine, water evaporating into the gas, molecular diffusion, gravity (1D, horizontal). This is the deliberate "Path B" simplification while the Søreide–Whitson route (Step 4A) is validated.
- Brine properties set explicitly. JutulDarcy's defaults for the immiscible water phase are a reference density of 1 kg/m³ and a viscosity of 1 mPa·s. We use 990 kg/m³ at 150 bar, a compressibility of 4.5e-10 1/Pa and a viscosity of 0.35 mPa·s.
- PR CO2 density vs JutulDarcy's NIST-based CO2 table at 95 °C, 40–150 bar: within ±1.7% (no volume shift needed for the gas).

## Case
- 1 km × 100 m × 20 m; k = 100 mD; φ = 0.2; connate water S_wc = 0.2 (immobile); initially pure CH4 at 150 bar, 95 °C.
- Schedule:
  - 10 yr production at 2,000 t/yr CH4 (bottom-hole pressure floor 30 bar, not reached);
  - 1 yr shut-in;
  - 10 yr CO2 injection at 8,000 t/yr at the opposite end;
  - 10 yr shut-in.
- Grids: 100 and 400 cells.

## Checks
| Check | Result |
|---|---|
| Component mass conservation (well flows integrated over every ministep) | CO2, CH4 and water all within 1–3 × 10⁻⁶ of initial + injected |
| p/Z material balance during depletion | Within 0.012% (year 1) to 0.122% (year 10) of the closed-tank line p/Z = (p/Z)_i (1 − G_p/G) |
| Explanation of the p/Z drift | Expansion of connate water as pressure falls takes up gas volume. The simple p/Z equation ignores this; rock compressibility is zero in this model, so no further correction is needed. Predicted size c_w Δp S_w/(1 − S_w) with c_w = 4.5e-10 1/Pa, Δp = 150 − 42.1 = 107.9 bar, S_w = 0.2: 4.5e-10 × 1.079e7 × 0.25 = 0.121%, vs 0.122% simulated. The same estimate matches every year of the decline (e.g. year 5: Δp = 53.7 bar → 0.060%, simulated 0.060%). |
| Phase state | Single gas phase throughout |
| Production rate | 20,007.5 t over 10 years vs 20,000 t target |

The p/Z test is an independent check: the line comes from classical gas reservoir engineering and the PR compressibility factor, not from the simulator.

## Findings about JutulDarcy 0.3.10 well controls
- `TotalMassRateTarget` cannot be used for producers in this version: the automatic "mrat" limit is rejected by the limit checker ("Unsupported well injector constraint/limit mrat"). Surface-volume targets (`TotalRateTarget`) were used instead.
- `TotalRateTarget` means different things for the two well types in a compositional model:
  - **injector:** surface volume × the `density` given to `InjectorControl` = mass (80,000 t injected exactly);
  - **producer:** volume at standard conditions (1 atm, 15 °C) from a flash of the produced fluid. A target of "2,000" first produced 1,361 t/yr; the ratio 1.47 m³/kg is the standard volume of 1 kg of methane. The producer target is now converted with the PR standard-condition density.

## Results
- Average pressure falls from 150 to 42 bar over 10 years of production (20,000 t CH4). In this closed 1D numerical experiment, injecting 80,000 t of CO2 afterwards restored the pressure to about its original value (150.4 / 150.2 bar on the 100- / 400-cell grids).
  - The ~4:1 mass ratio is specific to this model and its operating conditions (residual CH4, pressure, temperature, water saturation, compressibilities). It is consistent with a first-order voidage-replacement estimate (CO2 at 344 kg/m³ vs CH4 at 85 kg/m³ at 150 bar, 95 °C), but it is not a general storage-capacity relationship.
- CO2 displaces CH4 as a front, with a mixing zone between them (numbers below).
- Nothing changes during the final shut-in: there is no diffusion and no gravity in this 1D model.

### Mixing zone and numerical dispersion (end of injection, year 21)
| Grid / max time step | Ministeps | Front (z_CO2 = 0.5) | 90%→10% width |
|---|---|---|---|
| 100 cells (10 m), 30 d | 407 | 741.5 m | 138.3 m |
| 400 cells (2.5 m), 30 d | 407 | 738.2 m | 101.8 m |
| 400 cells (2.5 m), 3 d | 3,784 | 734.9 m | 60.9 m |

- The front position is robust (735–742 m).
- The width is **numerical dispersion from both space and time discretisation.** Refining cells 4× narrows it by only 26% at a fixed time step; cutting the time step 10× narrows it a further 40%. Fully implicit time stepping smears fronts just as coarse cells do.
- The model contains no physical mixing process: neither molecular diffusion nor mechanical dispersion. The simulated width is therefore controlled by numerical dispersion, and **the physically realistic mixing width remains uncertain**.
  - Molecular diffusion alone would give a few metres over 10 years (√(Dt), D ~ 1e-7 m²/s).
  - Field-scale mechanical dispersion (dispersivity × velocity) is not modelled and could be comparable to the simulated widths.
  - *(Earlier draft said "10–30× wider than physical"; withdrawn after review, because it compared against diffusion only.)*
- Effect on pressure: when the mixing zone narrowed from 138 m to 102 m (grid refinement), the end-of-injection average pressure changed from 150.43 to 150.15 bar (0.2%). Small in this case, but not zero: mixing changes mixture density and compressibility.
- Larger effects are expected on quantities controlled by mixing: CO2 breakthrough at a producer, contamination of residual CH4, and the CO2 content of any back-produced gas.

## What this does and does not show
- Shows: the three-phase compositional set-up works in JutulDarcy 0.3.10, conserves mass, and reproduces classical gas material balance to within a physically explained 0.1%.
- Does not show: dissolution effects; gravity segregation of CO2 under CH4 (CO2 is about 4× denser and will slump in 2D/3D); behaviour near a fault. Those are the next stages (2D section, then the faulted model).
