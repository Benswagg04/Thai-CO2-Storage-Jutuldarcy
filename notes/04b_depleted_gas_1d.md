# Step 4B, part 1: one-dimensional depleted gas reservoir

Scripts: `src/gas.jl`, `scripts/04b_depletion_1d.jl`. Figure: `figures/04b_depletion_1d.png`. Uses the same environment as before (JutulDarcy 0.3.10).

## Fluid model and what it leaves out
There are three phases: brine as a separate, immiscible phase, and CO2 and CH4 described by Peng–Robinson with k_CO2–CH4 = 0.10 (a typical value, assumed). At 95 °C both gases are supercritical, so the mixture stays in one gas-like phase; I checked that the hydrocarbon-liquid saturation is exactly zero throughout.

This leaves out CO2 and CH4 dissolving in brine, water evaporating into the gas, molecular diffusion and (in 1D) gravity. It is a deliberate simplification while the Søreide–Whitson route from Step 4A is not yet usable.

JutulDarcy's defaults for the immiscible water phase are a reference density of 1 kg/m³ and a viscosity of 1 mPa·s, so I set the brine properties explicitly: 990 kg/m³ at 150 bar, compressibility 4.5e-10 1/Pa, viscosity 0.35 mPa·s. The Peng–Robinson CO2 density agrees with JutulDarcy's NIST-based CO2 table to within ±1.7% at 95 °C and 40–150 bar, so the gas needs no volume shift.

## Case
A 1 km × 100 m × 20 m reservoir with k = 100 mD, φ = 0.2 and immobile connate water S_wc = 0.2, initially full of CH4 at 150 bar and 95 °C. The schedule is 10 years of production at 2,000 t/yr CH4 (bottom-hole pressure floor 30 bar, never reached), 1 year shut in, 10 years of CO2 injection at 8,000 t/yr at the opposite end, and 10 more years shut in. I ran it on 100 and 400 cells.

## Checks

| Check | Result |
|---|---|
| Mass conservation per component (well flows integrated over every ministep) | CO2, CH4 and water all within 1–3 × 10⁻⁶ of initial + injected |
| p/Z material balance during depletion | Within 0.012% (year 1) to 0.122% (year 10) of the closed-tank line p/Z = (p/Z)_i (1 − G_p/G) |
| Reason for the p/Z difference | Connate water expands as the pressure falls and takes up gas volume, which the simple p/Z equation ignores (rock compressibility is zero here). The expected size is c_w Δp S_w/(1 − S_w) with c_w = 4.5e-10 1/Pa, Δp = 150 − 42.1 = 107.9 bar and S_w = 0.2: 4.5e-10 × 1.079e7 × 0.25 = 0.121%, against 0.122% simulated. The same estimate matches every year (e.g. year 5: Δp = 53.7 bar → 0.060%, simulated 0.060%). |
| Phase state | Single gas phase throughout |
| Production | 20,007.5 t over 10 years against a 20,000 t target |

The p/Z test is independent of the simulator: the line comes from classical gas reservoir engineering and the Peng–Robinson compressibility factor.

## Well controls in JutulDarcy 0.3.10
`TotalMassRateTarget` cannot be used for producers in this version: the automatic "mrat" limit is rejected ("Unsupported well injector constraint/limit mrat"), so I used surface-volume targets (`TotalRateTarget`). In a compositional model these mean different things for the two well types. For the injector, surface volume × the `density` given to `InjectorControl` is the mass, and exactly 80,000 t went in. For the producer, it is the volume at standard conditions (1 atm, 15 °C) from a flash of the produced fluid: my first target of "2,000" produced 1,361 t/yr, and the ratio, 1.47 m³/kg, is the standard volume of 1 kg of methane. The producer target is now converted with the Peng–Robinson standard-condition density.

## Results
The average pressure falls from 150 to 42 bar over the 10 years of production (20,000 t of CH4). In this closed 1D model, injecting 80,000 t of CO2 afterwards brought the pressure back to about its original value (150.4 bar on 100 cells, 150.2 bar on 400). The roughly 4:1 mass ratio is specific to this model and its conditions (residual CH4, pressure, temperature, water saturation, compressibilities). It fits a first-order voidage-replacement estimate (CO2 at 344 kg/m³ against CH4 at 85 kg/m³ at 150 bar and 95 °C), but it is not a general storage-capacity rule. CO2 displaces CH4 as a front with a mixing zone between them, and nothing changes in the final shut-in because there is neither diffusion nor gravity in 1D.

### Mixing zone and numerical dispersion (end of injection, year 21)

| Grid / max time step | Ministeps | Front (z_CO2 = 0.5) | 90%→10% width |
|---|---|---|---|
| 100 cells (10 m), 30 d | 407 | 741.5 m | 138.3 m |
| 400 cells (2.5 m), 30 d | 407 | 738.2 m | 101.8 m |
| 400 cells (2.5 m), 3 d | 3,784 | 734.9 m | 60.9 m |

The front position is robust (735–742 m), but the width comes from numerical dispersion in both space and time. Refining the cells 4 times narrows it by only 26% at a fixed time step, and cutting the time step 10 times narrows it by a further 40%: fully implicit time stepping smears fronts just as coarse cells do. The model has no physical mixing process (no molecular diffusion or mechanical dispersion), so the realistic mixing width is uncertain. Molecular diffusion alone would give a few metres over 10 years (√(Dt), D ~ 1e-7 m²/s), but field-scale mechanical dispersion (dispersivity × velocity) is not modelled and could be comparable to the simulated widths. (An earlier draft called the simulated zone "10–30 times wider than physical"; I withdrew that because it compared against diffusion only.)

When the mixing zone narrowed from 138 m to 102 m with grid refinement, the average pressure at the end of injection changed from 150.43 to 150.15 bar (0.2%). That is small here but not zero, because mixing changes the density and compressibility of the gas. I would expect larger effects on quantities controlled by mixing: CO2 breakthrough at a producer, contamination of the remaining CH4, and the CO2 content of any gas produced back.

## What this does and does not show
It shows that the three-phase compositional set-up works in JutulDarcy 0.3.10, conserves mass, and reproduces classical gas material balance to within an explained 0.1%. It says nothing yet about dissolution, about gravity segregation of CO2 under CH4 (CO2 is about 4 times denser and will slump in 2D or 3D), or about behaviour near a fault. Those come next: a 2D section (part 2), then the faulted model (Step 4D).
