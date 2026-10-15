# Step 4B, part 2: small 2D gas reservoir with a gas–water contact

Script: `scripts/04b2_gas_2d.jl` (`stability`, `gravity`). Same fluid model as part 1 (`src/gas.jl`): immiscible brine plus CO2–CH4 Peng–Robinson, with no dissolution, no diffusion and no capillary pressure at this stage.

## Model
A 600 m × 100 m × 60 m slab on 60 × 30 cells (10 m × 2 m, 1,800 cells), with its top at 1,500 m on the left, rising 3° to the right, so the crest is on the right. k_h = 100 mD, k_v = 10 mD and φ = 0.2, with closed boundaries (a volumetric compartment). The gas–water contact is at 1,520 m, with brine pressure 152 bar at the contact and 95 °C; 1,071 of the 1,800 cells are in the gas zone. One well at the crest is perforated over the top 8 m.

## Is the initial state at rest?
I equilibrated with `EquilibriumRegion`, putting both contacts (water–liquid and liquid–gas) at the gas–water contact so there is no hydrocarbon-liquid zone, and pure CH4, then ran 10 years without wells.

| Case | Initial pressure | max \|Δp\| after 10 yr | max \|ΔSw\| | max \|Δz\| |
|---|---|---|---|---|
| Gravity on (hydrostatic) | 151.6–155.8 bar | 4.3e-7 bar | 3.9e-12 | 1e-10 |
| Gravity off (uniform pressure) | 152.45 bar | 2.3e-3 bar | 7.9e-3 | 1e-10 |

This caught a problem in my set-up. JutulDarcy's equilibration takes the connate water saturation from tabulated relative-permeability curves (`ReservoirRelativePermeabilities`), but I had used the analytic Brooks–Corey object, which it does not read. It printed "Variable RelativePermeabilities is not a ReservoirRelativePermeabilities, skipping…" and left Sw = 0 in the gas zone, so there was 25% too much gas in place. That was my set-up choice, not a defect in JutulDarcy. The script now sets Sw = 0.20 in the gas zone after equilibration; that water is immobile (kr_w = 0 at S_wc) and there is no capillary pressure, so the pressure equilibrium is unchanged, as the rerun above confirms.

The small drift with gravity off is confined to cells either side of the contact (for example Sw 1.000 → 0.999 below and 0.000 → 0.0005 above), with pressure unchanged to 1e-4 bar. Without gravity the gas–water interface is neutrally stable: nothing holds the water below the gas, so pressure differences of about 1 Pa (within the solver tolerance) move water slowly across it. Measured with `scripts/04b2_gravityoff_drift.jl`: over 10 years, 11.9 m³ of water (3e-5 of the water pore volume) moved up across the contact in 39 of the 1,800 cells, total water mass was conserved to 2e-8, and the contact moved by at most 0.003 m. (I first called this round-off; it is a small real redistribution, negligible for the gravity comparison.) With gravity on, the same 1 Pa is nothing against the ~17,600 Pa the gas–brine density difference provides across one 2 m cell, which is why that case does not move.

To switch gravity off I set `parameters[:Reservoir][:TwoPointGravityDifference] .= 0`. This is the per-face g·Δz used by the two-point flux; I checked in the source that the standard TPFA scheme takes gravity only from this parameter (the `PhasePotentials` path is used by other discretisations). The gravity-off run starts from the same saturations and composition with uniform pressure, which is its own equilibrium.

## Production, CO2 injection, gravity on and off
Through the crest well: 10 years of CH4 production at 1,500 t/yr (bottom-hole floor 30 bar), 1 year shut in, 10 years of CO2 injection at 6,000 t/yr, and 30 years shut in. The two runs are identical apart from gravity (gdz = 0 and uniform initial pressure in the gravity-off case). Figures: `figures/04b2_gravity.png`, `figures/04b2_segregation.png`.

To measure segregation I use the mass-weighted mean depth of each component, Σ M_i z_i / Σ M_i, from the cell CO2 and CH4 masses. A positive "CO2 − CH4" means CO2 sits deeper on average.

| CO2 − CH4 mean depth | Year 21 (end of injection) | Year 31 | Year 51 |
|---|---|---|---|
| Gravity on | +9.8 m | +13.3 m | +13.3 m |
| Gravity off | −10.7 m | −10.7 m | −10.7 m |

Without gravity, the result shows the effect of the well location alone: CO2 injected at the crest stays there and pushes CH4 down-dip, so the CO2 ends up 10.7 m shallower on average. With gravity, CO2-rich gas moves beneath the CH4 along the whole structure, and by year 51 the gas zone is stratified, with CO2-rich gas next to the contact and CH4-rich gas at the crest. Segregation is essentially complete about 10 years after injection stops.

The two cases differ by 24 m (+13.3 − (−10.7)) in the relative depth of CO2 and CH4. That is strong evidence that gravity matters and that injecting at the crest works against segregation, but it is not the isolated gravity effect: subtracting the gravity-off case assumes the two effects add linearly, and without gravity the CO2 also spreads differently. (An earlier draft called this "the gravity-driven shift of about 24 m"; I corrected it.) The gravity-off reference starts from its own equilibrium and was checked with its own no-well test.

CO2 and CH4 here are components of one miscible gas phase, so there is no sharp interface between them, only a continuous vertical gradient in composition. How sharp it is depends on numerical dispersion (2 m cells, no diffusion or mechanical dispersion), so the shape of the gradient is not quantitatively reliable; the direction and approximate size of the shift are.

### Conservation and solver behaviour

| | Gravity on | Gravity off |
|---|---|---|
| Wall time | 33.8 min | 1.4 min |
| Ministeps accepted / cut | 2,288 / 3,462 | 666 / 0 |
| Max \|residual\| CO2 / CH4 / water | 1.9 / 2.3 / 10.8 t | 0.07 / 0.02 / 0.11 t |
| Relative to component inventory | ≤ 3e-5 | ≤ 3e-6 |

With gravity the solver cut more time steps than it accepted, and the run took about 24 times longer. Buoyant counter-current flow of CO2-rich and CH4-rich gas, and of gas next to brine at the contact, is a hard nonlinear problem, although conservation is still good (≤ 3e-5). This was the main practical risk for the faulted model in Step 4D, which has 4.5 times as many cells, so I looked at the solver settings first (below).

### Other observations
The average pressure at the end is 158.6 bar with gravity and 153.4 bar without. These are not directly comparable, because the averages include the water zone, whose initial pressure differs between hydrostatic and uniform initialisation; the initial CH4 in place differs slightly (29,640 against 29,759 t) for the same reason.

Injecting 60,000 t of CO2 after producing 15,000 t of CH4 pushed the pore-volume-averaged pressure above its starting value (158.6 bar against about 152.5 bar, gravity on). That is an average over different depths rather than a pressure at a reference depth; the faulted model should look at pressure next to and across the fault at matched depths. Pressure alone also says nothing about fault slip, which would need a mechanical model.

I first wrote that the direction and rough size of the segregation were "reliable", which was premature without any sensitivity check. The solver-tuning runs below, with different time-step sequences, give a first one.

Because CO2-rich gas collects at the gas–water contact, dissolution into the brine (left out here) would probably be enhanced there, which is one reason the Søreide–Whitson work in Step 4A still matters.

## Solver tuning (gravity on)
Script: `scripts/04b2_tuning.jl`. I kept the tolerances (tol_cnv 1e-3, tol_mb 1e-7) and varied only the time stepping and the nonlinear update controls.

| Run | Newton iterations | Ministeps ok / cut | Wall time | CO2−CH4 mean depth, yr 21 / 31 / 51 | p_avg yr 51 | Max rel. mass residual |
|---|---|---|---|---|---|---|
| Baseline (max Δt 30 d) | 61,053 | 2,288 / 3,462 | 16.2 min | +9.78 / +13.31 / +13.28 m | 158.60 bar | 7.8e-5 |
| A: max Δt 10 d | 75,591 | 3,751 / 4,171 | 20.8 min | +9.88 / +13.38 / +13.35 m | 158.55 bar | 2.5e-4 |
| B: relaxation | 52,987 | 2,084 / 2,974 | 15.1 min | +9.78 / +13.31 / +13.28 m | 158.60 bar | 6.8e-5 |
| C: step growth ≤ 2× | 49,693 | 3,231 / 2,635 | 15.6 min | +9.79 / +13.32 / +13.29 m | 158.59 bar | 9.4e-5 |
| D: dz_max 0.05 | 73,785 | 2,680 / 4,183 | 22.3 min | +9.80 / +13.32 / +13.29 m | 158.60 bar | 4.7e-5 |

No steps are cut during production and the first shut-in, about 460–1,200 during injection, and about 2,200–3,400 (80–90% of them) during the 30-year shut-in after injection.

The physical results do not depend on the time-step sequence: across the five runs the segregation measure varies by at most 0.07 m and the average pressure by at most 0.05 bar. That is a first time-step check for this case, not a grid check.

The gains from tuning are modest. Limiting time-step growth to 2× (C) cuts Newton iterations by 19% and relaxation (B) by 13%, while a smaller maximum step (A) or a tighter composition limit (D) increases the work, so smaller steps are not automatically cheaper. Wall time is not a reliable measure on this laptop (the same baseline took 33.8 minutes earlier and 16.2 minutes here, depending on what else was running), so I compare Newton iterations.

I then combined relaxation with step growth ≤ 2× (`scripts/04b2_diagnose.jl`): 48,294 Newton iterations (−21% against the baseline), 3,171 ministeps accepted and 2,551 cut, and the same physical results. My first guess for the cuts was the phase label switching, because dense CO2-rich gas (up to 306 kg/m³) could be labelled liquid, but there is no hydrocarbon-liquid saturation in any cell at any time. The cuts come in bursts in particular years (year 33: 934, year 27: 492, year 31: 289, year 18: 104) with none in most years between, so the time step grows, meets a difficult transient and collapses. I did not find the root cause; that would need per-iteration convergence logs (which equation and which cells fail), which the default report level does not keep. I stopped there because the physical results are robust and the cost is acceptable for a few runs.

For the faulted model (Step 4D) I use relaxation with step growth ≤ 2×, keep the default tolerances, budget about an hour per run (4.5 times the cells) and avoid sweeps on that grid. (In Step 4D I later had to tighten tol_mb; see that note.)
