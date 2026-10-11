# Step 4B, stage 2: small 2D gas reservoir with a gas–water contact

Script: `scripts/04b2_gas_2d.jl` (`stability`, `gravity`). Same fluid model as stage 1 (`src/gas.jl`): immiscible brine plus CO2–CH4 Peng–Robinson, no dissolution, no diffusion, **no capillary pressure** in this stage.

## Model
- 600 m × 100 m × 60 m slab, 60 × 30 cells (10 m × 2 m), 1,800 cells. Top at 1,500 m on the left, rising 3° to the right (crest on the right).
- k_h = 100 mD, k_v = 10 mD, φ = 0.2. Closed boundaries (a volumetric compartment).
- Gas–water contact (GWC) at 1,520 m; brine pressure 152 bar at the contact; 95 °C. 1,071 of 1,800 cells are in the gas zone.
- One well at the crest, perforated over the top 8 m.

## Stage B2a: is the initial state at rest?
Equilibration: `EquilibriumRegion` with both contacts (water–liquid and liquid–gas) at the GWC, so there is no hydrocarbon-liquid zone, and composition pure CH4. Then 10 years without wells.

| Case | Initial pressure | max \|Δp\| after 10 yr | max \|ΔSw\| | max \|Δz\| |
|---|---|---|---|---|
| Gravity on (hydrostatic) | 151.6–155.8 bar | 4.3e-7 bar | 3.9e-12 | 1e-10 |
| Gravity off (uniform pressure) | 152.45 bar | 2.3e-3 bar | 7.9e-3 | 1e-10 |

**Set-up limitation found and corrected: connate water.** JutulDarcy's equilibration takes the connate water saturation from tabulated relative-permeability curves (`ReservoirRelativePermeabilities`). Our set-up used the analytic Brooks–Corey object instead, which the equilibration does not read. It printed *"Variable RelativePermeabilities is not a ReservoirRelativePermeabilities, skipping…"* and left Sw = 0 in the gas zone, i.e. 25% too much gas in place. This is a consequence of our set-up choice, not a defect in JutulDarcy. The script now sets Sw = 0.20 in the gas zone after equilibration. That water is immobile (kr_w = 0 at S_wc) and there is no capillary pressure, so the pressure equilibrium is unchanged (confirmed by the rerun above).

**Gravity-off drift explained.** The changes are confined to cells on either side of the contact (e.g. Sw 1.000 → 0.999 below, 0.000 → 0.0005 above), and pressure is unchanged to 1e-4 bar. Without gravity, the gas–water interface is neutrally stable: there is nothing to hold water below gas, so pressure differences of ~1 Pa from solver round-off move a few litres of water a year. With gravity on, the same 1 Pa is negligible against the ~17,600 Pa that the gas–brine density difference provides across one 2 m cell, which is why that case shows no movement.

**How gravity is switched off:** `parameters[:Reservoir][:TwoPointGravityDifference] .= 0`. This is the per-face g·Δz used by the two-point flux; the source was checked to confirm that the standard TPFA scheme takes gravity only from this parameter (the alternative `PhasePotentials` path is used by non-TPFA discretisations). The gravity-off initial state uses the same saturations and composition with uniform pressure, which is its equilibrium.

## Stage B2b: production, CO2 injection, gravity on vs off
- Schedule through the crest well: 10 yr CH4 production at 1,500 t/yr (bottom-hole pressure floor 30 bar), 1 yr shut-in, 10 yr CO2 injection at 6,000 t/yr, 30 yr shut-in.
- The two runs are identical except gravity: gravity-off has gdz = 0 and a uniform initial pressure.
- Figures: `figures/04b2_gravity.png`, `figures/04b2_segregation.png`.

### Segregation metric
Mass-weighted mean depth of each component, Σ M_i z_i / Σ M_i, using cell CO2 and CH4 masses. Positive "CO2 − CH4" means CO2 sits deeper on average.

| CO2 − CH4 mean depth | Year 21 (end of injection) | Year 31 | Year 51 |
|---|---|---|---|
| Gravity on | +9.8 m | +13.3 m | +13.3 m |
| Gravity off | −10.7 m | −10.7 m | −10.7 m |

- **Gravity off isolates the effect of the well location.** CO2 injected at the crest stays there and pushes CH4 down-dip, so CO2 ends up 10.7 m *shallower* on average.
- **Gravity on:** CO2-rich gas moves beneath the CH4 along the whole structure. By year 51 the gas zone is stratified: CO2-rich gas next to the gas–water contact, CH4-rich gas at the crest.
- **The two simulated cases differ by 24 m (+13.3 − (−10.7)) in the relative mean depth of CO2 and CH4.** This is strong evidence that gravity matters, and that the crest well location works against segregation. It should not be read as the isolated gravity contribution: subtracting the gravity-off case assumes the well-location and gravity effects add linearly, but without gravity the CO2 also spreads in a different pattern. *(Earlier draft: "the gravity-driven shift is about 24 m"; revised after review.)*
- The gravity-off reference starts from its own equilibrium (same saturations, uniform pressure), not from the hydrostatic pressure field, and was checked with its own no-well test (stage B2a).
- Segregation is essentially complete about 10 years after injection stops (+13.3 m from year 31 onwards).
- These are components of one miscible gas phase, so there is no sharp CO2–CH4 interface; the maps show a continuous vertical composition gradient. Its sharpness is controlled by numerical dispersion (2 m cells, no diffusion or mechanical dispersion modelled), so the shape of the gradient is not quantitatively reliable; the direction and the approximate size of the shift are.

### Conservation and solver behaviour
| | Gravity on | Gravity off |
|---|---|---|
| Wall time | 33.8 min | 1.4 min |
| Ministeps accepted / cut | 2,288 / 3,462 | 666 / 0 |
| Max \|residual\| CO2 / CH4 / water | 1.9 / 2.3 / 10.8 t | 0.07 / 0.02 / 0.11 t |
| Relative to component inventory | ≤ 3e-5 | ≤ 3e-6 |

- With gravity the solver cut more steps than it accepted, and the run took about 24× longer. Buoyancy-driven counter-current flow of CO2-rich and CH4-rich gas (and gas next to brine at the contact) is a hard nonlinear problem. Conservation is still good (≤ 3e-5).
- **This is the main practical risk for stage B2c:** the faulted section has 4.5× as many cells. Before scaling up, test whether a smaller maximum time step reduces the cutting.

### Other observations
- Average pressure at the end: 158.6 bar (gravity on) vs 153.4 bar (gravity off). The two are not directly comparable: the averages include the water zone, whose initial pressure differs between hydrostatic and uniform initialisation. The initial CH4 in place also differs slightly (29,640 vs 29,759 t) for the same reason.
- Injecting 60,000 t CO2 after producing 15,000 t CH4 raised the pore-volume-averaged pressure above its initial value (158.6 bar vs about 152.5 bar initially, gravity on). This is an average over different depths, not a pressure at a reference depth. Stage B2c should measure the pressure next to the fault and across it, at matched depths. Pressure alone says nothing about fault slip; that would need a mechanical model.
- The statement in an earlier draft that the direction and rough size of the segregation are "reliable" was premature: it had no time-step or grid sensitivity behind it. The solver-tuning runs (different time-step sequences) provide a first time-step check.
- Because CO2-rich gas collects at the gas–water contact, dissolution into brine (excluded here) would likely be enhanced there. This is one reason the Søreide–Whitson work (Path A) still matters.

## Solver tuning (gravity-on case)
Script: `scripts/04b2_tuning.jl`. Tolerances unchanged (tol_cnv 1e-3, tol_mb 1e-7); only time stepping and nonlinear update controls varied.

| Run | Newton iterations | Ministeps ok / cut | Wall time | CO2−CH4 mean depth, yr 21 / 31 / 51 | p_avg yr 51 | Max rel. mass residual |
|---|---|---|---|---|---|---|
| Baseline (max Δt 30 d) | 61,053 | 2,288 / 3,462 | 16.2 min | +9.78 / +13.31 / +13.28 m | 158.60 bar | 7.8e-5 |
| A: max Δt 10 d | 75,591 | 3,751 / 4,171 | 20.8 min | +9.88 / +13.38 / +13.35 m | 158.55 bar | 2.5e-4 |
| B: relaxation | 52,987 | 2,084 / 2,974 | 15.1 min | +9.78 / +13.31 / +13.28 m | 158.60 bar | 6.8e-5 |
| C: step growth ≤ 2× | 49,693 | 3,231 / 2,635 | 15.6 min | +9.79 / +13.32 / +13.29 m | 158.59 bar | 9.4e-5 |
| D: dz_max 0.05 | 73,785 | 2,680 / 4,183 | 22.3 min | +9.80 / +13.32 / +13.29 m | 158.60 bar | 4.7e-5 |

- **Where steps are cut:** 0 during production and the first shut-in, about 460–1,200 during injection, and about 2,200–3,400 (about 80–90%) during the 30-year post-injection shut-in.
- **Physical results are insensitive to the time-step sequence:** the segregation metric varies by ≤ 0.07 m and the average pressure by ≤ 0.05 bar across five different sequences. This is the first time-step sensitivity check for stage B2b. It is not a grid check.
- **Tuning gives modest gains.** Limiting time-step growth to 2× (C) reduces Newton iterations by 19%, relaxation (B) by 13%. A smaller maximum step (A) and a tighter composition limit (D) increase the work. Smaller steps are not automatically cheaper.
- **Wall time is not a reliable metric on this laptop.** The same baseline took 33.8 min in the earlier run and 16.2 min here, depending on machine load. Newton iterations are used for comparison.
- **Diagnosis** (`scripts/04b2_diagnose.jl`, relaxation + step growth ≤ 2× combined):
  - Newton iterations 48,294 (−21% vs baseline); ministeps 3,171 ok / 2,551 cut; physical results unchanged.
  - **Phase-label switching ruled out:** no hydrocarbon-liquid saturation in any cell at any time (max S_l = 0). That was the first hypothesis, because dense CO2-rich gas (up to 306 kg/m³) could be labelled liquid.
  - **Cuts are episodic, not steady.** Bursts in particular years (e.g. year 33: 934, year 27: 492, year 31: 289, year 18: 104), with no cuts in most years between. The time step grows, meets a difficult transient, and collapses.
  - **Root cause not identified.** Finding it would need per-iteration convergence logs (which equation and which cells fail), which are not kept at the default report level. Timeboxed rather than pursued, because the physical results are robust and the cost is acceptable for a small number of runs.
- **Decision for stage B2c:** use relaxation + step growth ≤ 2×; keep default tolerances; budget about 1 hour per run on the faulted model (4.5× the cells); avoid sweeps on that grid.
