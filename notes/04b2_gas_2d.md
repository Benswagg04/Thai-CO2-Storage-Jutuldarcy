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

**Issue found and fixed: connate water.** JutulDarcy's equilibration takes the connate water saturation from tabulated relative-permeability curves (`ReservoirRelativePermeabilities`) only. With the Brooks–Corey object it printed *"Variable RelativePermeabilities is not a ReservoirRelativePermeabilities, skipping…"* and left Sw = 0 in the gas zone, i.e. 25% too much gas in place. The script now sets Sw = 0.20 in the gas zone after equilibration. That water is immobile (kr_w = 0 at S_wc) and there is no capillary pressure, so the pressure equilibrium is unchanged (confirmed by the rerun above).

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
- **The gravity-driven shift is therefore about 24 m (+13.3 − (−10.7)), not the 13 m that the gravity-on run alone would suggest.** The crest well location works against segregation.
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
- Injecting 60,000 t CO2 after producing 15,000 t CH4 raised the pressure above its initial value (to 158.6 bar vs about 152 bar initially, gravity on). In a closed compartment, that over-pressure is what would load a bounding fault in stage B2c.
- Because CO2-rich gas collects at the gas–water contact, dissolution into brine (excluded here) would likely be enhanced there. This is one reason the Søreide–Whitson work (Path A) still matters.
