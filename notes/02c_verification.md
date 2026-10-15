# Step 2c: Checks before the fault-property sweep

Before running the Step 3 sweep I went back over Steps 2 and 2b and checked the parts I was least sure about. Steps 2 and 2b used a first configuration (no padding, closed top, Pc capped at 10 Pe); their figures are kept in `figures/superseded/`. The current scripts (`scripts/02_run_cases.jl`, then `scripts/make_figures.jl 02`) produce the figures described here.

## 1. Geometry in metres
Fault position, throw, layer boundaries and the well are now given in metres instead of grid indices, so a finer grid keeps the same geology (the fault stays 10 m wide). The rewritten code reproduces Step 1 (juxtaposition, pore volumes, depths) and the Step 2 baseline exactly (to 0.1 t at year 100; mass-balance error ~1e-10 t).

## 2. Capillary-pressure curve range
The Brooks–Corey curve was capped at Se = 0.01 (10 × Pe), but the simulation reached S_CO2 = 0.795, beyond the cap (0.792), so the top of the plume sat on the flat, capped part of the sand curve (0.5 bar). I moved the cap to Se = 1e-3 (≈ 32 × Pe) and made the table denser at high saturation. In the sand-rich case, with nothing else changed, every region moves by at most 0.4 t. The cap did not affect the conclusions, because leakage is set by the CO2 pressure itself rather than by the shape of the sand curve.

## 3. Leakage accounting
The regions (footwall and hanging-wall storage sand, fault zone, caprock, overburden, underburden) cover every cell exactly once, so their sum is the CO2 inside the model, and CO2 that left the model is injected minus inside. Free and dissolved CO2 in each region come from each cell's phase masses: free ∝ S_g ρ_g Y_CO2, dissolved ∝ S_l ρ_l X_CO2. "Outside the storage formation" means everything except the two storage sands, plus what left the model.

There is no molecular diffusion, so dissolved CO2 moves only with flowing brine, and "zero CO2 in the caprock" with capillary pressure is partly true by construction. Diffusion over 100 years would only reach a few metres anyway (√(D t) ≈ 2 m for D = 1e-9 m²/s).

## 4. Boundaries

| Boundary | Before | Now | Meaning |
|---|---|---|---|
| Left/right | Constant pressure at x = 0 and 1,200 m | Constant pressure 2.54 km beyond each side of the core (7 padding columns, 20 → 1,280 m wide) | The formations continue sideways; pressure dissipates far from the fault |
| Top | No flow (a lid) | Constant hydrostatic pressure, brine inflow (z-direction) | Connection to shallower aquifers. CO2 that crosses it has left the modelled interval and counts as leaked |
| Bottom | No flow | No flow | Underburden shale (Pe 12 bar, 1 µD) |

Extending the model upwards would only move the artificial lid higher; an open top gives a clear measure instead, since anything leaving through it has escaped the modelled interval. With the open top the model is still in equilibrium, and I track CO2 in the side columns to make sure it only leaves through the top.

## 5. Pressure at the fault
Instead of estimating buoyancy from the column height, I now compute the leakage criterion directly: p_CO2 in the sand next to the fault minus p_brine in the fault cell at the same level, compared with the fault entry pressure. This includes buoyancy and injection overpressure together.

## 6. Grid refinement
I repeated the sand-rich case on a 5 m × 1 m grid (4 times as many cells). The fault stays 10 m wide (2 columns), and the well moves by 2.5 m to the nearest cell centre, which is negligible.

## Results

### Checks

| Check | Result |
|---|---|
| Rewritten code, first configuration | Identical to the Step 2 baseline at year 100 (all regions to 0.1 t) |
| Pc cap 10 Pe → 32 Pe (sand-rich fault) | All regions change by ≤ 0.4 t |
| Equilibrium with open top and padding | 0.0000 bar drift in all cases |
| CO2 reaching the side columns | ≤ 0.003 t in all cases, so CO2 only leaves through the top |
| Mass balance | "Injected = in model + left model" is not a real test, because "left model" is defined as the difference. The sealed case shows the real error: 0.6–2.2 t "left model" (≤ 0.006% of injected), changing sign, with no CO2 near any boundary. That is solver-tolerance noise, not leakage. |

### Main cases (2.54 km padding, open top), year 100, % of 37,500 t injected

| Case | Fault | Caprock | Overburden | Left via top | Outside storage formation |
|---|---|---|---|---|---|
| No Pc, clay-rich fault | 2.4 | 3.6 | 2.6 | 4.3 | **12.9%** |
| Pc, clay-rich fault (Pe 5.5 bar) | 0 (1 t dissolved) | 0 (2 t dissolved) | 0 | noise | **~0%** |
| Pc, sand-rich fault (Pe 0.55 bar) | 1.4 | 0 | 21.4 | 27.3 | **50.2%** |

Without capillary pressure, the overburden CO2 at year 100 is 88% dissolved (112 t free, 848 t dissolved). With the sand-rich fault, the overburden holds 3,503 t free and 4,535 t dissolved. Free CO2 that reaches the overburden keeps rising and leaves through the top; the dissolved part stays behind.

### Pressure at the fault (Pc, clay-rich fault)
The entry drive is 0.88 bar at year 5, 1.74 bar at year 25 and 1.09 bar at year 100, against an entry pressure of 5.48 bar, so injection adds about 0.65 bar to the ~1.1 bar of buoyancy at the end of injection. With the first configuration's nearby boundaries the injection added only ~0.1 bar; those boundaries were relieving pressure artificially.

The brine pressure inside the fault rises by at most 0.2 bar. What adds to the entry drive is the pressure difference across the low-permeability fault, because the fault's pressure lags behind the sand's. A leaky fault relieves pressure: the maximum build-up is 0.74 bar with the sand-rich fault against 1.51 bar with the clay-rich one.

After shut-in the "maximum pressure rise anywhere" stays at 0.5–1 bar. Part of that is the brine reference pressure inside the buoyant CO2 column, not residual overpressure, so it should not be read as overpressure alone.

### Grid refinement (sand-rich fault; 10 × 2 m vs 5 × 1 m, 8,040 vs 30,480 cells; 1.6 vs 29.6 min)

| Quantity | Coarse | Fine |
|---|---|---|
| Leakage onset | ~year 5 | ~year 5 |
| Entry drive, year 25 / 100 (bar) | 1.01 / 0.64 | 0.96 / 0.64 |
| Outside storage formation, year 25 | 31.3% | 36.1% |
| Outside storage formation, year 100 | 50.2% | 55.0% |
| Retained in overburden, year 100 | 8,038 t | 3,240 t |
| Left via top, year 100 | 10,248 t | 16,877 t |
| Residual / dissolved trapping (whole model) | 10,056 / 10,484 t | 7,019 / 7,186 t |

The threshold behaviour (onset time, entry drive) and the total leaving the storage formation hold up on the finer grid, to within about 5 percentage points. Where the leaked CO2 ends up, and how much is trapped, do not: the coarse grid smears CO2 over more brine (numerical dispersion), which gives about 40% more residual and dissolved trapping and more CO2 retained in the overburden. So for Step 3 I use the coarse grid and report leakage out of the storage formation and the threshold behaviour, but not trapping fractions or overburden retention.
