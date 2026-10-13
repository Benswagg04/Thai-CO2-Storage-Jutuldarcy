# Step 2c: Verification before the fault-property sweep

Prompted by a critical review of Steps 2 and 2b. Each item: what was checked, why, and the result.

Steps 2 and 2b used a first configuration (no padding, closed top, Pc cap at 10 Pe). Their figures are kept in `figures/superseded/` for the record. The current scripts (`scripts/02_run_cases.jl`, then `scripts/make_figures.jl 02`) produce the figures described here.

## 1. Geometry defined in metres (prerequisite for grid refinement)
- Fault position, throw, layer boundaries and well are now specified in metres, not grid indices, so a finer grid keeps the same geology (the fault stays 10 m wide).
- **Check:** the rewritten code reproduces Step 1 (juxtaposition, pore volumes, depths) and the Step 2 baseline exactly (to 0.1 t at year 100; mass-balance error ~1e-10 t).

## 2. Capillary-pressure curve range
- **Issue found:** the Brooks-Corey curve was capped at Se = 0.01 (10 × Pe). The simulation reached S_CO2 = 0.795, beyond the cap (0.792), so the top of the plume sat on the flat, capped part of the sand curve (0.5 bar).
- **Fix:** cap at Se = 1e-3 (≈ 32 × Pe), with a denser table at high saturation.
- **Check:** sand-rich-fault case, nothing else changed. Every region changes by ≤ 0.4 t. The cap did not affect the conclusions, because leakage is decided by the CO2 pressure itself, not by the shape of the sand curve.

## 3. Leakage accounting
- Regions partition every cell exactly once: footwall storage sand, hanging-wall storage sand, fault zone, caprock, overburden, underburden. Their sum is the CO2 inside the model.
- CO2 that left the model = injected − inside model (mass balance).
- Free vs dissolved in each region is split with each cell's phase masses: free ∝ S_g ρ_g Y_CO2, dissolved ∝ S_l ρ_l X_CO2.
- "Outside the storage formation" = everything except the footwall and hanging-wall storage sand, plus CO2 that left the model.
- **No molecular diffusion is modelled.** Dissolved CO2 moves only with flowing brine. So "zero CO2 in the caprock" with capillary pressure is partly true by construction. Diffusion over 100 years would reach only a few metres (√(D t) ≈ 2 m for D = 1e-9 m²/s).

## 4. Boundaries and their physical meaning
| Boundary | Before | Now | Meaning |
|---|---|---|---|
| Left/right | Constant pressure at x = 0 and 1,200 m | Constant pressure 2.54 km beyond each side of the core (7 padding columns, 20 → 1,280 m wide) | Formations continue laterally; pressure dissipates far from the fault |
| Top | No flow (a lid) | Constant hydrostatic pressure, brine inflow (z-direction) | Connection to shallower aquifers. CO2 that crosses it has left the modelled interval and is counted as leaked |
| Bottom | No flow | No flow | Underburden shale (Pe 12 bar, 1 µD) |

- Why not just extend the model upward? A thicker overburden would only move the artificial lid higher. Making the top open defines a clear metric: anything leaving through it has escaped the modelled interval.
- **Check:** with the open top, the equilibrium drift is still 0.0000 bar. CO2 reaching the side columns is tracked to confirm that exits are through the top only.

## 5. Pressure at the fault
- The leakage criterion is evaluated directly as p_CO2(sand next to fault) − p_brine(fault cell at the same level), compared with the fault entry pressure. This includes buoyancy and injection overpressure together, rather than estimating buoyancy from column height alone.

## 6. Grid refinement
- Sand-rich-fault case repeated on a 5 m × 1 m grid (4× as many cells). The fault stays 10 m wide (2 columns).
- The well moves 2.5 m (cell-centre position on the finer grid). Negligible.

## Results

### Checks
| Check | Result |
|---|---|
| Regression (rewritten code, first configuration) | Identical to the Step 2 baseline at year 100 (all regions to 0.1 t) |
| Pc cap 10 Pe → 32 Pe (sand-rich fault, first configuration) | All regions change by ≤ 0.4 t |
| Equilibrium with open top and padding | 0.0000 bar drift in all cases |
| CO2 reaching side-boundary columns | ≤ 0.003 t in all cases, so exits are through the top only |
| Mass-balance closure | Not tested by "injected = in model + left model": that identity is circular, because "left model" is defined as the difference. The real closure error shows in the sealed case: 0.6–2.2 t "left model" (≤ 0.006% of injected), fluctuating in sign, with no CO2 near any boundary. This is solver-tolerance noise, not leakage. |

### Main cases (2.54 km padding, open top), year 100, % of 37,500 t injected
| Case | Fault | Caprock | Overburden | Left via top | Outside storage formation |
|---|---|---|---|---|---|
| No Pc, clay-rich fault | 2.4 | 3.6 | 2.6 | 4.3 | **12.9%** |
| Pc, clay-rich fault (Pe 5.5 bar) | 0 (1 t dissolved) | 0 (2 t dissolved) | 0 | noise | **~0%** |
| Pc, sand-rich fault (Pe 0.55 bar) | 1.4 | 0 | 21.4 | 27.3 | **50.2%** |

Free vs dissolved, year 100:
- No Pc: the overburden CO2 is 88% dissolved (112 t free, 848 t dissolved).
- Sand-rich fault: the overburden holds 3,503 t free and 4,535 t dissolved.
- Free CO2 that reaches the overburden keeps rising and leaves through the top. The dissolved part stays behind.

### Pressure at the fault (Pc, clay-rich fault)
- Entry drive p_CO2(sand) − p_brine(fault): 0.88 bar at year 5, **1.74 bar at year 25**, 1.09 bar at year 100, against Pe = 5.48 bar.
- Injection therefore adds about 0.65 bar to the ~1.1 bar of buoyancy at the end of injection.
- With the first configuration's nearby boundaries, the injection contribution was only ~0.1 bar. Those boundaries were relieving pressure artificially.
- The brine pressure inside the fault rises by at most 0.2 bar. It is the pressure *difference across* the low-permeability fault that adds to the entry drive, because the fault's pressure lags behind the sand's.
- A leaky fault relieves pressure: maximum build-up is 0.74 bar with the sand-rich fault vs 1.51 bar with the clay-rich one.
- After shut-in, the "max pressure rise anywhere" stays at 0.5–1 bar. This partly reflects the brine reference pressure inside the buoyant CO2 column, not only residual overpressure. Do not read it as overpressure alone.

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

Interpretation:
- **Robust to grid resolution:** the threshold behaviour (onset time, entry drive) and the total leaving the storage formation (within ~5 percentage points).
- **Not robust:** where the leaked CO2 ends up, and how much is trapped. The coarse grid smears CO2 over more brine (numerical dispersion), giving ~40% more residual and dissolved trapping and more retention in the overburden.
- Conclusion for Step 3: use the coarse grid for the sweep, and report leakage out of the storage formation and threshold behaviour. Do **not** report trapping fractions or overburden retention as results.
