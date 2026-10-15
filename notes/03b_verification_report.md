# Step 3b: Checking Step 3

Before moving to the compositional depleted-gas model, I wanted to know how far the Step 3 conclusions could be trusted, so I went back and tested the weakest points: whether CO2 is really conserved, how much leaves through the boundaries, and how sensitive the results are to time step, grid and injection rate. All numbers here are in `results/03b_tables.md`, produced by `scripts/03b_report.jl` from the runs in `scripts/03b_verify.jl`.

## Summary

| Check | Result | Verdict |
|---|---|---|
| CO2 conservation (independent residual) | ≤ 3.2 t (8e-5 of injected) at default tolerances on the base grid; ≤ 0.34 t with a 100× tighter mass-balance tolerance; up to 14 t on the finest grid | Passed: small, and controlled by the solver tolerance |
| Boundary transport (integrated face flux) | ≤ 0.7 t through the top (dissolved only) in the cases tested; nothing through the sides | Passed. The earlier ±2.5 t was mostly conservation error, not flow through the boundary |
| Integration method (hook vs stored ministep states) | Agree to about 0.01 t | Passed |
| Time step (max Δt 90 → 30 days) | Leakage changes by ≤ 0.8% | Passed |
| Grid (4 grids at Pe 1.2 bar, 2 at 0.8 bar) | Leakage changes by up to ~20–27%, in different directions in different regimes; the regimes are the same on all grids | Partly: the qualitative results hold, the magnitudes are not converged |
| Injection rate at fixed total mass | Rate matters even more than I expected (§5) | Passed |
| Phase-resolved flux across internal interfaces | Not available through a public JutulDarcy function; I used net region balances instead | Not done |

## 1. Definitions
- Inventory: the CO2 held in a region at a given time, split into free (CO2-rich phase) and dissolved (in brine).
- Boundary transport: the cumulative CO2 that crossed the model boundaries, integrated from JutulDarcy's `compute_bc_mass_fluxes` after every accepted ministep (flux × ministep length, consistent with the solver's backward-Euler scheme).
- Conservation residual: injected − (reservoir inventory + well inventory + boundary transport). Each term is measured separately, so unlike the "injected − in model" figure I used before, this is a real check.
- Crossed the seal: the cumulative net CO2 moved into the region above the seal (the overburden sand plus the fault cells above the top of the footwall caprock), computed from that region's inventory plus what has left it through the boundaries. It does not depend on the conservation residual.
- Out of the storage formation: injected − storage-sand inventory − boundary transport from storage cells. It includes CO2 sitting in the fault and caprock, and it absorbs the conservation residual.

## 2. Conservation and boundary transport

| Case | Tolerances | Residual yr 25 / 100 (t) | Out via top, free / dissolved, yr 100 (t) | Crossed seal yr 100 (t) |
|---|---|---|---|---|
| Pe 5.5 bar, k 1/10 mD | default | −0.57 / −3.17 | 0 / 0.686 | 57.3 |
| Pe 5.5 bar, k 1/10 mD | tight | −0.02 / −0.34 | 0 / 0.678 | 57.2 |
| Pe 1.2 bar, k 0.1/1 mD | default | −0.03 / −1.30 | 0 / 0 | 162.5 |
| Pe 1.2 bar, k 0.1/1 mD | tight | −0.00 / −0.32 | 0 / 0 | 162.3 |

Default tolerances: tol_cnv = 1e-3, tol_mb = 1e-7, tol_mb_well = 1e-3. Tight: tol_cnv = 1e-4, tol_mb = 1e-9, tol_mb_well = 1e-9.

The residual is negative, so the model gains a little CO2. It grows after injection stops and with the number of ministeps, reaching 14.4 t on the 5 × 1 m grid at Pe 0.8 bar (about 0.4% of the CO2 that crossed the seal there). Tightening the mass-balance tolerance 100 times reduces it about 9–10 times. That is not proportional, but it clearly marks the residual as solver error, and the physical quantities change by less than 0.5% between the two settings.

What crosses the boundaries is small and entirely dissolved. In JutulDarcy 0.3.10 the flux through the open top depends only on the brine-pressure difference, without gravity across the half cell or capillary pressure (the source has a `# TODO: Capillary pressure` comment in `multicomponent/bc.jl`), so buoyant free CO2 only leaves if the top cell is overpressured. This affects where CO2 above the seal ends up, but not how much crosses the seal.

## 3. Fault invasion vs crossing the seal (Pe 1.2 bar, k 0.1/1 mD)

| | Year 25 | Year 100 |
|---|---|---|
| Free CO2 in fault | 127.8 t | 103.2 t |
| Dissolved CO2 in fault | 105.2 t | 118.6 t |
| Free CO2 in overburden | 3.9 t | 0 |
| Dissolved CO2 in overburden | 121.5 t | 139.6 t |
| Cumulative net crossing of the seal | 145.0 t | 162.5 t (peak 168.6 t at year 33) |

This corrects what I wrote in Step 3: invasion of the fault stops when injection stops, but crossing of the seal does not. CO2 already in the fault keeps moving up for about 8 years after injection, so crossing peaks at year 33, and then declines slowly as brine flowing back down carries some dissolved CO2 with it. The same pattern appears on all four grids (peaks at years 33–34). Most of the CO2 that crosses the seal does so in solution; the free CO2 that briefly reached the overburden (3.9 t at year 25) had dissolved by year 100. A few tonnes of dissolved CO2 in the overburden is a model result, not evidence of harm; its significance depends on the receiving formation and on what is being assessed.

## 4. Time-step and grid sensitivity

| Case (k 0.1/1 mD) | Cells | Out of storage, yr 100 | Δ vs 10×2 | Crossed seal, yr 100 | Δ vs 10×2 | Free CO2 in fault, yr 100 | Wall time |
|---|---|---|---|---|---|---|---|
| Pe 1.2, 10×2 m, max Δt 90 d | 8,040 | 361.8 t | – | 162.5 t | – | 103.2 t | 2.9 min |
| Pe 1.2, 10×2 m, max Δt 30 d | 8,040 | 363.0 t | +0.3% | 161.2 t | −0.8% | 103.1 t | 2.5 min |
| Pe 1.2, 10×1 m | 16,080 | 369.0 t | +2.0% | 176.8 t | +8.8% | 105.9 t | 11.1 min |
| Pe 1.2, 5×2 m | 15,240 | 284.0 t | −21.5% | 118.0 t | −27.4% | 88.0 t | 4.7 min |
| Pe 1.2, 5×1 m | 30,480 | 297.8 t | −17.7% | 134.0 t | −17.5% | 92.1 t | 46.5 min |
| Pe 0.8, 10×2 m | 8,040 | 2,988.5 t | – | 2,714 t | – | 282.8 t | 2.1 min |
| Pe 0.8, 5×1 m | 30,480 | 3,600.3 t | +20.5% | 3,350 t | +23.4% | 265.9 t | 109.4 min |

For comparison, Step 2c (Pe 0.55 bar, k 1/10 mD) went from 50.2% to 55.0% of injected leaving storage on the finer grid (+9.6%).

The time step makes no meaningful difference. Horizontal resolution matters most near the threshold: halving dx changes the Pe 1.2 leakage by about −20%, while halving dz changes it by only +2–9%. A plausible reason is that the entry drive is evaluated in the sand cell next to the fault, whose pressure and saturation depend on how finely that region is resolved. (Step 3c later tested the fault itself.)

The direction of the change depends on the regime: refinement reduces leakage near the threshold (Pe 1.2) and increases it when leakage is persistent (Pe 0.8 and 0.55), so a coarse grid is neither systematically conservative nor optimistic. Across all grids, whether and when the fault is invaded, the post-injection behaviour (invasion stops, crossing of the seal peaks about 8 years later) and the free/dissolved character are the same. The magnitudes are not converged, and I treat the Step 3 percentages as carrying about ±25% numerical uncertainty. A formal convergence study would need the next level (2.5 × 0.5 m, about 16 times the base cells). On the 5 × 1 m grid a single run already took 46–109 minutes, with more ministeps cut than accepted (1,750 cut, 1,662 accepted), which is impractical on this laptop.

## 5. Injection rate at fixed total mass (37,500 t; Pe 1.4 bar; k 0.1/1 mD)

| Rate × duration | Injected (simulator) | Peak entry drive | Peak Δp | Out of storage: end of inj. / +25 yr / yr 100 | Crossed seal: end / +25 / 100 | Max free CO2 in fault | Max free CO2 in overburden |
|---|---|---|---|---|---|---|---|
| 750 t/yr × 50 yr | 37,500.0 t | 1.42 bar | 1.04 bar | 16.5 / 16.9 / 15.7 t | 2.4 / 2.4 / 2.4 t | 0 | 0 |
| 1,500 t/yr × 25 yr | 37,500.0 t | 1.65 bar | 1.46 bar | 112.7 / 113.1 / 112.0 t | 16.4 / 17.6 / 17.8 t | 39.6 t | 0 |
| 2,500 t/yr × 15 yr | 37,500.0 t | 1.92 bar | 1.93 bar | 391.5 / 395.1 / 393.0 t | 158.4 / 186.0 / 179.8 t | 137.1 t | 12.6 t |

The injected mass is the simulator's own well output; the well has no pressure limit, so the rate was always achieved. Doubling the rate from 750 to 1,500 t/yr increases the CO2 out of storage about 7 times, and at 2,500 t/yr about 24 times. At 750 t/yr the peak entry drive (1.42 bar) barely exceeds the entry pressure (1.4 bar): no free CO2 entered the fault, and the 16 t out of storage is dissolved transport. Only the highest rate moves free CO2 past the seal.

The entry drive after injection returns to 1.08–1.10 bar in all three cases. So the plume shapes differ, but the buoyancy load on the fault after injection is essentially the same; the difference in outcome comes from the injection overpressure. For a fault whose entry pressure lies inside the injection window, the injection rate decides whether its seal is breached, not just the fault's properties.

## 6. What the Step 3 conclusions can and cannot claim

They support:
- three approximate regimes, controlled by the fault entry pressure relative to the available drive (buoyancy plus injection overpressure), with boundaries specific to this scenario;
- entry pressure controlling whether free CO2 enters the fault, and permeability how much moves once it does;
- invasion within the injection window being triggered by injection and controlled by rate, with crossing of the seal continuing for several years afterwards, mostly as dissolved CO2;
- a capillary seal stopping free CO2 but not dissolved CO2 carried by brine.

They do not support:
- leakage percentages to better than about ±25% (grid sensitivity);
- trapping fractions, or the split between "in the overburden" and "left through the top" (both depend on the grid and the boundary treatment);
- anything about fault reactivation or slip, since there is no geomechanics;
- transferring the thresholds to a depleted gas field, where gas–gas mixing has no capillary barrier and the CO2–brine interfacial tension differs from CH4–brine.

## 7. Known gaps
- No phase-resolved flux across internal interfaces. JutulDarcy's interior flux functions are internal and need solver objects, so I computed net transport from region balances and report the free/dissolved split for inventories only.
- No molecular diffusion. JutulDarcy's compositional flux supports it when a `Diffusivities` field is present, but I have not used it yet.
- One capillary-pressure curve per cell, so no directional entry pressure in the fault.
- Fault properties chosen rather than derived from geology.
- No formal grid convergence.

## 8. Environment and reproduction
Julia 1.12.6, JutulDarcy 0.3.10, Jutul 0.4.30, MultiComponentFlash 1.1.19, CairoMakie 0.15.15 (exact versions in `Manifest.toml`), on a Windows 11 laptop with 8 GB RAM. Report steps are 1 year, the maximum internal time step is 90 days unless stated, and the tolerances are as in §2.

```
julia --project=. scripts/03b_verify.jl flux
julia --project=. scripts/03b_verify.jl hookcheck
julia --project=. scripts/03b_verify.jl rate
julia --project=. scripts/03b_verify.jl dt
julia --project=. scripts/03b_verify.jl grid
julia --project=. scripts/03b_report.jl
```
