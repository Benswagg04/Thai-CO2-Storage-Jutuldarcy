# Step 3b: Verification report

Purpose: test whether the Step 3 conclusions survive numerical and physical scrutiny before moving to the compositional depleted-gas-field model (Step 4). This step was prompted by a critical review of Step 3.

All numbers come from `results/03b_tables.md`, produced by `scripts/03b_report.jl` from the runs in `scripts/03b_verify.jl`.

## Summary

| Check | Result | Verdict |
|---|---|---|
| CO2 conservation (independent residual) | Residual ≤ 3.2 t (8e-5 of injected) at default tolerances on the base grid; ≤ 0.34 t with 100× tighter mass-balance tolerance. Up to 14 t on the finest grid. | **Passed.** Numerical, tolerance-controlled, small. |
| Boundary transport (integrated face flux) | Measured directly: ≤ 0.7 t through the top (dissolved only) in the cases tested; 0 through the sides. | **Passed.** The earlier ±2.5 t was mostly conservation error, not boundary flow. |
| Integration method (hook vs stored ministep states) | Agree to ≈ 0.01 t. | **Passed** |
| Time-step sensitivity (max Δt 90 → 30 days) | Leakage changes ≤ 0.8%. | **Passed** |
| Grid sensitivity (4 grids at Pe 1.2 bar; 2 at 0.8 bar) | Leakage magnitudes change by up to ~20–27%, in different directions in different regimes. Regime behaviour is the same on all grids. | **Partially passed.** Qualitative results are robust; magnitudes are not converged. |
| Injection-rate hypothesis (fixed total mass) | Confirmed, and more strongly than expected (see §5). | **Passed** |
| Phase-resolved flux across internal interfaces | Not available through a public JutulDarcy function. Net region balances used instead. | **Not done (known gap)** |

## 1. Definitions

- **Inventory:** CO2 held in a region at a given time, split into free (CO2-rich phase) and dissolved (in brine).
- **Boundary transport:** cumulative CO2 that crossed the model boundaries. It is integrated from JutulDarcy's `compute_bc_mass_fluxes` after every accepted ministep (flux × ministep length, consistent with the solver's backward-Euler scheme).
- **Conservation residual:** injected − (reservoir inventory + well inventory + boundary transport). Each term is measured independently, so this is a genuine check, unlike the earlier "injected − in model" figure, which was circular.
- **Crossed the seal:** cumulative net CO2 transport into the region above the seal (the overburden sand plus the fault cells above the footwall caprock top). It is computed from that region's balance: its inventory plus what has left it through the boundaries. It does not depend on the conservation residual.
- **Out of the storage formation:** injected − storage-sand inventory − boundary transport from storage cells. This includes CO2 sitting inside the fault and caprock, and it absorbs the conservation residual.

## 2. Conservation and boundary transport

| Case | Tolerances | Residual yr 25 / 100 (t) | Out via top, free / dissolved, yr 100 (t) | Crossed seal yr 100 (t) |
|---|---|---|---|---|
| Pe 5.5 bar, k 1/10 mD | default | −0.57 / −3.17 | 0 / 0.686 | 57.3 |
| Pe 5.5 bar, k 1/10 mD | tight | −0.02 / −0.34 | 0 / 0.678 | 57.2 |
| Pe 1.2 bar, k 0.1/1 mD | default | −0.03 / −1.30 | 0 / 0 | 162.5 |
| Pe 1.2 bar, k 0.1/1 mD | tight | −0.00 / −0.32 | 0 / 0 | 162.3 |

Default tolerances: tol_cnv = 1e-3, tol_mb = 1e-7, tol_mb_well = 1e-3. Tight: tol_cnv = 1e-4, tol_mb = 1e-9, tol_mb_well = 1e-9.

- The residual is negative: the model gains a little CO2. It grows after injection stops and with the number of ministeps (up to 14.4 t on the 5 × 1 m grid at Pe 0.8 bar, about 0.4% of the CO2 that crossed the seal there).
- Tightening the mass-balance tolerance 100× reduces the residual about 9–10×. That is not proportional, but it clearly identifies the residual as solver error.
- The physical quantities change by < 0.5% between tolerance settings.
- Boundary transport is small and entirely dissolved CO2. The open-top boundary flux in JutulDarcy 0.3.10 is driven by the brine-pressure difference only, without gravity across the half cell or capillary pressure (source comment `# TODO: Capillary pressure` in `multicomponent/bc.jl`). Buoyant free CO2 therefore leaves only if the top cell is overpressured. This affects where CO2 above the seal ends up, but not the amount that crosses the seal.

## 3. Fault invasion vs crossing the seal (Pe 1.2 bar, k 0.1/1 mD)

| | Year 25 | Year 100 |
|---|---|---|
| Free CO2 in fault | 127.8 t | 103.2 t |
| Dissolved CO2 in fault | 105.2 t | 118.6 t |
| Free CO2 in overburden | 3.9 t | 0 |
| Dissolved CO2 in overburden | 121.5 t | 139.6 t |
| Cumulative net crossing of the seal | 145.0 t | 162.5 t (peak 168.6 t at year 33) |

- **Correction to Step 3:** fault invasion stops when injection stops, but crossing of the seal does not.
  - CO2 already in the fault keeps moving up for about 8 years after injection; crossing peaks at year 33.
  - It then declines slowly, as net brine back-flow carries some dissolved CO2 downward.
  - The same pattern appears on all four grids (peaks at years 33–34).
- Most CO2 that crosses the seal does so in solution. Free CO2 that briefly reached the overburden (3.9 t at year 25) had dissolved by year 100.
- A few tonnes of dissolved CO2 in the overburden is a model result, not evidence of environmental harm. Its significance depends on the receiving formation and on what is being assessed.

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

Step 2c, for comparison: Pe 0.55 bar, k 1/10 mD; leakage out of storage 50.2% → 55.0% of injected (+9.6%) on refinement.

Interpretation:
- **Time step:** no meaningful sensitivity.
- **Horizontal resolution matters most near the threshold.** Halving dx changes Pe 1.2 leakage by about −20%. Halving dz changes it by only +2–9%.
  - A plausible reason: the entry drive is evaluated in the sand cell next to the fault. Its pressure and saturation depend on how finely the region next to the fault is resolved.
- **The direction of the change depends on the regime.** Refinement reduces leakage near the threshold (Pe 1.2) and increases it when leakage is persistent (Pe 0.8, Pe 0.55). A coarse grid is therefore neither systematically conservative nor optimistic.
- **Robust across all grids:** whether the fault is invaded, when, the post-injection behaviour (invasion stops, seal crossing peaks about 8 years later), and the free/dissolved character.
- **Not converged:** leakage magnitudes. Treat Step 3 percentages as carrying roughly ±25% numerical uncertainty.
  - A formal convergence study would need the next level (2.5 × 0.5 m, ≈ 16× the base cells). On the 5 × 1 m grid one run already took 46–109 minutes, and more ministeps were cut than accepted (1,750 cut / 1,662 accepted). That is impractical on this laptop.

## 5. Injection rate at fixed total mass (37,500 t; Pe 1.4 bar; k 0.1/1 mD)

| Rate × duration | Injected (simulator) | Peak entry drive | Peak Δp | Out of storage: end of inj. / +25 yr / yr 100 | Crossed seal: end / +25 / 100 | Max free CO2 in fault | Max free CO2 in overburden |
|---|---|---|---|---|---|---|---|
| 750 t/yr × 50 yr | 37,500.0 t | 1.42 bar | 1.04 bar | 16.5 / 16.9 / 15.7 t | 2.4 / 2.4 / 2.4 t | **0** | 0 |
| 1,500 t/yr × 25 yr | 37,500.0 t | 1.65 bar | 1.46 bar | 112.7 / 113.1 / 112.0 t | 16.4 / 17.6 / 17.8 t | 39.6 t | 0 |
| 2,500 t/yr × 15 yr | 37,500.0 t | 1.92 bar | 1.93 bar | 391.5 / 395.1 / 393.0 t | 158.4 / 186.0 / 179.8 t | 137.1 t | **12.6 t** |

- The injected mass is the simulator's well output. The well has no pressure limit, so the prescribed rate was always achieved.
- **The hypothesis is confirmed.** Doubling the rate (750 → 1,500) increases CO2 out of storage about 7×. At 2,500 t/yr it is about 24×.
- At 750 t/yr, the peak entry drive (1.42 bar) barely exceeds Pe (1.4 bar). No free CO2 entered the fault; the 16 t out of storage is dissolved transport.
- Only the highest rate moves free CO2 past the seal.
- The post-injection entry drive returns to 1.08–1.10 bar in all three cases. So, as anticipated in the review, plume shapes differ, but the buoyancy loading on the fault after injection is essentially the same. The difference in outcome comes from injection overpressure.
- **Implication:** for a fault whose entry pressure lies within the injection window, the injection rate (not only the fault's properties) determines whether its seal is breached.

## 6. What the Step 3 conclusions can and cannot claim

Supported:
- Three approximate regimes controlled by fault entry pressure relative to the available drive (buoyancy + injection overpressure). The boundaries are scenario-specific.
- Entry pressure controls whether free CO2 enters the fault; permeability controls how much moves once it does.
- Within the injection window, fault invasion is injection-triggered and controlled by rate. Seal crossing continues for several years afterwards, mostly as dissolved CO2.
- A capillary seal stops free CO2 but not dissolved CO2 carried by brine.

Not supported:
- Leakage percentages to better than about ±25% (grid sensitivity).
- Trapping fractions, and the split between "in overburden" and "left through the top" (grid- and boundary-treatment-dependent).
- Any statement about fault reactivation or slip. This is a flow-only model, with no geomechanics.
- Transfer of the thresholds to a depleted gas field. Gas–gas mixing has no capillary barrier, and CO2–brine interfacial tension differs from CH4–brine.

## 7. Known gaps
- No phase-resolved flux across internal interfaces. JutulDarcy's interior flux functions are internal and need solver objects. Net transport is computed from region balances instead, and the free/dissolved split is reported for inventories only.
- No molecular diffusion in these runs. JutulDarcy's compositional flux supports it when a `Diffusivities` field is present; not yet used.
- The fault's capillary pressure is one scalar curve per cell, so there is no directional entry pressure.
- Prescribed fault properties, not derived from geology (e.g. clay content).
- No formal grid convergence.

## 8. Environment and reproduction
- Julia 1.12.6; JutulDarcy 0.3.10; Jutul 0.4.30; MultiComponentFlash 1.1.19; CairoMakie 0.15.15 (exact versions in `Manifest.toml`).
- Laptop: Windows 11, 8 GB RAM.
- Settings: 1-year report steps; max internal time step 90 days unless stated; tolerances as in §2.

```
julia --project=. scripts/03b_verify.jl flux
julia --project=. scripts/03b_verify.jl hookcheck
julia --project=. scripts/03b_verify.jl rate
julia --project=. scripts/03b_verify.jl dt
julia --project=. scripts/03b_verify.jl grid
julia --project=. scripts/03b_report.jl
```
