# Stage B2c: faulted depleted gas field (in progress)

Code: `src/gasfault.jl`, `scripts/04d_fault_gas.jl`. Fluid: immiscible brine + CO2–CH4 Peng–Robinson. **Dissolution of CO2 or CH4 in brine is not modelled**; any "dissolved" quantity is reported as *not modelled*, not as zero.

## Design (agreed after review)
- **Grid:** the exact Step 1–3 cross-section (8,040 cells, 10 m × 2 m, 10 m-wide fault zone with k_across 0.01 mD, k_along 0.1 mD), so the fault representation is identical to Steps 2–3. Other differences from Step 3 (fluid, history, schedule) are not controlled by this.
- **Outer boundaries closed.** The 2.54 km padding on each side acts as a finite aquifer.
- **Methane trap** in the footwall storage sand against the fault. Free-water level (FWL) at 1,534 m, giving a 13.9 m gas column at the fault.
- **Capillary pressure** (implementation verified in `notes/04c_pc_check.md`): Brooks–Corey per facies (λ = 2). Entry pressures are Leverett-scaled from the storage sand (0.05 bar) as in Steps 2–3; the fault value is set directly.
- **Fault entry pressure, methane–brine:** Pe,CH4 = 2.5 bar ≈ **2.0 ×** the column buoyancy at the fault (1.23 bar). This is a *consistency requirement* on the synthetic model (the initial trap must be capillary-stable), with a stated margin. It is not a measured fault property; for a real field, the methane-column argument is only a conditional lower bound (Naylor et al. 2011).
- **CO2-scaled case:** all entry pressures × 0.6 (Pe,fault 1.5 bar). The ratio Pe,CO2/Pe,CH4 ≈ (γ_CO2/γ_CH4)(cos θ_CO2/cos θ_CH4) is **provisional**. Naylor et al.'s reference values (25 vs 50 mN/m at 10 MPa, wettability factor 0.7–1.0) suggest 0.35–0.5, so 0.6 may be generous to the CO2 seal. Interfacial tension also varies with pressure, so a single ratio is a simplification. Composition-dependent entry pressure (e.g. from Liu, Li & Okuno 2016 data) is left for later.
- **Margin choice:** with margin 2, the CO2-scaled fault (1.2 × buoyancy) can still hold the original methane column. With a margin below 1/0.6 ≈ 1.67 it could not, and migration would be caused by the property change rather than by CO2.

## Initial state
JutulDarcy's `EquilibriumRegion` would put gas everywhere above the contact, including the hanging-wall sand opposite the fault, so the state is built by hand:
- brine pressure hydrostatic, integrated with the compressible brine density;
- in the trap, Pc = (ρw − ρg) g (z_FWL − z) with Sw from each cell's Brooks–Corey curve;
- every other cell water-filled.

| | Methane curves | CO2-scaled curves |
|---|---|---|
| Fault Pe / column buoyancy | 2.50 / 1.23 bar (2.03) | 1.50 / 1.23 bar (1.22) |
| Initial CH4 in place | 4,248 t | 4,530 t (thinner transition zone) |
| 10 yr without wells: CH4 into fault / beyond the seal | 0 / 0 t | 0 / 0 t |
| 10 yr without wells: max \|Δp\|, max \|ΔSw\| | 0.033 bar, 0.010 | 0.26 bar, 0.005 |
| Where the changes occur | 11 cells, all in the trap (pressure near the crest by the fault; Sw at the transition-zone base) | 4 cells, all in the trap |
| **Second 10 yr from the relaxed state** | 0.0028 bar, 0.0023 | 0.015 bar, 0.0004 |

- The hand-built profile is approximate (local densities instead of integrated phase pressures; the solver uses a tabulated curve where it is steep). A 10-year no-well relaxation is used as **numerical equilibration**, and the relaxed state is the starting point for the runs.
- The residual drift (≤ 0.015 bar, ≤ 0.002 in Sw per 10 years) is about 1/1000 of the pressure changes caused by production. It persists because transition-zone water near residual saturation is almost immobile, so it equilibrates very slowly. Recorded as a limitation.
- **The CO2-scaled curves hold the initial methane column** (no gas enters the fault in 20 years), so the comparison case is not invalid from the start.

## Planned runs
From one depleted state (baseline after 10 yr production + 1 yr shut-in), three branches of 30 years:
- **(a) baseline:** methane curves; CO2 injection 10 yr, then shut-in.
- **(b) switched + injection:** CO2-scaled curves; same injection.
- **(c) switched, no injection (control):** CO2-scaled curves, no wells. Measures migration caused by the property change alone.

Migration attributable to injection under the weakened seal ≈ (b) − (c).
