# Stage B2c: CO2 injection into a faulted, depleted methane field

Code: `src/gasfault.jl`, `scripts/04d_fault_gas.jl` (runs), `scripts/04d_report.jl` (tables and figures), `scripts/04d_conservation.jl` and `scripts/04d_tolerance.jl` (checks). Results: `results/04d_branches.csv`; figures `figures/04d_*.png`. Fluid: immiscible brine + CO2–CH4 Peng–Robinson. **Dissolution of CO2 or CH4 in brine is not modelled**; any "dissolved" quantity is *not modelled*, not zero.

## Question
After a methane field next to a fault has been produced, does re-injected CO2 cross the fault seal, and how much of that depends on the switch from methane–brine to CO2–brine capillary properties (as opposed to injection itself)?

## Design (agreed after review)
- **Grid:** the exact Step 1–3 cross-section (8,040 cells, 10 m × 2 m; 10 m-wide fault zone with k_across 0.01 mD, k_along 0.1 mD), so the fault representation is identical to Steps 2–3. A fault-resolution test (below) splits the fault into 5 and 10 columns.
- **Outer boundaries closed.** The 2.54 km padding on each side acts as a finite aquifer. The model top is also closed, so gas that leaks into the overburden pools under it (as in Step 2b); "beyond the seal" is a correct tally, but the gas's onward migration is not represented.
- **Methane trap** in the footwall storage sand against the fault. Free-water level (FWL) at 1,534 m, giving a 13.9 m gas column at the fault.
- **Capillary pressure** (implementation verified in `notes/04c_pc_check.md`): Brooks–Corey per facies (λ = 2). Entry pressures are Leverett-scaled from the storage sand (0.05 bar) as in Steps 2–3; the fault value is set directly.
- **Fault entry pressure, methane–brine:** Pe,CH4 = 2.5 bar ≈ **2.0 ×** the column buoyancy at the fault (1.23 bar). This is a *consistency requirement* on the synthetic model (the initial trap must be capillary-stable), with a stated margin. It is not a measured fault property; for a real field, the methane-column argument is only a conditional lower bound (Naylor et al. 2011).
- **CO2-scaled case:** all entry pressures × 0.6 (Pe,fault 1.5 bar). The ratio Pe,CO2/Pe,CH4 ≈ (γ_CO2/γ_CH4)(cos θ_CO2/cos θ_CH4) is **provisional**. Naylor et al.'s reference values (25 vs 50 mN/m at 10 MPa, wettability factor 0.7–1.0) suggest 0.35–0.5, so 0.6 may be generous to the CO2 seal. Interfacial tension also varies with pressure and composition, so a single ratio is a simplification.
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
| Where the changes occur | 11 cells, all in the trap | 4 cells, all in the trap |
| **Second 10 yr from the relaxed state** | 0.0028 bar, 0.0023 | 0.015 bar, 0.0004 |

The hand-built profile is approximate, so a 10-year no-well relaxation is used as **numerical equilibration**, and the relaxed state is the starting point. The residual drift (≤ 0.015 bar, ≤ 0.002 in Sw per 10 years) is about 1/1000 of the pressure changes caused by production; transition-zone water near residual saturation is almost immobile and equilibrates very slowly.

## Schedule
One history, then branches restarted from the **same depleted state** (year 11):
- **Depletion (methane curves):** 10 yr production at a 300 t/yr target (bottom-hole pressure floor 30 bar, not reached), then 1 yr shut-in. Producer near the crest (x = 665 m, perforated 44–50 m below the top surface).
- **(a)** methane curves; CO2 1,500 t/yr for 10 yr through the same well, then 20 yr shut-in.
- **(b)** CO2-scaled curves; same injection.
- **(c)** CO2-scaled curves; no injection (**control** for the property switch).
- **Extension, added after (a):** (a′) and (b′), the same injection with the injector's bottom-hole pressure capped at the original brine pressure at the top perforation (153.1 bar). Reason: 15,000 t over-pressured the crest to 168 bar, 15 bar above its original 153 bar, which an operator would not do in a depleted field.

## Verification
| Check | Result |
|---|---|
| **Component conservation**: (reservoir + well contents) − start − cumulative well flow, with the well flow integrated over every time step in a `post_ministep_hook` | Depletion: CH4 ≤ 0.008 t (3,095 t produced), water ≤ 0.02 t. Branches (a), (b), (c), (a′): CO2 ≤ 0.8 t of 11,000–15,000 t injected, CH4 ≤ 0.2 t. (b′): CO2 12.8 t (0.12%), CH4 1.4 t; this run predates the final well-flow bookkeeping, and the after-the-fact correction assumes all flow was injection, which is probably not exact for a pressure-capped well |
| **Solver tolerance** (see below) | Year-41 leakage changes by ≤ 1.2 t between `tol_mb` 1e-9 and 1e-7 (e.g. (a): 2,425.6 vs 2,424.8 t CO2) |
| **Fault resolution** (see below) | Leakage is grid-sensitive: −31% for (a) with 10 fault columns |
| **Control (c)** | 0.00 t of CO2 or CH4 enters the fault in 30 years |

**Mass-balance tolerance.** JutulDarcy's mass-balance convergence check is relative to the mass of the *whole model* (~1.27 × 10¹⁰ kg, almost all brine). With the default `tol_mb = 1e-7`, each time step may leave ~1.3 t of error in *any* component, including the ~1,000 t of CH4: the control lost 1.95 t of CH4 (0.17%) while re-equilibrating, with nothing leaving the model. `1e-10` brings this to 7 kg but stalls injection once gas reaches the fault (median step 0.14 day; 9.6 of 30 years in 40 min). `1e-9` completes branch (a) in 24–39 min (12 min at 1e-7) and is used for all B2c results. The physics is insensitive to this choice; conservation is not.

**Producer rate (JutulDarcy 0.3.10 behaviour).** The 300 t/yr production target is a standard-condition volume. The producer delivered 300 → 323 t/yr instead of a constant 300: the conversion weights *surface* densities by the *reservoir-condition* volume fractions in the well, so a trace of water in the well (~6 × 10⁻⁵ by volume) raises the mixture density and the produced mass. Mass balance is unaffected; the injector (fixed density) is exact. Total produced: 3,095 t.

## Results

### Depletion (years 0–11)
- 3,095 t of the 4,248 t of CH4 produced (73%).
- Crest pressure fell from 153.1 (original) to 90.4 bar (94.5 bar after shut-in), far less than for a closed gas tank, because **18,300 m³ of water flowed into the footwall storage sand** (38% of the initial gas pore volume of 48,500 m³). Sources: overburden on both sides 10,300 m³ (mostly through the fault's sand–sand window), hanging-wall storage sand 4,900, underburden 1,700, caprock 1,300 m³. The 0.001 mD caprock and underburden conduct measurable water over 3.7 km of section.
- The fault became a water conduit *into* the trap: the fault-zone brine pressure stayed up to 18 bar above the crest (entry drive −17.5 bar at year 10). No methane entered the fault.

### Branches (1 fault column, tol_mb 1e-9; year-41 values unless stated)
| | CO2 injected (t) | Peak crest p (bar) | Peak entry drive (bar) / fault Pe | First year drive > Pe | CO2 beyond seal, yr 21 / 41 (t) | % of injected | CH4 beyond seal (t) |
|---|---|---|---|---|---|---|---|
| (a) methane curves | 15,000 | 168.4 | 8.8 / 2.5 | 13 | 1,770 / 2,426 | 16.2% | 213 |
| (b) CO2-scaled | 15,000 | 168.3 | 7.7 / 1.5 | 13 | 1,986 / 2,784 | 18.6% | 201 |
| (c) control | 0 | 103.4 | 0.7 / 1.5 | never | 0 / 0 | – | 0 |
| (a′) capped | 10,966 | 153.1 | 8.8 / 2.5 | 13 | 1,338 / 1,581 | 14.4% | 163 |
| (b′) capped | 11,019 | 153.1 | 7.7 / 1.5 | 13 | 1,499 / 1,847 | 16.8% | 163 |

"Beyond the seal" = overburden sand plus fault cells above the footwall caprock top. The caprock itself is never invaded. CH4 beyond the seal is 14–18% of the 1,153 t left after depletion: injected CO2 pushes residual methane out first.

**Answer to the planned question.**
- **(b) − (c):** the control moves no gas at all, so all migration in (b) is caused by injection. The curve switch on its own causes none.
- **Effect of the switch under injection:** (b) − (a) = +358 t CO2 (+15%) at constant rate; (b′) − (a′) = +266 t (+17%) when capped, with the one-cell fault. **With a resolved fault the effect is larger: +676 t (+41%) with 10 fault columns** (see Fault resolution), because the spurious part of the pressure lag shrinks and the drive sits closer to the entry pressures, so their difference (2.5 vs 1.5 bar) matters more.
- **Pressure capping** reduces leakage (−35% for (a)) but does not prevent it.

### Mechanism: the fault's pore pressure lags the gas pressure
The entry drive is p_gas(crest) − p_brine(first fault cell). Split into p_brine(crest) − p_brine(fault) + Pc(crest):
- Pc(crest) is small, 0.3–0.9 bar, below both entry pressures throughout.
- The brine-pressure difference dominates: 4–8 bar during injection.
- In the gas-filled crest the water is at Sw = 0.2015–0.2036 (residual 0.20), relative permeability 3 × 10⁻⁶ to 2 × 10⁻⁵. It is essentially immobile, so its pressure rises with the injected gas instead of equalising with the fault.
- Inside the fault, brine can only drain slowly through 0.01 mD rock toward the hanging wall, which lags the trap by up to 16 bar (year 13: crest 125.3, hanging wall 114.9 bar).
- Gas therefore enters the fault at year 13, **while the crest is still ~28 bar below its original pressure (125.3 vs 153.1 bar)**. Depletion makes this worse: it drew the fault zone down to ~106 bar.
- After shut-in the pressures equalise; the drive falls below Pe 6–9 years after injection stops, and crossing stops.

### Fault resolution
| | 1 column (10 m) | 5 columns (2 m) | 10 columns (1 m) |
|---|---|---|---|
| (a) brine-pressure part of drive, year 15 (bar) | 7.69 | 4.92 | 4.78 |
| (a) CO2 beyond seal, year 41 (t) | 2,426 | 1,813 | 1,664 |
| (b) brine-pressure part of drive, year 15 (bar) | 6.92 | 3.89 | 3.66 |
| (b) CO2 beyond seal, year 41 (t) | 2,784 | 2,442 | 2,340 |

- I predicted that refining the fault would shrink the brine-pressure part of the drive to ~1/5 and ~1/10 and remove most of the leakage. **That was wrong:** it levels off at ~4.8 bar, because the immobile water in the gas zone, not the fault's internal gradient, sets most of the difference.
- The single 10 m fault cell **overestimates leakage**: by ~45% for (a) and ~19% for (b), relative to 10 columns. The 5 → 10 change is 8% for (a) and 4% for (b), so the converged values are probably ~1,600 t and ~2,300 t.
- **The curve-switch effect grows with resolution:** (b) − (a) = +358 t (1 column), +629 t (5), +676 t (10). The one-cell fault understates the importance of the CO2-vs-methane entry pressure.
- Cost: the 10-column runs took 50 min (a) and 187 min (b), against 24 and 13 min with one column; the narrow fault cells force small time steps while gas moves through them.
- **Consequence for Steps 2–3:** they used the same single-cell fault. Checked in Step 3c (`notes/03c_fault_resolution.md`): the bias goes both ways there, injection-triggered leakage overestimated (up to ~2.5×) but persistent leakage underestimated (~20%), because one wide cell also spreads the gas at low saturation and slows its flow up the fault.

## Predictions vs outcomes
| Prediction (made before the runs) | Outcome |
|---|---|
| ~3,000 t produced; crest falls to 50–60 bar (little aquifer support) | 3,095 t produced; **crest 90 bar**: 18,300 m³ of water inflow, from regions I had not counted |
| No CH4 enters the fault during depletion | Correct (drive −17.5 bar) |
| 15,000 t returns the trap to about its original pressure | **Wrong: 168 bar** (15 bar over the original 153 bar) |
| (a) holds unless injection pushes the drive above 2.5 bar | Drive reached 8.8 bar; **the fault failed**, via the pressure-lag mechanism, not the gas column |
| (c) shows essentially zero migration | Correct (0.00 t) |
| A one-cell fault exaggerates the drive; refinement removes most leakage | **Partly wrong**: leakage falls 31% (a) and 16% (b), not to zero; and the curve-switch effect grows instead of disappearing |

## Limitations
- One synthetic fault and one set of rock properties; 2D section; no rock compressibility; no relative-permeability or capillary hysteresis.
- Dissolution not modelled (it would remove some CO2 from the gas phase over decades).
- Closed model top: leaked gas pools in the overburden instead of migrating further.
- Pe,CO2/Pe,CH4 = 0.6 is provisional and may be generous.
- Leakage magnitudes depend on fault resolution (one cell overestimates by 19–45% and understates the curve-switch effect); the branch table uses the one-cell fault for comparability with Steps 2–3. Future fault work should use ≥ 5 fault columns.
- Producer rate quirk (above) and the (b′) residual of 12.8 t.
