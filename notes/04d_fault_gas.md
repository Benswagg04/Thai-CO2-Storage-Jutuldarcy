# Step 4D: CO2 injection into a faulted, depleted methane field

Code: `src/gasfault.jl`, `scripts/04d_fault_gas.jl` (runs), `scripts/04d_report.jl` (tables), `scripts/make_figures.jl 04d` (figures), and `scripts/04d_conservation.jl` and `scripts/04d_tolerance.jl` (checks). Results: `results/04d_branches.csv`, `results/04d_timeseries.csv`; figures `figures/04d_*.png`. Fluid: immiscible brine with CO2–CH4 Peng–Robinson. Dissolution of CO2 or CH4 in brine is not modelled, so any "dissolved" amount is unknown, not zero.

## Question
After a methane field next to a fault has been produced, does re-injected CO2 cross the fault seal? And how much of that is due to switching from methane–brine to CO2–brine capillary properties, as opposed to the injection itself?

## Design
I kept the exact cross-section from Steps 1–3 (8,040 cells, 10 m × 2 m, a 10 m fault zone with k_across 0.01 mD and k_along 0.1 mD), so the fault is represented the same way as before. Later in this note a resolution test splits the fault into 5 and 10 columns.

The outer boundaries are closed, so the 2.54 km of padding on each side acts as a finite aquifer. The top is closed too, so gas that leaks into the overburden pools under it, as in Step 2b. The "beyond the seal" tally is correct, but the gas's onward migration is not represented.

The methane trap sits in the footwall storage sand against the fault, with the free-water level at 1,534 m, giving a 13.9 m gas column at the fault. Capillary pressure is Brooks–Corey per facies (λ = 2), checked in Step 4C; entry pressures are Leverett-scaled from the storage sand (0.05 bar) as in Steps 2–3, and the fault value is set directly.

For the methane–brine case I set the fault entry pressure to Pe,CH4 = 2.5 bar, about 2.0 times the buoyancy of the column at the fault (1.23 bar). This is a consistency requirement on a synthetic model (the initial trap has to be capillary-stable), with a margin. It is not a measured property, and for a real field the methane-column argument only gives a conditional lower bound (Naylor et al. 2011).

For the CO2 case I multiplied all entry pressures by 0.6 (fault Pe 1.5 bar), using Pe,CO2/Pe,CH4 ≈ (γ_CO2/γ_CH4)(cos θ_CO2/cos θ_CH4). This ratio is provisional. Naylor et al.'s reference values (25 against 50 mN/m at 10 MPa, wettability factor 0.7–1.0) suggest 0.35–0.5, so 0.6 may be generous to the CO2 seal, and interfacial tension also varies with pressure and composition, so a single ratio is a simplification.

The margin of 2 matters: with it, the CO2-scaled fault (1.2 times the buoyancy) can still hold the original methane column. With a margin below 1/0.6 ≈ 1.67 it could not, and any migration would come from the property change rather than from the CO2.

## Initial state
JutulDarcy's `EquilibriumRegion` would put gas everywhere above the contact, including the hanging-wall sand across the fault, so I built the state by hand: the brine pressure is hydrostatic, integrated with the compressible brine density; in the trap, Pc = (ρw − ρg) g (z_FWL − z), with Sw from each cell's Brooks–Corey curve; and every other cell is full of water.

| | Methane curves | CO2-scaled curves |
|---|---|---|
| Fault Pe / column buoyancy | 2.50 / 1.23 bar (2.03) | 1.50 / 1.23 bar (1.22) |
| Initial CH4 in place | 4,248 t | 4,530 t (thinner transition zone) |
| 10 yr without wells: CH4 into fault / beyond the seal | 0 / 0 t | 0 / 0 t |
| 10 yr without wells: max \|Δp\|, max \|ΔSw\| | 0.033 bar, 0.010 | 0.26 bar, 0.005 |
| Where the changes occur | 11 cells, all in the trap | 4 cells, all in the trap |
| Second 10 yr from the relaxed state | 0.0028 bar, 0.0023 | 0.015 bar, 0.0004 |

The hand-built profile is approximate, so I let it relax for 10 years without wells and start from the relaxed state. The remaining drift (at most 0.015 bar and 0.002 in Sw per 10 years) is about 1/1000 of the pressure changes caused by production; water in the transition zone near residual saturation is almost immobile and equilibrates very slowly.

## Schedule
One history, then branches restarted from the same depleted state at year 11:
- Depletion (methane curves): 10 years of production at a 300 t/yr target (bottom-hole floor 30 bar, never reached), then 1 year shut in. The producer is near the crest (x = 665 m, perforated 44–50 m below the top surface).
- (a) methane curves; CO2 at 1,500 t/yr for 10 years through the same well, then 20 years shut in.
- (b) CO2-scaled curves, same injection.
- (c) CO2-scaled curves, no injection: the control for the property switch.
- (a′) and (b′), added after (a): the same injection with the injector's bottom-hole pressure capped at the original brine pressure at the top perforation (153.1 bar). In (a), 15,000 t pushed the crest to 168 bar, 15 bar above its original 153 bar, which no operator would do in a depleted field.

## Checks

| Check | Result |
|---|---|
| Component conservation: (reservoir + well contents) − start − cumulative well flow, with the well flow integrated over every time step in a `post_ministep_hook` | Depletion: CH4 ≤ 0.008 t (3,095 t produced), water ≤ 0.02 t. Branches (a), (b), (c), (a′): CO2 ≤ 0.8 t of 11,000–15,000 t injected, CH4 ≤ 0.2 t. (b′): CO2 12.8 t (0.12%), CH4 1.4 t; this run predates the final well-flow bookkeeping, and the correction applied afterwards assumes all flow was injection, which is probably not exact for a pressure-capped well |
| Solver tolerance (below) | Year-41 leakage changes by ≤ 1.2 t between `tol_mb` 1e-9 and 1e-7 (e.g. (a): 2,425.6 against 2,424.8 t CO2) |
| Fault resolution (below) | Leakage depends on it: −31% for (a) with 10 fault columns |
| Control (c) | 0.00 t of CO2 or CH4 enters the fault in 30 years |

The mass-balance tolerance needed attention. JutulDarcy's convergence check on mass balance is relative to the mass of the whole model (~1.27 × 10¹⁰ kg, almost all brine). With the default `tol_mb = 1e-7`, each time step can leave ~1.3 t of error in any component, including the ~1,000 t of CH4: the control lost 1.95 t of CH4 (0.17%) while re-equilibrating, with nothing leaving the model. `1e-10` brings that down to 7 kg but stalls the injection once gas reaches the fault (median step 0.14 day, 9.6 of 30 years in 40 minutes). `1e-9` completes branch (a) in 24–39 minutes (12 minutes at 1e-7), and I use it for all results here. The physics is insensitive to this choice; the conservation is not.

The producer also behaved unexpectedly (JutulDarcy 0.3.10). The 300 t/yr production target is a standard-condition volume, and the producer delivered 300 rising to 323 t/yr instead of a constant 300. The conversion weights surface densities by the reservoir-condition volume fractions in the well, so a trace of water in the well (~6 × 10⁻⁵ by volume) raises the mixture density and the produced mass. The mass balance is unaffected, and the injector (fixed density) is exact. In total 3,095 t was produced.

## Results

### Depletion (years 0–11)
3,095 t of the 4,248 t of CH4 was produced (73%). The crest pressure fell from 153.1 bar to 90.4 bar (94.5 bar after shut-in), much less than in a closed gas tank, because 18,300 m³ of water flowed into the footwall storage sand, 38% of the initial gas pore volume of 48,500 m³. It came from the overburden on both sides (10,300 m³, mostly through the fault's sand-to-sand window), the hanging-wall storage sand (4,900 m³), the underburden (1,700 m³) and the caprock (1,300 m³): even 0.001 mD rock conducts measurable water over 3.7 km of section.

The fault became a conduit for water into the trap: the brine pressure in the fault zone stayed up to 18 bar above the crest (entry drive −17.5 bar at year 10), and no methane entered the fault.

### Branches (one fault column, tol_mb 1e-9; year-41 values unless stated)

| | CO2 injected (t) | Peak crest p (bar) | Peak entry drive (bar) / fault Pe | First year drive > Pe | CO2 beyond seal, yr 21 / 41 (t) | % of injected | CH4 beyond seal (t) |
|---|---|---|---|---|---|---|---|
| (a) methane curves | 15,000 | 168.4 | 8.8 / 2.5 | 13 | 1,770 / 2,426 | 16.2% | 213 |
| (b) CO2-scaled | 15,000 | 168.3 | 7.7 / 1.5 | 13 | 1,986 / 2,784 | 18.6% | 201 |
| (c) control | 0 | 103.4 | 0.7 / 1.5 | never | 0 / 0 | – | 0 |
| (a′) capped | 10,966 | 153.1 | 8.8 / 2.5 | 13 | 1,338 / 1,581 | 14.4% | 163 |
| (b′) capped | 11,019 | 153.1 | 7.7 / 1.5 | 13 | 1,499 / 1,847 | 16.8% | 163 |

"Beyond the seal" means the overburden sand plus the fault cells above the top of the footwall caprock. The caprock itself is hardly invaded (at most 1.2 t of CO2 in (b), where its entry pressure is also scaled by 0.6, and none in (a)). The CH4 beyond the seal is 14–18% of the 1,153 t left after depletion: the injected CO2 pushes the remaining methane out first.

To answer the question: the control (c) moves no gas at all, so all the migration in (b) is caused by injection, and the switch of capillary curves on its own causes none. Under injection the switch does matter: (b) − (a) = +358 t of CO2 (+15%) at constant rate, and (b′) − (a′) = +266 t (+17%) when capped, both with the one-cell fault. With the fault resolved the effect is larger, +676 t (+41%) with 10 columns (see below): the exaggerated part of the pressure lag shrinks, the drive sits closer to the entry pressures, and their difference (2.5 against 1.5 bar) matters more. Capping the injection pressure reduces leakage (by 35% for (a)) but does not prevent it.

### Mechanism: the fault's pore pressure lags behind the gas pressure
The entry drive is p_gas(crest) − p_brine(first fault cell), which splits into p_brine(crest) − p_brine(fault) plus Pc(crest). The capillary part is small, 0.3–0.9 bar, below both entry pressures throughout; the brine-pressure difference dominates at 4–8 bar during injection.

In the gas-filled crest the water is at Sw = 0.2015–0.2036 (residual 0.20), with relative permeability 3 × 10⁻⁶ to 2 × 10⁻⁵. It is essentially immobile, so its pressure rises with the injected gas instead of equalising with the fault. Inside the fault, brine can only drain slowly through 0.01 mD rock towards the hanging wall, which lags behind the trap by up to 16 bar (at year 13 the crest is at 125.3 bar and the hanging wall at 114.9 bar). Gas therefore enters the fault at year 13, while the crest is still about 28 bar below its original pressure (125.3 against 153.1 bar). Depletion makes this worse, because it drew the fault zone down to ~106 bar. After shut-in the pressures equalise, the drive falls below Pe 6–9 years after injection stops, and crossing stops.

### Fault resolution

| | 1 column (10 m) | 5 columns (2 m) | 10 columns (1 m) |
|---|---|---|---|
| (a) brine-pressure part of drive, year 15 (bar) | 7.69 | 4.92 | 4.78 |
| (a) CO2 beyond seal, year 41 (t) | 2,426 | 1,813 | 1,664 |
| (b) brine-pressure part of drive, year 15 (bar) | 6.92 | 3.89 | 3.66 |
| (b) CO2 beyond seal, year 41 (t) | 2,784 | 2,442 | 2,340 |

I predicted that refining the fault would shrink the brine-pressure part of the drive to about 1/5 and 1/10 and remove most of the leakage. That was wrong: it levels off at about 4.8 bar, because the immobile water in the gas zone, not the gradient inside the fault, sets most of the difference.

The single 10 m fault cell overestimates leakage by about 45% for (a) and 19% for (b), relative to 10 columns. The change from 5 to 10 columns is 8% for (a) and 4% for (b), so the converged values are probably about 1,600 t and 2,300 t. The effect of the curve switch grows with resolution: (b) − (a) = +358 t with 1 column, +629 t with 5 and +676 t with 10, so the single-cell fault understates how much the CO2 entry pressure matters. The 10-column runs took 50 minutes (a) and 187 minutes (b), against 24 and 13 minutes with one column, because the narrow fault cells force small time steps while gas moves through them.

Could the entry be an artefact of where the fault pressure is sampled? The entry drive above uses the centre of the first fault cell, which moves closer to the face as the fault is refined. A more direct measure is the CO2 that has left the storage sand, which equals the cumulative crossing of the sand–fault interface (the caprock takes at most 1.2 t):

| | First year with CO2 in the fault: 1 / 5 / 10 columns | CO2 crossed by year 21: 1 / 5 / 10 columns (t) | by year 41 (t) |
|---|---|---|---|
| (a) methane curves | 13 / 12 / 12 | 2,091 / 1,616 / 1,484 | 2,460 / 1,850 / 1,703 |
| (b) CO2-scaled | 13 / 12 / 12 | 2,372 / 2,139 / 2,053 | 2,831 / 2,488 / 2,388 |

Gas enters at the same time, or a year earlier, when the fault is resolved. If the entry came only from sampling the fault pressure 5 m inside the fault, refinement would delay or prevent it, so the pressure-lag mechanism stands; refinement changes how much crosses (−16 to −31% by year 41), not whether or when.

The refinement is across the fault only (1, 5 and 10 columns, with the fault's pore volume and rock properties identical by construction). The 2 m vertical cells and the 10 m sand cells next to the fault were not refined, so this is approximate convergence across the fault, not a full grid-convergence study. Steps 2–3 used the same single-cell fault; Step 3c (`notes/03c_fault_resolution.md`) found the bias goes both ways there, with injection-triggered leakage overestimated (up to about 2.5 times) and persistent leakage underestimated (about 20%), because one wide cell also spreads the gas at low saturation and slows its flow up the fault.

## Predictions against outcomes

| What I predicted before the runs | What happened |
|---|---|
| ~3,000 t produced; crest falls to 50–60 bar (little aquifer support) | 3,095 t produced, but the crest only fell to 90 bar: 18,300 m³ of water flowed in from regions I had not counted |
| No CH4 enters the fault during depletion | Correct (drive −17.5 bar) |
| 15,000 t brings the trap back to about its original pressure | Wrong: 168 bar, 15 bar above the original 153 bar |
| (a) holds unless injection pushes the drive above 2.5 bar | The drive reached 8.8 bar and the fault failed, through the pressure lag rather than the gas column |
| (c) shows essentially no migration | Correct (0.00 t) |
| A one-cell fault exaggerates the drive, and refinement removes most of the leakage | Partly wrong: leakage falls by 31% (a) and 16% (b), not to zero, and the effect of the curve switch grows instead of disappearing |

## Limitations
- One synthetic fault and one set of rock properties, in 2D, with no rock compressibility and no hysteresis in relative permeability or capillary pressure.
- No dissolution, which would take some CO2 out of the gas phase over decades.
- A closed model top, so leaked gas pools in the overburden instead of migrating further.
- The CO2/CH4 entry-pressure ratio of 0.6 is provisional and may be generous.
- Leakage depends on fault resolution (the single cell overestimates it by 19–45% and understates the curve-switch effect). The branch table uses the single-cell fault for comparison with Steps 2–3; future fault work should use at least 5 columns.
- The producer-rate behaviour described above, and the 12.8 t residual in (b′).
