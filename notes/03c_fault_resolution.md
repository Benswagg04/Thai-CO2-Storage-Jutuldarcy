# Step 3c: fault-resolution check of the Step 3 sweep

Script: `scripts/03c_fault_resolution.jl`. Results: `results/03c_fault_resolution.csv`.

## Why
Stage 4B-2c (`notes/04d_fault_gas.md`) showed that one cell across the 10 m fault zone overestimates leakage in the depleted-field case. Steps 2–3 use the same one-cell fault (the 5 m grid of Step 3b had two cells). Rather than repeat the 36-run sweep, six representative scenarios covering all three Step 3 regimes were rerun with the fault split into **5 columns of 2 m** (`section_spec(fault_subcells = 5)`). Everything else (grid, fluid, schedule, solver settings) is as in Step 3, so the only difference is fault resolution. Each run took 2–4.5 min.

## Prediction (before the runs)
The regimes survive; persistent leakage changes little (it is limited by fault permeability); injection-triggered leakage falls the most, because its entry drive includes the pressure-lag part that one wide cell exaggerates; the seal case is unchanged.

## Results
| Scenario (Step 3 regime) | Fault k across / along (mD) | Outside storage, yr 25 / 100 (% of injected): 1 → 5 columns | Onset year | Above the seal, yr 100 (%) | Max entry drive (bar) |
|---|---|---|---|---|---|
| Pe 0.6 (persistent) | 0.1 / 1 | 7.71 / 18.10 → 8.53 / **22.11** | 6 → 8 | 16.62 → 20.75 | 1.43 → 1.27 |
| Pe 1.0 (near threshold) | 0.1 / 1 | 2.35 / 2.40 → 1.72 / 1.78 | 11 → 14 | 1.42 → 0.93 | 1.56 → 1.44 |
| Pe 1.2 (injection-triggered) | 0.1 / 1 | 0.96 / 0.97 → 0.52 / 0.52 | 15 → 19 | 0.37 → 0.13 | 1.61 → 1.49 |
| Pe 1.4 (injection-triggered, weak) | 0.1 / 1 | 0.30 / 0.30 → 0.12 / 0.12 | 20 → 25 | 0.051 → 0.005 | 1.65 → 1.55 |
| Pe 2.0 (capillary seal) | 0.1 / 1 | 0.042 / 0.041 → 0.070 / 0.069 | none → none | 0.007 → 0.004 | 1.71 → 1.59 |
| Pe 1.2, tight fault (injection-triggered) | 0.01 / 0.1 | 0.150 / 0.149 → 0.076 / 0.073 | 22 → **none** | 0.006 → ~0 | 1.60 → 1.48 |

"Outside storage" includes CO2 dissolved in brine in the fault and caprock; "onset" is the first year with > 0.1% outside; "above the seal" is overburden plus CO2 that left through the top.

## What it shows
- **The three-regime structure survives**, but the regime boundaries move. Near the thresholds (Pe 1.0–1.4) leakage falls by 26–60% and starts 3–5 years later; the tight-fault case (Pe 1.2, 0.01 mD) moves from injection-triggered into the seal regime. The maximum entry drive falls by 0.1–0.16 bar in every case.
- **Persistent leakage increases** (+22% at year 100, Pe 0.6). The prediction was wrong here.
- **Two grid effects act in opposite directions:**
  1. *Entry (whether gas gets in):* one wide fault cell places the brine pressure the gas must overcome in the middle of the fault, exaggerating the pressure-lag part of the drive, so gas enters too easily. This dominates near the thresholds.
  2. *Flow (how much moves up the fault once it is in):* one 10 m cell spreads the gas across the whole fault width at low saturation. Resolved, the fault carries gas at higher local saturation (year 100, Pe 0.6: peak 0.36–0.41 in the 2 m columns vs 0.29 in the single cell), and gas relative permeability rises steeply with saturation (Brooks–Corey exponent 2), so up-fault flow is faster. This dominates when leakage is persistent. The saturation comparison supports this explanation; it is not a separate proof.
- In the injection-triggered case the gas stays in the three fault columns nearest the sand (peak Sg 0.15, 0.12, 0.07; none in the hanging-wall half), consistent with entry-limited leakage.
- The capillary-seal case changes only in the dissolved transport (0.04 → 0.07%, still negligible), which is dominated by dissolved CO2 moving into the fault zone.

## Consequences
- Step 3's qualitative conclusion stands: entry pressure decides *whether* free CO2 crosses the fault, permeability decides *how much*.
- Step 3 magnitudes are biased in both directions: injection-triggered leakage is overestimated (by up to a factor ~2.5 at Pe 1.4), persistent leakage underestimated (~20%), and the regime boundary lies at a somewhat lower entry pressure than the sweep shows.
- Only 1 vs 5 columns was tested here; Stage 4B-2c found the 5 → 10 change to be 4–8%, so 5 columns is a reasonable working resolution.
- For new fault studies in this project: use ≥ 5 columns across the fault zone.
