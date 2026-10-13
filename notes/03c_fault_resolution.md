# Step 3c: fault-resolution check of the Step 3 sweep

Script: `scripts/03c_fault_resolution.jl`. Results: `results/03c_fault_resolution.csv`.

## Why
Stage 4B-2c (`notes/04d_fault_gas.md`) showed that one cell across the 10 m fault zone changes leakage in the depleted-field case. Steps 2–3 use the same one-cell fault (the 5 m grid of Step 3b had two cells). Rather than repeat the 36-run sweep, six representative scenarios covering all three Step 3 regimes were rerun with the fault split into **5 columns of 2 m** (`section_spec(fault_subcells = 5)`). Everything else (grid, fluid, schedule, solver settings) is as in Step 3, so the only difference is fault resolution. Fault pore volume and per-rock-type properties are identical by construction. Each run took 2–4.5 min.

**Scope:** this is a mesh-sensitivity test with two levels, refined only across the fault. The vertical cell size (2 m) and the sand cells (10 m) are unchanged. Stage 4B-2c tested three levels (1, 5, 10 columns) and found a 4–8% change from 5 to 10.

## Prediction (before the runs)
The regimes survive; persistent leakage changes little (it is limited by fault permeability); injection-triggered leakage falls the most, because its entry drive includes the pressure-lag part that one wide cell exaggerates; the seal case is unchanged.

## Metrics
- **Outside storage:** CO2 not in the storage sand, as % of injected: free + dissolved in the fault, caprock, overburden and underburden, plus CO2 that left through the open top.
- **0.1% threshold year:** first year with > 0.1% outside storage. A *reporting threshold*, called "onset" in earlier versions; it includes dissolved CO2 and is not the time of first entry.
- **First invasion:** first year with > 0.01 t of *free* (CO2-rich phase) CO2 in the fault.
- **Above the seal:** overburden plus CO2 that left through the top.

## Results
| Scenario (Step 3 regime) | Fault k across / along (mD) | Columns | Outside storage, yr 25 / 100 (%) | 0.1% threshold year | First invasion year | Free / dissolved outside, yr 100 (%) | Left through top, yr 100 (%) | Above the seal, yr 100 (%) |
|---|---|---|---|---|---|---|---|---|
| Pe 0.6 (persistent) | 0.1 / 1 | 1 | 7.71 / 18.10 | 6 | 5 | 6.30 / 8.65 | 3.15 | 16.62 |
| | | 5 | 8.53 / **22.11** | 8 | 6 | 3.66 / 4.40 | **14.05** | 20.75 |
| Pe 1.0 (near threshold) | 0.1 / 1 | 1 | 2.35 / 2.40 | 11 | 9 | 0.60 / 1.80 | ~0 | 1.42 |
| | | 5 | 1.72 / 1.78 | 14 | 10 | 0.51 / 1.25 | 0.03 | 0.93 |
| Pe 1.2 (injection-triggered) | 0.1 / 1 | 1 | 0.96 / 0.97 | 15 | 13 | 0.28 / 0.69 | ~0 | 0.37 |
| | | 5 | 0.52 / 0.52 | 19 | 14 | 0.16 / 0.36 | ~0 | 0.13 |
| Pe 1.4 (injection-triggered, weak) | 0.1 / 1 | 1 | 0.30 / 0.30 | 20 | 18 | 0.09 / 0.21 | ~0 | 0.051 |
| | | 5 | 0.12 / 0.12 | 25 | 20 | 0.02 / 0.10 | ~0 | 0.005 |
| Pe 2.0 (capillary seal) | 0.1 / 1 | 1 | 0.042 / 0.041 | none | none | 0 / 0.041 | ~0 | 0.007 |
| | | 5 | 0.070 / 0.069 | none | none | 0 / 0.068 | ~0 | 0.004 |
| Pe 1.2, tight fault (injection-triggered) | 0.01 / 0.1 | 1 | 0.150 / 0.149 | 22 | 16 | 0.05 / 0.10 | ~0 | 0.006 |
| | | 5 | 0.076 / 0.073 | none | **14** | 0.03 / 0.05 | ~0 | ~0 |

The maximum entry drive falls by 0.1–0.16 bar in every case (e.g. 1.61 → 1.49 bar at Pe 1.2).

## What it shows
- **The three-regime structure survives**, but magnitudes change, in both directions.
- **Near the thresholds (Pe 1.0–1.4), leakage falls by 26–60%.** The 0.1% threshold year moves 3–5 years later, but **first invasion moves only 1–2 years** (and in the tight-fault case it comes 2 years *earlier*). So refinement mainly changes *how much* crosses, not *when* free CO2 first enters. The tight-fault case stays below the 0.1% threshold, i.e. it joins the seal regime by that measure, although some free CO2 still enters the fault.
- **Persistent leakage increases** (Pe 0.6: 18.1 → 22.1% outside at year 100). The prediction was wrong here. The increase is CO2 reaching the model top faster (14.1% vs 3.1% left through the top), while free + dissolved CO2 held in the fault and overburden fell (15.0 → 8.1%).
- **Two grid effects act in opposite directions:**
  1. *Entry:* one wide fault cell places the brine pressure the gas must overcome in the middle of the fault, exaggerating the pressure-lag part of the entry drive. This matters most near the thresholds.
  2. *Flow:* one 10 m cell spreads the gas across the whole fault width at low saturation. Resolved, the fault carries gas at higher local saturation (year 100, Pe 0.6: peak 0.36–0.41 in the 2 m columns vs 0.29 in the single cell), and gas relative permeability rises steeply with saturation (Brooks–Corey exponent 2), so up-fault transport is faster. This dominates when leakage is persistent. The saturation comparison supports this explanation; it is not a separate proof.
- In the injection-triggered case the gas stays in the three fault columns nearest the sand (peak Sg 0.15, 0.12, 0.07; none in the hanging-wall half), consistent with entry-limited leakage.
- The capillary-seal case changes only in dissolved transport (0.04 → 0.07%, still negligible).

## Consequences
- Step 3's qualitative picture stands, with a refinement of wording: entry pressure mainly controls *whether* free CO2 enters the fault and permeability mainly controls *how much* moves, but the two interact: permeability also shapes the pressure field and therefore whether the entry condition is reached (the tight-fault case is the clearest example).
- Step 3 magnitudes are biased in both directions: injection-triggered leakage is overestimated (by up to ~2.5× at Pe 1.4), persistent leakage underestimated (~20%).
- For new fault studies in this project: use ≥ 5 columns across the fault zone, and report first invasion separately from reporting thresholds.
