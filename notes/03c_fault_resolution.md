# Step 3c: Fault resolution in the Step 3 sweep

Script: `scripts/03c_fault_resolution.jl`. Results: `results/03c_fault_resolution.csv`.

## Why
In the depleted-field model (Step 4D, `notes/04d_fault_gas.md`) I found that using a single cell across the 10 m fault zone changes the leakage. Steps 2–3 use the same single-cell fault (the 5 m grid in Step 3b had two cells), so I needed to know whether the Step 3 conclusions survive a resolved fault. Rather than repeat all 36 runs, I reran six scenarios covering the three regimes with the fault split into five 2 m columns (`section_spec(fault_subcells = 5)`). Everything else (grid, fluid, schedule, solver settings) is as in Step 3, and the fault's pore volume and rock properties are identical by construction, so fault resolution is the only difference. Each run took 2–4.5 minutes.

This is a two-level sensitivity test, refined only across the fault; the 2 m vertical cells and the 10 m sand cells are unchanged. In Step 4D I tested three levels (1, 5 and 10 columns) and the change from 5 to 10 was 4–8%.

## Prediction before the runs
I expected the regimes to survive, persistent leakage to change little (it is limited by fault permeability), injection-triggered leakage to fall the most (its entry drive includes the pressure-lag part that a single wide cell exaggerates), and the sealed case to be unchanged.

## What is measured
- Outside storage: CO2 not in the storage sand, as % of injected (free and dissolved in the fault, caprock, overburden and underburden, plus what left through the open top).
- 0.1% threshold year: the first year with more than 0.1% outside storage. This is a reporting threshold, called "onset" in earlier versions. It includes dissolved CO2, so it is not the time of first entry.
- First invasion: the first year with more than 0.01 t of free (CO2-rich phase) CO2 in the fault.
- Above the seal: the overburden plus what left through the top.

## Results

| Scenario (Step 3 regime) | Fault k across / along (mD) | Columns | Outside storage, yr 25 / 100 (%) | 0.1% threshold year | First invasion year | Free / dissolved outside, yr 100 (%) | Left through top, yr 100 (%) | Above the seal, yr 100 (%) |
|---|---|---|---|---|---|---|---|---|
| Pe 0.6 (persistent) | 0.1 / 1 | 1 | 7.71 / 18.10 | 6 | 5 | 6.30 / 8.65 | 3.15 | 16.62 |
| | | 5 | 8.53 / 22.11 | 8 | 6 | 3.66 / 4.40 | 14.05 | 20.75 |
| Pe 1.0 (near threshold) | 0.1 / 1 | 1 | 2.35 / 2.40 | 11 | 9 | 0.60 / 1.80 | ~0 | 1.42 |
| | | 5 | 1.72 / 1.78 | 14 | 10 | 0.51 / 1.25 | 0.03 | 0.93 |
| Pe 1.2 (injection-triggered) | 0.1 / 1 | 1 | 0.96 / 0.97 | 15 | 13 | 0.28 / 0.69 | ~0 | 0.37 |
| | | 5 | 0.52 / 0.52 | 19 | 14 | 0.16 / 0.36 | ~0 | 0.13 |
| Pe 1.4 (injection-triggered, weak) | 0.1 / 1 | 1 | 0.30 / 0.30 | 20 | 18 | 0.09 / 0.21 | ~0 | 0.051 |
| | | 5 | 0.12 / 0.12 | 25 | 20 | 0.02 / 0.10 | ~0 | 0.005 |
| Pe 2.0 (capillary seal) | 0.1 / 1 | 1 | 0.042 / 0.041 | none | none | 0 / 0.041 | ~0 | 0.007 |
| | | 5 | 0.070 / 0.069 | none | none | 0 / 0.068 | ~0 | 0.004 |
| Pe 1.2, tight fault (injection-triggered) | 0.01 / 0.1 | 1 | 0.150 / 0.149 | 22 | 16 | 0.05 / 0.10 | ~0 | 0.006 |
| | | 5 | 0.076 / 0.073 | none | 14 | 0.03 / 0.05 | ~0 | ~0 |

The maximum entry drive falls by 0.1–0.16 bar in every case (for example 1.61 → 1.49 bar at Pe 1.2).

## What it shows
The three regimes survive, but the magnitudes change, and not all in the same direction.

Near the thresholds (Pe 1.0–1.4) leakage falls by 26–60%. The 0.1% threshold year moves 3–5 years later, but the first free CO2 still reaches the fault within 1–2 years of the single-cell timing, and in the tight-fault case it actually arrives 2 years earlier. So refinement mainly changes how much crosses, not when free CO2 first enters. The tight-fault case now stays below the 0.1% threshold, so by that measure it joins the sealed regime, although some free CO2 still enters the fault.

Persistent leakage increases (Pe 0.6: 18.1% → 22.1% outside storage at year 100), so my prediction was wrong here. The increase comes from CO2 reaching the top of the model faster (14.1% against 3.1% leaving through the top), while the free and dissolved CO2 held in the fault and overburden fell (15.0% → 8.1%).

I think two grid effects act in opposite directions. The first is entry: a single wide fault cell puts the brine pressure the gas has to overcome in the middle of the fault, which exaggerates the pressure-lag part of the entry drive. This matters most near the thresholds. The second is flow: a single 10 m cell spreads the gas across the whole width of the fault at low saturation. With the fault resolved, the gas travels at higher local saturation (at year 100, Pe 0.6: a peak of 0.36–0.41 in the 2 m columns against 0.29 in the single cell), and since the gas relative permeability rises steeply with saturation (Brooks–Corey exponent 2), it moves up the fault faster. This dominates when leakage is persistent. The saturation comparison supports this explanation but does not prove it.

In the injection-triggered case the gas stays in the three fault columns nearest the sand (peak gas saturation 0.15, 0.12 and 0.07, none in the hanging-wall half), which fits leakage limited by entry. In the sealed case only the dissolved transport changes (0.04% → 0.07%, still negligible).

## Consequences
Step 3's qualitative picture stands, with more careful wording: entry pressure mostly decides whether free CO2 enters the fault and permeability mostly decides how much moves, but the two interact, because permeability also shapes the pressure field and so whether the entry condition is reached. The tight-fault case is the clearest example.

The Step 3 magnitudes are biased in both directions: injection-triggered leakage is overestimated (by up to about 2.5 times at Pe 1.4) and persistent leakage underestimated (by about 20%). For new fault studies in this project I will use at least five columns across the fault zone, and report the first invasion separately from any reporting threshold.
