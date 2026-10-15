# Step 3: Fault-property sweep

## Design
In Steps 2–2c the fault's entry pressure and permeability were tied together by Leverett scaling, so the sand-rich versus clay-rich comparison changed both at once. Here I varied them independently:

- Entry pressure Pe (bar): 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0, 5.5, spaced closely around the entry drive measured in Step 2c (≈1.1 bar from buoyancy alone, ≈1.7 bar at the end of injection).
- Fault permeability, across / along (mD): 0.01/0.1, 0.1/1, 1/10.

That makes 36 runs on the 10 m × 2 m grid, which Step 2c showed is good enough for the threshold behaviour and for the total leaving the storage formation. Everything else is as in Step 2c (capillary pressure on, open top, 2.54 km padding, 1,500 t/yr for 25 years, 100 years in total). Each run took 2.2–3 minutes, about 1.5 hours altogether. The script is `scripts/03_sweep.jl` (it can be restarted), the summary is `results/03_sweep_summary.csv`, and the figures are `figures/03_*.png`.

The measure throughout is the CO2 outside the storage formation: injected minus what is in the footwall and hanging-wall storage sand. It includes the fault zone, the caprock, the overburden and anything that left through the top.

## Results (% of 37,500 t injected outside the storage formation)

| Pe (bar) | 0.01/0.1 mD, yr 25 → 100 | 0.1/1 mD, yr 25 → 100 | 1/10 mD, yr 25 → 100 |
|---|---|---|---|
| 0.4 | 1.0 → 4.8 | 11.6 → 31.0 | 41.0 → 62.3 |
| 0.6 | 0.7 → 2.5 | 7.7 → 18.1 | 28.0 → 45.5 |
| 0.8 | 0.5 → 1.0 | 4.6 → 8.0 | 16.7 → 25.7 |
| 0.9 | 0.4 → 0.6 | 3.4 → 4.6 | 12.0 → 15.6 |
| 1.0 | 0.28 → 0.29 | 2.35 → 2.40 | 8.0 → 8.0 |
| 1.1 | 0.21 → 0.21 | 1.56 → 1.58 | 4.9 → 4.9 |
| 1.2 | 0.15 → 0.15 | 0.96 → 0.97 | 2.6 → 2.6 |
| 1.4 | 0.07 → 0.06 | 0.30 → 0.30 | 0.49 → 0.47 |
| ≥ 1.7 | 0.01 | 0.04 | 0.22 |

With a sealing fault (Pe ≥ 1.7), the entry drive after injection is 1.09–1.11 bar for all three permeabilities, set by the ~18 m CO2 column. The peak at the end of injection is 1.74, 1.71 and 1.54 bar for the three permeability levels: a more permeable fault lets brine through and relieves pressure.

## Interpretation
The runs fall into three rough regimes, depending on where the fault's entry pressure sits relative to the CO2 pressure acting on it. The boundaries are not universal thresholds; they depend on the fault properties, the injection and the model, and the transitions are gradual.

**Pe below the buoyancy drive (≲ 0.9 bar): persistent leakage.** CO2 enters the fault during injection and keeps leaking afterwards on buoyancy alone; leakage grows 1.5–5 times between years 25 and 100. The column bleeds down towards the height the fault can hold, so the drive falls with time.

**Pe between the buoyancy drive and the injection peak (≈ 1.0–1.7 bar): the fault is invaded only during injection**, while the injection overpressure adds to buoyancy. Afterwards the CO2 outside the storage formation hardly increases: the year-25 and year-100 values agree to within a few hundredths of a percent. That does not mean migration stops, though. Step 3b showed that at Pe = 1.2 bar the CO2 already in the fault keeps moving up past the top of the seal: the net amount crossing the seal rises from 145 t at year 25 to 162 t at year 100. Entry into the fault stops; what entered keeps redistributing. The width of this window equals the injection overpressure at the fault, so it depends on how the site is operated (rate, duration) as well as on the fault. Pe = 1.0–1.1 bar is marginal: the post-injection drive (~1.07–1.09 bar) is only slightly above Pe, and flow nearly stops because the fault's own capillary pressure rises as CO2 enters it, and CO2 cannot flow until it exceeds the critical saturation (0.05).

**Pe above the injection peak (≳ 1.7 bar): the capillary seal holds.** No free CO2 enters the fault. What remains (0.01–0.22%) is dissolved CO2 carried by brine flowing through the fault, and it increases with fault permeability: for the most permeable fault (1/10 mD), 0.13% reaches the overburden in solution even though the fault is a perfect capillary seal.

In short, entry pressure decides *whether* the fault leaks and permeability decides *how much*. At Pe = 0.4 bar the leakage at year 100 is 4.8%, 31% and 62% for the three permeabilities. This is only approximately true, because permeability also shapes the pressure field and therefore whether the entry condition is reached; the tight-fault case in Step 3c shows this.

## Caveats
- The thresholds (1.1 and ~1.7 bar) are specific to this geometry (an ~18 m column at the fault), this injection rate and these boundaries. What should carry over is the three-regime structure and what controls it: the column height sets the lower threshold, and the injection overpressure sets the width of the window.
- There is no molecular diffusion, so the dissolved leakage in the third regime is advective only.
- Capillary pressure is one curve per cell, so the fault's entry pressure is the same across and along it.
- The fault properties are chosen values, not derived from geology (for example from clay smear or PREDICT-type estimates).
- Fault resolution, added after Step 3c: the sweep uses one cell across the 10 m fault. Rerunning six of the scenarios with five 2 m columns (`notes/03c_fault_resolution.md`) keeps the three regimes but shifts the magnitudes both ways. Leakage near the threshold falls by 26–60% (although the first free CO2 still reaches the fault within 1–2 years of the one-cell timing, and the tight-fault case at Pe 1.2 stays below the 0.1% threshold), and persistent leakage rises by about 20%. The "onset year" in `results/03_sweep_summary.csv` is now called `threshold_0p1pct_year` (first year with more than 0.1% of the injected CO2 outside the storage sand, free plus dissolved), and the first year with free CO2 in the fault is reported separately.
