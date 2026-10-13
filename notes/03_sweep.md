# Step 3: Fault-property sweep

## Design
- The fault's capillary entry pressure (Pe) and permeability were varied **independently**. In Steps 2–2c they were tied together by Leverett scaling, so the sand-rich vs clay-rich comparison changed both at once.
- Pe values (bar): 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0, 5.5. They are spaced closely around the entry drive measured in Step 2c (≈1.1 bar from buoyancy alone, ≈1.7 bar at the end of injection).
- Fault permeability (across / along, mD): 0.01/0.1, 0.1/1, 1/10.
- 36 runs. Grid 10 m × 2 m, which Step 2c showed is adequate for threshold behaviour and for total leakage out of the storage formation. All other settings as in Step 2c (Pc on, open top, 2.54 km padding, 1,500 t/yr for 25 years, 100 years total).
- Script: `scripts/03_sweep.jl` (restartable). Summary: `results/03_sweep_summary.csv`. Figures: `figures/03_*.png`.
- Wall time: 2.2–3 min per run, about 1.5 h in total.
- Metric: CO2 outside the storage formation = injected − CO2 in the footwall and hanging-wall storage sand. This includes the fault zone, caprock, overburden, and CO2 that left through the top.

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

Entry drive at the fault with a sealing fault (Pe ≥ 1.7):
- After injection: 1.09–1.11 bar, the same for all permeabilities. This is set by the ~18 m CO2 column.
- End-of-injection peak: 1.74, 1.71 and 1.54 bar for the three permeability levels. A more permeable fault lets brine through and relieves pressure.

## Interpretation
Three approximate regimes, set by where the fault's entry pressure sits relative to the CO2 pressure acting on it. The boundaries between them are not universal thresholds: they depend on the fault properties, the injection conditions and the numerical model, and the transitions are gradual rather than sharp.

1. **Pe below the buoyancy drive (≲ 0.9 bar): persistent leakage.**
   - CO2 enters during injection and keeps leaking afterwards, driven by buoyancy alone. Leakage grows 1.5–5× between year 25 and 100.
   - The column bleeds down toward the height the fault can hold, so the drive falls with time.
2. **Pe between the buoyancy drive and the injection peak (≈ 1.0–1.7 bar): injection-triggered fault invasion.**
   - The fault is invaded only while injection overpressure adds to buoyancy. After injection, the CO2 leaving the storage formation shows a negligible further increase: year-25 and year-100 values agree to within a few hundredths of a percent.
   - **Correction (Step 3b):** this does *not* mean migration stops. At Pe = 1.2 bar, the CO2 already inside the fault keeps moving up past the top of the seal after injection. Cumulative net crossing of the seal rises from 145 t (year 25) to 162 t (year 100). Entry into the fault stops; redistribution of what entered does not.
   - The width of this window equals the injection overpressure at the fault, so it depends on how the site is operated (rate, duration), not just on fault properties.
   - Pe = 1.0–1.1 bar is marginal: the post-injection drive (~1.07–1.09 bar) is only slightly above Pe. Flow nearly stops, because the fault's own Pc rises as CO2 enters it, and CO2 cannot flow until it exceeds the critical saturation (0.05).
3. **Pe above the injection peak (≳ 1.7 bar): the capillary seal holds.**
   - No free CO2 enters the fault. What remains (0.01–0.22%) is dissolved CO2 carried by brine flowing through the fault. It increases with fault permeability.
   - For the most permeable fault (1/10 mD), 0.13% reaches the overburden in solution even though the fault is a perfect capillary seal.

**Entry pressure decides *whether* the fault leaks; permeability decides *how much*.** (Approximately: the two interact, because permeability also shapes the pressure field and therefore whether the entry condition is reached; see Step 3c, tight-fault case.) At Pe = 0.4 bar, leakage at year 100 is 4.8%, 31% and 62% for the three permeability levels.

## Caveats
- The threshold values (1.1 and ~1.7 bar) are specific to this geometry (an ~18 m column at the fault), this injection rate and these boundaries. The transferable result is the three-regime structure and its controls: column height sets the lower threshold, injection overpressure sets the window.
- No molecular diffusion. The dissolved-CO2 leakage in regime 3 is advective only.
- Pc is a single scalar curve per cell, so the fault's entry pressure does not differ between across-fault and along-fault directions.
- The fault properties are prescribed values, not derived from geology (e.g. clay smear / PREDICT-type estimates).
- **Fault resolution (added after Step 3c):** the sweep uses one cell across the 10 m fault. Rerunning six scenarios with five 2 m columns (`notes/03c_fault_resolution.md`) keeps the three regimes but shifts magnitudes both ways: injection-triggered leakage falls by 26–60% (first free-CO2 invasion moves by only 1–2 years; the tight-fault case at Pe 1.2 stays below the 0.1% threshold), persistent leakage rises by ~20%. The "onset year" in `results/03_sweep_summary.csv` is now named `threshold_0p1pct_year` (first year with > 0.1% of the injected CO2 outside the storage sand, free + dissolved); the first year with free CO2 in the fault is reported separately.
