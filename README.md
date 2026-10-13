# CO2 storage in a faulted reservoir: JutulDarcy.jl simulations with verification

**Author:** Paul Tivere Benjamin · **Status:** ongoing learning project (October 2026). All models are simplified 2D cross-sections with synthetic properties; results describe mechanisms and numerical behaviour, not predictions for a specific site.

A sequence of small, reproducible reservoir simulations in [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl) (SINTEF), built step by step to learn compositional reservoir simulation. The focus is CO2 migration across and along a fault, in a saline aquifer and in a depleted gas field. Every step pairs a model with checks of its numerical and physical consistency, and the notes record what each check found, including corrections.

## What is in the repository

| Step | Question | Main result | Notes |
|---|---|---|---|
| 1 | Build a faulted cross-section | Tilted fault block, 30 m throw > 20 m caprock → 10 m sand–sand leak window; fault as its own anisotropic zone | [01](notes/01_geometry.md) |
| 2, 2b, 2c | CO2–brine injection; effect of capillary pressure; verification | Capillary entry pressure turns fault seal into a threshold process. Boundaries, Pc range, accounting and grid resolution checked | [02](notes/02_baseline.md), [02b](notes/02b_capillary.md), [02c](notes/02c_verification.md) |
| 3 | Fault entry pressure × permeability sweep (36 runs) | Three approximate regimes: persistent leakage; injection-triggered fault invasion; capillary seal with dissolved-only transport. Entry pressure mainly controls *whether*, permeability mainly *how much*; the two interact, since permeability also shapes the pressure field | [03](notes/03_sweep.md) |
| 3b | Verification of Step 3 | Conservation residual ≤ 0.34 t (tight tolerances); boundary flux integrated per ministep; time-step insensitive; grid sensitivity ~±25% in magnitudes (not converged); injection rate controls free-CO2 seal breach at fixed total mass | [03b report](notes/03b_verification_report.md) |
| 3c | Fault-resolution check of Step 3 (6 scenarios, fault split into five 2 m columns) | The three regimes survive. Near the thresholds leakage falls 26–60% (one wide fault cell exaggerates the entry drive), while first free-CO2 invasion moves by only 1–2 years; persistent leakage rises ~20% (one wide cell also spreads the gas at low saturation, slowing its flow up the fault). Two-level, fault-normal mesh-sensitivity test | [03c](notes/03c_fault_resolution.md) |
| 4A | Søreide–Whitson equation of state for CO2–CH4–brine (standalone) | CO2 solubility vs JutulDarcy's CO2-brine correlations; weak salting-out; Peng–Robinson brine density 15–19% low; a CH4 term discrepancy in the installed version, already fixed upstream | [04a](notes/04a_flash_benchmark.md) |
| 4B stage 1 | 1D depleted gas reservoir (CO2–CH4 Peng–Robinson, immiscible brine) | Mass conserved to ~1e-6; p/Z material balance within 0.12%, explained by connate-water expansion; mixing zone controlled by numerical dispersion | [04b](notes/04b_depleted_gas_1d.md) |
| 4B stage 2 | 2D gas–water reservoir: equilibrium and gravity | Stable initial state; CO2-rich gas collects beneath CH4 with gravity on (CO2–CH4 mean-depth difference +13 m vs −11 m with gravity off) | [04b2](notes/04b2_gas_2d.md) |
| 4C | Capillary pressure in the three-phase compositional model | Sign convention, phase ordering and per-rock-type curves verified (equilibrium slope ratio 0.9996; barrier holds and breaks at the predicted column heights) | [04c](notes/04c_pc_check.md) |
| 4B stage 2c | Faulted depleted methane field: production, then CO2 injection, with methane vs CO2-scaled capillary curves and a no-injection control | Depletion draws ~18,000 m³ of water into the trap, largely through the fault. During injection the fault's pore pressure lags the gas pressure (residual water in the gas zone is immobile; the 0.01 mD fault drains slowly), so gas crosses the fault while the reservoir is still ~28 bar below its original pressure. The control moves no gas; CO2-scaled capillary curves increase leakage (+15% with a one-cell fault, +41% with the fault resolved by 10 columns); capping injection at the original pressure reduces but does not stop it. A single 10 m fault cell overestimates leakage by 19–45% | [04d](notes/04d_fault_gas.md) |

<p align="center"><img src="figures/03_leakage_vs_entry_pressure.png" width="85%"></p>

*CO2 leaving the storage formation vs fault capillary entry pressure, for three fault permeabilities (Step 3). The grey band is a scenario-specific reference range of the CO2 pressure acting on the fault.*

## Models and main assumptions
- **Saline aquifer (Steps 1–3):** JutulDarcy's `:co2brine` K-value model (H2O–CO2, CO2 dissolution and water evaporation; correlations from Saló-Salgado et al. 2024), 95 °C, 30,000 mg/L NaCl-equivalent (a *hypothetical* baseline), Brooks–Corey relative permeability with Killough hysteresis, Leverett-scaled Brooks–Corey capillary pressure per facies.
- **Depleted gas field (Step 4B):** three-phase compositional model with brine as an immiscible phase and CO2–CH4 described by Peng–Robinson. Excludes gas dissolution in brine, water evaporation and diffusion.
- 2D cross-sections; isothermal; no geomechanics. Pressure results say nothing about fault slip.
- Fault properties are prescribed, not derived from geology (e.g. clay content).

## Known limitations (see each step's notes)
- Leakage magnitudes carry roughly ±25% grid sensitivity; trapping fractions are grid-dependent and are not reported as results.
- Steps 2–4 represent the 10 m fault zone with one cell across. Resolving it (5–10 columns, `section_spec(fault_subcells = …)`) changes leakage magnitudes in both directions: injection-triggered leakage falls (by 19–45% in Stage 4B-2c, 26–60% near the Step 3 thresholds) while persistent leakage rises (~20%). The three-regime structure of Step 3 survives; see Step 3c.
- JutulDarcy's mass-balance tolerance is relative to the whole model's mass (mostly brine); minor components need a tighter `tol_mb` (1e-9 is used in stage 4B-2c).
- JutulDarcy 0.3.10's open boundary flux ignores gravity across the half cell and capillary pressure.
- No molecular diffusion or mechanical dispersion; mixing widths are numerical.
- Capillary pressure is one scalar curve per cell (no directional entry pressure in the fault).

## Reproducing
Requires Julia 1.12 (via [juliaup](https://github.com/JuliaLang/juliaup)). Exact package versions are pinned in `Manifest.toml` (JutulDarcy 0.3.10, Jutul 0.4.30, MultiComponentFlash 1.1.19, CairoMakie 0.15.15).

Reproducibility was checked on 2026-10-13 by cloning this repository into an empty folder, instantiating it from the committed `Manifest.toml` and rerunning Steps 1 and 4C, which reproduced the committed numbers exactly (e.g. equilibrium slope ratio 0.9996; barrier test 0.000 t and 51.52%). That check reused the local package cache; a run with an empty package depot has not been done yet.

Run times below are wall-clock on the author's Windows laptop (8 GB RAM, one simulation at a time); the first call in a fresh Julia session adds a few minutes of compilation.

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/01_build_section.jl              # seconds
julia --project=. scripts/02_run_cases.jl all              # ~40 min (the fine-grid case takes ~30)
julia --project=. scripts/03_sweep.jl                      # 36 runs, ~1.5–2 h
julia --project=. scripts/03b_verify.jl flux               # also: rate, dt, grid (grid: 46–109 min per run)
julia --project=. scripts/03c_fault_resolution.jl          # 6 runs with a 5-column fault
julia --project=. scripts/04a_flash_benchmark.jl           # minutes
julia --project=. scripts/04b_depletion_1d.jl              # minutes
julia --project=. scripts/04b2_gas_2d.jl stability
julia --project=. scripts/04c_pc_check.jl                  # minutes
julia --project=. scripts/04d_fault_gas.jl runs            # depletion + branches a, b, c: ~40 min
julia --project=. scripts/04d_fault_gas.jl runs_bhp        # pressure-capped branches: ~50 min
julia --project=. scripts/04d_fault_gas.jl runs_fault5     # 5-column fault: ~65 min
julia --project=. scripts/04d_fault_gas.jl runs_fault10    # 10-column fault: ~4 h
julia --project=. scripts/04d_report.jl                    # tables and figures from the runs above
```

Simulation outputs (`output/`, Julia binary files of up to several hundred MB) are not version-controlled. The numbers behind every table and figure are, as plain CSV/Markdown in `results/`:

| File | Contents |
|---|---|
| `03_sweep_summary.csv` | Step 3: leakage (free and dissolved), 0.1% threshold year, first free-CO2 invasion of the fault, and entry drive for all 36 fault scenarios |
| `03b_tables.md` | Step 3b verification tables: boundary transport and conservation, fault invasion vs seal crossing, injection rate at fixed mass, time step and grid |
| `03c_fault_resolution.csv` | Step 3c: six Step 3 scenarios with 1 vs 5 fault columns |
| `04a_co2_brine.csv`, `04a_ch4_water.csv` | Step 4A flash benchmark |
| `04d_branches.csv` | Stage 4B-2c: one summary row per branch (leakage, peak pressures, conservation residuals, solver statistics) |
| `04d_timeseries.csv` | Stage 4B-2c: year-by-year pressures, entry-drive components, gas mass by region and conservation residuals for every run (incl. 5- and 10-column fault runs) |

## Layout
- `src/`: model code (geometry, cases, run/verification helpers, gas-field fluid model)
- `scripts/`: one script per step
- `notes/`: reasoning, results, checks and corrections for each step; `notes/walkthrough_04b_1d.md` is a guided walkthrough for running and modifying the 1D depleted gas model
- `results/`, `figures/`: summary tables and figures (`figures/superseded/`: first-configuration figures from Steps 2–2b)

## Acknowledgements
JutulDarcy.jl and Jutul.jl (SINTEF Digital); MultiComponentFlash.jl.

## License
MIT (see `LICENSE`).
