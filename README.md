# CO2 storage in a faulted reservoir: JutulDarcy.jl simulations

Paul Tivere Benjamin · October 2026 · work in progress

I started this project to learn [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl) and compositional reservoir simulation, using a question I find interesting: when CO2 is stored next to a fault, what decides whether it leaks across, and how much? The models are small 2D cross-sections with synthetic rock properties, so the results are about mechanisms and numerical behaviour rather than any real site. Each step comes with checks (mass conservation, grid and time-step sensitivity, comparison with simple analytical estimates), and the notes record what those checks found, including the places where I got something wrong and corrected it.

The project has two parts. Steps 1–3 inject CO2 into a saline aquifer next to a fault. Step 4 moves to the case I am most interested in: a depleted methane field against a fault, where CO2 is injected after production.

<p align="center"><img src="figures/03_leakage_vs_entry_pressure.png" width="85%"></p>

*Step 3: CO2 leaving the storage formation against the fault's capillary entry pressure, for three fault permeabilities. The grey band is the range of entry drive the fault experiences in this scenario.*

## Steps

| Step | What it does | Main finding | Notes |
|---|---|---|---|
| 1 | Builds the faulted cross-section | 30 m of throw against a 20 m caprock leaves a 10 m sand-to-sand window across the fault | [01](notes/01_geometry.md) |
| 2 | CO2–brine injection, with and without capillary pressure | Capillary entry pressure turns fault sealing into a threshold process | [02](notes/02_baseline.md), [02b](notes/02b_capillary.md), [02c](notes/02c_verification.md) |
| 3 | 36-run sweep of fault entry pressure and permeability | Three regimes: persistent leakage, leakage only during injection, and a capillary seal. Entry pressure mostly decides *whether* CO2 crosses, permeability *how much* | [03](notes/03_sweep.md) |
| 3b | Checks on Step 3 | Mass conserved to within 0.34 t; insensitive to time step; leakage magnitudes change by about ±25% with grid size | [03b](notes/03b_verification_report.md) |
| 3c | Step 3 rerun with the fault resolved by five columns | The regimes hold, but leakage near the threshold drops by 26–60% and persistent leakage rises by about 20% | [03c](notes/03c_fault_resolution.md) |
| 4A | Søreide–Whitson flash calculations for CO2–CH4–brine | Peng–Robinson brine density is 15–19% low; found a CH4 solubility term that had already been fixed upstream | [04a](notes/04a_flash_benchmark.md) |
| 4B | 1D and 2D depleted gas reservoirs | p/Z material balance matches to 0.12%; with gravity, injected CO2 settles below the remaining methane | [04b](notes/04b_depleted_gas_1d.md), [04b2](notes/04b2_gas_2d.md) |
| 4C | Capillary pressure in the three-phase model | Checked the sign convention and that a capillary barrier holds and breaks at the expected column heights | [04c](notes/04c_pc_check.md) |
| 4D | Faulted depleted methane field: production, then CO2 injection | CO2 entered the fault while the reservoir was still about 28 bar below its original pressure, because pressure in the low-permeability fault lagged behind. A no-injection control moved no gas. Lowering entry pressures for CO2 increased leakage by 15–41% | [04d](notes/04d_fault_gas.md) |

## Assumptions and limitations

- The aquifer cases (Steps 1–3) use JutulDarcy's `:co2brine` model at 95 °C and 30,000 mg/L NaCl. The salinity is a hypothetical value, not taken from a specific basin.
- The depleted-field cases (Step 4) use a three-phase compositional model: CO2 and CH4 with Peng–Robinson, and brine as a separate phase. Gas dissolution in brine is not included yet.
- Everything is 2D and isothermal, with no geomechanics, so nothing here says anything about fault slip.
- Fault properties are chosen, not derived from geology.
- The fault resolution matters. Steps 2–4 use a single 10 m cell across the fault; resolving it changes leakage in both directions (Steps 3c and 4D), so magnitudes from the single-cell runs should be read with that in mind.
- JutulDarcy's mass-balance tolerance is relative to the mass of the whole model, which is mostly brine. For the methane runs I tightened it (`tol_mb = 1e-9`); see the Step 4D notes.

## Running it

Julia 1.12 is needed ([juliaup](https://github.com/JuliaLang/juliaup) is the easiest way to install it). Package versions are pinned in `Manifest.toml` (JutulDarcy 0.3.10, Jutul 0.4.30, MultiComponentFlash 1.1.19, CairoMakie 0.15.15). Times are from my Windows laptop with 8 GB RAM.

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/01_build_section.jl              # seconds
julia --project=. scripts/02_run_cases.jl all              # ~40 min
julia --project=. scripts/03_sweep.jl                      # 36 runs, ~1.5–2 h
julia --project=. scripts/03b_verify.jl flux               # also: rate, dt, grid (grid runs are slow)
julia --project=. scripts/03c_fault_resolution.jl          # 6 runs
julia --project=. scripts/04a_flash_benchmark.jl           # minutes
julia --project=. scripts/04b_depletion_1d.jl              # minutes
julia --project=. scripts/04b2_gas_2d.jl stability
julia --project=. scripts/04b2_gas_2d.jl gravity
julia --project=. scripts/04c_pc_check.jl                  # minutes
julia --project=. scripts/04d_fault_gas.jl runs            # ~40 min
julia --project=. scripts/04d_fault_gas.jl runs_bhp        # ~50 min
julia --project=. scripts/04d_fault_gas.jl runs_fault5     # ~65 min
julia --project=. scripts/04d_fault_gas.jl runs_fault10    # ~4 h
julia --project=. scripts/04d_report.jl
julia --project=. scripts/make_figures.jl                  # all figures, from saved results
```

I checked that a fresh clone reproduces the committed numbers exactly for Steps 1 and 4C. The raw simulation output is too large to include, but the numbers behind every table and figure are in `results/` as CSV or Markdown.

## Layout

- `src/` model code: geometry, cases, fluid models, verification helpers, plotting style
- `scripts/` one script per step, plus `make_figures.jl`
- `notes/` what each step does, what I found and what I corrected; `walkthrough_04b_1d.md` walks through the 1D gas model in detail
- `results/`, `figures/` summary tables and figures

## Acknowledgements and license

Built on JutulDarcy.jl and Jutul.jl (SINTEF Digital) and MultiComponentFlash.jl. MIT license (see `LICENSE`).
