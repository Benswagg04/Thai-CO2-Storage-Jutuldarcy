# Walkthrough: the 1D depleted gas reservoir (Step 4B, stage 1)

Goal: be able to run, explain and modify `scripts/04b_depletion_1d.jl` and `src/gas.jl` yourself. Work through it with the two files open side by side.

## 0. Run it
In a terminal in the project folder:
```
julia --project=.
```
Then, in the Julia REPL (keeping the session open avoids recompiling every time):
```julia
SKIP_MAIN = true                            # load the functions without running everything
include("scripts/04b_depletion_1d.jl")
r = run_1d(100)                             # one run, 100 cells (~1–3 min after compilation)
names, mb = mass_balance(r)
pz_check(r, mb)
```
`main()` runs the full script (two grids, figures). The first call in a fresh session is slow because Julia compiles the code.

## 1. The physics being solved
For each component i (CO2, CH4, water), mass conservation in every cell:

  ∂/∂t ( φ Σ_phases ρ_α S_α x_iα ) + ∇·( Σ_phases ρ_α x_iα v_α ) = q_i

- φ porosity; ρ_α, S_α density and saturation of phase α; x_iα mass fraction of component i in phase α; q_i well source/sink.
- Darcy velocity of each phase: v_α = −(k k_rα / μ_α)(∇p_α − ρ_α g ∇z). Gravity plays no role in this 1D horizontal model.
- Phase equilibrium: CO2 and CH4 split between hydrocarbon liquid and vapour so that their fugacities are equal (a "flash"), with Peng–Robinson giving fugacities, densities and the compressibility factor Z. At 95 °C both are supercritical, so there is one gas phase (checked: hydrocarbon-liquid saturation = 0).
- Water is a separate phase that never mixes with the gas; its density is ρ_w(p) = ρ_ref (1 + c_w (p − p_ref)) and its viscosity is constant.
- JutulDarcy discretises this with finite volumes (two-point fluxes between neighbouring cells), solves it fully implicitly in time, and uses Newton's method with automatic differentiation for the nonlinear equations.

## 2. `src/gas.jl`: the fluid
| Function | What it does | Why |
|---|---|---|
| `gas_spec()` | All fluid inputs in one place (95 °C, 150 bar, Swc 0.2, k_CO2–CH4 0.10, brine properties) | Change inputs in one spot |
| `gas_eos(gs)` | Peng–Robinson for CO2 and CH4 (critical properties from MultiComponentFlash's table), with a binary interaction coefficient | The equation of state behind densities, Z and the flash |
| `gas_system(gs)` | `MultiPhaseCompositionalSystemLV` with phases (aqueous, liquid, vapour). The aqueous phase is the "other", immiscible phase. | This is what makes brine separate from the EOS components |
| `set_brine_properties!` | Replaces the default water properties (1 kg/m³, 1 mPa·s) with brine values | **Without it the water would be 1,000× too light.** Try removing it and see what happens. |
| `co2_injector_control` | CO2 injection at a mass rate | For injectors, surface volume × `density` = mass |
| `gas_producer_control` | Production at a mass rate, converted to standard-condition volume (1 atm, 15 °C) with PR density | For compositional producers JutulDarcy flashes the produced gas to standard conditions. Getting this wrong produced 1,361 instead of 2,000 t/yr. |
| `gas_Z`, `gas_density` | Z and density from PR | Used in the p/Z check and the rate conversion |

## 3. `run_1d`: building and running the case, line by line
1. `reservoir_mesh((nx, 1, 1), (L, W, H))`: a 1D Cartesian grid, 1,000 m × 100 m × 20 m.
2. `reservoir_domain(...)`: attaches rock properties (k, φ) and temperature to the mesh.
3. `setup_well(domain, [(1,1,1)], ...)`: injector in the first cell, producer in the last. `simple_well = true` means a single well node, no wellbore flow model.
4. `setup_reservoir_model(domain, sys, wells = ...)`: builds the coupled reservoir and well model. It returns a *MultiModel* with parts `:Reservoir`, `:Injector`, `:Producer` and their controls.
5. `set_brine_properties!`, then `replace_variables!(..., RelativePermeabilities = kr)`: Brooks–Corey curves; residual water 0.2 means the connate water never flows.
6. `setup_reservoir_state(...)`: initial pressure, gas composition (pure CH4) and water saturation. These are the **primary variables** the solver iterates on. Everything else (densities, saturations, masses) is computed from them as **secondary variables**.
7. `setup_reservoir_forces(...)`: one set of well controls per period (production, shut-in, injection). `limits` adds the bottom-hole pressure floor.
8. The 1-second, wells-shut simulation: a trick to obtain the initial component masses (`TotalMasses` is a secondary variable, so it is not in `state0`).
9. `simulate_reservoir(state0, model, dt; forces, ...)`: runs the schedule. `dt` is a vector of report steps (1 year each); the solver takes smaller internal steps ("ministeps") as needed, up to `max_timestep`. `output_substates = true` keeps every ministep, so that well flows can be integrated exactly.

## 4. The checks, and why each one is independent
- **`mass_balance`:** for each component, initial mass + Σ(injected − produced) × ministep length should equal the mass in place. The well rates and the cell masses are computed separately, so agreement (1e-6 here) is a real test.
- **`pz_check`:** for a closed gas tank, p/Z = (p/Z)_i (1 − G_p/G). This comes from classical reservoir engineering, not from the simulator. The 0.12% offset equals the connate-water expansion c_w Δp S_w/(1 − S_w); try setting `brine_compressibility = 0` and see whether the offset disappears.
- **Phase check:** the maximum hydrocarbon-liquid saturation must be 0 (single gas phase).
- **PR CO2 density vs a NIST-based table:** checks the equation of state independently of the flow model.

## 5. Exercises: predict first, then run
1. Halve the production rate. How long does it take to reach 42 bar? Is the p/Z line still within 0.12%?
2. Set `swc = 0.1` in `gas_spec`. How do the initial CH4 in place and the p/Z offset change? (Predict both numbers from the formulas above.)
3. Lower `bhp_min_bar` to 50. When does the producer switch from rate to pressure control, and what does the pressure curve look like after that?
4. Remove `set_brine_properties!`. What breaks first, and why?
5. Change `max_dt` from 30 to 3 days (as in the mixing-zone test). Which outputs change and which do not?

## 6. Troubleshooting seen in this project
- *"Package X not found in current path"*: the package is not a direct dependency. Add it with `Pkg.add(...)` (pin the version).
- *"Unsupported well injector constraint/limit mrat"*: mass-rate producers are not supported in JutulDarcy 0.3.10. Use surface-volume rates.
- *A rate that comes out wrong by a constant factor*: check whether it is a surface or a standard-condition volume.
- *A "steady" state that moves*: check the initial state is in equilibrium (hydrostatic pressure; with capillary pressure, the transition zone; see `notes/04c_pc_check.md`).
- *Very slow runs with many cut steps*: see the tuning notes in `notes/04b2_gas_2d.md`. Compare Newton iterations, not wall time.
