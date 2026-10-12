# Step 4B, stage 2c-0: gas–water capillary pressure in the three-phase compositional model

Script: `scripts/04c_pc_check.jl`. Model: immiscible brine + CO2–CH4 Peng–Robinson (`src/gas.jl`), JutulDarcy 0.3.10.

## How capillary pressure enters this model
- Phases (aqueous, liquid, vapour); the **reference pressure is the hydrocarbon-liquid phase** (index 2, checked with `get_reference_phase_index`). That phase never exists in our runs.
- From `src/variables/capillary.jl` (three phases, two curves): p_water = p + Pc_ow(Sw), p_gas = p + Pc_og(Sg).
- For a gas–water system: **Pc_ow ≡ 0 and Pc_og(Sg) = Pc_gw(Sw = 1 − Sg)**, giving p_gas − p_water = Pc_gw. Brooks–Corey: Pc_gw = Pe·Se^(−1/λ), λ = 2, Se = (Sw − Swr)/(1 − Swr), capped at Se = 1e-3.
- Curves are assigned per region (rock type) with `SimpleCapillaryPressure((pcow_tables, pcog_tables), regions = …)`.

## Tests
| Test | Set-up | Expected | Result |
|---|---|---|---|
| 1. Equilibrium slope | 60 m column, Pe = 0.2 bar, sharp gas-over-brine contact, 100 yr | At capillary–gravity equilibrium, Pc_gw(Sw(z)) falls linearly with depth at slope (ρw − ρg)·g | Fitted −8,858 Pa/m vs expected −8,862 Pa/m (**ratio 0.9996**); ΔSw over the last 10 yr = 0 |
| 2a. Barrier holds | 4 m layer with Pe = 1.0 bar above a 10 m gas column (buoyancy ≈ 0.88 bar) | No gas crosses | **0.000 t** of 8.4 t above the barrier after 50 yr |
| 2b. Barrier breached | Same, 20 m column (≈ 1.77 bar) | Gas crosses | **51.5%** above the barrier after 50 yr |

This confirms the sign convention, phase ordering, saturation argument, per-region assignment and threshold behaviour.

## Lesson for initial states (important for stage B2c)
The first version of test 2a failed: 33% of the gas crossed a barrier it should not have. **The test set-up was at fault, not the model.**
- Starting gas cells at Sw = Swr means Pc at its maximum (the table cap, about 32 Pe = 6.3 bar for Pe = 0.2 bar). That 6.3 bar of artificial capillary overpressure forced gas through the 1 bar barrier at time zero.
- In a real reservoir, gas at Sw ≈ Swr exists only high above the free-water level. A column in equilibrium has a **transition zone**.
- The fix is to start the column from the equilibrium profile: at height h above the free-water level, Pc = Δρ g h, and Sw = Swr + (1 − Swr)(Pe/(Δρ g h))^λ (Sw = 1 where Δρ g h < Pe).

**Consequence for B2c:** with capillary pressure on, the gas-filled footwall trap must start from a capillary–gravity-consistent state, with each cell's saturation from its own rock type's curve at Pc = Δρ g (z_FWL − z). The stages B2a/B2b state (Sw = Swr throughout the gas zone) was consistent only because those runs had no capillary pressure.
