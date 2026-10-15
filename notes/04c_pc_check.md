# Step 4C: Gas–water capillary pressure in the three-phase compositional model

Script: `scripts/04c_pc_check.jl`. Model: immiscible brine with CO2–CH4 Peng–Robinson (`src/gas.jl`), JutulDarcy 0.3.10.

## How capillary pressure enters this model
The phases are aqueous, liquid and vapour, and the reference pressure belongs to the hydrocarbon-liquid phase (index 2, checked with `get_reference_phase_index`), which never exists in my runs. With three phases and two curves, `src/variables/capillary.jl` gives p_water = p + Pc_ow(Sw) and p_gas = p + Pc_og(Sg). For a gas–water system that means setting Pc_ow ≡ 0 and Pc_og(Sg) = Pc_gw(Sw = 1 − Sg), so that p_gas − p_water = Pc_gw. I use Brooks–Corey, Pc_gw = Pe·Se^(−1/λ) with λ = 2 and Se = (Sw − Swr)/(1 − Swr), capped at Se = 1e-3, and assign one curve per rock type with `SimpleCapillaryPressure((pcow_tables, pcog_tables), regions = …)`.

## Tests

| Test | Set-up | Expected | Result |
|---|---|---|---|
| 1. Equilibrium slope | 60 m column, Pe = 0.2 bar, sharp gas-over-brine contact, 100 yr | At capillary–gravity equilibrium, Pc_gw(Sw(z)) falls linearly with depth at slope (ρw − ρg)·g | Fitted −8,858 Pa/m against −8,862 Pa/m expected (ratio 0.9996); ΔSw over the last 10 yr = 0 |
| 2a. Barrier holds | 4 m layer with Pe = 1.0 bar above a 10 m gas column (buoyancy ≈ 0.88 bar) | No gas crosses | 0.000 t of 8.4 t above the barrier after 50 yr |
| 2b. Barrier breached | Same, 20 m column (≈ 1.77 bar) | Gas crosses | 51.5% above the barrier after 50 yr |

Together these confirm the sign convention, the phase ordering, the saturation argument, the assignment by rock type and the threshold behaviour.

## A lesson about initial states
The first version of test 2a failed: 33% of the gas crossed a barrier it should not have crossed. The fault was in my test set-up, not the model. I had started the gas cells at Sw = Swr, which puts the capillary pressure at its maximum (the table cap, about 32 Pe = 6.3 bar for Pe = 0.2 bar), and that 6.3 bar of artificial capillary overpressure pushed gas through the 1 bar barrier at time zero. In a real reservoir, gas at Sw ≈ Swr exists only well above the free-water level; a column in equilibrium has a transition zone. The fix is to start from the equilibrium profile: at height h above the free-water level, Pc = Δρ g h and Sw = Swr + (1 − Swr)(Pe/(Δρ g h))^λ, with Sw = 1 where Δρ g h < Pe.

This matters for the faulted model (Step 4D): with capillary pressure on, the gas-filled trap has to start from a state consistent with capillary–gravity equilibrium, with each cell's saturation taken from its own rock type's curve at Pc = Δρ g (z_FWL − z). The Step 4B initial state (Sw = Swr throughout the gas zone) was only consistent because those runs had no capillary pressure.
