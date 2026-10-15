# Step 2b: Capillary pressure

## Set-up
Each facies gets one Brooks–Corey drainage curve, Pc = Pe·Se^(−1/λ), with λ = 2 and Se = (1 − S_CO2 − Swr)/(1 − Swr), capped at Se = 0.01. Entry pressures are scaled from the storage sand (Pe = 0.05 bar at 300 mD, φ = 0.25) with Leverett J, Pe ∝ √(φ/k):

| Unit | Pe (bar) | CO2 column it can hold (Δρ ≈ 630 kg/m³) |
|---|---|---|
| Storage sand | 0.05 | – |
| Overburden sand | 0.057 | – |
| Caprock / underburden shale | 12.2 | ~200 m |
| Clay-rich fault (k_across 0.01 mD) | 5.48 | ~89 m |
| Sand-rich fault (k_across 1 mD) | 0.55 | ~9 m |

In JutulDarcy brine is the reference phase and p_CO2 = p_brine + Pc(S_CO2) (I checked this in the source), so the table starts at Pe when S_CO2 = 0. The fault's entry pressure uses the across-fault permeability, because CO2 has to enter it from the side. One scalar curve per cell cannot represent different pore-throat sizes in different directions; that is a limitation.

## Checks
The equilibrium drift is 0.0000 bar in every case, and the rewritten code reproduces the Step 2 baseline to 0.1 t. The sand-rich fault leaks, so capillary pressure is not simply blocking all flow. The threshold arithmetic also works out:

| Case | Column at fault | Buoyancy (Δρ g h) | Fault Pe | Outcome |
|---|---|---|---|---|
| Clay-rich fault | 18 m | 1.11 bar | 5.48 bar | Holds |
| Sand-rich fault, yr 25 | 16 m | 0.99 bar | 0.55 bar | Leaks |
| Sand-rich fault, yr 100 | 14 m | 0.87 bar | 0.55 bar | Still leaking; column draining towards ~9 m |

## Results at year 100 (% of 37,500 t injected outside the footwall storage sand)

| Case | Outside storage sand |
|---|---|
| No Pc, clay-rich fault (Step 2) | 11% (up the fault and into the caprock) |
| Pc, clay-rich fault | ~0% |
| Pc, sand-rich fault (k 1/10 mD) | 46% (overburden 9,730 t; 7,120 t left the model) |

Capillary pressure removes the caprock seepage completely, as expected with 1 bar of buoyancy against 12 bar of entry pressure. It also turns the fault seal into a threshold process: the fault holds completely until the CO2 column exceeds Pe/(Δρ g), and then it leaks. In the sand-rich case leakage starts at about year 8, once the column is tall enough. After injection the leaking fault bleeds the column down towards the height it can hold (~9 m), so the leak slows.

The sand-rich fault also has 100 times the permeability, so this comparison changes two things at once. Step 3 separates them.

## Limitations
The model top is a no-flow boundary, so CO2 that reaches the overburden pools against it, which is artificial. In the sand-rich case some CO2 also migrates up-dip along the top of the hanging wall and leaves through the right boundary. Leverett scaling and λ = 2 are assumptions; I have no measured capillary data.
