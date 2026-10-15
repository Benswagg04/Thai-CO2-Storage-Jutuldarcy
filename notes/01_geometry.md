# Step 1: Faulted cross-section

## Model
A 2D vertical section of 120 × 1 × 60 cells (10 m × 100 m × 2 m), 7,200 cells in total. The top is at 1,500 m depth and the beds dip at 2°, rising towards the fault, so the footwall forms a tilted fault-block trap. The fault is vertical, at column 70, and the hanging wall (right) is downthrown by 30 m (15 cells).

| Unit | k (mD) | kv/kh | φ |
|---|---|---|---|
| Overburden sand | 200 | 0.1 | 0.22 |
| Caprock shale (20 m) | 0.001 | 1 | 0.05 |
| Storage sand (40 m) | 300 | 0.1 | 0.25 |
| Underburden shale | 0.001 | 1 | 0.05 |
| Fault zone (1 column) | k_across = 0.01, k_along = 0.1 | | 0.10 |

## Choices
I used a 2D slice so that a sweep over fault properties runs in minutes on an 8 GB laptop. The price is that CO2 cannot go around the fault along strike.

The beds rise towards the fault because buoyant CO2 migrates up-dip; this geometry delivers it to the fault, which is the case where fault properties decide containment.

The throw is represented by shifting the facies on the hanging-wall side rather than by a non-matching grid. This gives the right juxtaposition (which rock faces which across the fault), although the rock is not physically displaced in depth. At 30 m of throw and a 2° dip that is a secondary effect.

The throw (30 m) is larger than the caprock thickness (20 m), which leaves a 10 m window where storage sand faces overburden sand across the fault. That is the cross-fault leak path.

The fault is its own column with separate permeabilities across and along it. A transmissibility multiplier can only reduce flow across the fault; a fault column can also carry CO2 up the fault.

The fault zone is 10 m wide, much wider than a real fault core (often under 1 m), so its permeabilities are equivalent values. For a core of thickness t inside a cell of width w, flow across the fault in series gives k_eq = w / (t/k_core + (w − t)/k_damage).

## Simplifications to revisit
There is no capillary pressure yet, so the caprock and the fault seal only through low permeability. In reality capillary entry pressure is often what seals. My plan for Step 2b: Brooks–Corey capillary pressure per facies, with entry pressure scaled by Leverett J (Pe ∝ √(φ/k)) from Pe ≈ 0.05 bar in the storage sand. That gives roughly 12 bar for the shales and 5 bar for the fault rock (using k_across).

A seal holds a CO2 column of height h = Pe / (Δρ g). With Δρ ≈ 500 kg/m³, 1 bar holds about 20 m. The footwall trap here has roughly 24 m of structural relief (695 m × tan 2°), so the fault entry pressure matters.

The fault is also vertical and planar, with no damage-zone gradient, and the model is isothermal.
