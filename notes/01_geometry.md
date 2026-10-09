# Step 1: Faulted cross-section, design choices

## Model
- 2D vertical section, 120 × 1 × 60 cells (10 m × 100 m × 2 m), 7,200 cells.
- Top at 1,500 m depth, beds tilted 2° and rising toward the fault (tilted fault-block trap).
- Vertical fault at column 70; hanging wall (right) downthrown 30 m (15 cells).

| Unit | k (mD) | kv/kh | φ |
|---|---|---|---|
| Overburden sand | 200 | 0.1 | 0.22 |
| Caprock shale (20 m) | 0.001 | 1 | 0.05 |
| Storage sand (40 m) | 300 | 0.1 | 0.25 |
| Underburden shale | 0.001 | 1 | 0.05 |
| Fault zone (1 column) | k_across = 0.01, k_along = 0.1 | | 0.10 |

## Why these choices
- **2D slice:** thousands of cells instead of millions, so a fault-property sweep runs in minutes on an 8 GB laptop. The cost is that CO2 cannot go around the fault along strike.
- **Beds rising toward the fault:** buoyant CO2 migrates up-dip. This geometry delivers it to the fault, which is the case where fault properties decide containment.
- **Throw by shifting facies index:** gives the correct juxtaposition (which rock faces which across the fault) without a non-matching grid. The rock itself is not physically displaced in depth across the fault; at 30 m throw and 2° dip this is a secondary effect.
- **Throw (30 m) > caprock thickness (20 m):** opens a 10 m window where storage sand faces overburden sand across the fault. That is the cross-fault leak path.
- **Fault as its own column with anisotropic k:** a single transmissibility multiplier can only reduce cross-fault flow. A fault column with separate k_across and k_along can also carry CO2 *up* the fault.
- **Fault zone 10 m wide:** wider than a real fault core (often < 1 m). The cell permeabilities are *equivalent* values. For a core of thickness t inside a cell of width w, across-fault flow in series gives k_eq = w / (t/k_core + (w − t)/k_damage).

## Known simplifications (to revisit)
- No capillary pressure yet: the caprock and fault seal only through low permeability. In reality, capillary entry pressure is often what seals. Plan (Step 2b):
  - Brooks-Corey Pc per facies, with entry pressure scaled by Leverett J: Pe ∝ sqrt(φ/k), anchored at Pe ≈ 0.05 bar for the storage sand.
  - This gives roughly 12 bar for the shales and roughly 5 bar for the fault rock (using k_across).
  - Seal check: a seal holds a CO2 column of height h = Pe / (Δρ g). With Δρ ≈ 500 kg/m³, 1 bar holds about 20 m.
  - The footwall trap here has roughly 24 m of structural relief (695 m × tan 2°).
- Fault is vertical and planar; no damage-zone gradient.
- Isothermal.
