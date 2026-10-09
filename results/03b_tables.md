## Boundary transport and CO2 conservation (tonnes)

| Case | Year | Injected | In reservoir | In well | Out via top (free / dissolved) | Out via sides | Residual | Max abs residual | Ministeps (cut) | Newton its | Wall (min) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| flux_pe5.5_leaky | 25 | 37500.00 | 37499.63 | 0.275 | 0.000 / 0.672 | 0.000 | -0.574 | 3.169 | 639 (188) | 5306 | 3.2 |
| flux_pe5.5_leaky | 100 | 37500.00 | 37502.45 | 0.028 | 0.000 / 0.686 | 0.000 | -3.169 | 3.169 | 639 (188) | 5306 | 3.2 |
| flux_pe5.5_leaky_tight | 25 | 37500.00 | 37499.08 | 0.275 | 0.000 / 0.664 | 0.000 | -0.022 | 0.342 | 669 (195) | 5983 | 2.6 |
| flux_pe5.5_leaky_tight | 100 | 37500.00 | 37499.63 | 0.028 | 0.000 / 0.678 | 0.000 | -0.336 | 0.342 | 669 (195) | 5983 | 2.6 |
| flux_pe1.2_mid | 25 | 37500.00 | 37499.75 | 0.276 | 0.000 / 0.000 | 0.000 | -0.027 | 1.299 | 641 (193) | 5421 | 2.3 |
| flux_pe1.2_mid | 100 | 37500.00 | 37501.27 | 0.028 | 0.000 / 0.000 | 0.000 | -1.299 | 1.299 | 641 (193) | 5421 | 2.3 |
| flux_pe1.2_mid_tight | 25 | 37500.00 | 37499.73 | 0.276 | 0.000 / 0.000 | 0.000 | -0.004 | 0.322 | 676 (204) | 6093 | 2.3 |
| flux_pe1.2_mid_tight | 100 | 37500.00 | 37500.29 | 0.028 | 0.000 / 0.000 | 0.000 | -0.322 | 0.322 | 676 (204) | 6093 | 2.3 |

## Fault invasion vs crossing the seal (tonnes)

| Case | Year | Free CO2 in fault | Dissolved CO2 in fault | Free in overburden | Dissolved in overburden | Cumulative net crossing of the seal | Net out of storage formation |
|---|---|---|---|---|---|---|---|
| flux_pe5.5_leaky | 25 | 0.0 | 32.7 | 0.0 | 50.3 | 60.2 | 84.6 |
| flux_pe5.5_leaky | 100 | 0.0 | 31.6 | 0.0 | 50.2 | 57.3 | 80.5 |
| flux_pe1.2_mid | 25 | 127.8 | 105.2 | 3.9 | 121.5 | 145.0 | 360.3 |
| flux_pe1.2_mid | 100 | 103.2 | 118.6 | 0.0 | 139.6 | 162.5 | 361.8 |

## Injection rate at fixed total mass (37,500 t), fault Pe 1.4 bar, k 0.1/1 mD

| Rate (t/yr) | Injected by simulator (t) | Peak entry drive (bar) | Peak Δp (bar) | Out of storage: end of injection / +25 yr / yr 100 (t) | Crossed seal: end / +25 / 100 (t) |
|---|---|---|---|---|---|
| 750 × 50 yr | 37500.0 | 1.42 | 1.04 | 16.5 / 16.9 / 15.7 | 2.4 / 2.4 / 2.4 |
| 1500 × 25 yr | 37500.0 | 1.65 | 1.46 | 112.7 / 113.1 / 112.0 | 16.4 / 17.6 / 17.8 |
| 2500 × 15 yr | 37500.0 | 1.92 | 1.93 | 391.5 / 395.1 / 393.0 | 158.4 / 186.0 / 179.8 |

## Time-step and grid sensitivity (fault k 0.1/1 mD)

| Case | Cells | Out of storage, yr 25 / 100 (t) | Crossed seal, yr 25 / 100 (t) | Free CO2 in fault, yr 100 (t) | Max abs residual (t) | Ministeps (cut) | Wall (min) |
|---|---|---|---|---|---|---|---|
| dt90_pe1.2_mid | 8040 | 360.3 / 361.8 | 145.0 / 162.5 | 103.2 | 1.299 | 641 (193) | 2.9 |
| dt30_pe1.2_mid | 8040 | 359.5 / 363.0 | 143.5 / 161.2 | 103.1 | 2.609 | 1398 (183) | 2.5 |
| grid10x2_pe1.2_mid | 8040 | 360.3 / 361.8 | 145.0 / 162.5 | 103.2 | 1.299 | 641 (193) | 2.9 |
| grid10x1_pe1.2_mid | 16080 | 369.4 / 369.0 | 158.2 / 176.8 | 105.9 | 5.240 | 1036 (749) | 11.1 |
| grid5x2_pe1.2_mid | 15240 | 281.5 / 284.0 | 100.8 / 118.0 | 88.0 | 2.187 | 680 (243) | 4.7 |
| grid5x1_pe1.2_mid | 30480 | 296.5 / 297.8 | 115.2 / 134.0 | 92.1 | 2.886 | 1662 (1750) | 46.5 |
| grid10x2_pe0.8_mid | 8040 | 1721.5 / 2988.5 | 1403.5 / 2714.0 | 282.8 | 5.969 | 636 (179) | 2.1 |
| grid5x1_pe0.8_mid | 30480 | 1812.2 / 3600.3 | 1501.2 / 3350.0 | 265.9 | 14.366 | 1675 (1774) | 109.4 |
