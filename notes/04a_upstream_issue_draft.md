# DRAFT: possible issue for sintefmath/MultiComponentFlash.jl (not submitted)

To be posted only after confirming against Søreide & Whitson (1992), eq. 12, and checking the latest MultiComponentFlash release in case it is already fixed.

---

**Title:** Søreide–Whitson hydrocarbon–water aqueous BIP: third term uses T_r instead of T_r²?

**Version:** MultiComponentFlash 1.1.19 (via JutulDarcy 0.3.10), Julia 1.12.6

In `src/eos/soreide_whitson.jl`, `soreide_whitson_hc_aqueous_bic` returns

```julia
A_0*(1 + α_0*c_sw) + A_1*T_ri*(1 + α_1*c_sw) + A_2*T_ri*(1 + α_2*c_sw)
```

Søreide & Whitson (1992), eq. 12 (as also written in generalised form in arXiv:2506.02582, eq. 18) has the third term proportional to T_r²:

  k_ij^AQ = A0(1 + α0 c_sw) + A1 T_ri (1 + α1 c_sw) + A2 T_ri² (1 + α2 c_sw)

For CH4 at 368.15 K (T_r ≈ 1.93), fresh water, 200 bar, with vapour BIP H2O–CH4 = 0.485:
- as implemented: k_aq ≈ +0.211, x_CH4 → 0.0106 mol/kg
- with T_r²: k_aq ≈ −0.095, x_CH4 → 0.108 mol/kg (closer to typical literature values of ~0.1 mol/kg)

Separately (minor): `mixture_compressibility_factor(eos, cond)` with a `SoreideWhitson` EOS errors because the default `forces = force_coefficients(eos, cond)` builds a single coefficient set, while `get_force_coefficients` expects `forces.liquid`/`forces.vapor` for per-phase EOS.

A minimal reproduction script is available on request (it is in the project as `scripts/04a_flash_benchmark.jl`).
