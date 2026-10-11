# Upstream issue draft: NOT TO BE FILED

**Status (11 Oct 2026):** the discrepancy is already fixed upstream. Do not file an issue.

| MultiComponentFlash version | Third term of `soreide_whitson_hc_aqueous_bic` |
|---|---|
| 1.1.19 (this project, required by JutulDarcy 0.3.10) | `A_2*T_ri*(1 + α_2*c_sw)` |
| 1.1.20 | `A_2*T_ri*(1 + α_2*c_sw)` |
| 2.0.1 | `A_2*T_ri^2*(1 + α_2*c_sw)` |

The change came in pull request #31 ("Modernize docs + EOS fixes + Clapeyron comparisons"), merged 13 Sep 2026, with the commit message "Fixes from audit". Checked by reading the source of each tagged version on GitHub (moyner/MultiComponentFlash.jl).

MultiComponentFlash 2.x requires JutulDarcy 0.4 and Jutul 0.5, so the fix is not reachable from this project's JutulDarcy 0.3.10 environment without a major upgrade.

**Accurate description of what this project found:** in Step 4A we independently detected a discrepancy between the installed (older) release and the published Søreide–Whitson equation (T_r instead of T_r²), quantified its effect (about 10× lower CH4 solubility at 95 °C, 200 bar), and then identified that it had already been corrected upstream. We did not discover a new bug in the current release.

Lesson: check the latest release and changelog before suggesting a defect.

The separate minor point (`mixture_compressibility_factor(eos, cond)` failing for SW because its default argument builds a single coefficient set) has not been checked against 2.0.1.
