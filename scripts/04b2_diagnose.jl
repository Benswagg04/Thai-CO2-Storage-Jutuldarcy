# Step 4B, stage 2: why are steps cut in the post-injection shut-in?
# Runs the gravity-on case with relaxation + step growth ≤ 2× (the two settings
# that helped), then reports cuts per year and checks whether a hydrocarbon
# "liquid" phase ever appears (phase-label switching would destabilise Newton).
# Run from the project folder:  julia --project=. scripts/04b2_diagnose.jl

include(joinpath(@__DIR__, "04b2_gas_2d.jl"))

r = gravity_case(true; sim_kw = (relaxation = true, timestep_max_increase = 2.0))
@printf("Relaxation + growth ≤ 2×: %.1f min; ministeps %d ok / %d cut; Newton %d\n", r.wall, r.n_ok, r.n_cut, r.newton)
println("Cuts per year: ", join(["y$y:$c" for (y, c) in enumerate(r.cuts_per_year) if c > 0], " "))

sl = [maximum(s[:Saturations][2, :]) for s in r.states]
n_liq = [count(>(1e-8), s[:Saturations][2, :]) for s in r.states]
@printf("Hydrocarbon-liquid saturation: max over run %.3e; max number of cells with S_l > 1e-8: %d\n",
    maximum(sl), maximum(n_liq))
if maximum(n_liq) > 0
    first_i = findfirst(>(0), n_liq)
    @printf("First liquid appears at t = %.2f yr (ministep %d)\n", r.t[first_i] / yr, first_i)
end
# Where is the gas densest? (CO2-rich gas near the contact)
s = r.states[end]
@printf("Final gas density range: %.1f–%.1f kg/m³; brine %.1f kg/m³\n",
    minimum(s[:PhaseMassDensities][3, s[:Saturations][3, :] .> 0.01]),
    maximum(s[:PhaseMassDensities][3, s[:Saturations][3, :] .> 0.01]),
    maximum(s[:PhaseMassDensities][1, :]))
