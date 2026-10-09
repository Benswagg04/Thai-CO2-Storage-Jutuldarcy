# Step 2: baseline, capillary pressure, and verification runs.
# Run from the project folder:
#   julia --project=. scripts/02_run_cases.jl checks    # regression + Pc-cap check
#   julia --project=. scripts/02_run_cases.jl main      # three cases, current set-up
#   julia --project=. scripts/02_run_cases.jl fine      # grid-refinement case
#   julia --project=. scripts/02_run_cases.jl all

include(joinpath(@__DIR__, "..", "src", "run.jl"))

group = isempty(ARGS) ? "all" : ARGS[1]
outdir = joinpath(@__DIR__, "..", "output")
sandy = (k_across = 1.0, k_along = 10.0)   # sand-rich fault

# A. Checks against the first configuration (no padding, closed top), whose
#    results are saved in output/02_baseline.jls and 02_capillary_sandyfault.jls.
if group in ("all", "checks")
    first_config(; kw...) = build_section(section_spec(; pad_cells = 0, kw...))

    # Regression: rewritten geometry code, same configuration → same answer.
    old = deserialize(joinpath(outdir, "02_baseline.jls"))
    new = run_case("02_check_regression", first_config(), case_spec(open_top = false))
    for k in keys(old.masses)
        @printf("Regression %-12s old %9.1f  new %9.1f t\n", k, old.masses[k][end], new.acc.total[k][end])
    end

    # Capillary cap: only the Pc table changes (cap 10 Pe → 32 Pe, denser table).
    old = deserialize(joinpath(outdir, "02_capillary_sandyfault.jls"))
    new = run_case("02_check_pccap", first_config(; sandy...), case_spec(capillary = true, open_top = false))
    for k in keys(old.masses)
        @printf("Pc cap     %-12s old %9.1f  new %9.1f t\n", k, old.masses[k][end], new.acc.total[k][end])
    end
    @printf("Pc cap     %-12s old %9.1f  new %9.1f t\n", "left model",
        old.inventory[end][:outside_domain] / 1000, new.exited[end])
end

# B. Current set-up: 2.5 km padding each side, open top.
if group in ("all", "main")
    run_case("02_nopc_clay", build_section(), case_spec())
    run_case("02_pc_clay", build_section(), case_spec(capillary = true))
    run_case("02_pc_sandy", build_section(section_spec(; sandy...)), case_spec(capillary = true))
end

# C. Grid refinement: same geology, cells halved in x and z (4× as many).
if group in ("all", "fine")
    run_case("02_pc_sandy_fine", build_section(section_spec(; dx = 5.0, dz = 1.0, sandy...)),
        case_spec(capillary = true))
end
