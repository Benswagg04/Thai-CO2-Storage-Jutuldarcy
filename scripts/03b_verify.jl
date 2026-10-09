# Step 3b: verification of the Step 3 results.
# Run from the project folder:
#   julia --project=. scripts/03b_verify.jl flux    # boundary flux + conservation residual (+ tight tolerances)
#   julia --project=. scripts/03b_verify.jl rate    # injection rate at fixed total mass
#   julia --project=. scripts/03b_verify.jl dt      # time-step sensitivity
#   julia --project=. scripts/03b_verify.jl grid    # directional grid refinement (slow)
# Restartable: finished runs are skipped.

include(joinpath(@__DIR__, "..", "src", "verify.jl"))

group = isempty(ARGS) ? "flux" : ARGS[1]
outdir = joinpath(@__DIR__, "..", "output", "03b")
mkpath(outdir)
done(name) = isfile(joinpath(outdir, "$name.jls"))
kmid = (k_across = 0.1, k_along = 1.0)
kleaky = (k_across = 1.0, k_along = 10.0)
tight = (tol_mb = 1e-9, tol_cnv = 1e-4, tol_mb_well = 1e-9)

function go(f, name)
    done(name) && return println("skip $name (done)")
    try
        f()
    catch err
        println("FAILED $name: ", sprint(showerror, err)[1:min(end, 400)])
    end
    GC.gc()
end

if group == "hookcheck"
    # The flux runs were made by storing every intermediate state; later runs
    # accumulate fluxes in a post-ministep hook instead (needed for memory on
    # the fine grids). Rerun one flux case with the hook to confirm both
    # methods give the same answer.
    go("hookcheck_pe5.5_leaky") do
        r = verify_case("hookcheck_pe5.5_leaky", build_section(section_spec(; kleaky...)),
            case_spec(capillary = true, pe_fault_bar = 5.5), output_dir = outdir)
        old = deserialize(joinpath(outdir, "flux_pe5.5_leaky.jls"))
        for (k, lab) in [(:top_diss, "out via top, dissolved"), (:residual, "residual"),
                         (:crossed_seal, "crossed seal"), (:out_of_storage, "out of storage")]
            @printf("Hook vs stored states, year 100, %-24s %10.4f vs %10.4f t\n", lab, getfield(r, k)[end], getfield(old, k)[end])
        end
    end
elseif group == "flux"
    # Sealing fault with the largest dissolved-only leakage, and a case in the
    # injection-triggered window. Default and tight solver tolerances.
    for (name, k, pe, kw) in [
            ("flux_pe5.5_leaky", kleaky, 5.5, (;)),
            ("flux_pe5.5_leaky_tight", kleaky, 5.5, tight),
            ("flux_pe1.2_mid", kmid, 1.2, (;)),
            ("flux_pe1.2_mid_tight", kmid, 1.2, tight)]
        go(name) do
            verify_case(name, build_section(section_spec(; k...)),
                case_spec(capillary = true, pe_fault_bar = pe), sim_kw = kw, output_dir = outdir)
        end
    end
elseif group == "rate"
    # Same total mass (37,500 t), different rate and duration, 100 years total.
    for (rate, years) in [(750.0, 50), (1500.0, 25), (2500.0, 15)]
        name = "rate$(Int(rate))_pe1.4_mid"
        go(name) do
            verify_case(name, build_section(section_spec(; kmid...)),
                case_spec(capillary = true, pe_fault_bar = 1.4, rate_t_per_yr = rate,
                    inject_years = years, post_years = 100 - years), output_dir = outdir)
        end
    end
elseif group == "dt"
    # Maximum time step 90 days (default) vs 30 days.
    for (name, mx) in [("dt90_pe1.2_mid", 90day), ("dt30_pe1.2_mid", 30day)]
        go(name) do
            verify_case(name, build_section(section_spec(; kmid...)),
                case_spec(capillary = true, pe_fault_bar = 1.2), sim_kw = (max_timestep = mx,), output_dir = outdir)
        end
    end
elseif group == "grid"
    # Fault stays 10 m wide in metres; only the cell size changes.
    for (name, pe, dx, dz) in [
            ("grid10x2_pe1.2_mid", 1.2, 10.0, 2.0),
            ("grid10x1_pe1.2_mid", 1.2, 10.0, 1.0),
            ("grid5x2_pe1.2_mid", 1.2, 5.0, 2.0),
            ("grid5x1_pe1.2_mid", 1.2, 5.0, 1.0),
            ("grid10x2_pe0.8_mid", 0.8, 10.0, 2.0),
            ("grid5x1_pe0.8_mid", 0.8, 5.0, 1.0)]
        go(name) do
            verify_case(name, build_section(section_spec(; dx = dx, dz = dz, kmid...)),
                case_spec(capillary = true, pe_fault_bar = pe), output_dir = outdir)
        end
    end
end
