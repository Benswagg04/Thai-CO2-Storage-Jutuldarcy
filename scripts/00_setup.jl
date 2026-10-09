# One-time environment setup: adds the plotting package to this project.
# Run from the project folder:  julia --project=. scripts/00_setup.jl
using Pkg
Pkg.add("CairoMakie")
Pkg.precompile()
using CairoMakie
println("CairoMakie OK ", pkgversion(CairoMakie))
