"This program should install everything required to run SeitComp-1D in your system"
## Lets begin making your computer ready
# In linux, remember to run "sudo julia" for this code to work!
using Pkg
## Install all required Packages
Pkg.add("Plots")
Pkg.add("DelimitedFiles")
Pkg.add("Interpolations")
Pkg.add("JLD2")
Pkg.add("ScatteredInterpolation")
Pkg.add("LinearAlgebra")
Pkg.add("StructArrays")
Pkg.add("FFTW")
Pkg.add("Statistics")
Pkg.add("QuadGK")
Pkg.add("PolynomialRoots")
Pkg.add("Arpack")
Pkg.add("SparseArrays")
Pkg.add(Pkg.PackageSpec(;name="Interpolations", version="0.13.1"));
ENV["GRDIR"] = ""
Pkg.build("GR")
ENV["GKSwstype"] = "nul"


println(" ")
printstyled("    All should be installed by now!           ",color=:light_green)
println(" ")
printstyled("    Unfortunatly I can't install GMT ",color=:light_red)
println(" ")
printstyled("    Go to julia and hit ] to access the package manager          ",color=:light_blue)
println(" ")
printstyled("    write: add GMT#master and hit RETURN          ",color=:light_blue)
println(" ")
printstyled("    is should install GMT.jl and update your gmt           ",color=:light_blue)
println(" ")
