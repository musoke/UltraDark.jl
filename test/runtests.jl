using JultraDark
using Test, Documenter
using SafeTestsets
using MPI

MPI.Init()

doctest(JultraDark)

@safetestset "grids.jl" begin
    include("grids.jl")
end

@testset "Phase grad Grids" begin
    include("phase_gradient.jl")
end

@safetestset "Evolution" begin
    include("evolution.jl")
end

# Run examples
@safetestset "full sim" begin
    include("sims/full.jl")
end

# Put notebook in module to emulate SafeTestsets
# https://github.com/YingboMa/SafeTestsets.jl/issues/3
module GrowthNotebook
    using NBInclude
    @nbinclude("../examples/growth.ipynb")
end
