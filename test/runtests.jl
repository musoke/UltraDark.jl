using JultraDark
using Test, Documenter
using SafeTestsets
using MPI

MPI.Init()

doctest(JultraDark)

@safetestset "grids.jl" begin
    include("grids.jl")
end

@safetestset "Evolution" begin
    include("evolution.jl")
end

@safetestset "full sim" begin
    include("sims/full.jl")
    include("../examples/growth.jl")
end
