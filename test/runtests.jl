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
    using JultraDark
    
    resol = 16

    grids = JultraDark.Grids(zeros(Complex{Float64}, resol, resol, resol), 1)

    output_dir = "output"
    output_times = 0:10

    output_config = OutputConfig(output_dir, output_times)
    options = Config.SimulationConfig(10, t->1)

    @test simulate(grids, options, output_config) == nothing
end