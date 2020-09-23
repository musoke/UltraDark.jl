using JultraDark
using Test, Documenter
using SafeTestsets

doctest(JultraDark)

@safetestset "grids.jl" begin
    using JultraDark
    grids = Grids(zeros(Complex{Float64}, 16, 16, 16), 1)

    @test typeof(grids) == Grids
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
