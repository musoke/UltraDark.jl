using JultraDark
using Test, Documenter

doctest(JultraDark)

@testset "grids.jl" begin
    grids = Grids(zeros(Complex{Float64}, 16, 16, 16), 1)

    @test typeof(grids) == Grids
end

@testset "Evolution" begin
    include("evolution.jl")
end
@testset "full sim" begin
    @test simulate() == nothing
end