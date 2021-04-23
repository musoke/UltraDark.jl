using UltraDark
using Test, Documenter
using SafeTestsets
using MPI

doctest(UltraDark)

@safetestset "grids.jl" begin
    include("grids.jl")
end

@safetestset "output.jl" begin
    include("output.jl")
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

@safetestset "soliton" begin
    include("sims/soliton_static.jl")
end

@safetestset "soliton" begin
    include("../examples/soliton_velocity.jl")
end

# Put notebook in module to emulate SafeTestsets
# https://github.com/YingboMa/SafeTestsets.jl/issues/3
module GrowthNotebook
    using NBInclude
    @nbinclude("../examples/growth.ipynb")
end

module TwoDNotebook
    using NBInclude
    @nbinclude("../examples/2d.ipynb")
end
