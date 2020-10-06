using JultraDark

using JultraDark

@testset "Initialise PencilGrids" begin
    grids = PencilGrids(zeros(Complex{Float64}, 16, 16, 16), 1)
    @test typeof(grids) == PencilGrids

    grids = PencilGrids(1.0, 4)
    @test typeof(grids) == PencilGrids

end
