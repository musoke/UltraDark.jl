using JultraDark

using JultraDark

@testset "Initialise Grids" begin
    grids = Grids(zeros(Complex{Float64}, 16, 16, 16), 1)
    @test typeof(grids) == Grids

    grids = Grids(1.0, 4)
    @test typeof(grids) == Grids

end