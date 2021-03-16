using JultraDark

for grids_type in [Grids, PencilGrids]
    grids = grids_type(zeros(Complex{Float64}, 16, 16, 16), 1)
    @test typeof(grids) == grids_type

    grids = grids_type(1.0, 4)
    @test typeof(grids) == grids_type
end
