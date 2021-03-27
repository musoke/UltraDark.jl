using Test
using JultraDark
using NPZ
import MPI
import PencilFFTs

@testset "Can initialize" begin
    for grids_type in [Grids, PencilGrids]
        grids = grids_type(zeros(Complex{Float64}, 16, 16, 16), 1)
        @test typeof(grids) == grids_type

        grids = grids_type(1.0, 4)
        @test typeof(grids) == grids_type
    end
end

@testset "Consistent initialisation" begin
    resol = 16
    dir = mktempdir()

    grids = Grids(1.0, resol)
    pencil_grids = PencilGrids(1.0, resol)

    gathered_pg = PencilFFTs.gather(pencil_grids.dist)

    if MPI.Comm_rank(pencil_grids.MPI_COMM) == 0
        @test all(grids.dist .== gathered_pg)

        npzwrite(dir * "/dist_grids.npy", grids.dist)

        npzwrite(dir * "/dist_pencilgrids.npy", gathered_pg)

        @test all(npzread(dir * "/dist_grids.npy") .== npzread(dir * "/dist_pencilgrids.npy"))
    else
        @test size(gathered_pg) == 0
    end

end
