using Test
using JultraDark
using NPZ
import MPI
import PencilFFTs

for grids_type in [Grids, PencilGrids]
    @testset "Can initialize $grids_type" begin
        resol = 8

        grids = grids_type(1.0, resol)
        @test typeof(grids) == grids_type
        @test all(grids.ψx .== 0)
        @test reshape(grids.x, resol) == reshape(grids.y, resol)
        @test reshape(grids.x, resol) == reshape(grids.z, resol)

        grids = grids_type((1.0, 1.0, 1.0), (resol, resol, resol))
        @test typeof(grids) == grids_type
        @test all(grids.ψx .== 0)
        @test reshape(grids.x, resol) == reshape(grids.y, resol)
        @test reshape(grids.x, resol) == reshape(grids.z, resol)

        grids = grids_type((4.0, 2.0, 1.0), (4resol, 2resol, 1resol))
        @test typeof(grids) == grids_type
        @test all(grids.ψx .== 0)
        @test size(grids.x) == (4resol, 1, 1)
        @test size(grids.y) == (1, 2resol, 1)
        @test size(grids.z) == (1, 1, 1resol)
        @test reshape(grids.x, 4resol)[1 + 4resol÷2 - 4resol÷8:4resol÷2 + 4resol÷8] == reshape(grids.z, resol)
        @test reshape(grids.y, 2resol)[1 + 2resol÷2 - 2resol÷4:2resol÷2 + 2resol÷4] == reshape(grids.z, resol)
    end
end

@testset "Consistent initialisation" begin
    resol = 16
    dir = mktempdir()

    grids = Grids(1.0, resol)
    pencil_grids = PencilGrids(1.0, resol)

    gathered_pg = PencilFFTs.gather(pencil_grids.k)

    if MPI.Comm_rank(pencil_grids.MPI_COMM) == 0
        @test all(grids.k .== gathered_pg)

        npzwrite(joinpath(dir, "k_grids.npy"), grids.k)

        npzwrite(joinpath(dir, "k_pencilgrids.npy"), gathered_pg)

        @test all(npzread(joinpath(dir, "k_grids.npy")) .== npzread(joinpath(dir, "k_pencilgrids.npy")))
    else
        @test size(gathered_pg) == 0
    end

end
