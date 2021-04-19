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
