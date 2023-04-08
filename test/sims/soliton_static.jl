using UltraDark
using Test
using HDF5
using MPI
using NPZ
using PencilFFTs

resol = 64

include("../../examples/init_soliton.jl")

output_dir = mktempdir()

for grid_type in [Grids, PencilGrids]

    this_output_dir = joinpath(output_dir, "$grid_type")
    output_times = [0, 1]

    output_config = OutputConfig(
        this_output_dir,
        output_times;
        box = true,
        slice = true,
        npy = true,
        h5 = grid_type == Grids,
    )
    options = Config.SimulationConfig()

    grids = grid_type(10.0, resol)

    mass = 10
    position = [0, 0, 0]
    velocity = [0, 0, 0]
    phase = 0
    t0 = 0

    add_soliton(grids, mass, position, velocity, phase, t0)

    @test simulate!(grids, options, output_config) == nothing
end

# for grid_type in [Grids, PencilGrids]
#     output_initial = npzread(joinpath(output_dir, "$grid_type", "psi_1.npy"))
#     output_final = npzread(joinpath(output_dir, "$grid_type", "psi_2.npy"))

#     @test all(output_initial .≈ output_final)
# end

for i in [1, 2]
    ## Boxes
    grids_output_npy = npzread(joinpath(output_dir, "Grids", "psi_$i.npy"))
    grids_output_h5 = h5read(joinpath(output_dir, "Grids", "psi_$i.h5"), "psi")

    @test grids_output_npy == grids_output_h5

    pencil_grids_output_npy = npzread(joinpath(output_dir, "PencilGrids", "psi_$i.npy"))
    @test_throws ErrorException pencil_grids_output_h5 =
        h5read(joinpath(output_dir, "PencilGrids", "psi_$i.h5"), "psi")

    @test_broken all(pencil_grids_output_npy .≈ pencil_grids_output_h5)

    # Fail for multiprocess - there are order 1e-12 differences
    @test all(grids_output_npy .≈ pencil_grids_output_npy)

    ## Slices
    grids_slice_npy = npzread(joinpath(output_dir, "Grids", "psi_$i.npy"))
    grids_slice_h5 = h5read(joinpath(output_dir, "Grids", "psi_$i.h5"), "psi")

    @test grids_slice_npy == grids_slice_h5

    pencil_grids_slice_npy = npzread(joinpath(output_dir, "PencilGrids", "psi_$i.npy"))
    @test_throws ErrorException pencil_slice_output_h5 =
        h5read(joinpath(output_dir, "PencilGrids", "psi_$i.h5"), "psi")

    @test_broken all(pencil_grids_output_npy .≈ pencil_grids_output_h5)

    # Fail for multiprocess - there are order 1e-12 differences
    @test all(grids_slice_npy .≈ pencil_grids_slice_npy)

end
