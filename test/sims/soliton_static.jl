using UltraDark
using Test
using MPI
using NPZ
using PencilFFTs

resol = 64

include("../../examples/init_soliton.jl")

output_dir = mktempdir()

for grid_type in [Grids, PencilGrids]

    this_output_dir = joinpath(output_dir, "$grid_type")
    output_times = [0, 1]

    output_config = OutputConfig(this_output_dir, output_times; box=true, slice=false)
    options = Config.SimulationConfig(10)

    grids = grid_type(10.0, resol)

    mass = 10
    position = [0, 0, 0]
    velocity = [0, 0, 0]
    phase = 0
    t0 = 0

    add_soliton(grids, mass, position, velocity, phase, t0)

    @test simulate(grids, options, output_config) == nothing
end

# for grid_type in [Grids, PencilGrids]
#     output_initial = npzread(joinpath(output_dir, "$grid_type", "psi_1.npy"))
#     output_final = npzread(joinpath(output_dir, "$grid_type", "psi_2.npy"))

#     @test all(output_initial .≈ output_final)
# end

for i in [1, 2]
    grids_output = npzread(joinpath(output_dir, "Grids", "psi_$i.npy"))
    pencil_grids_output = npzread(joinpath(output_dir, "PencilGrids", "psi_$i.npy"))

    # Fail for multiprocess - there are order 1e-12 differences
    @test all(grids_output .≈ pencil_grids_output)
end
