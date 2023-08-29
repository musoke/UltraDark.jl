using UltraDark
using Test
using CSV
using HDF5
using MPI
using NPZ
using PencilFFTs

resol = 64

include("../../examples/init_soliton.jl")

output_dir = mktempdir()

for a in [UltraDark.constant_scale_factor, t -> t]

    this_output_dir = joinpath(output_dir, "$a")
    output_times = [1.0, 2.0]

    output_config = OutputConfig(
        this_output_dir,
        output_times;
        box = true,
        slice = true,
        npy = true,
        summary_statistics = (Summary.SimulationTime, Summary.ScaleFactor),
    )
    options = Config.SimulationConfig(a = a)

    grids = Grids(10.0, resol)

    mass = 10
    position = [0, 0, 0]
    velocity = [0, 0, 0]
    phase = 0
    t0 = 0

    add_soliton(grids, mass, position, velocity, phase, t0)

    @test simulate!(grids, options, output_config) == nothing

    summary = CSV.File("$(output_config.directory)/summary.csv")

    @test summary.a == a.(summary.t)
end
