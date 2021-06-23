using MPI
using UltraDark

MPI.Init()
comm = MPI.COMM_WORLD

println("Process ID $(MPI.Comm_rank(comm)+1) of $(MPI.Comm_size(comm))\n")

include(joinpath(@__DIR__, "..", "init_soliton.jl"))

output_dir = "output/vel"

num_snapshots = 20
output_times = LinRange(0, 5, num_snapshots)

output_config = OutputConfig(output_dir, output_times; box=true, slice=false)
options = Config.SimulationConfig(10)

mass = 10
position0 = [-4, 0, 0]
velocity = [1, 0, 0]
phase = 0
t0 = 0

resol = 64
grids = PencilGrids(10.0, resol)
add_soliton(grids, mass, position0, velocity, phase, t0)

simulate(grids, options, output_config) == nothing
