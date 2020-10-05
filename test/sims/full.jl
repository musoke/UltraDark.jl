using JultraDark
using Test
using MPI

if ~MPI.Initialized()
    MPI.Init()
end

comm = MPI.COMM_WORLD
# Disable output on all but one process.
const DEV_NULL = @static Sys.iswindows() ? "nul" : "/dev/null"
rank = MPI.Comm_rank(comm)
rank == 0 || redirect_stdout(open(DEV_NULL, "w"))

resol = 16

grids = JultraDark.Grids(zeros(Complex{Float64}, resol, resol, resol), 1)

output_dir = "output"
output_times = 0:10

output_config = OutputConfig(output_dir, output_times; box=false)
options = Config.SimulationConfig(10, t->1)

@test simulate(grids, options, output_config) == nothing
