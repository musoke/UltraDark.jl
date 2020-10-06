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

grids = Grids(1.0, 16)
grids.ψx .= grids.dist ./ 1e9 # Set ψx to something non-zero

output_dir = "output"
output_times = 0.0:0.01:0.1

output_config = OutputConfig(output_dir, output_times; box=false)
options = Config.SimulationConfig(10, t->1)

@test simulate(grids, options, output_config) == nothing
