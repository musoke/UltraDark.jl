using BenchmarkTools
using MPI
using Statistics
using UltraDark
using UltraDark: take_steps!

const DEV_NULL = @static Sys.iswindows() ? "nul" : "/dev/null"

nthreads = Threads.nthreads()

resol = try
    parse(Int64, ENV["ULTRADARK_RESOL"])
catch
    @warn "Resolution not supplied, using default"
    64
end

grids_type = try
    eval(Meta.parse(ENV["ULTRADARK_GRIDS_TYPE"]))
catch
    @warn "ULTRADARK_GRIDS_TYPE not supplied, using default"
    Grids
end

ntasks = if grids_type <: PencilGrids
    if ~MPI.Initialized()
        MPI.Init()
    end
    MPI.Comm_size(MPI.COMM_WORLD)
else
    1
end

Δt = 0.1
n_steps = 10

grids = grids_type(1.0, resol)

output_config = OutputConfig(mktempdir(), 1:2)

# Run once to ensure functions are precompiled
take_steps!(grids, 1.0, Δt, n_steps, output_config, Config.constant_scale_factor, nothing)

# Collect data
res = @timed take_steps!(
    grids,
    1.0,
    Δt,
    n_steps,
    output_config,
    Config.constant_scale_factor,
    nothing,
)

time_mean = mean(res[2] / n_steps)

if grids_type <: PencilGrids
    # Find maximum time across MPI jobs
    time_mean = MPI.Allreduce(time_mean, MPI.MAX, grids.MPI_COMM)

    # Disable output on all but one process.
    MPI.Comm_rank(grids.MPI_COMM) == 0 || redirect_stdout(open(DEV_NULL, "w"))
end

println("$ntasks,$nthreads,$resol,$time_mean,$grids_type")
