using BenchmarkTools
using Statistics
using JultraDark
using JultraDark: take_steps!

nthreads = Threads.nthreads()

resol = try
    parse(Int64, ARGS[1])
catch
    @warn "Resolution not supplied, using default"
    64
end

Δt = 0.1
n_steps = 10

grids = Grids(1.0, resol)

# b = @benchmark take_steps!($grids, $1.0, $Δt, 10, $(OutputConfig("output", 1:2)), $(t -> 1))
# time_mean = mean(b.times / 1e9 / n_steps)

res = @timed take_steps!(grids, 1.0, Δt, n_steps, OutputConfig("output", 1:2), t -> 1)

time_mean = mean(res[2] / n_steps)

println("$nthreads, $resol, $time_mean")
