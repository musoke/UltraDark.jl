using BenchmarkTools
using JultraDark
using JultraDark: psi_half_step!, psi_whole_step!, phi_whole_step!, take_steps!
using MPI

if ~MPI.Initialized()
    MPI.Init()
end

n = 64
grids = PencilGrids(1.0, n)

Δt = 0.1

suite = BenchmarkGroup()

suite["psi_half_step"] = @benchmarkable psi_half_step!($Δt, $grids)
suite["psi_whole_step"] = @benchmarkable psi_whole_step!($Δt, $grids)
suite["phi_whole_step"] = @benchmarkable phi_whole_step!($Δt, $grids)

suite["take_steps"] = @benchmarkable take_steps!($grids, $1.0, $Δt, 1, $(OutputConfig("output", 1:2)), $(t -> 1))

tune!(suite)

results = run(suite)

@show res_min = minimum(results)
@show res_med = median(results)
