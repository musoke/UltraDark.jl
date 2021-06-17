using Base.Threads
using BenchmarkTools
using LinearAlgebra
using LoopVectorization
using MPI
using Strided
using UltraDark: phi_whole_step!, Grids, PencilGrids

if ~MPI.Initialized()
    MPI.Init()
end

n = 128

Δt = 1.0

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
suite = BenchmarkGroup()

for (label, grids) in [("Grids", Grids(1.0, n)), ("PencilGrids", PencilGrids(1.0, n))]

    grids.ψx .= randn(size(grids.ψx))
    grids.ψk .= randn(size(grids.ψk))

    suite[label] = BenchmarkGroup()

    for f in [
              phi_whole_step!,
             ]
        suite[label]["$f"] = @benchmarkable $f($Δt, $grids, $nothing)
    end

end

@info "Tuning"
@show tune!(suite)

@info "Running benchmarks"
results = run(suite)

@info "Results"
@show res_min = minimum(results)
# @show res_med = median(results)

@show judge(res_min["Grids"], res_min["PencilGrids"])
