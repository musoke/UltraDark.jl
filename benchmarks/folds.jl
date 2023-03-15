using BenchmarkTools
import Folds
using Random
using Strided
using UltraDark

Random.seed!(0)

n = 128

grids = Grids(1.0, n)
grids.ρx .= rand(size(grids.ρx)...)

suite = BenchmarkGroup()

suite["Base"] = BenchmarkGroup([])
suite["Folds.jl"] = BenchmarkGroup([])

suite["Base"]["maximum"] = @benchmarkable maximum($(grids.ρx))
suite["Folds.jl"]["maximum"] = @benchmarkable Folds.maximum($(grids.ρx))

suite["Base"]["sum"] = @benchmarkable sum($(grids.ρx))
suite["Folds.jl"]["sum"] = @benchmarkable Folds.sum($(grids.ρx))


tune!(suite)

results = run(suite)

@show res_min = minimum(results)
@show res_med = median(results)

judge(res_med["Folds.jl"], res_med["Base"])
