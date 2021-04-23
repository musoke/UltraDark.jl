using BenchmarkTools
using UltraDark
using LinearAlgebra
using MPI

if ~MPI.Initialized()
    MPI.Init()
end

n = 64

suite = BenchmarkGroup()

function fft_broadcast!(grids)
    grids.ψk .= grids.fft_plan * grids.ψx
end

function fft_inplace!(grids)
    mul!(grids.ψk, grids.fft_plan, grids.ψx)
end

function fft_inv_broadcast!(grids)
    grids.ψx .= grids.fft_plan \ grids.ψk
end

function fft_inv_inplace!(grids)
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)
end

for (label, grids) in [("Grids", Grids(1.0, n)), ("PencilGrids", PencilGrids(1.0, n))]

    suite[label] = BenchmarkGroup()

    suite[label]["fft_broadcast!"] = @benchmarkable fft_broadcast!($grids)
    suite[label]["fft_inplace!"] = @benchmarkable fft_inplace!($grids)
    suite[label]["fft_inv_broadcast!"] = @benchmarkable fft_inv_broadcast!($grids)
    suite[label]["fft_inv_inplace!"] = @benchmarkable fft_inv_inplace!($grids)

end

@show tune!(suite)

results = run(suite)

@show res_min = minimum(results)
@show res_med = median(results)

judge(res_med["Grids"], res_med["PencilGrids"])
