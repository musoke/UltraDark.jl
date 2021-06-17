using BenchmarkTools
using UltraDark
using LinearAlgebra
using MPI

if ~MPI.Initialized()
    MPI.Init()
end

n = 128

function fft_broadcast!(grids)
    grids.ψk .= grids.fft_plan * grids.ψx
end

function fft_inplace!(grids)
    mul!(grids.ψk, grids.fft_plan, grids.ψx)
end

function fft_inplace_sep!(psi_k, psi_x, fft_plan)
    mul!(psi_k, fft_plan, psi_x)
end

function fft_inv_broadcast!(grids)
    grids.ψx .= grids.fft_plan \ grids.ψk
end

function fft_inv_inplace!(grids)
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)
end

suite = BenchmarkGroup()

for (label, grids) in [("Grids", Grids(1.0, n)), ("PencilGrids", PencilGrids(1.0, n))]

    grids.ψx .= randn(size(grids.ψx))
    grids.ψk .= randn(size(grids.ψk))

    suite[label] = BenchmarkGroup()

    suite[label]["fft_broadcast!"] = @benchmarkable fft_broadcast!($grids)
    suite[label]["fft_inplace!"] = @benchmarkable fft_inplace!($grids)
    suite[label]["fft_inv_broadcast!"] = @benchmarkable fft_inv_broadcast!($grids)
    suite[label]["fft_inv_inplace!"] = @benchmarkable fft_inv_inplace!($grids)

    suite[label]["fft_inplace_sep!"] = @benchmarkable fft_inplace_sep!($(grids.ψk), $(grids.ψx), $(grids.fft_plan))
end

@show tune!(suite)

results = run(suite)

@show res_min = minimum(results)
@show res_med = median(results)

judge(res_med["Grids"], res_med["PencilGrids"])
