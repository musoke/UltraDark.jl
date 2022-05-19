using Base.Threads
using BenchmarkTools
using LinearAlgebra
using LoopVectorization
using MPI
using Strided
using UltraDark: psi_whole_step!, Grids, PencilGrids

if ~MPI.Initialized()
    MPI.Init()
end

n = 128

Δt = 1.0

function step_prim_for!(Δt, psi, phi, nothing)
    for i in eachindex(psi)
        psi[i] *= exp(-im * Δt / 1 * phi[i])
    end
end

function step_prim_for_threads!(Δt, psi, phi, nothing)
    @threads for i in eachindex(psi)
        psi[i] *= exp(-im * Δt / 1 * phi[i])
    end
end

function step_for!(Δt, grids, nothing)
    for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_cartesian!(Δt, grids, nothing)
    for i in CartesianIndices(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_threads!(Δt, grids, nothing)
    @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_inbounds_threads!(Δt, grids, nothing)
    @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_fastmath_inbounds_threads!(Δt, grids, nothing)
    @fastmath @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_avx!(Δt, grids, nothing)
    @avx for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_avxt!(Δt, grids, nothing)
    @avxt for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_turbo!(Δt, grids, nothing)
    @turbo for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_for_tturbo!(Δt, grids, nothing)
    @tturbo for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function step_broadcast!(Δt, grids, nothing)
    @. grids.ψx *= exp(-im * Δt / 1 * grids.Φx)
end

function step_broadcast_strided!(Δt, grids, nothing)
    @strided @. grids.ψx *= exp(-im * Δt / 1 * grids.Φx)
end

# function step_broadcast_turbo!(Δt, grids, nothing)
#     @turbo @. grids.ψx *= exp(- im * Δt / 1 * grids.Φx)
# end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
suite = BenchmarkGroup()

for (label, grids) in [("Grids", Grids(1.0, n)), ("PencilGrids", PencilGrids(1.0, n))]

    grids.ψx .= randn(size(grids.ψx))
    grids.ψk .= randn(size(grids.ψk))

    suite[label] = BenchmarkGroup()

    for f in [step_prim_for!, step_prim_for_threads!]
        suite[label]["$f"] = @benchmarkable $f($Δt, $(grids.ψx), $(grids.Φx), $nothing)
    end

    for f in [
        psi_whole_step!,
        step_for!,
        step_for_cartesian!,
        step_for_threads!,
        step_for_inbounds_threads!,
        step_for_fastmath_inbounds_threads!,
        step_for_avx!,
        step_for_avxt!,
        step_for_turbo!,
        step_broadcast!,
        # step_broadcast_strided!,
        # step_broadcast_turbo!,
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
