module UltraDark

using Base.Threads: @threads
using LinearAlgebra
using Statistics
using AbstractFFTs: fftfreq, rfftfreq
using PencilFFTs, MPI

using FFTW

export simulate
export Grids, PencilGrids
export Config, OutputConfig

const PHASE_GRAD_LIMIT = π / 4

include("grids.jl")
include("pencil_grids.jl")
include("phase_diff.jl")
include("time_step.jl")
include("output.jl")
include("config.jl")

function psi_half_step!(Δt::Real, grids)
    @fastmath @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 2 * grids.Φx[i])
    end
end

function psi_whole_step!(Δt::Real, grids)
    @fastmath @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 1 * grids.Φx[i])
    end
end

function phi_whole_step!(Δt::Real, grids; a::Real=1.0)
    # TODO: not all part of Φ update

    mul!(grids.ψk, grids.fft_plan, grids.ψx)
    @fastmath @inbounds @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt/2 * grids.k[i]^2 / a^2)
    end
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)

    @fastmath @inbounds @threads for i in eachindex(grids.ρx)
        grids.ρx[i] = abs2(grids.ψx[i])
    end

    mul!(grids.Φk, grids.rfft_plan, grids.ρx)
    @fastmath @inbounds @threads for i in eachindex(grids.Φk)
        grids.Φk[i] *= -4 * π / (a * grids.rk[i]^2)
    end
    grids.Φk[1, 1, 1] = 0
    ldiv!(grids.Φx, grids.rfft_plan, grids.Φk)
end

"""
    add_external_potential!(t, grids)

Add a gravitational potential to the grid.
By default this does nothing, but can be overridden in multiple dispatch.
"""
function add_external_potential!(t, grids)
end

"""
Take `n` steps with time step `Δt`

# Examples

```jldoctest
julia> using UltraDark: take_steps!, Grids, OutputConfig

julia> take_steps!(Grids(1.0, 16), 0, 0.5, 10, OutputConfig(mktempdir(), []), t->1)
5.0

```
"""
function take_steps!(grids, t_start, Δt, n, output_config, a)

    t = t_start

    half_step = true

    for step in 1:n
        if half_step
            psi_half_step!(Δt, grids)
            t += Δt / 2
            half_step = false
        else
            psi_whole_step!(Δt, grids)
            t += Δt
        end

        phi_whole_step!(Δt, grids, a=a(t))
        add_external_potential!(t, grids)

        output_summary_row(grids, output_config, t, a(t), Δt)
    end

    psi_half_step!(Δt, grids)
    t += Δt / 2

    output_summary_row(grids, output_config, t, a(t), Δt)

    t
end

"""
Evolve `grids` forward from `t_start` to `t_end`
"""
function evolve_to!(t_start, t_end, grids, output_config, config::Config.SimulationConfig)

    @assert t_start < t_end

    t = t_start

    while (t < t_end) && ~(t ≈ t_end)

        if max_normed_phase_grad(grids) > PHASE_GRAD_LIMIT
            output_grids(grids, output_config, -1)
            throw("Phase gradient is too large to continue")
        end

        Δt, n_steps = actual_time_step(
            max_time_step(grids, config.a(t)),
            t_end - t,
            config.time_step_update_period,
        )
        @debug "t = $t, max_time_step = $(max_time_step(grids, config.a(t))), Δt = $Δt"

        t = take_steps!(grids, t, Δt, n_steps, output_config, config.a)
    end

    t
end

function simulate(grids, options::Config.SimulationConfig, output_config::OutputConfig)

    # Setup output
    mkpath(output_config.directory)
    output_summary_header(output_config)

    t_begin = output_config.output_times[1]

    # Initialise vars other than ψx
    grids.ψk .= (grids.fft_plan * grids.ψx)
    grids.ρx .= abs2.(grids.ψx)
    grids.Φk .= -4 * π * (grids.rfft_plan * grids.ρx) ./ (options.a(t_begin) * grids.rk.^2)
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= grids.rfft_plan \ grids.Φk

    output_grids(grids, output_config, 1)

    if max_normed_phase_grad(grids) > PHASE_GRAD_LIMIT
        throw("Phase gradient is too large to start")
    end

    for (index, t_end) in enumerate(output_config.output_times[2:end])
        t_begin = evolve_to!(
            t_begin,
            t_end,
            grids,
            output_config,
            options,
        )
        @info "Reached time $t_begin"
        output_grids(grids, output_config, index + 1)
    end

    if max_normed_phase_grad(grids) > PHASE_GRAD_LIMIT
        throw("Phase gradient is too large to end")
    end

end

end # module
