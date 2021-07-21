module UltraDark

using Base.Threads: @threads
using LinearAlgebra
using AbstractFFTs: fftfreq, rfftfreq
using PencilFFTs, MPI

using FFTW

export simulate
export Grids, PencilGrids
export Config, SimulationConfig, constant_scale_factor
export OutputConfig

const PHASE_GRAD_LIMIT = π / 4

include("grids.jl")
include("pencil_grids.jl")
include("phase_diff.jl")
include("time_step.jl")
include("output.jl")
include("config.jl")

import .Output: OutputConfig
import .Output: output_summary_row, output_summary_header, output_grids
import .Config: SimulationConfig, constant_scale_factor

function psi_half_step!(Δt, grids, constants)
    @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 2 * grids.Φx[i])
    end
end

function psi_whole_step!(Δt, grids, constants)
    @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 1 * grids.Φx[i])
    end
end

function phi_whole_step!(Δt, grids, constants; a=1.0)
    # TODO: not all part of Φ update

    mul!(grids.ψk, grids.fft_plan, grids.ψx)
    @inbounds @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt/2 * grids.k[i]^2 / a^2)
    end
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)

end

function update_gravitational_potential!(grids, constants; a=1.0)
    @inbounds @threads for i in eachindex(grids.ρx)
        grids.ρx[i] = abs2(grids.ψx[i])
    end

    mul!(grids.Φk, grids.rfft_plan, grids.ρx)
    @inbounds @threads for i in eachindex(grids.Φk)
        grids.Φk[i] *= -4 * π / (a * grids.rk[i]^2)
    end
    grids.Φk[grids.rk_vanish_indices] .= 0
    ldiv!(grids.Φx, grids.rfft_plan, grids.Φk)
end

"""
    add_external_potential!(t, grids, constants)

Add a gravitational potential to the grid.
By default this does nothing, but can be overridden in multiple dispatch.
"""
function add_external_potential!(t, grids, constants)
end

"""
Take `n` steps with time step `Δt`

# Examples

```jldoctest
julia> using UltraDark: take_steps!, Grids, OutputConfig, Config

julia> take_steps!(Grids(1.0, 16), 0, 0.5, 10, OutputConfig(mktempdir(), []), Config.constant_scale_factor, nothing)
5.0

```
"""
function take_steps!(grids, t_start, Δt, n, output_config, a, constants)

    t::Float64 = t_start

    half_step = true

    for step in 1:n
        if half_step
            psi_half_step!(Δt, grids, constants)
            t += Δt / 2
            half_step = false
        else
            psi_whole_step!(Δt, grids, constants)
            t += Δt
        end

        phi_whole_step!(Δt, grids, constants; a=a(t))
        update_gravitational_potential!(grids, constants; a=a(t))
        add_external_potential!(t, grids, constants)

        output_summary_row(grids, output_config, t, a(t), Δt)
    end

    psi_half_step!(Δt, grids, constants)
    t += Δt / 2

    output_summary_row(grids, output_config, t, a(t), Δt)

    t
end

"""
Evolve `grids` forward from `t_start` to `t_end`
"""
function evolve_to!(t_start, t_end, grids, output_config, config::Config.SimulationConfig; constants=nothing)

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

        t = take_steps!(grids, t, Δt, n_steps, output_config, config.a, constants)
    end

    t
end

function simulate(grids, options::Config.SimulationConfig, output_config::OutputConfig; constants=nothing)

    # Setup output
    mkpath(output_config.directory)
    output_summary_header(output_config)

    t_begin = output_config.output_times[1]

    # Initialise vars other than ψx
    # This is required so the initial time step can be calculated
    t_initial = output_config.output_times[1]
    update_gravitational_potential!(grids, constants; a=options.a(t_initial))

    # Output initial conditions
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
            options;
            constants=constants,
        )
        @info "Reached time $t_begin"
        output_grids(grids, output_config, index + 1)
    end

    if max_normed_phase_grad(grids) > PHASE_GRAD_LIMIT
        throw("Phase gradient is too large to end")
    end

end

end # module
