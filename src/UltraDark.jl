module UltraDark

using Base.Threads: @threads
using LinearAlgebra
using AbstractFFTs: fftfreq, rfftfreq
using PencilFFTs, MPI

using FFTW

export simulate!
export Grids, PencilGrids
export dV
export Config, SimulationConfig, constant_scale_factor, TimeStepOptions
export OutputConfig
export Summary

include("grids.jl")
include("pencil_grids.jl")
include("phase_diff.jl")
include("time_step.jl")
include("summary.jl")
include("output.jl")
include("config.jl")

import .Output: OutputConfig
import .Output: output_state, output_xyz, output_external_states_headers
import .Output: output_summary_row, output_summary_header
import .Config: SimulationConfig, constant_scale_factor, TimeStepOptions

"""
    outer_step!(Δt, grids, constants; a=1.0)
    outer_step!(Δt, grids, constants, s; a=1.0)

Perform the "outer" time step in the symmetrized split-step Fourier method.

This step only updates the phase of ψ applying accelerations due to gravity,
the amplitude is not changed.
"""
function outer_step!(Δt, grids, constants; a = 1.0)
    @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(-im * Δt / 1 * grids.Φx[i])
    end
end

function outer_step!(Δt, grids, constants, s; a = 1.0) end

"""
    inner_step!(Δt, grids, constants; a=1.0)
    inner_step!(Δt, grids, constants, s; a=1.0)

Perform the "inner" time step in the symmetrized split-step Fourier method.

This step applies the diffusion terms and updates the gravitational potential.
"""
function inner_step!(Δt, grids, constants; a = 1.0)

    mul!(grids.ψk, grids.fft_plan, grids.ψx)
    @inbounds @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt / 2 * grids.k[i]^2 / a^2)
    end
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)

end

function inner_step!(Δt, grids, constants, s; a = 1.0) end

function update_gravitational_potential!(grids, constants; a = 1.0)
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
    auxiliary_step!(Δt, grids, t, constants)
    auxiliary_step!(Δt, grids, t, constants, s; a = 1.0)

Do an auxiliary inner step.
By default this does nothing, but can be overridden in multiple dispatch.
"""
function auxiliary_step!(Δt, grids, t, constants) end

function auxiliary_step!(Δt, grids, t, constants, s; a = 1.0) end

"""
    add_external_potential!(t, grids, constants)

Add a gravitational potential to the grid.
By default this does nothing, but can be overridden in multiple dispatch.
"""
function add_external_potential!(t, grids, s; a = 1.0)
    grids.Φx .+= gravitational_potential(s, grids; a = a)
end

function gravitational_potential(s, grids; a = 1.0)
    0.0
end

"""
Take `n` steps with time step `Δt`

# Examples

```jldoctest
julia> using UltraDark: take_steps!, Grids, OutputConfig, Config

julia> take_steps!(
           Grids(1.0, 16),
           0,
           0.5,
           10,
           OutputConfig(mktempdir(), []),
           Config.constant_scale_factor,
           nothing,
           (),
       )
5.0
```
"""
function take_steps!(grids, t_start, Δt, n, output_config, a, constants, external_states)

    t::Float64 = t_start

    half_step = true

    for step in 1:n
        if half_step
            outer_step!(Δt / 2, grids, constants; a = a(t))
            for s in external_states
                outer_step!(Δt / 2, grids, constants, s; a = a(t))
            end

            t += Δt / 2
            half_step = false
        else
            outer_step!(Δt, grids, constants; a = a(t))
            for s in external_states
                outer_step!(Δt, grids, constants, s; a = a(t))
            end

            t += Δt
        end

        inner_step!(Δt, grids, constants; a = a(t))
        for s in external_states
            inner_step!(Δt, grids, constants, s; a = a(t))
        end

        update_gravitational_potential!(grids, constants; a = a(t))
        for s in external_states
            add_external_potential!(t, grids, s; a = a(t))
        end

        auxiliary_step!(Δt, grids, t, constants)
        for s in external_states
            auxiliary_step!(Δt, grids, t, constants, s; a = a(t))
        end

        output_summary_row(grids, output_config, t, a(t), Δt, constants, external_states)
    end

    outer_step!(Δt / 2, grids, constants)
    t += Δt / 2

    t
end

"""
Evolve `grids` forward from `t_start` to `t_end`
"""
function evolve_to!(
    t_start,
    t_end,
    grids,
    output_config,
    sim_config;
    constants = nothing,
    external_states = (),
)

    @assert t_start < t_end

    t = t_start

    while (t < t_end) && ~(t ≈ t_end)

        if ~good_phase_diff(grids, sim_config)
            output_state(grids, external_states, output_config, -1)
            throw("Phase gradient is too large to continue")
        end

        Δt, n_steps = actual_time_step(
            max_time_step(grids, sim_config.a(t), external_states),
            t_end - t,
            sim_config.time_step_options,
        )

        t = take_steps!(
            grids,
            t,
            Δt,
            n_steps,
            output_config,
            sim_config.a,
            constants,
            external_states,
        )
    end

    t
end

function simulate!(
    grids,
    sim_config,
    output_config::OutputConfig;
    constants = nothing,
    external_states = (),
)

    # Setup output
    mkpath(output_config.directory)
    output_summary_header(output_config)
    output_xyz(grids, output_config)
    Output.output_external_states_headers(external_states, output_config)

    t_begin = output_config.output_times[1]

    # Initialise vars other than ψx
    # This is required so the initial time step can be calculated
    t_initial = output_config.output_times[1]

    update_gravitational_potential!(grids, constants; a = sim_config.a(t_initial))
    for s in external_states
        add_external_potential!(t_initial, grids, s, a = sim_config.a(t_initial))
    end

    # Output initial conditions
    output_state(grids, external_states, output_config, 1)
    if ~good_phase_diff(grids, sim_config)
        throw("Initial phase gradient is already in untrusted regime")
    end

    for (index, t_end) in enumerate(output_config.output_times[2:end])
        t_begin = evolve_to!(
            t_begin,
            t_end,
            grids,
            output_config,
            sim_config;
            constants = constants,
            external_states = external_states,
        )
        @info "Reached time $t_begin"
        output_state(grids, external_states, output_config, index + 1)
    end

    if ~good_phase_diff(grids, sim_config)
        throw("Phase gradient is too large to end")
    end

end

end # module
