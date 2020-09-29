module JultraDark

using Base.Threads: @threads
using Statistics
using FFTW

export simulate
export Grids
export Config, OutputConfig

include("grids.jl")
include("output.jl")
include("config.jl")

function psi_half_step!(Δt::Real, grids)
    @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 2 * grids.Φx[i])
    end
end

function psi_whole_step!(Δt::Real, grids)
    @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 1 * grids.Φx[i])
    end
end

function phi_whole_step!(Δt::Real, grids; a::Real=1.0)
    # TODO: not all part of Φ update

    grids.ψk .= grids.fft_plan * grids.ψx
    @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt/2 * grids.k[i]^2 / a^2)
    end
    grids.ψx .= grids.fft_plan \ grids.ψk

    @threads for i in eachindex(grids.ρx)
        grids.ρx[i] = abs2(grids.ψx[i])
    end

    grids.Φk .= grids.rfft_plan * grids.ρx
    @threads for i in eachindex(grids.Φk)
        grids.Φk[i] *= -4 * π / (a * grids.rk[i]^2)
    end
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= grids.rfft_plan \ grids.Φk
end

"""
    max_time_step(grids, a)

Calculate an upper bound on the time step
"""
function max_time_step(grids, a)
    max_time_step_gravity = 2π / maximum(grids.Φx)
    max_time_step_pressure = 2π * 2 / maximum(grids.k)^2 * a^2  # TODO: cache k_max

    @assert isfinite(max_time_step_gravity)
    @assert isfinite(max_time_step_pressure)

    time_step = min(max_time_step_gravity, max_time_step_pressure)

    time_step
end

"""
    actual_time_step(max_timestep::Real, time_interval::Real, n::Integer)

Actual size and number of time steps that should be taken if the maximum 
is `max_timestep`, no more than `n` steps should be taken, and they should
fit in `time_interval`.

# Examples

```jldoctest
julia> using JultraDark: actual_time_step

julia> actual_time_step(0.11, 1, 20)
(0.1, 10)
```
"""
function actual_time_step(max_timestep, time_interval, n)::Tuple{Real, Integer}
    if max_timestep * n > time_interval
        num_steps = ceil(time_interval / max_timestep)
        time_interval / num_steps, num_steps
    else
        max_timestep, n
    end
end

"""
Take `n` steps with time step `Δt`

# Examples

```jldoctest
julia> using JultraDark: take_steps!, Grids, OutputConfig

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

        output_summary(grids, output_config, t, a(t), Δt)
    end

    psi_half_step!(Δt, grids)
    t += Δt / 2

    output_summary(grids, output_config, t, a(t), Δt)

    t
end

"""
Evolve `grids` forward from `t_start` to `t_end`
"""
function evolve_to!(t_start, t_end, grids, output_config, config::Config.SimulationConfig)

    @assert t_start < t_end

    t = t_start

    while (t < t_end) && ~(t ≈ t_end)
        @debug "t = $t"

        Δt, n_steps = actual_time_step(
            max_time_step(grids, config.a(t)),
            t_end - t,
            config.time_step_update_period,
        )

        t = take_steps!(grids, t, Δt, n_steps, output_config, config.a)
    end

    t
end

function simulate(grids::Grids, options::Config.SimulationConfig, output_config::OutputConfig)

    mkpath(output_config.directory)

    t_begin = output_config.output_times[1]
    output_grids(grids, output_config, 0)

    for (index, t_end) in enumerate(output_config.output_times[2:end])
        t_begin = evolve_to!(
            t_begin,
            t_end,
            grids,
            output_config,
            options,
        )
        @info "Reached time $t_begin"
        output_grids(grids, output_config, index)
    end

end

end # module
