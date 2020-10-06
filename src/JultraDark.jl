module JultraDark

using Base.Threads: @threads
using Statistics
using AbstractFFTs: fftfreq, rfftfreq
using PencilFFTs, MPI

using FFTW
FFTW.set_num_threads(Threads.nthreads())

export simulate
export PencilGrids
export Config, OutputConfig

include("mpi_grids.jl")
include("output.jl")
include("config.jl")

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
julia> using JultraDark: take_steps!, PencilGrids, OutputConfig

julia> take_steps!(PencilGrids(1.0, 16), 0, 0.5, 10, OutputConfig(mktempdir(), []), t->1)
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

    mkpath(output_config.directory)

    t_begin = output_config.output_times[1]

    # Initialise vars other than ψx
    grids.ψk .= (grids.fft_plan * grids.ψx)
    grids.ρx .= abs2.(grids.ψx)
    grids.Φk .= -4 * π * (grids.rfft_plan * grids.ρx) ./ (options.a(t_begin) * grids.rk.^2)
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= grids.rfft_plan \ grids.Φk

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
