"""
    actual_time_step(max_time_step, time_interval, time_step_options)

Actual size and number of time steps that should be taken if the maximum
is `max_time_step`. No more than `time_step_options.update_period` steps should
be taken, and they should fit in `time_interval`.

# Examples

```jldoctest
julia> using UltraDark: actual_time_step, TimeStepOptions

julia> actual_time_step(0.11, 1, TimeStepOptions())
(0.1, 10)
```
"""
function actual_time_step(max_time_step, time_interval, time_step_options)::Tuple{Float64, Integer}
    max_time_step *= time_step_options.multiplier
    if max_time_step * time_step_options.update_period > time_interval
        num_steps = ceil(time_interval / max_time_step)
        time_interval / num_steps, num_steps
    else
        max_time_step, time_step_options.update_period
    end
end

"""
    max_time_step(grids, a)

Calculate an upper bound on the time step
"""
function max_time_step(grids, a, time_step_options)
    max_time_step_gravity = 2π / maximum(abs.(grids.Φx))
    max_time_step_pressure = 2π * 2 / maximum(grids.k)^2 * a^2  # TODO: cache k_max

    @assert isfinite(max_time_step_gravity)
    @assert isfinite(max_time_step_pressure)

    time_step = min(max_time_step_gravity, max_time_step_pressure)

    time_step
end

"""
    max_time_step(grids, a)

Calculate an upper bound on the time step for `PencilGrids`
"""
function max_time_step(grids::PencilGrids, a, time_step_options)
    # Find maximum on local grid
    local_max_time_step_gravity = 2π / maximum(abs.(grids.Φx))
    local_max_time_step_pressure = 2π * 2 / maximum(grids.k)^2 * a^2  # TODO: cache k_max

    # Maximize over other grids
    max_time_step_gravity = MPI.Allreduce(local_max_time_step_gravity, MPI.MIN, grids.MPI_COMM)
    max_time_step_pressure = MPI.Allreduce(local_max_time_step_pressure, MPI.MIN, grids.MPI_COMM)

    @assert isfinite(max_time_step_gravity)
    @assert isfinite(max_time_step_pressure)

    time_step = min(max_time_step_gravity, max_time_step_pressure)

    time_step
end
