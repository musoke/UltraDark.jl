module JultraDark

using Statistics
using FFTW

export simulate
export Grids
export OutputConfig

include("grids.jl")
include("output.jl")

function psi_half_step!(Δt::Real, grids)
    grids.ψx .*= exp.(- im * Δt / 2 * grids.Φx)
end

function psi_whole_step!(Δt::Real, grids)
    grids.ψx .*= exp.(- im * Δt / 1 * grids.Φx)
end

function phi_whole_step!(Δt::Real, grids; a::Real=1.0)
    # TODO: not all part of Φ update

    grids.ψk .= fft(grids.ψx) .* exp.(-im * Δt/2 * grids.k.^2 / a^2)
    grids.ψx .= ifft(grids.ψk)

    grids.ρx .= abs2.(grids.ψx)
    grids.Φk .= -4 * π * rfft(grids.ρx) ./ (a * grids.rk.^2)
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= irfft(grids.Φk, size(grids.Φx, 1))
end

function max_timestep(grids)
    1e-1  # TODO
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

julia> take_steps!(Grids(1.0, 16), 0, 0.5, 10, OutputConfig(mktempdir()), t->1)
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

        output_summary(grids, output_config, t, a(t))
    end

    psi_half_step!(Δt, grids)
    t += Δt / 2

    t
end

"""
Evolve `grids` forward from `t_start` to `t_end`
"""
function evolve_to!(t_start, t_end, grids, output_config, a)

    @assert t_start < t_end

    t = t_start

    max_num_steps = 10  # Steps to take before recalculating step size
    while t < t_end
        Δt, n_steps = actual_time_step(max_timestep(grids), t_end - t, max_num_steps)

        t = take_steps!(grids, t, Δt, n_steps, output_config, a)
    end

    t
end

function simulate()
    resol = 16

    grids = Grids(zeros(Complex{Float64}, resol, resol, resol), 1)

    a = t -> 1

    output_times = 0:10

    output_config = OutputConfig("output")
    mkpath(output_config.directory)

    t_begin = output_times[1]

    for (index, t_end) in enumerate(output_times[2:end])
        t_begin = evolve_to!(
            t_begin,
            t_end,
            grids,
            output_config,
            a,
        )
        output_grids(grids, output_config, index)
    end

end

end # module
