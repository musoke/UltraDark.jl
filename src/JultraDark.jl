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

function actual_timestep(max_timestep, time_interval)
    time_interval / ceil(time_interval / max_timestep)
end

function evolve_to!(t_start, t_end, grids, output_config)

    half_step = true

    t = t_start
    Δt = actual_timestep(max_timestep(grids), t_end - t_start)

    steps = Int((t_end - t_start)/Δt)

    for step in 1:steps
        if half_step
            psi_half_step!(Δt, grids)
            t += Δt / 2
            half_step = false
        else
            psi_whole_step!(Δt, grids)
            t += Δt
        end
        a = 1  # TODO a
        phi_whole_step!(Δt, grids, a=a)

        output_summary(grids, output_config, t)
    end

    psi_half_step!(Δt, grids)
    t += Δt / 2

    t
end

function simulate()
    resol = 16

    grids = Grids(zeros(Complex{Float64}, resol, resol, resol), 1)

    output_times = 0:10

    t_begin = output_times[1]

    output_config = OutputConfig("output")
    mkpath(output_config.directory)

    for (index, t_end) in enumerate(output_times[2:end])
        t_begin = evolve_to!(
            t_begin,
            t_end,
            grids,
            output_config,
        )
        output_grids(grids, output_config, index)
    end

end

end # module
