using MPI
using NPZ
using Plots
using UltraDark

if ~MPI.Initialized()
    MPI.Init()
end
comm = MPI.COMM_WORLD
@info "Process ID $(MPI.Comm_rank(comm)+1) of $(MPI.Comm_size(comm)) with $(Threads.nthreads()) threads\n"

include(joinpath(@__DIR__, "init_soliton.jl"))

const resol = 64

const num_snapshots = 20
const output_times = LinRange(0, 5, num_snapshots)

const mass = 10
const position0 = [-4, 0, 0]
const velocity = [1, 0, 0]
const phase = 0
const t0 = 0

const output_dir = "output/vel"

function run_sim()

    output_config = OutputConfig(output_dir, output_times; box = true, slice = false)
    options = Config.SimulationConfig()

    grids = PencilGrids(10.0, resol)
    @info "Grids created" MPI.Comm_rank(comm) size(grids.ψk)

    mkpath(output_config.directory)
    npzwrite(joinpath(output_config.directory, "x.npy"), grids.x)
    npzwrite(joinpath(output_config.directory, "y.npy"), grids.y)
    npzwrite(joinpath(output_config.directory, "z.npy"), grids.z)

    add_soliton(grids, mass, position0, velocity, phase, t0)
    @info "Initialized" MPI.Comm_rank(comm)

    simulate!(grids, options, output_config) == nothing
    @info "Simulation complete" MPI.Comm_rank(comm)

end

function plot_results()
    @info "Starting plots" MPI.Comm_rank(comm)
    pos_the =
        reshape(position0, (1, :)) .+
        (reshape(velocity, (1, :)) .* reshape(output_times, (:, 1)))

    pos_sim = Array{Float64}(undef, (num_snapshots, 3))
    maxs = Array{Float64}(undef, num_snapshots)

    x = npzread(joinpath(output_dir, "x.npy"))
    y = npzread(joinpath(output_dir, "y.npy"))
    z = npzread(joinpath(output_dir, "z.npy"))

    for i in 1:num_snapshots
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))

        rhomax, maxindex = findmax(rho)

        maxs[i] = rhomax
        pos_sim[i, :] .= [x[maxindex[1]], y[maxindex[2]], z[maxindex[3]]]

    end

    yerr = abs(x[2] - x[1]) / 2
    plot_pos = plot(legend = :bottomright)
    scatter!(
        output_times,
        pos_sim[:, 1],
        yerr = yerr,
        markershape = :auto,
        label = "simulation",
    )
    plot!(output_times, pos_the[:, 1], label = "theory")
    plot!(xlabel = raw"$t$", ylabel = raw"$x$")

    savefig(plot_pos, "soliton_position.pdf")

    plot_max = plot(legend = :bottomright)
    plot!(output_times, maxs, label = "theory")

    anim = @animate for i in 1:num_snapshots
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))
        psi = npzread(joinpath(output_dir, "psi_$i.npy"))

        density_plot = contourf(
            rho[:, :, resol÷2],
            title = raw"$\rho$",
            clims = (0, maximum(maxs)),
            showaxis = false,
            ticks = false,
            aspect_ratio = :equal,
        )
        phase_plot = contourf(
            angle.(psi[:, :, resol÷2]),
            title = "phase",
            clims = (-π, +π),
            showaxis = false,
            ticks = false,
            aspect_ratio = :equal,
        )

        plot(density_plot, phase_plot)
    end

    gif(anim, "soliton_velocity.gif", fps = 15)
end

run_sim()

if MPI.Comm_rank(comm) == 0
    plot_results()
end
