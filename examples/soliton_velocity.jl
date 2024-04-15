using MPI
using NPZ
using CairoMakie
using UltraDark

if ~MPI.Initialized()
    MPI.Init()
end
comm = MPI.COMM_WORLD
@info "Process ID $(MPI.Comm_rank(comm)+1) of $(MPI.Comm_size(comm)) with $(Threads.nthreads()) threads\n"

const resol = 64

const num_snapshots = 20
const output_times = LinRange(0, 5, num_snapshots)

const mass = 10
const position0 = [-4, 0, 0]
const velocity = [1, 0, 0]
const phase = 0
const t0 = 0

const output_dir = joinpath(mktempdir(), "output/vel/centered_on_soliton")

function run_sim()

    output_config = OutputConfig(output_dir, output_times; box = true, slice = false)
    options = Config.SimulationConfig()

    grids = PencilGrids(10.0, resol)
    @info "Grids created" MPI.Comm_rank(comm) size(grids.ψk)

    mkpath(output_config.directory)
    npzwrite(joinpath(output_config.directory, "x.npy"), grids.x)
    npzwrite(joinpath(output_config.directory, "y.npy"), grids.y)
    npzwrite(joinpath(output_config.directory, "z.npy"), grids.z)

    UltraDark.Initialise.add_fdm_soliton!(grids, mass, position0, velocity, phase, t0)
    @info "Initialized" MPI.Comm_rank(comm)

    simulate!(grids, options, output_config) == nothing
    @info "Simulation complete" MPI.Comm_rank(comm)

end

function plot_results()
    @info "Starting plots" MPI.Comm_rank(comm)
    pos_the =
        reshape(position0, (1, :)) .+
        (reshape(velocity, (1, :)) .* reshape(output_times, (:, 1)))

    x = npzread(joinpath(output_dir, "x.npy"))
    y = npzread(joinpath(output_dir, "y.npy"))
    z = npzread(joinpath(output_dir, "z.npy"))

    pos_sim = Array{Float64}(undef, (num_snapshots, 3))
    maxs = Array{Float64}(undef, num_snapshots)

    for i in 1:num_snapshots
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))

        rhomax, maxindex = findmax(rho)

        maxs[i] = rhomax
        pos_sim[i, :] .= [x[maxindex[1]], y[maxindex[2]], z[maxindex[3]]]
    end

    yerr = abs(x[2] - x[1]) / 2

    fig_pos = Figure()
    ax_pos = Axis(fig_pos[1, 1], ylabel = "position", xlabel = L"t")

    lines!(ax_pos, output_times, pos_the[:, 1], label = "theory")
    errorbars!(
        ax_pos,
        output_times,
        pos_sim[:, 1],
        yerr,
        whiskerwidth = 10,
        label = "simulation",
    )
    axislegend(ax_pos, position = :rc)

    save("soliton_position.pdf", fig_pos)

    fig_max, ax_max, p_max = scatter(
        output_times,
        maxs,
        label = "theory",
        axis = (limits = (nothing, nothing, 0, maximum(maxs) * 1.1),),
    )

    save("soliton_max.pdf", fig_max)

    fig_rho = Figure()
    ax_rho = Axis(fig_rho[1, 1], aspect = DataAspect())
    ax_psi = Axis(fig_rho[2, 1], aspect = DataAspect())
    colorrange = (0, maximum(maxs))
    Colorbar(fig_rho[1, 2], limits = colorrange, label = L"\rho")
    Colorbar(fig_rho[2, 2], limits = (-π, π), label = L"arg(\psi)")

    record(fig_rho, "soliton_velocity.mp4", 1:num_snapshots) do i
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))
        psi = npzread(joinpath(output_dir, "psi_$i.npy"))

        heatmap!(fig_rho[1, 1], x[:], y[:], rho[:, :, resol÷2])

        heatmap!(fig_rho[2, 1], x[:], y[:], angle.(psi[:, :, resol÷2]))
    end
end

run_sim()

if MPI.Comm_rank(comm) == 0
    plot_results()
end
