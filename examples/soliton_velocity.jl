# # PencilGrids example
#
# This example uses [`PencilGrids`](@ref UltraDark.PencilGrids) to demonstrate potential
# use of [`MPI`](@extref MPI MPI.jl)-based parallelism.
#
# Note, however, that the online documentation is rendered by a single process and so using
# `PencilGrids` has little practical benefit other than demonstrating how it is done.
# If you run this example with multiple MPI processes, many of the `@info` statements will
# print once for each process.
#
#
# ## Initialise a simulation
#
using MPI
using NPZ
using CairoMakie
using UltraDark

# Initialise [`MPI`](@extref MPI MPI.jl) and set up [`MPI.Comm`](@extref) for inter-process communication.
if ~MPI.Initialized()
    MPI.Init()
end
comm = MPI.COMM_WORLD
@info "Process ID $(MPI.Comm_rank(comm)+1) of $(MPI.Comm_size(comm)) with $(Threads.nthreads()) threads\n"

# Define grids on which simulation will be performed.
# We use a relatively low resolution since this is just an example.
const resol = 64
const box_length = 10.0

grids = PencilGrids(box_length, resol)
@info "Grids created" MPI.Comm_rank(comm) size(grids.ψk)

# Define initial conditions consisting of a soliton of mass `mass` at position `position0`
# with velocity `velocity`.
const mass = 10
const position0 = [-4, 0, 0]
const velocity = [1, 0, 0]
const phase = 0
const t0 = 0

UltraDark.Initialise.add_fdm_soliton!(grids, mass, position0, velocity, phase, t0)
UltraDark.update_gravitational_potential!(grids) # ensure the density is up to date
@info "Initialized" MPI.Comm_rank(comm) UltraDark.mass(grids)

# Define output configuration
const output_dir = joinpath(mktempdir(), "output")

const num_snapshots = 20
const output_times = LinRange(0, 5, num_snapshots)

output_config = OutputConfig(output_dir, output_times; box = true, slice = false);

# ## Run the simulation
simulate!(grids, output_config) == nothing
@info "Simulation complete" MPI.Comm_rank(comm)

# ## Plot the results
#
# Plotting and other analysis should only be done once, not repeated in every MPI job.
# Here we do it on the root process by checking against the output of
# [`MPI.Comm_rank`](@extref).
#
# First, plot the theoretical and numerical positions of the soliton's peak against time.
# The error bars represent uncertainty due to the grid's finite resolution.
#
function plot_position()
    x = npzread(joinpath(output_dir, "x.npy"))
    y = npzread(joinpath(output_dir, "y.npy"))
    z = npzread(joinpath(output_dir, "z.npy"))

    position_sim = Array{Float64}(undef, (num_snapshots, 3))
    position_theory =
        reshape(position0, (1, :)) .+
        (reshape(velocity, (1, :)) .* reshape(output_times, (:, 1)))

    for i in 1:num_snapshots
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))
        rhomax, maxindex = findmax(rho)
        position_sim[i, :] .= [x[maxindex[1]], y[maxindex[2]], z[maxindex[3]]]
    end

    yerr = abs(x[2] - x[1]) / 2

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"$t / \mathcal{T}$", ylabel = L"$x / \mathcal{L}$")
    lines!(
        ax,
        output_times,
        position_theory[:, 1],
        color = :gray,
        linestyle = :dash,
        label = "theory",
    )
    errorbars!(
        ax,
        output_times,
        position_sim[:, 1],
        yerr,
        whiskerwidth = 10,
        label = "simulation",
    )
    axislegend(ax, position = :rb)

    fig
end

if MPI.Comm_rank(comm) == 0
    @info "Plotting" MPI.Comm_rank(comm)
    plot_position()
end

# Next, plot the maximum density as a function of time.
# This would be approximately constant if the resolution was sufficiently high.
function plot_maximum()
    max_sim = Array{Float64}(undef, num_snapshots)

    for i in 1:num_snapshots
        rho = npzread(joinpath(output_dir, "rho_$i.npy"))
        rhomax, maxindex = findmax(rho)
        max_sim[i] = rhomax
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = L"$t / \mathcal{T}$",
        ylabel = L"$\rho / \rho_c$",
        limits = (nothing, (0, nothing)),
    )
    scatter!(ax, output_times, max_sim, label = "simulation")

    fig
end

if MPI.Comm_rank(comm) == 0
    @info "Plotting" MPI.Comm_rank(comm)
    plot_maximum()
end
