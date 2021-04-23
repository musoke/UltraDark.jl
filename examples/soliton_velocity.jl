using UltraDark
using NPZ
using Plots

include(joinpath(@__DIR__, "init_soliton.jl"))

output_dir = mktempdir()
output_dir = "output/vel"

num_snapshots = 20
output_times = LinRange(0, 5, num_snapshots)

output_config = OutputConfig(output_dir, output_times; box=true, slice=false)
options = Config.SimulationConfig(10, t->1)

mass = 10
position0 = [-4, 0, 0]
velocity = [1, 0, 0]
phase = 0
t0 = 0

resol = 64
grids = Grids(10.0, resol)
add_soliton(grids, mass, position0, velocity, phase, t0)

simulate(grids, options, output_config) == nothing

pos_the = reshape(position0, (1, :)) .+ (reshape(velocity, (1, :)) .* reshape(output_times, (:, 1)))

pos_sim = Array{Float64}(undef, (num_snapshots, 3))
maxs = Array{Float64}(undef, num_snapshots)

for i in 1:num_snapshots
    rho = npzread(joinpath(output_dir, "rho_$i.npy"))

    rhomax, maxindex = findmax(rho)

    maxs[i] = rhomax
    pos_sim[i, :] .= [grids.x[maxindex[1]], grids.y[maxindex[2]], grids.z[maxindex[3]]]

end

yerr = abs(grids.x[2]-grids.x[1])/2
plot_pos = plot(legend=:bottomright)
scatter!(output_times, pos_sim[:, 1], yerr=yerr, markershape=:auto, label="simulation")
plot!(output_times, pos_the[:, 1], label="theory")
plot!(xlabel=raw"$t$", ylabel=raw"$x$")

savefig(plot_pos, "soliton_position.pdf")

plot_max = plot(legend=:bottomright)
plot!(output_times, maxs, label="theory")

anim = @animate for i ∈ 1:num_snapshots
    rho = npzread(joinpath(output_dir, "rho_$i.npy"))
    psi = npzread(joinpath(output_dir, "psi_$i.npy"))

    density_plot = contourf(
                            rho[:, :, resol÷2],
                            title=raw"$\rho$",
                            clims=(0, maximum(maxs)),
                            showaxis=false,
                            ticks=false,
                            aspect_ratio=:equal,
                           )
    phase_plot = contourf(
                          angle.(psi[:, :, resol÷2]),
                          title="phase",
                          clims=(-π, +π),
                          showaxis=false,
                          ticks=false,
                          aspect_ratio=:equal,
                         )

    plot(density_plot, phase_plot, )
end

gif(anim, "soliton_velocity.gif", fps = 15)
