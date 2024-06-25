using UltraDark
using Test
using NPZ
using CairoMakie
using CSV

Threads.nthreads()

resol = 128
len = 10.0

grids = Grids((len, len, len / resol), (resol, resol, 1));

mass = 10
position_1 = [-len / 5, -len / 5, 0]
position_2 = [-len / 5, +len / 5, 0]
velocity = [1, 0, 0]
phase_1 = 0
phase_2 = π
t0 = 0

UltraDark.Initialise.add_fdm_soliton!(grids, mass, position_1, velocity, phase_1, t0)
UltraDark.Initialise.add_fdm_soliton!(grids, mass, position_2, velocity, phase_2, t0)

output_dir = joinpath(mktempdir(), "output", "2D")

output_times = 0:0.1:5
output_config = OutputConfig(output_dir, output_times);

options = Config.SimulationConfig();

@time simulate!(grids, options, output_config)

summary = CSV.File(joinpath(output_config.directory, "summary.csv"));

rho_last = npzread("$(output_config.directory)/rho_$(length(output_times)).npy");
δ_lims = extrema(rho_last .- 1)

fig_anim = Figure()
ax_anim = Axis(fig_anim[1, 1], aspect = DataAspect())

hidedecorations!(ax_anim)

cb_density = Colorbar(fig_anim[1, 2], limits = δ_lims, label = L"$\rho/\rho_{\text{crit}}$")

frame = Observable(1)

rho = lift(i -> npzread(joinpath(output_config.directory, "rho_$(i).npy"))[:, :, 1], frame)
δ = @lift($rho .- 1)

contourf!(ax_anim, grids.x[:, 1, 1], grids.y[1, :, 1], δ, levels = range(-1, δ_lims[2], 10))

record(fig_anim, "2d.mkv", 1:length(output_times); framerate = 5) do f
    frame[] = f
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
