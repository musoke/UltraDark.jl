using Test
using UltraDark
using MPI

stats_list = [
              Summary.WallTime,
              Summary.SimulationTime,
              Summary.ScaleFactor,
              Summary.TimeStep,
              Summary.MeanDensity,
              Summary.MaxDensity,
              Summary.RmsDensityContrast,
              Summary.TotalMass,
             ]

sim_time = 0.0
a = 1.
Δt = 1e-1
constants = nothing

grids = Grids(1.0, 8)
pencilgrids = PencilGrids(1.0, 8)

for g in [grids, pencilgrids]
    g.ψx .= 1.
    @. g.ρx = abs2(g.ψx)
end

for s in stats_list
    stat_grids = s(sim_time, a, Δt, grids, constants)

    local_stat_pencilgrids = s(sim_time, a, Δt, pencilgrids, constants)
    global_stat_pencilgrids = MPI.Reduce(local_stat_pencilgrids, Summary.pool_summarystat, 0, pencilgrids.MPI_COMM)

    if s !== Summary.WallTime
        @assert stat_grids == global_stat_pencilgrids
    end
end
