using Test
using UltraDark

sim_time = 0.0
a = 1.
Δt = 1e-1
g = Grids(1.0, 8)
constants = nothing

for s in [
          Summary.WallTime,
          Summary.SimulationTime,
          Summary.ScaleFactor,
          Summary.TimeStep,
          Summary.MeanDensity,
          Summary.MaxDensity,
          Summary.RmsDensityContrast,
          Summary.TotalMass,
         ]

    s(sim_time, a, Δt, g, constants)
end
