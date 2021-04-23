using BenchmarkTools
using UltraDark: Grids, max_phase_grad

resol = 64

grids = Grids(1.0, resol)

@benchmark max_phase_grad(grids.ψx)
