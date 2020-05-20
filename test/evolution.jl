using JultraDark
import JultraDark

@test JultraDark.actual_timestep(0.5, 1) ≈ 0.5
@test JultraDark.actual_timestep(0.7, 1) ≈ 0.5
@test JultraDark.actual_timestep(0.7, 2) ≈ 2/3

# Test that evolve to arrives at right time
@testset "Evolve from $t_begin to $t_end" for t_begin in [0, 0.123, 1], t_end in [1.23, 2]
    grids = Grids(1.0, 16)
    output_config = OutputConfig(mktempdir())
    @test JultraDark.evolve_to!(t_begin, t_end, grids, output_config, t -> 1) ≈ t_end
end