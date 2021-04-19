using JultraDark
using Test
using MPI

@testset "Actual time step" begin
    @test JultraDark.actual_time_step(0.5, 1, 2)[1] ≈ 0.5
    @test JultraDark.actual_time_step(0.5, 1, 2)[2] == 2

    @test JultraDark.actual_time_step(0.5, 1, 3)[1] ≈ 0.5
    @test JultraDark.actual_time_step(0.5, 1, 3)[2] == 2

    @test JultraDark.actual_time_step(0.01, 1, 10)[1] ≈ 0.01
    @test JultraDark.actual_time_step(0.01, 1, 10)[2] == 10

    @test JultraDark.actual_time_step(0.7, 1, 2)[1] ≈ 0.5
    @test JultraDark.actual_time_step(0.7, 1, 2)[2] == 2

    @test JultraDark.actual_time_step(0.7, 2, 3)[1] ≈ 2/3
    @test JultraDark.actual_time_step(0.7, 2, 3)[2] == 3
end

for grid_type in [Grids, PencilGrids]
    @testset "Take n steps, $grid_type" begin
        grids = grid_type(1.0, 16)
        output_config = OutputConfig(mktempdir(), [])
        @test JultraDark.take_steps!(grids, 0, 1, 10, output_config, t -> 1) == 10.0
    end
end

# Test that evolve to arrives at right time
for grid_type in [Grids, PencilGrids]
    @testset "Evolve from $t_begin to $t_end, $grid_type" for t_begin in [0, 0.123, 1], t_end in [1.23, 2]
        grids = grid_type(1.0, 16)

        a = 1
        grids.ψx .= (grids.x.^2 .+ grids.y.^2 .+ grids.z.^2).^0.5 ./ 1e9 # Set ψx to something non-zero
        grids.ψk .= (grids.fft_plan * grids.ψx)
        grids.ρx .= abs2.(grids.ψx)
        grids.Φk .= -4 * π * (grids.rfft_plan * grids.ρx) ./ (a * grids.rk.^2)
        grids.Φk[1, 1, 1] = 0
        grids.Φx .= grids.rfft_plan \ grids.Φk

        output_config = OutputConfig(mktempdir(), [])
        options = JultraDark.Config.SimulationConfig(10, t->1)
        @test JultraDark.evolve_to!(t_begin, t_end, grids, output_config, options) ≈ t_end
    end
end
