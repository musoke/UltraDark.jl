using Test
using UltraDark: Grids, PencilGrids, output_grids, OutputConfig

for grid_type in [Grids, PencilGrids]
    for stat in [SummaryStatistics, SummaryStatisticsMeanMaxRms]
        @test stat(grid_type())
    end
end
