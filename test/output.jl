using Test
using UltraDark: Grids, PencilGrids, OutputConfig
using UltraDark: Output.output_grids

for grid_type in [Grids]
    grids = grid_type(1.0, 16)

    for (name, output_config) in [
        ("default", OutputConfig(mktempdir(), [])),
        (
            "nothing",
            OutputConfig(
                mktempdir(),
                [];
                box = false,
                slice = false,
                psi = true,
                rho = true,
            ),
        ),
        (
            "box",
            OutputConfig(
                mktempdir(),
                [];
                box = true,
                slice = false,
                psi = true,
                rho = true,
            ),
        ),
        (
            "slice",
            OutputConfig(mktempdir(), []; box = true, slice = true, psi = true, rho = true),
        ),
        (
            "box & slice",
            OutputConfig(mktempdir(), []; box = true, slice = true, psi = true, rho = true),
        ),
        (
            "box & slice, npy & h5",
            OutputConfig(
                mktempdir(),
                [];
                box = true,
                slice = true,
                psi = true,
                rho = true,
                npy = true,
                h5 = true,
            ),
        ),
    ]
        mkpath(output_config.directory)
        @test output_grids(grids, output_config, 0) == nothing
    end
end

for grid_type in [PencilGrids]
    grids = grid_type(1.0, 16)

    for (name, output_config) in [
        ("default", OutputConfig(mktempdir(), [])),
        (
            "nothing",
            OutputConfig(
                mktempdir(),
                [];
                box = false,
                slice = false,
                psi = true,
                rho = true,
            ),
        ),
        (
            "box",
            OutputConfig(
                mktempdir(),
                [];
                box = true,
                slice = false,
                psi = true,
                rho = true,
            ),
        ),
    ]
        mkpath(output_config.directory)
        @test output_grids(grids, output_config, 0) == nothing
    end

    for (name, output_config) in [
        (
            "slice",
            OutputConfig(mktempdir(), []; box = true, slice = true, psi = true, rho = true),
        ),
        (
            "box & slice",
            OutputConfig(mktempdir(), []; box = true, slice = true, psi = true, rho = true),
        ),
        (
            "box & h5",
            OutputConfig(mktempdir(), []; box = true, psi = true, npy = false, h5 = true),
        ),
    ]
        mkpath(output_config.directory)
        @test_broken output_grids(grids, output_config, 0) == nothing
    end
end
