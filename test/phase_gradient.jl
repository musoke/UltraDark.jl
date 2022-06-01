@safetestset "correct phase diff" begin
    using UltraDark: phase_diff

    n = 4

    @testset "correct size" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))

        @test size(phase_diff(a, [1, 0, 0])) == (n, n, n)
        @test size(phase_diff(a, [0, 1, 0])) == (n, n, n)
        @test size(phase_diff(a, [0, 0, 1])) == (n, n, n)
    end

    @testset "purely real" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        @test all(phase_diff(a, [1, 0, 0]) .== 0)
        @test all(phase_diff(a, [0, 1, 0]) .== 0)
        @test all(phase_diff(a, [0, 0, 1]) .== 0)
    end

    @testset "purely imaginary" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n)) * im
        @test all(phase_diff(a, [1, 0, 0]) .== 0)
        @test all(phase_diff(a, [0, 1, 0]) .== 0)
        @test all(phase_diff(a, [0, 0, 1]) .== 0)
    end

    @testset "real with 1 row negated" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[1, :, :] *= -1

        diff = zeros(Real, size(a))
        diff .= 0
        diff[1, :, :] .= -π
        diff[2, :, :] .= +π

        @test all(phase_diff(a, [1, 0, 0]) .≈ diff)
        @test all(phase_diff(a, [0, 1, 0]) .≈ 0)
        @test all(phase_diff(a, [0, 0, 1]) .≈ 0)
    end

    @testset "real with 1 row negated" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[:, 3, :] *= -1

        diff = zeros(Real, size(a))
        diff .= 0
        diff[:, 3, :] .= -π
        diff[:, 4, :] .= +π

        @test all(phase_diff(a, [1, 0, 0]) .≈ 0)
        @test all(phase_diff(a, [0, 1, 0]) .≈ diff)
        @test all(phase_diff(a, [0, 0, 1]) .≈ 0)
    end

    @testset "real with 1 row imaginary" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[:, :, 2] *= im

        diff = zeros(Real, size(a))
        diff .= 0
        diff[:, :, 2] .= -π / 2
        diff[:, :, 3] .= +π / 2

        @test all(phase_diff(a, [1, 0, 0]) .≈ 0)
        @test all(phase_diff(a, [0, 1, 0]) .≈ 0)
        @test all(phase_diff(a, [0, 0, 1]) .≈ diff)
    end

end

@safetestset "terminate if grad too large" begin
    using UltraDark
    using UltraDark: evolve_to!

    n = 16

    # Set up grid with discontinuous phase
    grids = Grids(1.0, n)
    grids.ψx .= 1
    grids.ψx[1, 1, 1] = -1

    grids_orig = deepcopy(grids)

    options = UltraDark.Config.SimulationConfig()
    output_config = OutputConfig(mktempdir(), [1, 2])

    @test_throws "Phase gradient is too large to start" simulate!(
        grids,
        options,
        output_config,
    )
    @test grids.ψx == grids_orig.ψx

    @test_throws "Phase gradient is too large to continue" evolve_to!(
        0,
        1,
        grids,
        output_config,
        options,
    )
    @test grids.ψx == grids_orig.ψx

end
