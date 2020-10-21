@safetestset "correct max phase grad" begin
    using JultraDark:max_phase_grad

    n = 16

    @testset "purely real" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        @test max_phase_grad(a) == 0
    end

    @testset "purely imaginary" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n)) * im
        @test max_phase_grad(a) == 0
    end

    @testset "real with 1 row negated" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[1, :, :] *= -1
        @test max_phase_grad(a) ≈ π
    end

    @testset "real with 1 row negated" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[:, 3, :] *= -1
        @test max_phase_grad(a) ≈ π
    end

    @testset "real with 1 row negated" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[:, :, 5] *= -1
        @test max_phase_grad(a) ≈ π
    end

    @testset "real with 1 row imaginary" begin
        a = reshape(Vector{Complex}(1:n^3), (n, n, n))
        a[1, :, :] *= im
        @test max_phase_grad(a) ≈ π / 2
    end

end

@safetestset "terminate if grad too large" begin
    using JultraDark
    using JultraDark: evolve_to!, max_phase_grad

    n = 16

    # Set up grid with discontinuous phase
    grids = Grids(1.0, n)
    grids.ψx .= 1
    grids.ψx[1, 1, 1] = -1

    grids_orig = deepcopy(grids)

    options = JultraDark.Config.SimulationConfig(10, t->1)
    output_config = OutputConfig(mktempdir(), [1, 2])

    @test_throws "Phase gradient is too large to start" simulate(grids, options, output_config)
    @test grids.ψx == grids_orig.ψx

    @test_throws "Phase gradient is too large to continue" evolve_to!(0, 1, grids, output_config, options)
    @test grids.ψx == grids_orig.ψx

end
