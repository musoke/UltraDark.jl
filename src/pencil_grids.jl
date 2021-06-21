"""
struct containing grids used in a simulation

# Examples


```jldoctest
julia> using UltraDark

julia> len = 1;

julia> resol = 16;

julia> PencilGrids(len, resol);

```
"""
struct PencilGrids{K, RX, RK, CX, CK, FFT, RFFT, M}
    "Array of x positions"
    x::Array{Float64,3}
    "Array of y positions"
    y::Array{Float64,3}
    "Array of z positions"
    z::Array{Float64,3}
    "Fourier space postition array"
    k::K
    "Fourier space postition array for use with `rfft`"
    rk::K
    "ψ field"
    ψx::CX
    "ψ field in Fourier space"
    ψk::CK
    "density field ρ"
    ρx::RX
    "density field ρ in Fourier space"
    ρk::RK
    "gravitational potential field Φ"
    Φx::RX
    "gravitational potential field Φ in fourier space"
    Φk::RK
    fft_plan::FFT
    rfft_plan::RFFT
    MPI_COMM::M

    function PencilGrids(x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan, MPI_COMM)
        n_dims = 3
        resol_tuple = (size(x)[1], size(y)[2], size(z)[3])
        resol_tuple_realfft = (size(x)[1] ÷ 2 + 1, size(y)[2], size(z)[3])

        for var in [x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk]
            @assert(ndims(var) == n_dims)
        end

        for var in [k, ψx, ψk, ρx, Φx]
            @assert(size_global(var) == resol_tuple)
        end

        for var in [rk, ρk, Φk]
            @assert(size_global(var) == resol_tuple_realfft)
        end

        @assert(size(x) == (resol_tuple[1], 1, 1))
        @assert(size(y) == (1, resol_tuple[2], 1))
        @assert(size(z) == (1, 1, resol_tuple[3]))

        K = typeof(k)
        RX = typeof(Φx)
        RK = typeof(Φk)
        CX = typeof(ψx)
        CK = typeof(ψk)
        FFT = typeof(fft_plan)
        RFFT = typeof(rfft_plan)
        M = typeof(MPI_COMM)

        new{K, RX, RK, CX, CK, FFT, RFFT, M}(x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan, MPI_COMM)
    end
end

"""
    PencilGrids(length, resol)

Constructor for `PencilGrids`

Create an empty grid with length `length` and resolution `resol`.  Uses `PencilFFTs` to create `PencilArrays`.

# Examples

```jldoctest
julia> using UltraDark

julia> PencilGrids(1.0, 64);

```
"""
function PencilGrids(length, resol::Int)::PencilGrids

    PencilGrids((length, length, length), (resol, resol, resol))

end

"""
    PencilGrids(length_tuple, resol_tuple::Tuple{Int, Int, Int})

Constructor for `PencilGrids`

Create an empty `length[1]`x`length[2]`x`length[3]` grid with resolution
`resol[1]`x`resol[2]`x`resol[3]`.
Each grid is a `PencilArray`, allowing multiprocess FFTs.

# Examples

```jldoctest
julia> using UltraDark

julia> PencilGrids((1.0, 1.0, 0.5), (64, 64, 32));

```
"""
function PencilGrids(length_tuple, resol_tuple::Tuple{Int, Int, Int})::PencilGrids

    resol_tuple_realfft = (resol_tuple[1] ÷ 2 + 1, resol_tuple[2], resol_tuple[3])

    @assert length_tuple[1] / resol_tuple[1] ≈ length_tuple[2] / resol_tuple[2]
    @assert length_tuple[1] / resol_tuple[1] ≈ length_tuple[3] / resol_tuple[3]

    if ~MPI.Initialized()
        MPI.Init()
    end

    # MPI topology information
    comm = MPI.COMM_WORLD
    Nproc = MPI.Comm_size(comm)

    # Let MPI_Dims_create choose the decomposition.
    proc_dims = let pdims = zeros(Int, 2)
        MPI.Dims_create!(Nproc, pdims)
        pdims[1], pdims[2]
    end

    FFTW.set_num_threads(Threads.nthreads())
    # Plan a 3D complex-to-complex (c2c) FFT.
    fft_plan = PencilFFTPlan(resol_tuple, Transforms.FFT(), proc_dims, comm)

    # Allocate ψx and ψk arrays
    ψx = allocate_input(fft_plan)
    ψx .= 0
    ψk = allocate_output(fft_plan)

    # Plan a 3D real-to-complex (r2c) FFT.
    rfft_plan = PencilFFTPlan(resol_tuple, Transforms.RFFT(), proc_dims, comm)

    ρx = allocate_input(rfft_plan)
    ρk = allocate_output(rfft_plan)

    Φx = allocate_input(rfft_plan)
    Φk = allocate_output(rfft_plan)

    x = reshape(
                range(
                      -length_tuple[1] / 2 + length_tuple[1] / 2resol_tuple[1],
                      +length_tuple[1] / 2 - length_tuple[1] / 2resol_tuple[1],
                      length=resol_tuple[1]
                     ),
                resol_tuple[1], 1, 1
    )

    y = reshape(
                range(
                      -length_tuple[2] / 2 + length_tuple[2] / 2resol_tuple[2],
                      +length_tuple[2] / 2 - length_tuple[2] / 2resol_tuple[2],
                      length=resol_tuple[2]
                     ),
                1, resol_tuple[2], 1
    )

    z = reshape(
                range(
                      -length_tuple[3] / 2 + length_tuple[3] / 2resol_tuple[3],
                      +length_tuple[3] / 2 - length_tuple[3] / 2resol_tuple[3],
                      length=resol_tuple[3]
                     ),
                1, 1, resol_tuple[3]
    )

    k = similar(ψk, Float64)
    k_glob = global_view(k)
    _k = k_norm(length_tuple, resol_tuple)
    for I in CartesianIndices(k_glob)
        k_glob[I] = _k[I]
    end

    rk = similar(ρk, Float64)
    rk_glob = global_view(rk)
    _rk = rk_norm(length_tuple, resol_tuple)
    for I in CartesianIndices(rk_glob)
        rk_glob[I] = _rk[I]
    end

    PencilGrids(
        x, y, z,
        k, rk,
        ψx, ψk,
        ρx, ρk,
        Φx, Φk,
        fft_plan, rfft_plan,
        comm
    )
end
