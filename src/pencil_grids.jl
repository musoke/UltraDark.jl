"""
struct containing grids used in a simulation

# Examples


```jldoctest
julia> using JultraDark

julia> len = 1;

julia> resol = 16;

julia> PencilGrids(len, resol);

```
"""
struct PencilGrids
    "Array of x positions"
    x::Array{Float64,3}
    "Array of y positions"
    y::Array{Float64,3}
    "Array of z positions"
    z::Array{Float64,3}
    "Fourier space postition array"
    k
    "Fourier space postition array for use with `rfft`"
    rk
    "ψ field"
    ψx
    "ψ field in Fourier space"
    ψk
    "density field ρ"
    ρx
    "density field ρ in Fourier space"
    ρk
    "gravitational potential field Φ"
    Φx
    "gravitational potential field Φ in fourier space"
    Φk
    fft_plan
    rfft_plan
    MPI_COMM

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

        new(x, y, x, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan, MPI_COMM)
    end
end

"""
    PencilGrids(length::Real, resol::Integer)

Constructor for `PencilGrids`

Create an empty grid with length `length` and resolution `resol`.  Uses `PencilFFTs` to create `PencilArrays`.

# Examples

```jldoctest
julia> using JultraDark

julia> PencilGrids(1.0, 64);

```
"""
function PencilGrids(length::Real, resol::Integer)::PencilGrids

    resol_tuple = (resol, resol, resol)
    resol_tuple_realfft = (resol ÷ 2 + 1, resol, resol)

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

    gridvec = range(
        -length / 2 + length / 2resol,
        +length / 2 - length / 2resol,
        length=resol
    )

    x = reshape(gridvec, resol, 1, 1)
    y = reshape(gridvec, 1, resol, 1)
    z = reshape(gridvec, 1, 1, resol)

    k = similar(ψk, Float64)
    k_glob = global_view(k)
    _k = k_norm((length, length, length), (resol, resol, resol))
    for I in CartesianIndices(k_glob)
        k_glob[I] = _k[I]
    end

    rk = similar(ρk, Float64)
    rk_glob = global_view(rk)
    _rk = rk_norm((length, length, length), (resol, resol, resol))
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

"""
    PencilGrids(ψx::Array{Complex{Float64}}, length::Real)

Constructor for `PencilGrids`

Create a grid with given ψ field, length `length` and resolution inferred
from `ψx`

# Examples

Can be contructed from a ψ field and box length,
```jldoctest
julia> using JultraDark

julia> ψ = zeros(Complex{Float64}, 16, 16, 16);

julia> len = 1;

julia> PencilGrids(ψ, len);

```
"""
function PencilGrids(ψx::Array{Complex{Float64}}, length::Real)::PencilGrids
    @assert(
        ndims(ψx) == 3,
        "Invalid ψ: only three dimensions supported"
    )
    @assert(
        size(ψx, 1) == size(ψx, 2) == size(ψx, 3),
        "Invalid ψ: heterogenous resolutions"
    )

    resol = size(ψx, 1)
    grids = PencilGrids(length, resol)

    # Assign each ψx element desired value
    ψx_glob = global_view(grids.ψx)
    for I in CartesianIndices(ψx_glob)
        ψx_glob[I] = ψx[I]
    end

    grids

end

"""
    max_time_step(grids, a)

Calculate an upper bound on the time step
"""
function max_time_step(grids::PencilGrids, a)
    # Find maximum on local grid
    local_max_time_step_gravity = 2π / maximum(abs.(grids.Φx))
    local_max_time_step_pressure = 2π * 2 / maximum(grids.k)^2 * a^2  # TODO: cache k_max

    # Maximize over other grids
    max_time_step_gravity = MPI.Allreduce(local_max_time_step_gravity, MPI.MIN, grids.MPI_COMM)
    max_time_step_pressure = MPI.Allreduce(local_max_time_step_pressure, MPI.MIN, grids.MPI_COMM)

    @assert isfinite(max_time_step_gravity)
    @assert isfinite(max_time_step_pressure)

    time_step = min(max_time_step_gravity, max_time_step_pressure)

    time_step
end

"""
    phase_diff(field::PencilArray, dir)

Compute point-to-point difference of phase on a grid along a direction

Returns an array of size `(size(field)[1], size(field)[2], size(field)[2])`
containing gradients in direction `dir`.
"""
function phase_diff(field::PencilArray, dir)
    out = similar(field, Real)

    # TODO diff at boundaries
    out .= 0
    out[1:end-dir[1], 1:end-dir[2], 1:end-dir[3]] .= angle.(circshift(field[1:end-dir[1], 1:end-dir[2], 1:end-dir[3]], dir) ./ field[1+dir[1]:end, 1+dir[2]:end, 1+dir[3]:end])

    out
end
