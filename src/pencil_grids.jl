function psi_half_step!(Δt::Real, grids)
    @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 2 * grids.Φx[i])
    end
end

function psi_whole_step!(Δt::Real, grids)
    @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 1 * grids.Φx[i])
    end
end

function phi_whole_step!(Δt::Real, grids; a::Real=1.0)
    # TODO: not all part of Φ update

    grids.ψk .= grids.fft_plan * grids.ψx
    @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt/2 * grids.k[i]^2 / a^2)
    end
    grids.ψx .= grids.fft_plan \ grids.ψk

    @threads for i in eachindex(grids.ρx)
        grids.ρx[i] = abs2(grids.ψx[i])
    end

    grids.Φk .= grids.rfft_plan * grids.ρx
    @threads for i in eachindex(grids.Φk)
        grids.Φk[i] *= -4 * π / (a * grids.rk[i]^2)
    end
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= grids.rfft_plan \ grids.Φk
end

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
    "Real space distance array"
    dist
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

    function PencilGrids(dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan)
        n_dims = 3
        resol_tuple = size_global(dist)
        resol_tuple_realfft = (size_global(dist)[1] ÷ 2 + 1, size_global(dist)[2], size_global(dist)[3])

        for var in [dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk]
            @assert(ndims(var) == n_dims)
        end

        for var in [dist, k, ψx, ψk, ρx, Φx]
            @assert(size_global(var) == resol_tuple)
        end

        for var in [rk, ρk, Φk]
            @assert(size_global(var) == resol_tuple_realfft)
        end

        new(dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan)
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

    # MPI topology information
    comm = MPI.COMM_WORLD  # we assume MPI.Comm_size(comm) == 12
    Nproc = MPI.Comm_size(comm)

    # Let MPI_Dims_create choose the decomposition.
    proc_dims = let pdims = zeros(Int, 2)
        MPI.Dims_create!(Nproc, pdims)
        pdims[1], pdims[2]
    end

    # Plan a 3D complex-to-complex (c2c) FFT.
    fft_plan = PencilFFTPlan(resol_tuple, Transforms.FFT(), proc_dims, comm)

    # Allocate ψx and ψk arrays
    ψx = allocate_input(fft_plan)
    ψk = allocate_output(fft_plan)

    # Plan a 3D real-to-complex (r2c) FFT.
    rfft_plan = PencilFFTPlan(resol_tuple, Transforms.RFFT(), proc_dims, comm)

    ρx = allocate_input(rfft_plan)
    ρk = allocate_output(rfft_plan)

    Φx = allocate_input(rfft_plan)
    Φk = allocate_output(rfft_plan)

    dist = allocate_input(rfft_plan)
    dist_glob = global_view(dist)
    _dist = dist_array(length, resol)
    for I in CartesianIndices(dist_glob)
        dist_glob[I] = _dist[I]
    end

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
        dist,
        k, rk,
        ψx, ψk,
        ρx, ρk,
        Φx, Φk,
        fft_plan, rfft_plan,
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

    ψx_glob = global_view(grids.ψx)
    for I in CartesianIndices(ψx_glob)
        ψx_glob[I] = ψx[I]
    end

    grids.ρx .= abs2.(grids.ψx)
    a_init = 1  # TODO

    grids.ρk .= grids.rfft_plan * grids.ρx
    grids.Φk .= -4 * π * grids.ρk ./ (a_init * grids.rk.^2)
    grids.Φk[1, 1, 1] = 0
    grids.Φx .= grids.rfft_plan \ grids.Φk

    grids

end

function dist_array(length, resol::Integer)
    gridvec = range(
        -length / 2 + length / 2resol,
        +length / 2 - length / 2resol,
        length=resol
    )

    x = reshape(gridvec, resol, 1, 1)
    y = reshape(gridvec, 1, resol, 1)
    z = reshape(gridvec, 1, 1, resol)

    (x.^2 .+ y.^2 .+ z.^2).^0.5
end

"""
    k_vec(lengths, resols)

Calculate the Fourier frequencies of a box with side lengths `lengths` and resolutions `resols`

# Examples

```jldoctest
julia> using JultraDark: k_vec

julia> kvec = k_vec((2π, 2π, 2π), (4, 4, 4));

julia> kvec[1]
4-element AbstractFFTs.Frequencies{Float64}:
  0.0
  1.0
 -2.0
 -1.0

```
"""
function k_vec(lengths, resols)
    sample_rate = 2π .* resols ./ lengths

    kx = fftfreq(resols[1], sample_rate[1])

    ky = fftfreq(resols[2], sample_rate[2])
    kz = fftfreq(resols[3], sample_rate[3])

    kvec = (kx, ky, kz)
end

function rk_vec(lengths, resols)
    sample_rate = 2π .* resols ./ lengths

    kx = rfftfreq(resols[1], sample_rate[1])

    ky = fftfreq(resols[2], sample_rate[2])
    kz = fftfreq(resols[3], sample_rate[3])

    rkvec = (kx, ky, kz)
end

function k_norm(lengths, resols)
    kvec = k_vec(lengths, resols)

    kx = reshape(kvec[1], size(kvec[1])[1], 1, 1)
    ky = reshape(kvec[2], 1, size(kvec[2])[1], 1)
    kz = reshape(kvec[3], 1, 1, size(kvec[3])[1])

    (kx.^2 .+ ky.^2 .+ kz.^2).^0.5
end

function rk_norm(lengths, resols)
    kvec = rk_vec(lengths, resols)

    kx = reshape(kvec[1], size(kvec[1])[1], 1, 1)
    ky = reshape(kvec[2], 1, size(kvec[2])[1], 1)
    kz = reshape(kvec[3], 1, 1, size(kvec[3])[1])

    (kx.^2 .+ ky.^2 .+ kz.^2).^0.5
end
