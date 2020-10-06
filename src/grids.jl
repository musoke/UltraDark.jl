"""
struct containing grids used in a simulation

# Examples


```jldoctest
julia> using JultraDark

julia> len = 1;

julia> resol = 16;

julia> Grids(len, resol);

```
"""
struct Grids
    "Real space distance array"
    dist::Array{Float64,3}
    "Fourier space postition array"
    k::Array{Float64,3}
    "Fourier space postition array for use with `rfft`"
    rk::Array{Float64,3}
    "ψ field"
    ψx::Array{Complex{Float64},3}
    "ψ field in Fourier space"
    ψk::Array{Complex{Float64},3}
    "density field ρ"
    ρx::Array{Float64,3}
    "density field ρ in Fourier space"
    ρk::Array{Complex{Float64},3}
    "gravitational potential field Φ"
    Φx::Array{Float64,3}
    "gravitational potential field Φ in fourier space"
    Φk::Array{Complex{Float64},3}
    fft_plan
    rfft_plan

    function Grids(dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk)
        n_dims = 3
        resol_tuple = size(dist)
        resol_tuple_realfft = (size(dist, 1) ÷ 2 + 1, size(dist, 2), size(dist, 3))

        for var in [dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk]
            @assert(ndims(var) == n_dims)
        end
        for var in [dist, k, ψx, ψk, ρx, Φx]
            @assert(size(var) == resol_tuple)
        end
        for var in [rk, ρk, Φk]
            @assert(size(var) == resol_tuple_realfft)
        end

        fft_plan = plan_fft(
            zeros(Complex{Float64}, resol_tuple),
            flags=FFTW.MEASURE
        )
        inv(fft_plan)

        rfft_plan = plan_rfft(
            zeros(Float64, resol_tuple),
            flags=FFTW.MEASURE
        )
        inv(rfft_plan)

        new(dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan)
    end
end

"""
    Grids(length::Real, resol::Integer)

Constructor for `Grids`

Create an empty grid with length `length` and resolution `resol`

# Examples

```jldoctest
julia> using JultraDark

julia> Grids(1.0, 64);

```
"""
function Grids(length::Real, resol::Integer)::Grids

    resol_tuple = (resol, resol, resol)
    resol_tuple_realfft = (resol ÷ 2 + 1, resol, resol)

    ψx = zeros(Complex{Float64}, resol_tuple)
    ψk = zeros(Complex{Float64}, resol_tuple)

    ρx = zeros(Float64, resol_tuple)
    ρk = zeros(Complex{Float64}, resol_tuple_realfft)

    Φx = zeros(Float64, resol_tuple)
    Φk = zeros(Complex{Float64}, resol_tuple_realfft)

    Grids(
        dist_array(length, resol),
        k_norm((length, length, length), (resol, resol, resol)),
        rk_norm((length, length, length), (resol, resol, resol)),
        ψx,
        ψk,
        ρx,
        ρk,
        Φx,
        Φk,
    )
end

"""
    Grids(ψx::Array{Complex{Float64}}, length::Real)

Constructor for `Grids`

Create a grid with given ψ field, length `length` and resolution inferred
from `ψx`

# Examples

Can be contructed from a ψ field and box length,
```jldoctest
julia> using JultraDark

julia> ψ = zeros(Complex{Float64}, 16, 16, 16);

julia> len = 1;

julia> Grids(ψ, len);

```
"""
function Grids(ψx::Array{Complex{Float64}}, length::Real)::Grids
    @assert(
        ndims(ψx) == 3,
        "Invalid ψ: only three dimensions supported"
    )
    @assert(
        size(ψx, 1) == size(ψx, 2) == size(ψx, 3),
        "Invalid ψ: heterogenous resolutions"
    )

    resol = size(ψx, 1)
    grids = Grids(length, resol)

    grids.ψx .= ψx

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

function psi_half_step!(Δt::Real, grids)
    @fastmath @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 2 * grids.Φx[i])
    end
end

function psi_whole_step!(Δt::Real, grids)
    @fastmath @inbounds @threads for i in eachindex(grids.ψx)
        grids.ψx[i] *= exp(- im * Δt / 1 * grids.Φx[i])
    end
end

function phi_whole_step!(Δt::Real, grids; a::Real=1.0)
    # TODO: not all part of Φ update

    mul!(grids.ψk, grids.fft_plan, grids.ψx)
    @fastmath @inbounds @threads for i in eachindex(grids.ψk)
        grids.ψk[i] *= exp(-im * Δt/2 * grids.k[i]^2 / a^2)
    end
    ldiv!(grids.ψx, grids.fft_plan, grids.ψk)

    @fastmath @inbounds @threads for i in eachindex(grids.ρx)
        grids.ρx[i] = abs2(grids.ψx[i])
    end

    mul!(grids.Φk, grids.rfft_plan, grids.ρx)
    @fastmath @inbounds @threads for i in eachindex(grids.Φk)
        grids.Φk[i] *= -4 * π / (a * grids.rk[i]^2)
    end
    grids.Φk[1, 1, 1] = 0
    ldiv!(grids.Φx, grids.rfft_plan, grids.Φk)
end
