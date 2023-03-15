"""
    AbstractGrids

Abstract type for grids containing simulation data
"""

abstract type AbstractGrids end

"""
struct containing grids used in a simulation

# Examples

```jldoctest
julia> using UltraDark

julia> len = 1;

julia> resol = 16;

julia> Grids(len, resol);

```
"""
struct Grids <: AbstractGrids
    "Array of x positions"
    x::Array{Float64,3}
    "Array of y positions"
    y::Array{Float64,3}
    "Array of z positions"
    z::Array{Float64,3}
    "Array of x Fourier modes"
    kx::Array{Float64,3}
    "Array of y Fourier modes"
    ky::Array{Float64,3}
    "Array of z Fourier modes"
    kz::Array{Float64,3}
    "Fourier space postition array"
    k::Array{Float64,3}
    "Array of x Fourier modes for use with `rfft`"
    rkx::Array{Float64,3}
    "Array of y Fourier modes for use with `rfft`"
    rky::Array{Float64,3}
    "Array of z Fourier modes for use with `rfft`"
    rkz::Array{Float64,3}
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

    "FFT plan for complex-to-complex transforms"
    fft_plan::Any
    "FFT plan for real-to-complex transforms"
    rfft_plan::Any

    "Indices at which k==0.0"
    k_vanish_indices::Vector{CartesianIndex{3}}
    "Indices at which rk==0.0"
    rk_vanish_indices::Any#::Vector{CartesianIndex{3}}

    function Grids(x, y, z, kx, ky, kz, k, rkx, rky, rkz, rk, ψx, ψk, ρx, ρk, Φx, Φk)
        n_dims = 3
        resol_tuple = (size(x)[1], size(y)[2], size(z)[3])
        resol_tuple_realfft = (size(x)[1] ÷ 2 + 1, size(y)[2], size(z)[3])

        for var in [x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk]
            @assert(ndims(var) == n_dims)
        end
        for var in [k, ψx, ψk, ρx, Φx]
            @assert(size(var) == resol_tuple)
        end
        for var in [rk, ρk, Φk]
            @assert(size(var) == resol_tuple_realfft)
        end
        @assert(size(x) == (resol_tuple[1], 1, 1))
        @assert(size(y) == (1, resol_tuple[2], 1))
        @assert(size(z) == (1, 1, resol_tuple[3]))

        FFTW.set_num_threads(Threads.nthreads())
        fft_plan = plan_fft(zeros(Complex{Float64}, resol_tuple), flags = FFTW.MEASURE)
        inv(fft_plan)

        rfft_plan = plan_rfft(zeros(Float64, resol_tuple), flags = FFTW.MEASURE)
        inv(rfft_plan)

        k_vanish_indices = findall(x -> x == 0.0, k)
        rk_vanish_indices = findall(x -> x == 0.0, rk)

        new(
            x,
            y,
            z,
            kx,
            ky,
            kz,
            k,
            rkx,
            rky,
            rkz,
            rk,
            ψx,
            ψk,
            ρx,
            ρk,
            Φx,
            Φk,
            fft_plan,
            rfft_plan,
            k_vanish_indices,
            rk_vanish_indices,
        )
    end
end

"""
    Grids(length, resol::Int)

Constructor for `Grids`

Create an empty grid with length `length` and resolution `resol`

# Examples

```jldoctest
julia> using UltraDark

julia> Grids(1.0, 64);

```
"""
function Grids(length, resol::Int)::Grids

    Grids((length, length, length), (resol, resol, resol))

end

"""
    Grids(length_tuple, resol_tuple::Tuple{Int, Int, Int})

Constructor for `Grids`

Create an empty `length[1]`x`length[2]`x`length[3]` grid with resolution
`resol[1]`x`resol[2]`x`resol[3]`.

# Examples

```jldoctest
julia> using UltraDark

julia> Grids((1.0, 1.0, 0.5), (64, 64, 32));

```
"""
function Grids(length_tuple, resol_tuple::Tuple{Int,Int,Int})::Grids

    resol_tuple_realfft = (resol_tuple[1] ÷ 2 + 1, resol_tuple[2], resol_tuple[3])

    @assert length_tuple[1] / resol_tuple[1] ≈ length_tuple[2] / resol_tuple[2]
    @assert length_tuple[1] / resol_tuple[1] ≈ length_tuple[3] / resol_tuple[3]


    ψx = zeros(Complex{Float64}, resol_tuple)
    ψk = zeros(Complex{Float64}, resol_tuple)

    ρx = zeros(Float64, resol_tuple)
    ρk = zeros(Complex{Float64}, resol_tuple_realfft)

    Φx = zeros(Float64, resol_tuple)
    Φk = zeros(Complex{Float64}, resol_tuple_realfft)

    x = reshape(
        range(
            -length_tuple[1] / 2 + length_tuple[1] / 2resol_tuple[1],
            +length_tuple[1] / 2 - length_tuple[1] / 2resol_tuple[1],
            length = resol_tuple[1],
        ),
        :,
        1,
        1,
    )

    y = reshape(
        range(
            -length_tuple[2] / 2 + length_tuple[2] / 2resol_tuple[2],
            +length_tuple[2] / 2 - length_tuple[2] / 2resol_tuple[2],
            length = resol_tuple[2],
        ),
        1,
        :,
        1,
    )

    z = reshape(
        range(
            -length_tuple[3] / 2 + length_tuple[3] / 2resol_tuple[3],
            +length_tuple[3] / 2 - length_tuple[3] / 2resol_tuple[3],
            length = resol_tuple[3],
        ),
        1,
        1,
        :,
    )

    kvec = k_vec(length_tuple, resol_tuple)

    kx = reshape(kvec[1], :, 1, 1)
    ky = reshape(kvec[2], 1, :, 1)
    kz = reshape(kvec[3], 1, 1, :)
    k_norm = (kx .^ 2 .+ ky .^ 2 .+ kz .^ 2) .^ 0.5

    rkvec = rk_vec(length_tuple, resol_tuple)

    rkx = reshape(rkvec[1], :, 1, 1)
    rky = reshape(rkvec[2], 1, :, 1)
    rkz = reshape(rkvec[3], 1, 1, :)
    rk_norm = (rkx .^ 2 .+ rky .^ 2 .+ rkz .^ 2) .^ 0.5

    Grids(x, y, z, kx, ky, kz, k_norm, rkx, rky, rkz, rk_norm, ψx, ψk, ρx, ρk, Φx, Φk)
end

"""
    k_vec(lengths, resols)

Calculate the Fourier frequencies of a box with side lengths `lengths` and resolutions `resols`

# Examples

```jldoctest
julia> using UltraDark: k_vec

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

    kx = reshape(kvec[1], :, 1, 1)
    ky = reshape(kvec[2], 1, :, 1)
    kz = reshape(kvec[3], 1, 1, :)

    (kx .^ 2 .+ ky .^ 2 .+ kz .^ 2) .^ 0.5
end

function rk_norm(lengths, resols)
    kvec = rk_vec(lengths, resols)

    kx = reshape(kvec[1], :, 1, 1)
    ky = reshape(kvec[2], 1, :, 1)
    kz = reshape(kvec[3], 1, 1, :)

    (kx .^ 2 .+ ky .^ 2 .+ kz .^ 2) .^ 0.5
end

"""
    dV(grids)

Calculate the volume of each grid cell

# Examples

```jldoctest
julia> using UltraDark

julia> box_length = 1.0;

julia> resol = 16;

julia> g = Grids(box_length, resol);

julia> dV(g) * resol^3 == box_length^3
true
```
"""
function dV(grids)
    diff(grids.x, dims = 1)[1] * diff(grids.y, dims = 2)[1] * diff(grids.z, dims = 3)[1]
end

"""
    radius_spherical(grids, r0)

Calculate the radial coordinate in a spherical coordinate system

# Examples

```jldoctest
julia> using UltraDark

julia> import UltraDark: radius_spherical, polar_angle, azimuthal_angle

julia> box_length = 1.0;

julia> resol = 16;

julia> g = Grids(box_length, resol);

julia> all(radius_spherical(g) .* sin.(polar_angle(g)) .* cos.(azimuthal_angle(g)) .≈ g.x)
true

julia> all(
           radius_spherical(g, (1.0, 0.0, 0.0)) .* sin.(polar_angle(g, (1.0, 0.0, 0.0))) .*
           cos.(azimuthal_angle(g, (1.0, 0.0, 0.0))) .+ 1.0 .≈ g.x,
       )
true
```
"""
function radius_spherical(x, y, z)
    rs = (x^2 + y^2 + z^2)^(1 // 2)
end

function radius_spherical(grids)
    radius_spherical.(grids.x, grids.y, grids.z)
end

function radius_spherical(grids, r0)
    radius_spherical.(grids.x .- r0[1], grids.y .- r0[2], grids.z .- r0[3])
end

"""
    polar_angle(grids, r0)

Calculate the polar angle in spherical coordinates

This is \\theta in conventional physics notation.
"""
function polar_angle(x, y, z)
    rs = radius_spherical(x, y, z)
    acos(z / rs)
end

function polar_angle(grids)
    polar_angle.(grids.x, grids.y, grids.z)
end

function polar_angle(grids, r0)
    polar_angle.(grids.x .- r0[1], grids.y .- r0[2], grids.z .- r0[3])
end


"""
    polar_angle(grids, r0)

Calculate the azimuthal angle in spherical or cylindrical coordinates

This is \\phi in conventional physics notation.
"""
function azimuthal_angle(x, y, z)
    atan(y, x)
end

function azimuthal_angle(grids)
    azimuthal_angle.(grids.x, grids.y, grids.z)
end

function azimuthal_angle(grids, r0)
    azimuthal_angle.(grids.x .- r0[1], grids.y .- r0[2], grids.z .- r0[3])
end

"""
    radius_cylindrical(grids, r0)

Calculate the radial coordinate in cylindrical coordinates

# Examples

```jldoctest
julia> using UltraDark

julia> import UltraDark: radius_cylindrical, azimuthal_angle

julia> box_length = 1.0;

julia> resol = 16;

julia> g = Grids(box_length, resol);

julia> all(radius_cylindrical(g) .* cos.(azimuthal_angle(g)) .≈ g.x)
true

julia> all(
           radius_cylindrical(g, (0.0, 0.0, 1.0)) .* cos.(azimuthal_angle(g, (0.0, 0.0, 1.0))) .≈
           g.x,
       )
true
```
"""
function radius_cylindrical(x, y, z)
    (x^2 + y^2)^(1 // 2)
end

function radius_cylindrical(grids)
    radius_cylindrical.(grids.x, grids.y, grids.z)
end

function radius_cylindrical(grids, r0)
    radius_cylindrical.(grids.x .- r0[1], grids.y .- r0[2], grids.z .- r0[3])
end

"""
    mass(grids)
    mass(grids, rho)

Calculate total mass of a density field

# Examples

```jldoctest
julia> using UltraDark

julia> g = Grids(1.0, 16);

julia> g.ρx .= 0.0;

julia> g.ρx[1, 1, 1] = 1.0;

julia> UltraDark.mass(g) == 1.0 * (1.0 / 16)^3
true
```
"""
function mass(grids, rho)
    Folds.sum(rho) * dV(grids)
end

function mass(grids)
    mass(grids, grids.ρx)
end

"""
    E_grav(grids)
    E_grav(grids, psi)

Gravitational potential energy
"""
function E_grav(grids, psi)
    Folds.sum(grids.Φx .* abs2.(psi)) * dV(grids) / 2.0
end

function E_grav(grids::AbstractGrids)
    E_grav(grids, grids.ψx)
end

"""
    E_kq(grids)
    E_kq(grids, psi)

Sum of kinetic and quantum energies
"""
function E_kq(grids)
    E_kq(grids, grids.ψx)
end

function E_kq(grids, psi)

    f = grids.fft_plan * psi
    @inbounds @threads for i in eachindex(f)
        f[i] *= -grids.k[i]^2
    end
    g = grids.fft_plan \ f
    f = nothing

    out = Folds.sum(@. real(-0.5 * conj(psi) * g)) * dV(grids)
end

"""
    E_total(grids; constants=nothing)

Total energy of the scalar field: the sum of the kinetic and quantum energies.
"""
function E_total(grids; constants = nothing)
    E_grav(grids) + E_kq(grids)
end

"""
    angular_momentum_density(grids)
    angular_momentum_density(grids, ψx, ρx)

Calculate angular momentum density at each grid point

## Returns

L: AbstractArray with dimensions 3 x resol_x x resol_y x resol_z
"""
function angular_momentum_density(grids)
    angular_momentum_density(grids, grids.ψx, grids.ρx)
end

function angular_momentum_density(grids, ψx, ρx)

    x = grids.x
    y = grids.y
    z = grids.z

    momentum = zeros(3, size(ρx)...)
    momentum[1, :, :, :] .= phase_diff(ψx, (1, 0, 0))
    momentum[2, :, :, :] .= phase_diff(ψx, (0, 1, 0))
    momentum[3, :, :, :] .= phase_diff(ψx, (0, 0, 1))

    momentum .*= reshape(ρx, 1, size(ρx)...)

    angular_momentum = similar(momentum)

    @inbounds @threads for I in CartesianIndices(angular_momentum)
        _, i, j, k = Tuple(I)
        # angular_momentum[:, i, j, k] = cross([grids.x[i], grids.y[j], grids.z[k]], momentum[:, i, j, k])
        angular_momentum[1, i, j, k] =
            y[j] * momentum[3, i, j, k] - z[k] * momentum[2, i, j, k]
        angular_momentum[2, i, j, k] =
            z[k] * momentum[1, i, j, k] - x[i] * momentum[3, i, j, k]
        angular_momentum[3, i, j, k] =
            x[i] * momentum[2, i, j, k] - y[j] * momentum[1, i, j, k]
    end

    angular_momentum
end

"""
    angular_momentum(grids)
    angular_momentum(grids, ψx, ρx)

Calculate total angular momentum

## Returns

L: AbstractArray with length 3
"""
function angular_momentum(grids)
    angular_momentum(grids, grids.ψx, grids.ρx)
end

function angular_momentum(grids, ψx, ρx)
    L = angular_momentum_density(grids, ψx, ρx)
    reshape(sum(L, dims = [2, 3, 4]), 3)
end
