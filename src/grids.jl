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
        
        new(dist, k, rk, ψx, ψk, ρx, ρk, Φx, Φk)
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
        k_array(length, resol),
        rk_array(length, resol),
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

Create a grid with given ψ field, length `length` and resolution `resol` inferred
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
    resol_tuple = size(ψx)
    resol_tuple_realfft = (resol ÷ 2 + 1, resol, resol)
    
    ρx = abs2.(ψx)
    a_init = 1  # TODO

    ρk = rfft(ρx)
    Φk = -4 * π * ρk ./ (a_init * rk_array(length, resol).^2)
    Φk[1, 1, 1] = 0
    Φx = irfft(Φk, resol)

    Grids(
        dist_array(length, resol),
        k_array(length, resol),
        rk_array(length, resol),
        ψx,
        fft(ψx),
        ρx,
        ρk,
        Φx,
        Φk
    )
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

function k_array(length, resol::Integer)
    gridvec = fftfreq(resol, length/resol)

    kx = reshape(gridvec, resol, 1, 1)
    ky = reshape(gridvec, 1, resol, 1)
    kz = reshape(gridvec, 1, 1, resol)

    (kx.^2 .+ ky.^2 .+ kz.^2).^0.5
end

function rk_array(length, resol::Integer)
    resol_r = resol ÷ 2 + 1

    gridvec = fftfreq(resol, length/resol)
    r_gridvec = rfftfreq(resol, length/resol)

    kx = reshape(r_gridvec, resol_r, 1, 1)
    ky = reshape(gridvec, 1, resol, 1)
    kz = reshape(gridvec, 1, 1, resol)

    (kx.^2 .+ ky.^2 .+ kz.^2).^0.5
end