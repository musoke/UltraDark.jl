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
struct Grids
    "Array of x positions"
    x::Array{Float64,3}
    "Array of y positions"
    y::Array{Float64,3}
    "Array of z positions"
    z::Array{Float64,3}
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

    "FFT plan for complex-to-complex transforms"
    fft_plan
    "FFT plan for real-to-complex transforms"
    rfft_plan

    "Indices at which k==0.0"
    k_vanish_indices::Vector{CartesianIndex{3}}
    "Indices at which rk==0.0"
    rk_vanish_indices#::Vector{CartesianIndex{3}}

    function Grids(x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk)
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

        k_vanish_indices = findall(x -> x==0.0, k)
        rk_vanish_indices = findall(x -> x==0.0, rk)

        new(x, y, z, k, rk, ψx, ψk, ρx, ρk, Φx, Φk, fft_plan, rfft_plan, k_vanish_indices, rk_vanish_indices)
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
function Grids(length_tuple, resol_tuple::Tuple{Int, Int, Int})::Grids

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

    Grids(
        x, y, z,
        k_norm(length_tuple, resol_tuple),
        rk_norm(length_tuple, resol_tuple),
        ψx,
        ψk,
        ρx,
        ρk,
        Φx,
        Φk,
    )
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
