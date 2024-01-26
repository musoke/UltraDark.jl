module Initialise

using PencilFFTs
using NPZ

using ..UltraDark: AbstractGrids

export add_fdm_soliton!

"""
    add_fdm_soliton!(grids::AbstractGrids, mass, position, velocity, phase, t0)
    add_fdm_soliton!(grids::AbstractGrids, psi, mass, position, velocity, phase, t0)

Add a fuzzy dark matter soliton `grids`.  The
density profile is rescaled to the desired mass and the phase is set to the desired
velocity.

If the argument `psi` is passed, the soliton is added to the array-like `psi` rather than
than `grids.ψx`.

Note that due to coarse grid effects, the mass of the added soliton may not match the input
mass exactly.

The included density profile from comes [from PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight/blob/28d942cb3100b6df9976a6e6f392eaeeb2ce9689/Soliton%20Profile%20Files/initial_f.npy).

# Examples

The following  script adds a soliton to a grid and checks that the resulting total mass is
correct.

```jldoctest
using UltraDark
using UltraDark.Initialise

resol = 64
grids = Grids(10.0, resol)

mass = 10.0
position = [0.0, 0.0, 0.0]
velocity = [0.0, 0.0, 0.0]
phase = 0.0
t0 = 0.0

add_fdm_soliton!(grids, mass, position, velocity, phase, t0)

# Ensure density is correct
UltraDark.update_gravitational_potential!(grids, ())

actual_mass = UltraDark.mass(grids)

isapprox(mass, actual_mass, rtol = 1e-3)

# output

true
```
"""
function add_fdm_soliton!(grids::AbstractGrids, mass, position, velocity, phase, t0)
    add_fdm_soliton!(grids, grids.ψx, mass, position, velocity, phase, t0)
end

function add_fdm_soliton!(grids::AbstractGrids, psi, mass, position, velocity, phase, t0)

    profile::Vector{Float64} =
        npzread(joinpath(@__DIR__, "..", "examples", "initial_f.npy"))
    delta_x = 0.00001
    alpha = (mass / 3.883)^2
    beta = 2.454

    if typeof(psi) <: PencilArray
        ψ_glob = global_view(psi)
    elseif typeof(psi) <: Array
        ψ_glob = psi
    else
        throw("Unrecognised array type")
    end

    for I in CartesianIndices(ψ_glob)

        i, j, k = Tuple(I)  # unpack indices

        x = grids.x[i]
        y = grids.y[j]
        z = grids.z[k]

        vx = velocity[1]
        vy = velocity[2]
        vz = velocity[3]

        distfromcentre =
            ((x - position[1])^2 + (y - position[2])^2 + (z - position[3])^2)^0.5

        if alpha^0.5 * distfromcentre <= 5.6
            f = alpha * profile[trunc(Int, alpha^0.5 * (distfromcentre / delta_x + 1))]
            f *= exp(
                im * (
                    alpha * beta * t0 + (x * vx + y * vy + z * vz) -
                    0.5 * (vx^2 + vy^2 + vz^2) * t0
                ),
            )
            f *= exp(im * phase)
            ψ_glob[I] += f
        else
            ψ_glob[I] += 0
        end
    end
end

end # module
