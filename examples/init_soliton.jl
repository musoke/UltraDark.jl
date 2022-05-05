using UltraDark
using PencilFFTs
using NPZ

function add_soliton(grids, mass, position, velocity, phase, t0)
    add_soliton(grids, grids.ψx, mass, position, velocity, phase, t0)
end

function add_soliton(grids, psi, mass, position, velocity, phase, t0)

    profile = npzread(joinpath(@__DIR__, "initial_f.npy"))
    delta_x = 0.00001
    alpha = (mass / 3.883)^2
    beta = 2.454

    if typeof(grids) <: PencilGrids
        ψ_glob = global_view(psi)
    elseif typeof(grids) <: Grids
        ψ_glob = psi
    else
        throw("Unrecognised grids type")
    end

    for I in CartesianIndices(ψ_glob)

        i, j, k = Tuple(I)  # unpack indices

        x = grids.x[i]
        y = grids.y[j]
        z = grids.z[k]

        vx = velocity[1]
        vy = velocity[2]
        vz = velocity[3]

        distfromcentre = ((x - position[1])^2 + (y - position[2])^2 + (z - position[3])^2)^0.5

        if alpha^0.5 * distfromcentre <= 5.6
            f = alpha * profile[trunc(Int, alpha^0.5 * (distfromcentre / delta_x + 1))]
            f *= exp(im * (alpha * beta * t0 + (x * vx + y * vy + z * vz) - 0.5 * (vx^2 + vy^2 + vz^2) * t0 ))
            f *= exp(im * phase)
            ψ_glob[I] += f
        else
            ψ_glob[I] += 0
        end
    end
end
