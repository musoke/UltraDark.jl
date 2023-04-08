module Output

using ..UltraDark
import HDF5
import NPZ
import PencilFFTs, MPI
using Statistics

import ..Summary: output_summary_header, output_summary_row

export Summary
export OutputConfig

"""
    OutputConfig

struct containing information about what to output.

summary_statistics should be another struct whose constructor generates summary
statistics from `t`, `a`, `Δt` and a `grids` object.
If it is to be used with a `PencilGrids` object, each field must be concrete
and have a binary representation that MPI can handle.
"""
struct OutputConfig
    "where to write output"
    directory::String
    "times at which to output"
    output_times::Array{Float64}

    "whether to output boxes"
    box::Bool
    "whether to output slices"
    slice::Bool
    "whether to output ψ"
    psi::Bool
    "whether to output ρ"
    rho::Bool

    "write .npy files"
    npy::Bool
    "write HDF5 files"
    h5::Bool

    "Type of summary statistics to collect"
    summary_statistics::Tuple
end

function OutputConfig(
    directory,
    output_times;
    box = true,
    slice = false,
    psi = true,
    rho = true,
    npy = true,
    h5 = false,
    summary_statistics = (Summary.WallTime, Summary.SimulationTime),
)

    OutputConfig(directory, output_times, box, slice, psi, rho, npy, h5, summary_statistics)
end

"""
    output_state(grids, external_states, output_config, step)

Write out the grids and possible external states
"""
function output_state(grids, external_states, output_config, step)
    output_grids(grids, output_config, step)

    for (index, s) in enumerate(external_states)
        output_external_state(s, output_config, step, index)
    end
end

"""
    output_grids(grids, output_config, step)

Write output from `grids` as specified in `output_config`
"""
function output_grids(grids, output_config, step)

    if output_config.box
        if output_config.psi
            if output_config.npy
                NPZ.npzwrite(joinpath(output_config.directory, "psi_$step.npy"), grids.ψx)
            end
            if output_config.h5
                HDF5.h5open(joinpath(output_config.directory, "psi_$step.h5"), "w") do file
                    write(file, "psi", grids.ψx)
                end
            end
        end
        if output_config.rho
            if output_config.npy
                NPZ.npzwrite(joinpath(output_config.directory, "rho_$step.npy"), grids.ρx)
            end
            if output_config.h5
                HDF5.h5open(joinpath(output_config.directory, "rho_$step.h5"), "w") do file
                    write(file, "rho", grids.ρx)
                end
            end
        end
    end

    if output_config.slice
        if output_config.psi
            if output_config.npy
                NPZ.npzwrite(
                    joinpath(output_config.directory, "psi_slice_$step.npy"),
                    grids.ψx[1, :, :],
                )
            end
            if output_config.h5
                HDF5.h5open(
                    joinpath(output_config.directory, "rho_slice_$step.h5"),
                    "w",
                ) do file
                    write(file, "psi", grids.ψx[1, :, :])
                end
            end
        end
        if output_config.rho
            if output_config.npy
                NPZ.npzwrite(
                    joinpath(output_config.directory, "rho_slice_$step.npy"),
                    grids.ρx[1, :, :],
                )
            end
            if output_config.h5
                HDF5.h5open(
                    joinpath(output_config.directory, "rho_slice_$step.h5"),
                    "w",
                ) do file
                    write(file, "rho", grids.ρx[1, :, :])
                end
            end
        end
    end
end

"""
    output_external_state(external_state, output_config, step, index)

Output states other than the ψ field.

By default this does nothing, but can be overloaded.

# Arguments

index::Integer index of state in external states
"""
function output_external_state(external_state, output_config, step, index) end

function output_grids(grids::PencilGrids, output_config, step)

    #TODO: don't use gather.  This sends all data to one node.  Should instead use multithreaded HDF5 output
    if output_config.box
        if output_config.psi
            if output_config.npy
                output = PencilFFTs.gather(grids.ψx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(joinpath(output_config.directory, "psi_$step.npy"), output)
                end
            end
            if output_config.h5
                throw(Core.ArgumentError("h5 output unimplemented for PencilGrids"))
            end
        end

        if output_config.rho
            if output_config.npy
                output = PencilFFTs.gather(grids.ρx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(joinpath(output_config.directory, "rho_$step.npy"), output)
                end
            end
            if output_config.h5
                throw(Core.ArgumentError("h5 output unimplemented for PencilGrids"))
            end
        end
    end

    if output_config.slice
        if output_config.psi
            if output_config.npy
                output = PencilFFTs.gather(grids.ψx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "psi_slice_$step.npy"),
                        output[1, :, :],
                    )
                end
            end
            if output_config.h5
                throw(Core.ArgumentError("h5 output unimplemented for PencilGrids"))
            end
        end

        if output_config.rho
            if output_config.npy
                output = PencilFFTs.gather(grids.ρx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "rho_slice_$step.npy"),
                        output[1, :, :],
                    )
                end
            end
            if output_config.h5
                throw(Core.ArgumentError("h5 output unimplemented for PencilGrids"))
            end
        end
    end
end

"""
    output_xyz(grids, output_config)

Output the spatial coordinates defining the grid
"""
function output_xyz(grids, output_config)
    if output_config.npy
        NPZ.npzwrite(joinpath(output_config.directory, "x.npy"), grids.x)
        NPZ.npzwrite(joinpath(output_config.directory, "y.npy"), grids.y)
        NPZ.npzwrite(joinpath(output_config.directory, "z.npy"), grids.z)
    end
    if output_config.h5
        HDF5.h5open(joinpath(output_config.directory, "x.h5"), "w") do file
            write(file, "x", grids.x)
        end
        HDF5.h5open(joinpath(output_config.directory, "y.h5"), "w") do file
            write(file, "y", grids.y)
        end
        HDF5.h5open(joinpath(output_config.directory, "z.h5"), "w") do file
            write(file, "z", grids.z)
        end
    end
end

"""
    output_output_times(output_times, output_config)

Output the times corresponding to slices
"""
function output_output_times(output_times, output_config)
    if output_config.npy
        NPZ.npzwrite(joinpath(output_config.directory, "output_times.npy"), output_times)
    end
    if output_config.h5
        HDF5.h5open(joinpath(output_config.directory, "output_times.h5"), "w") do file
            write(file, "output_times", output_times)
        end
    end
end

"""
    output_external_states_headers(external_states, output_config)
"""
function output_external_states_headers(external_states, output_config)
    for (index, s) in enumerate(external_states)
        output_external_state_header(s, output_config, index)
    end
end

"""
    output_external_state_header(state, output_config)
"""
function output_external_state_header(state, output_config, index) end

end # module
