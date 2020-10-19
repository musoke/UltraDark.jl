import NPZ

struct OutputConfig
    "where to write output"
    directory
    "times at which to output"
    output_times::Array{Real}

    "whether to output boxes"
    box::Bool
    "whether to output slices"
    slice::Bool
    "whether to output ψ"
    psi::Bool
    "whether to output ρ"
    rho::Bool
end

function OutputConfig(
    directory, output_times;
    box=true, slice=false, psi=true, rho=true
)

    OutputConfig(directory, output_times, box, slice, psi, rho)
end

"""
    output_grids(grids, output_config, step)

Write output from `grids` as specified in `output_config`
"""
function output_grids(grids, output_config, step)

    if output_config.box
        if output_config.psi
            NPZ.npzwrite(
                joinpath(output_config.directory, "psi_$step.npy"),
                grids.ψx
            )
        end
        if output_config.rho
            NPZ.npzwrite(
                joinpath(output_config.directory, "rho_$step.npy"),
                grids.ρx
            )
        end
    end

    if output_config.slice
        if output_config.psi
            NPZ.npzwrite(
                joinpath(output_config.directory, "psi_slice_$step.npy"),
                grids.ψx[1, :, :]
            )
        end
        if output_config.rho
            NPZ.npzwrite(
                joinpath(output_config.directory, "rho_slice_$step.npy"),
                grids.ρx[1, :, :]
            )
        end
    end
end

function output_grids(grids::PencilGrids, output_config, step)

    #TODO: don't use gather.  This sends all data to one node.  Should instead use multithreaded HDF5 output
    if output_config.box
        if output_config.psi
            NPZ.npzwrite(
                joinpath(output_config.directory, "psi_$step.npy"),
                gather(grids.ψx)
            )
        end
        if output_config.rho
            NPZ.npzwrite(
                joinpath(output_config.directory, "rho_$step.npy"),
                gather(grids.ρx)
            )
        end
    end

    if output_config.slice
        if output_config.psi
            NPZ.npzwrite(
                joinpath(output_config.directory, "psi_slice_$step.npy"),
                gather(grids.ψx[1, :, :])
            )
        end
        if output_config.rho
            NPZ.npzwrite(
                joinpath(output_config.directory, "rho_slice_$step.npy"),
                gather(grids.ρx[1, :, :])
            )
        end
    end
end

"""
    output_summary_header(output_config)

Write a header for a summary file

The header contains labels for each column of the summary CSV file.
This function overwrites the current contents of the file.
"""
function output_summary_header(output_config)

    open(joinpath(output_config.directory, "summary.csv"), "w") do file
        write(file, "t,a,Δt,ρx_mean,δx_rms\n")
    end
end

"""
    output_summary_row(grids, output_config, t, a, Δt)

Write a new row to the summary file
"""
function output_summary_row(grids, output_config, t, a, Δt)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        ρx_mean = mean(grids.ρx)
        δx_rms = mean(((grids.ρx .- ρx_mean).^2))^0.5
        write(file, "$t, $a, $Δt, $ρx_mean, $δx_rms\n")
    end
end
