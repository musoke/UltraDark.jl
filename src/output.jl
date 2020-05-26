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

function output_summary(grids, output_config, t, a, Δt)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, "$t, $a, $Δt, $(mean(grids.ρx)), \n")
    end
end
