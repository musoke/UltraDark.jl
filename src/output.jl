import NPZ

struct OutputConfig
    directory
end

function output_grids(grids, output_config, step)
    NPZ.npzwrite(joinpath(output_config.directory, "psi_slice_$step.npy"), grids.ψx[1, :, :])
end

function output_summary(grids, output_config, t, a)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, "$t, $a, $(mean(grids.ρx)), \n")
    end
end
