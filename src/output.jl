import NPZ

struct OutputConfig
    "where to write output"
    directory
    "times at which to output"
    output_times::Array{Real}
end

function output_grids(grids, output_config, step)
    NPZ.npzwrite(joinpath(output_config.directory, "psi_slice_$step.npy"), grids.ψx[1, :, :])
end

function output_summary(grids, output_config, t, a)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, "$t, $a, $(mean(grids.ρx)), \n")
    end
end
