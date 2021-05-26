import NPZ

"""

Parallel HDF5 output is possible when using `PencilGrids`.
See the docs for `PencilArrays` for instructions to enable it.
https://jipolanco.github.io/PencilArrays.jl/v0.9/PencilIO/#setting_up_parallel_hdf5
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

    "whether to output HDF5; only has an effect with PencilGrids"
    hdf5::Bool
end

function OutputConfig(
    directory, output_times;
    box=true, slice=false, psi=true, rho=true,
    hdf5=false,
)

    OutputConfig(directory, output_times, box, slice, psi, rho, hdf5)
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

    if output_config.box
        if output_config.psi
            if output_config.hdf5
                open(PencilArrays.PencilIO.PHDF5Driver(), "psi_$step.hdf5", grids.MPI_COMM_WORLD; write=true) do ff
                    ff["psi_plane"] = ψx
                end
            else
                output = gather(grids.ψx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "psi_$step.npy"),
                        output,
                    )
                end
            end
        end

        if output_config.rho
            if output_config.hdf5
                open(PencilArrays.PencilIO.PHDF5Driver(), "rho_$step.hdf5", grids.MPI_COMM_WORLD; write=true) do ff
                    ff["psi_plane"] = ρx
                end
            else
                output = gather(grids.ρx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "rho_$step.npy"),
                        output,
                    )
                end
            end
        end
    end

    if output_config.slice
        if output_config.psi
            if output_config.hdf5
                open(PencilArrays.PencilIO.PHDF5Driver(), "psi_slice_$step.hdf5", grids.MPI_COMM_WORLD; write=true) do ff
                    ff["psi_plane"] = ψx[1, :, :]
                end
            else
                output = gather(grids.ψx[1, :, :])
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "psi_slice_$step.npy"),
                        output[1, :, :],
                    )
                end
            end
        end

        if output_config.rho
            if output_config.hdf5
                open(PencilArrays.PencilIO.PHDF5Driver(), "rho_slice_$step.hdf5", grids.MPI_COMM_WORLD; write=true) do ff
                    ff["rho_plane"] = ρx[1, :, :]
                end
            else
                output = gather(grids.ρx)
                if MPI.Comm_rank(grids.MPI_COMM) == 0
                    NPZ.npzwrite(
                        joinpath(output_config.directory, "rho_slice_$step.npy"),
                        output[1, :, :],
                    )
                end
            end
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
    s = SummaryStat(grids)

    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, "$t, $a, $Δt, $(s.ρx_mean), $(s.δx_rms)\n")
    end
end

"""
    output_summary_row(grids::PencilGrids, output_config, t, a, Δt)

Write a new row to the summary file
"""
function output_summary_row(grids::PencilGrids, output_config, t, a, Δt)
    root = 0
    s = MPI.Reduce(SummaryStat(grids), pool_summarystat, root, grids.MPI_COMM)

    if MPI.Comm_rank(grids.MPI_COMM) == 0
        open(joinpath(output_config.directory, "summary.csv"), "a") do file
            write(file, "$t, $a, $Δt, $(s.ρx_mean), $(s.δx_rms)\n")
        end
    end
end

"""
    SummaryStat
"""
struct SummaryStat
    "mean of density"
    ρx_mean::Float64
    "RMS of density contrast"
    δx_rms::Float64
    "number of grid points summarized"
    n::Float64
end

function SummaryStat(grids)
    ρx_mean = mean(grids.ρx)
    δx_rms = mean(((grids.ρx .- ρx_mean).^2))^0.5
    n = prod(size(grids.ρx))

    SummaryStat(ρx_mean, δx_rms, n)
end

"""
    pool_summarystat(S1::SummaryStat, S2::SummaryStat)

Custom MPI reduction operator for summary statistics.
"""
function pool_summarystat(S1::SummaryStat, S2::SummaryStat)
    n = S1.n + S2.n
    ρx_mean = (S1.ρx_mean * S1.n + S2.ρx_mean * S2.n) / n
    δx_rms = ((S1.n * S1.δx_rms^2 + S2.n * S2.δx_rms^2) / n)^0.5

    SummaryStat(ρx_mean, δx_rms, n)
end
