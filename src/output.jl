module Output

using ..UltraDark
import Dates
import PencilFFTs, MPI
import NPZ
using Statistics

export OutputConfig
export SummaryStatistics, SummaryStatisticsMaxRms

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

    "Type of summary statistics to collect"
    summary_statistics::DataType
end

function OutputConfig(
    directory, output_times;
    box=true, slice=false, psi=true, rho=true,
    summary_statistics=SummaryStatistics,
)

    OutputConfig(directory, output_times, box, slice, psi, rho, summary_statistics)
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
            output = PencilFFTs.gather(grids.ψx)
            if MPI.Comm_rank(grids.MPI_COMM) == 0
                NPZ.npzwrite(
                    joinpath(output_config.directory, "psi_$step.npy"),
                    output,
                )
            end
        end

        if output_config.rho
            output = PencilFFTs.gather(grids.ρx)
            if MPI.Comm_rank(grids.MPI_COMM) == 0
                NPZ.npzwrite(
                    joinpath(output_config.directory, "rho_$step.npy"),
                    output,
                )
            end
        end
    end

    if output_config.slice
        if output_config.psi
            output = PencilFFTs.gather(grids.ψx[1, :, :])
            if MPI.Comm_rank(grids.MPI_COMM) == 0
                NPZ.npzwrite(
                    joinpath(output_config.directory, "psi_slice_$step.npy"),
                    output[1, :, :],
                )
            end
        end

        if output_config.rho
            output = PencilFFTs.gather(grids.ρx)
            if MPI.Comm_rank(grids.MPI_COMM) == 0
                NPZ.npzwrite(
                    joinpath(output_config.directory, "rho_slice_$step.npy"),
                    output[1, :, :],
                )
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
        write(file, column_titles(output_config.summary_statistics))
    end
end

function generate_summary_row(summary)::String
    line = ""
    for i in 1:nfields(summary)
        line = line * "$(getfield(summary, i)),"
    end
    line = line * "\n"
end

"""
    output_summary_row(grids, output_config, t, a, Δt)

Write a new row to the summary file
"""
function output_summary_row(grids, output_config, t, a, Δt)
    summary = output_config.summary_statistics(Dates.now(), t, a, Δt, grids)
    line = generate_summary_row(summary)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, line)
    end
end

"""
    output_summary_row(grids::PencilGrids, output_config, t, a, Δt)

Write a new row to the summary file

If the grids is a PencilGrids, this uses `MPI.Reduce` to compute partial
summaries in each task and combine them.
"""
function output_summary_row(grids::PencilGrids, output_config, t, a, Δt)
    root = 0
    summary = MPI.Reduce(
                   output_config.summary_statistics(Dates.now(), t, a, Δt, grids),
                   pool_summarystat,
                   root,
                   grids.MPI_COMM
                  )

    line = generate_summary_row(summary)
    if MPI.Comm_rank(grids.MPI_COMM) == 0
        open(joinpath(output_config.directory, "summary.csv"), "a") do file
            write(file, line)
        end
    end
end

"""
    SummaryStatistics

Summary statistics including the time, scale factor and current time step.
"""
struct SummaryStatistics
    "wall time"
    date::Dates.DateTime
    "time"
    t::Float64
    "scale factor"
    a::Float64
    "time step"
    Δt::Float64
end

function SummaryStatistics(wall_time, sim_time, a, Δt, grids)
    SummaryStatistics(wall_time, sim_time, a, Δt)
end

"""
    SummaryStatisticsMaxRms

Summary statistics including the maximum and RMS values of the density grid ρx.
"""
struct SummaryStatisticsMaxRms
    "wall time"
    date::Dates.DateTime
    "time"
    t::Float64
    "scale factor"
    a::Float64
    "time step"
    Δt::Float64
    "mean of density"
    ρx_mean::Float64
    "RMS of density contrast"
    δx_rms::Float64
    "number of grid points summarized"
    n::Int64
end

function SummaryStatisticsMaxRms(wall_time, sim_time, a, Δt, grids)
    ρx_mean = mean(grids.ρx)
    δx_rms = mean(((grids.ρx .- ρx_mean).^2))^0.5
    n = prod(size(grids.ρx))

    SummaryStatisticsMaxRms(wall_time, sim_time, a, Δt, ρx_mean, δx_rms, n)
end

"""
    column_titles

Generate column titles for a CSV file from fields names of a struct
"""
function column_titles(stat_struct)
    reduce(replace, [":"=>"", "("=>"", ")"=>"", " "=>""], init="$(fieldnames(stat_struct))") * ",\n"
end

"""
    pool_summarystat(S1::SummaryStatistics, S2::SummaryStatistics)

MPI reduction operator for summary statistics.

Check that t, a, Δt are equal and return them.
"""
function pool_summarystat(S1::SummaryStatistics, S2::SummaryStatistics)

    if (S1.t != S2.t) || (S1.a != S2.a) || (S1.Δt != S2.Δt)
        @error "Summaries incompatible across nodes" S1 S2
    end

    S1
end

"""
    pool_summarystat(S1::SummaryStatisticsMaxRms, S2::SummaryStatisticsMaxRms)

MPI reduction operator for summary statistics.
"""
function pool_summarystat(S1::SummaryStatisticsMaxRms, S2::SummaryStatisticsMaxRms)

    if (S1.t != S2.t) || (S1.a != S2.a) || (S1.Δt != S2.Δt)
        @error "Summaries incompatible across nodes" S1 S2
    end

    n = S1.n + S2.n
    ρx_mean = (S1.ρx_mean * S1.n + S2.ρx_mean * S2.n) / n
    δx_rms = ((S1.n * S1.δx_rms^2 + S2.n * S2.δx_rms^2) / n)^0.5

    SummaryStatisticsMaxRms(S1.date, S1.t, S1.a, S1.Δt, ρx_mean, δx_rms, n)
end

end # module
