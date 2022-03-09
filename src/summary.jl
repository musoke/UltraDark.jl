module Summary

using ..UltraDark
import Dates
import MPI
using Statistics

export SummaryStatistics, SummaryStatisticsMeanMaxRms

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
    line = "$(getfield(summary, 1))"
    for i in 2:nfields(summary)
        line = line * ",$(getfield(summary, i))"
    end
    line = line * "\n"
end

"""
    output_summary_row(grids, output_config, t, a, Δt)

Write a new row to the summary file
"""
function output_summary_row(grids, output_config, t, a, Δt, constants)
    summary = output_config.summary_statistics(Dates.now(), t, a, Δt, grids, constants)
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
function output_summary_row(grids::PencilGrids, output_config, t, a, Δt, constants)
    root = 0
    summary = MPI.Reduce(
                   output_config.summary_statistics(Dates.now(), t, a, Δt, grids, constants),
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

function SummaryStatistics(wall_time, sim_time, a, Δt, grids, constants)
    SummaryStatistics(wall_time, sim_time, a, Δt)
end

"""
    SummaryStatisticsMeanMaxRms

Summary statistics including the mean density and RMS density contrast.
"""
struct SummaryStatisticsMeanMaxRms
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
    "max of density"
    ρx_max::Float64
    "RMS of density contrast"
    δx_rms::Float64
    "number of grid points summarized"
    n::Int64
end

function SummaryStatisticsMeanMaxRms(wall_time, sim_time, a, Δt, grids, constants)
    ρx_mean = mean(grids.ρx)
    ρx_max = maximum(grids.ρx)
    δx_rms = mean(((grids.ρx .- ρx_mean).^2))^0.5
    n = prod(size(grids.ρx))

    SummaryStatisticsMeanMaxRms(wall_time, sim_time, a, Δt, ρx_mean, ρx_max, δx_rms, n)
end

"""
    column_titles

Generate column titles for a CSV file from fields names of a struct
"""
function column_titles(stat_struct)
    reduce(replace, [":"=>"", "("=>"", ")"=>"", " "=>""], init="$(fieldnames(stat_struct))") * "\n"
end

"""
    pool_summarystat(S1::SummaryStatistics, S2::SummaryStatistics)

MPI reduction operator for summary statistics.

Check that t, a, Δt are equal and return them.
"""
function pool_summarystat(S1::SummaryStatistics, S2::SummaryStatistics)::SummaryStatistics

    if (S1.t != S2.t) || (S1.a != S2.a) || (S1.Δt != S2.Δt)
        @error "Summaries incompatible across nodes" S1 S2
    end

    S1
end

"""
    pool_summarystat(S1::SummaryStatisticsMeanMaxRms, S2::SummaryStatisticsMeanMaxRms)

MPI reduction operator for summary statistics.
"""
function pool_summarystat(S1::SummaryStatisticsMeanMaxRms, S2::SummaryStatisticsMeanMaxRms)::SummaryStatisticsMeanMaxRms

    if (S1.t != S2.t) || (S1.a != S2.a) || (S1.Δt != S2.Δt)
        @error "Summaries incompatible across nodes" S1 S2
    end

    n = S1.n + S2.n
    ρx_mean = (S1.ρx_mean * S1.n + S2.ρx_mean * S2.n) / n
    ρx_max = maximum(S1.ρx_max, S2.ρx_max)
    δx_rms = ((S1.n * S1.δx_rms^2 + S2.n * S2.δx_rms^2) / n)^0.5

    SummaryStatisticsMeanMaxRms(S1.date, S1.t, S1.a, S1.Δt, ρx_mean, ρx_max, δx_rms, n)
end

end # module
