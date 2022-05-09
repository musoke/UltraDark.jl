module Summary

using ..UltraDark
using ..UltraDark: AbstractGrids
import Dates
import MPI
using Statistics

export WallTime
export SimulationTime
export ScaleFactor
export TimeStep
export MeanDensity, MaxDensity, RmsDensityContrast

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

"""
    output_summary_row(grids, output_config, t, a, Δt)

Write a new row to the summary file
"""
function output_summary_row(grids, output_config, t, a, Δt, constants)
    summaries = map_summary_statistics(output_config.summary_statistics, t, a, Δt, grids, constants)
    line = generate_summary_row(summaries)
    open(joinpath(output_config.directory, "summary.csv"), "a") do file
        write(file, line)
    end
end

"""
    map_summary_statistics(summary_statistics, sim_time, a, Δt, grids, constants)

Calculate each summary statistic in summary_statistics

If the grids is a PencilGrids, this uses `MPI.Reduce` to compute partial
summaries in each task and combine them.
"""

function map_summary_statistics(summary_statistics, sim_time, a, Δt, grids, constants)
    summaries = map(x->x(sim_time, a, Δt, grids, constants), summary_statistics)
end

function map_summary_statistics(summary_statistics, sim_time, a, Δt, grids::PencilGrids, constants)
    root = 0

    local_summaries = map(x -> x(sim_time, a, Δt, grids, constants), summary_statistics)

    global_summaries = map(
                           x -> MPI.Reduce(x, pool_summarystat, root, grids.MPI_COMM,),
                           local_summaries
                          )

end

"""
    column_titles(summaries)

Generate column titles from an iterable of summaries
"""
function column_titles(summaries)
    mapreduce(column_title, (s1, s2)->s1 * "," * s2, summaries) * "\n"
end

"""
    column_title(summary_struct)

Get column title

Assumes that the desired column title is the name of the first field.  This can
be refined for other structs by defining a specialization for the given
datatype.
"""
function column_title(summary_struct)
    "$(fieldname(summary_struct, 1))"
end

"""
    generate_summary_row(summaries)
"""
function generate_summary_row(summaries)::String
    mapreduce(data, (s1, s2)-> "$s1,$s2", summaries) * "\n"
end


"""
    data(summary_struct)

Get column entry

Assumes that the desired data is in the first field of the struct.  This can
be refined for other structs by defining a specialization for the given
datatype.
"""
function data(summary_struct)
    getfield(summary_struct, 1)
end

"""
    WallTime

The current time in the real world.


# Examples

`WallTime` created after another contains a later time.

```jldoctest
julia> using UltraDark

julia> t1 = Summary.WallTime();

julia> t2 = Summary.WallTime(0., 1., 1e-1, Grids(1.0, 16), nothing);

julia> t1.date < t2.date
true
```
"""
struct WallTime
    "wall time"
    date::Dates.DateTime
end

function WallTime()
    WallTime(Dates.now())
end

function WallTime(sim_time, a, Δt, grids, constants)
    WallTime()
end

"""
    pool_summarystat(S1::WallTime, S2::WallTime)

MPI reduction operator for wall time

Return the time of the first argument.

"""
function pool_summarystat(S1::WallTime, S2::WallTime)::WallTime

    S1

end

struct SimulationTime
    "time"
    t::Float64
end

function SimulationTime(sim_time, a, Δt, grids, constants)
    SimulationTime(sim_time)
end

"""
    pool_summarystat(S1::SimulationTime, S2::SimulationTime)

MPI reduction operator for simulation time

Check that times are equal and return
"""
function pool_summarystat(S1::SimulationTime, S2::SimulationTime)::SimulationTime

    if (S1.t != S2.t)
        @error "Summaries incompatible across nodes" S1 S2
    end

    S1
end

struct ScaleFactor
    "scale factor"
    a::Float64
end

function ScaleFactor(sim_time, a, Δt, grids, constants)
    ScaleFactor(a)
end

"""
    pool_summarystat(S1::ScaleFactor, S2::ScaleFactor)

MPI reduction operator for scale factor

Check that scale factors are equal and return
"""
function pool_summarystat(S1::ScaleFactor, S2::ScaleFactor)::ScaleFactor

    if (S1.a != S2.a)
        @error "Summaries incompatible across nodes" S1 S2
    end

    S1
end

struct TimeStep
    "time step"
    Δt::Float64
end

function TimeStep(sim_time, a, Δt, grids, constants)
    TimeStep(Δt)
end

"""
    pool_summarystat(S1::TimeStep, S2::TimeStep)

MPI reduction operator for scale factor

Check that scale factors are equal and return
"""
function pool_summarystat(S1::TimeStep, S2::TimeStep)::TimeStep

    if (S1.Δt != S2.Δt)
        @error "Summaries incompatible across nodes" S1 S2
    end

    S1
end

"""
    MeanDensity

Summary statistic containing the mean density and the number of cells over
which it was calculated.

# Examples

```jldoctest
julia> using UltraDark

julia> g = Grids(1.0, 16);

julia> Summary.MeanDensity(g)
UltraDark.Summary.MeanDensity(0.0, 4096)

```
"""
struct MeanDensity
    "mean of density"
    ρx_mean::Float64
    "number of grid points summarized"
    n::Int64
end

function MeanDensity(grids)
    ρx_mean = mean(grids.ρx)
    n = prod(size(grids.ρx))

    MeanDensity(ρx_mean, n)
end

function MeanDensity(sim_time, a, Δt, grids, constants)
    ρx_mean = mean(grids.ρx)
    n = prod(size(grids.ρx))

    MeanDensity(ρx_mean, n)
end

"""
    pool_summarystat(S1::MeanDensity, S2::MeanDensity)

MPI reduction operator for mean density
"""
function pool_summarystat(S1::MeanDensity, S2::MeanDensity)::MeanDensity

    n = S1.n + S2.n
    ρx_mean = (S1.ρx_mean * S1.n + S2.ρx_mean * S2.n) / n
    ρx_max = maximum(S1.ρx_max, S2.ρx_max)

    MeanDensity(ρx_mean, n)
end

"""
    MaxDensity
"""
struct MaxDensity
    "max of density"
    ρx_max::Float64
end

function MaxDensity(sim_time, a, Δt, grids, constants)
    ρx_max = maximum(grids.ρx)

    MaxDensity(ρx_max)
end

"""
    pool_summarystat(S1::MaxDensity, S2::MaxDensity)

MPI reduction operator for max density
"""
function pool_summarystat(S1::MaxDensity, S2::MaxDensity)::MaxDensity

    n = S1.n + S2.n
    ρx_max = maximum(S1.ρx_max, S2.ρx_max)

    MaxDensity(ρx_max, n)
end

"""
    RmsDensityContrast
"""
struct RmsDensityContrast
    "RMS of density contrast"
    δx_rms::Float64
    "number of grid points summarized"
    n::Int64
end

function RmsDensityContrast(sim_time, a, Δt, grids, constants)
    ρx_mean = mean(grids.ρx)
    δx_rms = mean(((grids.ρx .- ρx_mean).^2))^0.5
    n = prod(size(grids.ρx))

    RmsDensityContrast(δx_rms, n)
end

"""
    pool_summarystat(S1::RmsDensityContrast, S2::RmsDensityContrast)

MPI reduction operator for summary statistics.
"""
function pool_summarystat(S1::RmsDensityContrast, S2::RmsDensityContrast)::RmsDensityContrast

    n = S1.n + S2.n
    δx_rms = ((S1.n * S1.δx_rms^2 + S2.n * S2.δx_rms^2) / n)^0.5

    RmsDensityContrast(δx_rms, n)
end

"""
    TotalMass

Total mass on a grid

# Examples

```jldoctest
julia> using UltraDark

julia> g = Grids(1.0, 16);

julia> g.ρx .= 0.;

julia> g.ρx[1, 1, 1] = 1.;

julia> getfield(Summary.TotalMass(g), 1) == 1.0 * (1.0/16)^3
true

```
"""
struct TotalMass
    mass::Float64
end

function TotalMass(grids, rho)
    mass = sum(rho * dV(grids))

    TotalMass(mass)
end

function TotalMass(grids::AbstractGrids)
    TotalMass(grids, grids.ρx)
end

function TotalMass(sim_time, a, Δt, grids, constants)
    TotalMass(grids)
end

function pool_summarystat(S1::TotalMass, S2::TotalMass)::TotalMass
    mass = S1.mass + S2.mass
    TotalMass(mass)
end

end # module
