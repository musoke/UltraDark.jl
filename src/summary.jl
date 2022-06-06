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
function output_summary_row(grids, output_config, t, a, Δt, constants, external_states)
    summaries = map_summary_statistics(
        output_config.summary_statistics,
        t,
        a,
        Δt,
        grids,
        constants,
        external_states,
    )
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

function map_summary_statistics(
    summary_statistics,
    sim_time,
    a,
    Δt,
    grids,
    constants,
    external_states,
)
    summaries =
        map(x -> x(sim_time, a, Δt, grids, constants, external_states), summary_statistics)
end

function map_summary_statistics(
    summary_statistics,
    sim_time,
    a,
    Δt,
    grids::PencilGrids,
    constants,
    external_states,
)
    root = 0

    local_summaries =
        map(x -> x(sim_time, a, Δt, grids, constants, external_states), summary_statistics)

    global_summaries =
        map(x -> MPI.Reduce(x, pool_summarystat, root, grids.MPI_COMM), local_summaries)

end

"""
    column_titles(summaries)

Generate column titles from an iterable of summaries
"""
function column_titles(summary_types)::String
    mapreduce(t -> Val(t) |> column_title, (s1, s2) -> s1 * "," * s2, summary_types) * "\n"
end

"""
    column_title(summary_struct)

Get column title

Assumes that the desired column title is the name of the first field.  This can
be refined for other structs by defining a specialization for the given
datatype.

Note that the input argument is of type `Val{T}`, so that one can dispatch on
`T`.
"""
function column_title(v::Val{T})::String where {T}
    "$(fieldname(T, 1))"
end

"""
    generate_summary_row(summaries)
"""
function generate_summary_row(summaries)::String
    mapreduce(get_relevant_data, (s1, s2) -> "$s1,$s2", summaries) * "\n"
end


"""
    get_relevant_data(summary_struct)

Format column entry as string

By default, assumes that the desired data is in the first field of the struct.
This can be refined for other structs by defining a specialization for the
given datatype and returning a comma separated string.
"""
function get_relevant_data(summary_data)::String
    "$(getfield(summary_data, 1))"
end

"""
    WallTime

The current time in the real world.

# Examples

`WallTime` created after another contains a later time.

```jldoctest
julia> using UltraDark


julia> t1 = Summary.WallTime();


julia> t2 = Summary.WallTime(0.0, 1.0, 1e-1, Grids(1.0, 16), nothing, ());


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

function WallTime(sim_time, a, Δt, grids, constants, external_states)
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

function SimulationTime(sim_time, a, Δt, grids, constants, external_states)
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

function ScaleFactor(sim_time, a, Δt, grids, constants, external_states)
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

function TimeStep(sim_time, a, Δt, grids, constants, external_states)
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

function MeanDensity(sim_time, a, Δt, grids, constants, external_states)
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

function MaxDensity(sim_time, a, Δt, grids, constants, external_states)
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
    MaxDensityIndex

This struct contains 4 useful pieces of information: the maxiumum value and
three indices.  Each is output to a summary.
"""
struct MaxDensityIndex
    "max of density"
    ρx_max::Float64
    "x index of max density"
    max_index_x::Int32
    "y index of max density"
    max_index_y::Int32
    "z index of max density"
    max_index_z::Int32
end

function MaxDensityIndex(sim_time, a, Δt, grids, constants, external_states)
    I = argmax(grids.ρx)

    ρx_max = grids.ρx[I]

    MaxDensityIndex(ρx_max, Tuple(I)...)
end

"""
    pool_summarystat(S1::MaxDensity, S2::MaxDensity)

MPI reduction operator for max density index
"""
function pool_summarystat(S1::MaxDensityIndex, S2::MaxDensityIndex)::MaxDensityIndex
    if S1.ρx_max > S2.ρx_max
        S1
    else
        S2
    end
end

function column_title(v::Val{T})::String where {T<:MaxDensityIndex}
    mapreduce(x -> "$x", (s1, s2) -> "$s1,$s2", fieldnames(T))
end

function get_relevant_data(summary_data::MaxDensityIndex)::String
    mapreduce(i -> getfield(summary_data, i), (s1, s2) -> "$s1,$s2", 1:4)
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

function RmsDensityContrast(sim_time, a, Δt, grids, constants, external_states)
    ρx_mean = mean(grids.ρx)
    δx_rms = mean(((grids.ρx .- ρx_mean) .^ 2))^0.5
    n = prod(size(grids.ρx))

    RmsDensityContrast(δx_rms, n)
end

"""
    pool_summarystat(S1::RmsDensityContrast, S2::RmsDensityContrast)

MPI reduction operator for summary statistics.
"""
function pool_summarystat(
    S1::RmsDensityContrast,
    S2::RmsDensityContrast,
)::RmsDensityContrast

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


julia> Summary.TotalMass(0.0, 1.0, 1e-1, g, nothing, ())
UltraDark.Summary.TotalMass(0.0)
```
"""
struct TotalMass
    mass::Float64
end

function TotalMass(sim_time, a, Δt, grids, constants, external_states)::TotalMass
    TotalMass(UltraDark.mass(grids))
end

function pool_summarystat(S1::TotalMass, S2::TotalMass)::TotalMass
    TotalMass(S1.mass + S2.mass)
end

"""
    EnergyGravity

Gravitational potential energy
"""
struct EnergyGravity
    E_grav::Float64
end

function EnergyGravity(sim_time, a, Δt, grids, constants, external_states)::EnergyGravity
    EnergyGravity(UltraDark.E_grav(grids))
end

function pool_summarystat(S1::EnergyGravity, S2::EnergyGravity)::EnergyGravity
    EnergyGravity(S1.E_grav + S2.E_grav)
end

"""
    EnergyKineticQuantum

Gravitational potential energy
"""
struct EnergyKineticQuantum
    E_kq::Float64
end

function EnergyKineticQuantum(
    sim_time,
    a,
    Δt,
    grids,
    constants,
    external_states,
)::EnergyKineticQuantum
    EnergyKineticQuantum(UltraDark.E_kq(grids))
end

function pool_summarystat(
    S1::EnergyKineticQuantum,
    S2::EnergyKineticQuantum,
)::EnergyKineticQuantum
    EnergyKineticQuantum(S1.E_kq + S2.E_kq)
end

"""
    AngularMomentum

Total angular momentum in each of 3 directions

# Examples

A stationary field has no angular momentum

```jldoctest
julia> using UltraDark

julia> g = Grids(1.0, 16);

julia> g.ψx .= 1;

julia> Summary.AngularMomentum(0., 1., 1e-1, g, nothing)
UltraDark.Summary.AngularMomentum(0.0, 0.0, 0.0)

"""

struct AngularMomentum
    Lx::Float64
    Ly::Float64
    Lz::Float64
end

function AngularMomentum(sim_time, a, Δt, grids, constants, external_states)
    Lx, Ly, Lz = UltraDark.angular_momentum(grids)

    AngularMomentum(Lx, Ly, Lz)
end

function column_title(v::Val{T})::String where {T<:AngularMomentum}
    mapreduce(x -> "$x", (s1, s2) -> "$s1,$s2", fieldnames(T))
end

function get_relevant_data(summary_data::AngularMomentum)::String
    mapreduce(i -> getfield(summary_data, i), (s1, s2) -> "$s1,$s2", 1:3)
end

function pool_summarystat(S1::AngularMomentum, S2::AngularMomentum)::AngularMomentum
    AngularMomentum(S1.Lx + S2.Lx, S1.Ly + S2.Ly, S1.Lz + S2.Lz)
end

end # module
