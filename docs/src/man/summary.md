# Summary statistics

It is often useful to collect summary statistics as a simulation runs and avoid storing complete simulation boxes.
The [`UltraDark.Summary`](@ref) module makes this easy.

## Example

Say you want to extract the minimum density at each time step.
Define a struct to contain this data:
```@example continued=true
struct MinDensity
    ρx_min::Float64
end
```
and a function to construct it
```@example continued=true
function MaxDensity(sim_time, a, Δt, grids, constants, external_states)
    ρx_min = minimum(grids.ρx)

    MinDensity(ρx_max)
end
```

Passing this object into a simulation will add a column to `$output_path/summary.csv` containing the single value contained by `MinDensity`.

More complex statistics are also possible.

## Docstrings

```@autodocs
Modules = [UltraDark.Summary]
```
