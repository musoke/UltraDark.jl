"""
The `Config` module contains structs and utility functions for controlling the simulation.

#Exports
$(EXPORTS)
"""
module Config

using DocStringExtensions

export SimulationConfig, constant_scale_factor
export TimeStepOptions

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct SimulationConfig
    "function defining scale factor a(t)"
    a::Any
    phase_grad_limit::Float64
    density_threshold::Float64
    time_step_options::Any
end

function SimulationConfig(;
    a = constant_scale_factor,
    phase_grad_limit = π / 4,
    density_threshold = 1e-6,
    time_step_options = TimeStepOptions(),
)
    SimulationConfig(a, phase_grad_limit, density_threshold, time_step_options)
end

"""
$(TYPEDSIGNATURES)

A function that always returns `1.0`, useful for simulations with no background expansion.
"""
function constant_scale_factor(t)
    1.0
end

"""
$(TYPEDEF)

struct containing options controlling the size and calculation of time steps.

See also: [`SimulationConfig`](@ref)

# Examples

```jldoctest
julia> using UltraDark

julia> TimeStepOptions()
TimeStepOptions(10, 1.0)
```

# Fields
$(TYPEDFIELDS)
"""
struct TimeStepOptions
    "controls how many steps are taken before the timestep is updated"
    update_period::Int64
    "multiplies the calculated maximum time step by a constant"
    multiplier::Float64
end

"""
    TimeStepOptions(; update_period=10, multiplier=1.0)
"""
function TimeStepOptions(; update_period = 10, multiplier = 1.0)
    TimeStepOptions(update_period, multiplier)
end

end # module
