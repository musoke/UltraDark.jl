module Config

export SimulationConfig, constant_scale_factor
export TimeStepOptions

struct SimulationConfig
    "function defining a(t)"
    a
    phase_grad_limit::Float64
    density_threshold::Float64
    time_step_options
end

function constant_scale_factor(t)
    1.0
end

function SimulationConfig(
                          ;
                          a=constant_scale_factor,
                          phase_grad_limit=π/4,
                          density_threshold=1e-6,
                          time_step_options=TimeStepOptions(),
                         )
    SimulationConfig(constant_scale_factor, phase_grad_limit, density_threshold, time_step_options)
end

"""
    TimeStepOptions(update_period=10, multiplier=1.0)

struct containing options controlling the size and calculation of time steps.

update_period::Int64 controls how many steps are taken before the timestep is
updated.

multiplier::Float64 multiplies the calculated maximum time step by a constant

See also: [`SimulationConfig`](@

# Examples
```jldoctest
julia> using UltraDark

julia> TimeStepOptions()
TimeStepOptions(10, 1.0)

```
"""
struct TimeStepOptions
    update_period::Int64
    multiplier::Float64
end

function TimeStepOptions(
                         ;
                         update_period=10,
                         multiplier=1.0,
                        )
    TimeStepOptions(update_period, multiplier)
end

end # module
