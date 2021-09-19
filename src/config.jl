module Config

export SimulationConfig, constant_scale_factor

struct SimulationConfig
    time_step_update_period::Int64
    time_step_multiplier::Float64
    "function defining a(t)"
    a
    phase_grad_limit::Float64
    density_threshold::Float64
end

function constant_scale_factor(t)
    1.0
end

function SimulationConfig(time_step_update_period)
    SimulationConfig(time_step_update_period=time_step_update_period)
end

function SimulationConfig(time_step_update_period, a)
    SimulationConfig(time_step_update_period=time_step_update_period, a=a)
end

function SimulationConfig(
                          ;
                          time_step_update_period=10,
                          time_step_multiplier=1,
                          a=constant_scale_factor,
                          phase_grad_limit=π/4,
                          density_threshold=1e-6
                         )
    SimulationConfig(time_step_update_period, time_step_multiplier, constant_scale_factor, phase_grad_limit, density_threshold)
end

end # module
