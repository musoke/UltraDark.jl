module Config

export SimulationConfig, constant_scale_factor

struct SimulationConfig
    time_step_update_period::Int64
    "function defining a(t)"
    a
    phase_grad_limit::Float64
    density_threshold::Float64
end

function constant_scale_factor(t)
    1.0
end

function SimulationConfig(time_step_update_period)
    SimulationConfig(time_step_update_period, constant_scale_factor, π/4, 1e-6)
end

function SimulationConfig(time_step_update_period, a)
    SimulationConfig(time_step_update_period, a, π/4, 1e-6)
end

function SimulationConfig(
                          ;
                          time_step_update_period=10,
                          a=constant_scale_factor,
                          phase_grad_limit=π/4,
                          density_threshold=1e-6
                         )
    SimulationConfig(time_step_update_period, constant_scale_factor, phase_grad_limit, density_threshold)
end

end # module
