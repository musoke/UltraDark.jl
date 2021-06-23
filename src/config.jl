module Config

export SimulationConfig, constant_scale_factor

struct SimulationConfig
    time_step_update_period::UInt8
    "function defining a(t)"
    a
end

function constant_scale_factor(t)
    1.0
end

function SimulationConfig(time_step_update_period)
    SimulationConfig(time_step_update_period, constant_scale_factor)
end

end # module
