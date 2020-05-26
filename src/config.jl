module Config

struct SimulationConfig
    time_step_update_period::UInt8
    "function defining a(t)"
    a
end

end # module