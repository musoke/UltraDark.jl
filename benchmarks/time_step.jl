using BenchmarkTools
using UltraDark
using UltraDark: psi_half_step!, psi_whole_step!, phi_whole_step!, take_steps!
using MPI

if ~MPI.Initialized()
    MPI.Init()
end

n = 64
Δt = 0.1

suite = BenchmarkGroup()

for (label, grids) in [("Grids", Grids(1.0, n)), ("PencilGrids", PencilGrids(1.0, n))]

    suite[label] = BenchmarkGroup()

    suite[label]["psi_half_step"] = @benchmarkable psi_half_step!($Δt, $grids)
    suite[label]["psi_whole_step"] = @benchmarkable psi_whole_step!($Δt, $grids)
    suite[label]["phi_whole_step"] = @benchmarkable phi_whole_step!($Δt, $grids)


    suite[label]["take_steps"] = @benchmarkable take_steps!(
        $grids,
        $1.0,
        $Δt,
        1,
        $(OutputConfig("output", 1:2)),
        $(Config.constant_scale_factor),
    )
    suite[label]["take_steps 10"] = @benchmarkable take_steps!(
        $grids,
        $1.0,
        $Δt,
        10,
        $(OutputConfig("output", 1:2)),
        $(Config.constant_scale_factor),
    )

end

@show tune!(suite)

results = run(suite)

@show res_min = minimum(results)
@show res_med = median(results)

judge(res_med["Grids"], res_med["PencilGrids"])
