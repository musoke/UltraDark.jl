using CSV
using Plots
using DataFrames

csv_file_jud = try
    ARGS[1]
catch
    "benchmark_output_jud.csv"
end

csv_file_pul = try
    ARGS[2]
catch
    "benchmark_output_pul.csv"
end

p_threads = plot(
     xlabel="threads",
     ylabel="time (s)",
     yscale=:log,
     legendtitle="resol",
     legend = :outerright,
)

df_jud = DataFrame(CSV.File(csv_file_jud, types=[Int64, Int64, Float64, Float64]))
df_pul = DataFrame(CSV.File(csv_file_pul, types=[Int64, Int64, Float64, Float64]))

for (resol, color) in zip(unique(df_jud.resol), [:blue, :red, :black])

    df_jud_filtered = df_jud[df_jud.resol .== resol, :]
    df_pul_filtered = df_pul[df_pul.resol .== resol, :]

    sort!(df_jud_filtered, :threads)

    df_jud_filtered[!, :theory] = df_jud_filtered[1, :time] .* df_jud_filtered[1, :threads] ./ df_jud_filtered[!, :threads]

    plot!(
        p_threads,
        df_jud_filtered.threads,
        df_jud_filtered.time,
        seriestype = :scatter,
        color=color,
        label = "UltraDark.jl $resol^3",
    )

    plot!(
        p_threads,
        df_pul_filtered.threads,
        df_pul_filtered.time,
        seriestype = :scatter,
        color=color,
        marker=:x,
        label = "PyUltraLight $resol^3",
    )

    plot!(
        p_threads,
        df_jud_filtered.threads,
        df_jud_filtered.theory,
        color=color,
        label = "theory $resol^3",
    )

end

savefig("threads.pdf")

p_resol = plot(
     xlabel="resol",
     ylabel="time (s)",
     yscale=:log,
     legendtitle="threads",
     legend = :outerright,
)

for threads in [1, 8]

    df_jud_filtered = df_jud[df_jud.threads .== threads, :]

    sort!(df_jud_filtered, :resol)

    df_jud_filtered[!, :theory] = df_jud_filtered[1, :time] ./ df_jud_filtered[1, :resol]^3 .* df_jud_filtered[!, :resol].^3

    plot!(
        p_resol,
        df_jud_filtered.resol,
        df_jud_filtered.time,
        seriestype = :scatter,
        label = "$threads",
    )

    plot!(
        p_resol,
        df_jud_filtered.resol,
        df_jud_filtered.theory,
        label = "theory, $threads",
    )

end

savefig("resol.pdf")

p = plot(p_threads, p_resol)
savefig("bench.pdf")
