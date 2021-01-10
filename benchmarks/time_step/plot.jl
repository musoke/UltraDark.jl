using CSV
using Plots
using DataFrames

csv_file_julia = try
    ARGS[1]
catch
    "jultradark_benchmark_output.csv"
end

p_threads = plot(
     xlabel="threads",
     ylabel="time (s)",
     yscale=:log,
     legendtitle="resol",
)

df_julia = DataFrame(CSV.File(csv_file_julia, types=[Int64, Int64, Float64, Float64]))

for resol in unique(df_julia.resol)

    df_filtered = df_julia[df_julia.resol .== resol, :]

    plot!(
        p_threads,
        df_filtered.threads,
        df_filtered.time,
        seriestype = :scatter,
        label = "$resol^3",
    )

end

savefig("threads.pdf")

p_resol = plot(
     xlabel="resol",
     ylabel="time (s)",
     yscale=:log,
     legendtitle="threads",
)

for threads in [1, 4]

    df_filtered = df_julia[df_julia.threads .== threads, :]

    plot!(
        p_resol,
        df_filtered.resol,
        df_filtered.time,
        seriestype = :scatter,
        label = "$threads",
    )

end

savefig("resol.pdf")

p = plot(p_threads, p_resol)
savefig("bench")
