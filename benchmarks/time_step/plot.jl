using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFrames

function main()

    df = load_generate_data()
    generate_data_resol(df) |> plot_resol
    generate_data_cpus(df) |> plot_cpus

end

function load_generate_data()

    clustername = "plasma"
    # clustername = "larb"

    csv_file_ud = try
        ENV["OUTPUT_FILE_UD"]
    catch
        "benchmark_output_ud-$clustername.csv"
    end

    csv_file_pul = try
        ENV["OUTPUT_FILE_PUL"]
    catch
        "benchmark_output_pul.csv"
    end

    df = CSV.read(
                  [csv_file_ud],
                  DataFrame,
                  types=[UInt32, UInt32, UInt32, Float64, String], strict=true,
                 )

    transform!(df, [:tasks, :threads] => ((tasks, threads) -> tasks.*threads) => :total_cpus)

end

function generate_data_cpus(df)

    filtered_df = copy(df)

    # Add column with expected time vs CPUs
    gdf = groupby(df, [:grids_type, :resol, ])
    maxdf = combine(gdf, :time => maximum => :max_time_cpu)
    leftjoin!(filtered_df, maxdf; on=[:grids_type, :resol, ])

    theory_time(cpus, max_time) = max_time ./ cpus
    transform!(filtered_df, [:total_cpus, :max_time_cpu] => theory_time => :theory_time)

    sort!(filtered_df, [:resol, :total_cpus])

    # Only consider the highest resolution
    resol_max(resol) = resol == maximum(filtered_df.resol)
    filtered_df = filter!(:resol => resol_max, filtered_df)

    # Only consider PencilGrids with single threads per job
    filter_pencilgrids_one_thread(grids_type, threads) = grids_type != "PencilGrids" || threads == 1
    filtered_df = filter!([:grids_type, :threads] => filter_pencilgrids_one_thread, filtered_df)

end

function plot_cpus(df)

    time_data = data(df) * mapping(:total_cpus => "number of CPUs", :time => "time (s)") * visual(Scatter) * mapping(marker=:grids_type)
    time_theory = data(df) * mapping(:total_cpus => "number of CPUs", :theory_time => "time (s)") * visual(Lines) * mapping(linestyle=:grids_type)

    plt = (time_data + time_theory) * mapping(color=:grids_type)

    fg = draw(plt, axis=(yscale=log10,))
    save("cpus.pdf", fg)
end

function generate_data_resol(df)

    filtered_df = copy(df)

    cpus_equal_1(cpus) = cpus == 1
    filtered_df = filter!(:total_cpus => cpus_equal_1, filtered_df)

    # Add column with expected time vs resolution
    gdf = groupby(filtered_df, [:grids_type, :tasks, :threads])
    maxdf = combine(gdf, [:time, :resol] .=> minimum .=> [:min_time_resol, :min_resol])
    leftjoin!(filtered_df, maxdf; on=[:grids_type, :tasks, :threads])

    fit_log_cube(resol, min_resol, min_time) = min_time .* (resol ./ min_resol).^3 .* log10.(resol ./ min_resol * 10)

    transform!(filtered_df, [:resol, :min_resol, :min_time_resol] => fit_log_cube => :theoretical_time)

    sort!(filtered_df)

end

function plot_resol(df)

    time_data = data(df) * mapping(:resol => "resolution", :time => "time (s)") * visual(Scatter) * mapping(marker=:grids_type)
    time_theory = data(df) * mapping(:resol => "resolution", :theoretical_time => "time (s)") * visual(Lines) * mapping(linestyle=:grids_type)

    plt = (time_data + time_theory) * mapping(color=:grids_type)

    fg = draw(plt, axis=(yscale=log10,),)
    save("resol.pdf", fg)
end
