using BenchmarkTools
using FFTW
using LinearAlgebra
using MPI
using PencilFFTs
using Random

function fft_inplace!(a, b, fft_plan)
    mul!(b, fft_plan, a)
end

function runtests()
    suite = BenchmarkGroup()

    if ~MPI.Initialized()
        MPI.Init()
    end

    n = 128
    dims = (n, n, n)

    suite["PencilFFTs"] = BenchmarkGroup()
    for (flag_name, flags) in [
        ("FFTW.ESTIMATE", FFTW.ESTIMATE),
        ("FFTW.MEASURE", FFTW.MEASURE),
        ("FFTW.PATIENT", FFTW.PATIENT),
    ]

        # Distribute 12 processes on a 1 × 1 grid.
        comm = MPI.COMM_WORLD
        Nproc = MPI.Comm_size(comm)

        # Let MPI_Dims_create choose the decomposition.
        proc_dims = let pdims = zeros(Int, 2)
            MPI.Dims_create!(Nproc, pdims)
            pdims[1], pdims[2]
        end

        transform = Transforms.FFT()

        fft_plan = PencilFFTPlan(dims, transform, proc_dims, comm; fftw_flags = flags)

        # Allocate and initialise input data, and apply transform.
        a_pen = allocate_input(fft_plan)
        b_pen = allocate_output(fft_plan)
        randn!(a_pen)
        randn!(b_pen)

        suite["PencilFFTs"][flag_name] =
            @benchmarkable fft_inplace!($a_pen, $b_pen, $fft_plan)
    end

    suite["FFTW"] = BenchmarkGroup()
    for (flag_name, flags) in [
        ("FFTW.ESTIMATE", FFTW.ESTIMATE),
        ("FFTW.MEASURE", FFTW.MEASURE),
        ("FFTW.PATIENT", FFTW.PATIENT),
    ]
        a_arr = Array{ComplexF64}(undef, dims)
        b_arr = Array{ComplexF64}(undef, dims)
        rand!(a_arr)
        rand!(b_arr)

        fft_plan = plan_fft(a_arr; flags = flags)

        suite["FFTW"][flag_name] = @benchmarkable fft_inplace!($a_arr, $b_arr, $fft_plan)
    end

    @show tune!(suite)

    results = run(suite)

    @show res_min = minimum(results)
    @show res_med = median(results)

    @show judge(res_med["FFTW"], res_med["PencilFFTs"])
end

runtests()
