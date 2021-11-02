#!/bin/bash -e

#SBATCH --job-name              "Benchmark UltraDark.jl threading"
#SBATCH --time                  00:20:00      # Walltime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1       # Can be overridden by sbatch args, e.g. in submit.bash
#SBATCH --mem                   30G
#SBATCH --output                log/benchmark.%j.out # Include the job ID in the names of
#SBATCH --error                 log/benchmark.%j.err # the output and error files

echo CLUSTER=$CLUSTER >&2
echo SLURM_CLUSTER_NAME=$SLURM_CLUSTER_NAME >&2

case $CLUSTER in
    "sdf")
        echo "Configuring for sdf" >&2
        module purge
        module load openmpi/4.0.4-gcc-4.8.5
        module load hdf5/1.12.0-openmpi-4.0.4-gcc-4.8.5
        # module load fftw3/3.3.8-openmpi-4.0.4-gcc-4.8.5
        SRUN_CMD=""
        ;;
    "plasma")
        echo "Configuring for plasma" >&2
        module load openmpi
        module load fftw3
        SRUN_CMD=""
        ;;
    "MAHUIKA")
        echo "Configuring for mahuika" >&2
        module purge
        module load Julia/1.5.1-GCC-9.2.0-VTune gimkl
        SRUN_CMD=srun
        ;;
    *)
        SRUN_CMD=srun
esac

julia --project -e 'using Pkg; Pkg.status()'

# MPI.jl warns that it is running with default MPI binary on a cluster
# Docs say to add something like this
echo "Building MPI.jl" >&2
julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'


resol_list="64 128 256 512"

#############
# UltraDark.jl
echo "Writing output to $OUTPUT_FILE_UD for UltraDark.jl" >&2

for resol in $resol_list
do
        export ULTRADARK_RESOL=$resol
        echo $(date) >&2
        echo "Starting Julia script with threads=$SLURM_CPUS_PER_TASK, resol=$ULTRADARK_RESOL" >&2
        $SRUN_CMD julia --project -t $SLURM_CPUS_PER_TASK time_step.jl $resol >> $OUTPUT_FILE_UD
done

##############
## PyUltraLight
#module purge
#module load Python/3.7.3-gimkl-2018b

#cd ../PyUltraLight

#echo "Writing output to $OUTPUT_FILE_PUL for PyUltraLight" >&2
#for resol in $resol_list
#do
#        echo $(date) >&2
#        echo "Starting PyUltraLight script with threads=$SLURM_CPUS_PER_TASK, resol=$resol" >&2
#        srun python benchmark_script.py $resol >> ../UltraDark.jl/$OUTPUT_FILE_PUL
#done
