#!/bin/bash -e

#SBATCH --job-name              "Benchmark UltraDark.jl threading"
#SBATCH --time                  00:20:00      # Walltime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1       # Can be overridden by sbatch args, e.g. in submit.bash
#SBATCH --mem                   4G
#SBATCH --output                log/benchmark.%j.out # Include the job ID in the names of
#SBATCH --error                 log/benchmark.%j.err # the output and error files

echo CLUSTER=$CLUSTER >&2
echo SLURM_CLUSTER_NAME=$SLURM_CLUSTER_NAME >&2
echo ULTRADARK_GRIDS_TYPE=$ULTRADARK_GRIDS_TYPE >&2
echo SLURM_NTASKS=$SLURM_NTASKS, SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK >&2
echo RESOL_LIST=$RESOL_LIST >&2

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
        module purge
        # module load openmpi
        # module load fftw3
        module load gcc/8.1.0
        module load craype-x86-rome
        module load craype-network-infiniband
        module load cray-libsci/20.03.1
        module load craype/2.6.4
        module load perftools-base/20.03.0
        module load cce/9.1.3
        module load PrgEnv-cray/1.0.6
        module load cray-impi/5
        module load cray-fftw/3.3.8.5
        # module load cray-fftw_impi/3.3.8.5
        SRUN_CMD="mpiexec"
        ;;
    "MAHUIKA")
        echo "Configuring for mahuika" >&2
        module purge
        module load Julia/1.5.1-GCC-9.2.0-VTune gimkl
        SRUN_CMD=srun
        ;;
    *)
        echo "Configuring for generic" >&2
        SRUN_CMD="mpiexec"
esac

if [ "$ULTRADARK_GRIDS_TYPE" == "Grids" ]; then
    SRUN_CMD=""
fi

echo SRUN_CMD=$SRUN_CMD >&2

julia --project -e 'using Pkg; Pkg.status()' >&2

# MPI.jl warns that it is running with default MPI binary on a cluster
# Docs say to add something like this
echo "Building MPI.jl" >&2
julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'


#############
# UltraDark.jl
echo "Writing output to $OUTPUT_FILE_UD for UltraDark.jl" >&2

for resol in $RESOL_LIST
do
        export ULTRADARK_RESOL=$resol
        echo $(date) >&2
        echo "Starting Julia script with tasks=$SLURM_NTASKS, threads=$SLURM_CPUS_PER_TASK, resol=$ULTRADARK_RESOL" >&2
        $SRUN_CMD julia --project -t $SLURM_CPUS_PER_TASK time_step.jl $resol >> $OUTPUT_FILE_UD
done

##############
## PyUltraLight
#module purge
#module load Python/3.7.3-gimkl-2018b

#cd ../PyUltraLight

#echo "Writing output to $OUTPUT_FILE_PUL for PyUltraLight" >&2
#for resol in $RESOL_LIST
#do
#        echo $(date) >&2
#        echo "Starting PyUltraLight script with threads=$SLURM_CPUS_PER_TASK, resol=$resol" >&2
#        srun python benchmark_script.py $resol >> ../UltraDark.jl/$OUTPUT_FILE_PUL
#done
