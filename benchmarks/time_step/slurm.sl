#!/bin/bash -e

#SBATCH --job-name              "Benchmark JultraDark.jl threading"
#SBATCH --account               uoa00492
#SBATCH --time                  00:10:00      # Walltime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8           # Can be overridden by sbatch args, e.g. in submit.bash
#SBATCH --mem                   16G
#SBATCH --output                log/MyJuliaJob.%j.out # Include the job ID in the names of
#SBATCH --error                 log/MyJuliaJob.%j.err # the output and error files

#############
# JultraDark.jl
module purge
module load Julia/1.5.1-GCC-9.2.0-VTune gimkl

# MPI.jl warns that it is running with default MPI binary on a cluster
# Docs say to add something like this
echo "Building MPI.jl"
#julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'
julia -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'


# SLURM_CPUS_PER_TASK=4
resol_list="64 128 256"

echo "Writing output to $OUTPUT_FILE_JUD for JultraDark.jl"
for resol in $resol_list
do
        echo $(date)
        echo "Starting Julia script with threads=$SLURM_CPUS_PER_TASK, resol=$resol"
        srun julia -t $SLURM_CPUS_PER_TASK benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE_JUD
        # julia -t $threads benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE_JUD
done

#############
# PyUltraLight
module purge
module load Python/3.7.3-gimkl-2018b

cd ../PyUltraLight

echo "Writing output to $OUTPUT_FILE_PUL for PyUltraLight"
for resol in $resol_list
do
        echo $(date)
        echo "Starting PyUltraLight script with threads=$SLURM_CPUS_PER_TASK, resol=$resol"
        srun python benchmark_script.py $resol >> ../JultraDark.jl/$OUTPUT_FILE_PUL
done
