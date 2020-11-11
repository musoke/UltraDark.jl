#!/bin/bash -e

#SBATCH --job-name              "Benchmark JultraDark.jl threading"
#SBATCH --account               uoa00492
#SBATCH --time                  00:30:00      # Walltime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem                   16G
#SBATCH --output                log/MyJuliaJob.%j.out # Include the job ID in the names of
#SBATCH --error                 log/MyJuliaJob.%j.err # the output and error files
#SBATCH --mail-type=ALL
#SBATCH --mail-user="n.musoke@auckland.ac.nz"

module purge
module load Julia/1.5.1-GCC-9.2.0-VTune gimkl

# MPI.jl warns that it is running with default MPI binary on a cluster
# Docs say to add something like this
echo "Building MPI.jl"
#julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'
julia -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'

OUTPUT_FILE=jultradark_benchmark_output.csv
rm $OUTPUT_FILE
echo "threads,resol,time" > $OUTPUT_FILE

# SLURM_CPUS_PER_TASK=4

for resol in 64 128 256 512
do
    for threads in $(seq 1 $SLURM_CPUS_PER_TASK)
    do
        echo $(date)
        echo "Starting Julia script with threads=$threads, resol=$resol"
        srun julia -t $threads benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE
        # julia -t $threads benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE
    done
done
