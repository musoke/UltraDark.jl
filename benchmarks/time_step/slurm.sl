#!/bin/bash -e

#SBATCH --job-name              "Benchmark JultraDark.jl threading"
#SBATCH --account               uoa00492
#SBATCH --time                  00:10:00      # Walltime
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8           # Can be overridden by sbatch args, e.g. in submit.bash
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


# SLURM_CPUS_PER_TASK=4

echo "Writing output to $OUTPUT_FILE"

for resol in 64 128 256
do
        echo $(date)
        echo "Starting Julia script with threads=$SLURM_CPUS_PER_TASK, resol=$resol"
        srun julia -t $SLURM_CPUS_PER_TASK benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE
        # julia -t $threads benchmarks/time_step/time_step.jl $resol >> $OUTPUT_FILE
done
