#!/bin/bash -e

#SBATCH --job-name       "Run an UltraDark test simulation"
#SBATCH --time           00:10:00   # Walltime
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem            4G
#SBATCH --output         log/%j.out
#SBATCH --error          log/%j.err

module load openmpi
module load fftw3

PATH=$PATH:~/julia-1.6.1/bin/

julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'

mpiexec julia --project ../soliton_velocity.jl
