#!/bin/bash -e

#SBATCH --job-name       "Run an UltraDark test simulation"
#SBATCH --time           00:10:00   # Walltime
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem            4G
#SBATCH --output         log/%j.out
#SBATCH --error          log/%j.err

module load openmpi/4.0.4-gcc-4.8.5
module load hdf5/1.12.0-openmpi-4.0.4-gcc-4.8.5
module load fftw3/3.3.8-openmpi-4.0.4-gcc-4.8.5

PATH=$PATH:~/julia-1.6.1/bin/

julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'

mpirun julia --project=. soliton_velocity.jl
