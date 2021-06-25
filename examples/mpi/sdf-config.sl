#!/bin/bash -e

#SBATCH --job-name       "Run an UltraDark test simulation"
#SBATCH --time           00:10:00   # Walltime
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem            100G
#SBATCH --output         log/%j.out
#SBATCH --error          log/%j.err

module load openmpi/4.0.4-gcc-4.8.5
module load hdf5/1.12.0-openmpi-4.0.4-gcc-4.8.5
module load fftw3/3.3.8-openmpi-4.0.4-gcc-4.8.5

PATH=$PATH:~/julia-1.6.1/bin/

echo "Building MPI.jl" >&2
julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'

echo "Running simulation" >&2
export JULIA_DEBUG=UltraDark
mpirun julia --project ../soliton_velocity.jl

echo "Job done" >&2
