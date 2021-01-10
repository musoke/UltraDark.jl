#! /bin/bash -e
# Submit jobs with varying numbers of cpus-per-task

export OUTPUT_FILE=jultradark_benchmark_output.csv
rm $OUTPUT_FILE
echo "threads,resol,time" > $OUTPUT_FILE

for threads in $(seq 1 16)
do
    echo "Submitting job with --cups-per-task=$threads"
    sbatch --cpus-per-task=$threads --export=ALL,OUTPUT_FILE benchmarks/time_step/slurm.sl
done
