#! /bin/bash -e
# Submit jobs with varying numbers of cpus-per-task

export OUTPUT_FILE_JUD=benchmark_output_jud.csv
export OUTPUT_FILE_PUL=benchmark_output_pul.csv
echo "threads,resol,time" > $OUTPUT_FILE_JUD
echo "threads,resol,time" > $OUTPUT_FILE_PUL

for threads in $(seq 1 16)
do
    echo "Submitting job with --cups-per-task=$threads"
    sbatch --cpus-per-task=$threads --export=ALL,OUTPUT_FILE_JUD,OUTPUT_FILE_PUL benchmarks/time_step/slurm.sl
done
