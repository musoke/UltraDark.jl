#! /bin/bash -e
# Submit jobs with varying numbers of cpus-per-task

case $(hostname) in
    "sdf-login01")
        export CLUSTER=sdf
        ;;
    "ln-0001")  # plasma
        export CLUSTER=plasma
        ;;
    "MAHUIKA")
        export CLUSTER=mahuika
        ;;
    *)
        export CLUSTER=$(hostname)
        ;;
esac

export OUTPUT_FILE_UD=benchmark_output_ud-$CLUSTER.csv
# export OUTPUT_FILE_PUL=benchmark_output_pul-$(CLUSTER).csv

echo "threads,resol,time" > $OUTPUT_FILE_UD
# echo "threads,resol,time" > $OUTPUT_FILE_PUL

export ULTRADARK_GRIDS_TYPE=Grids

for threads in $(seq 1 16)
do
    echo "Submitting job with --cups-per-task=$threads"
    sbatch --cpus-per-task=$threads --export=ALL,OUTPUT_FILE_UD,OUTPUT_FILE_PUL,ULTRADARK_GRIDS_TYPE,CLUSTER slurm.sh
done
