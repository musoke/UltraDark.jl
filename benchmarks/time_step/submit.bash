#! /bin/bash -e
# Submit jobs with varying numbers of cpus-per-task

case $(hostname) in
    "sdf-login01")
        export CLUSTER=sdf
        ;;
    "ln-0001")
        export CLUSTER=plasma
        ;;
    "ln-0002")
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

echo "tasks,threads,resol,time,grids_type" > $OUTPUT_FILE_UD
# echo "tasks,threads,resol,time,grids_type" > $OUTPUT_FILE_PUL

MAX_CPUS_TOTAL=16
MAX_THREADS=16
MAX_TASKS=16
declare -a grids_array=("Grids", "PencilGrids")

for nthreads in $(seq 1 $MAX_THREADS); do

    for ntasks in $(seq 1 $MAX_TASKS); do

        for ULTRADARK_GRIDS_TYPE in "Grids" "PencilGrids"; do

            export ULTRADARK_GRIDS_TYPE=$ULTRADARK_GRIDS_TYPE #So slurm job can inherit it

            if [ $nthreads = 1 ] && [ $ntasks = 1 ]; then
                export RESOL_LIST=$(seq 64 32 256)
            else
                export RESOL_LIST="64 128 256"
            fi

            if [ $ntasks -gt 1 ] && [ $ULTRADARK_GRIDS_TYPE != "PencilGrids" ]; then
                continue
            fi

            if [ $MAX_CPUS_TOTAL -lt $((nthreads*ntasks)) ]; then
                continue
            fi

            sbatch --cpus-per-task=$nthreads --ntasks=$ntasks --export=ALL,OUTPUT_FILE_UD,OUTPUT_FILE_PUL,ULTRADARK_GRIDS_TYPE,CLUSTER,RESOL_LIST slurm.sh
            echo "Submitted job with --cpus-per-task=$nthreads --ntasks=$ntasks, ULTRADARK_GRIDS_TYPE=$ULTRADARK_GRIDS_TYPE"
        done
    done
done
