#!/bin/bash

# launch snakemake to run jobs via SLURM
# ntasks
SM_PARAMS="job-name ntasks time mem mail-user mail-type error output"

#SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"
SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task}" #new

for P in ${SM_PARAMS}; do SM_ARGS="${SM_ARGS} --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs

conda config --set ssl_verify no

## activate snakemake conda env
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')
source ${activate_conda_path} snakemake_up

## Check if slurm_status.py is executable and change if not
if [[ ! -x "scripts/slurm_status_script.py" ]]
then
    chmod +x scripts/slurm_status_script.py;
    echo "Change 'scripts/slurm_status_script.py' to executable";
fi

## create tmp and conda storage directory in /ptmp/ (if not already present)
mkdir -p /ptmp/${USER}/tmp
mkdir -p /ptmp/${USER}/tools/conda-snakemake/

### run snakemake
# -j defines total number of jobs executed in parallel on cluster
# -n dryrun
# -p print command lines
# --use-conda allows to activate conda env necessary for rule
# --conda-prefix envs: builds environment in envs folder where it will be 
snakemake -p \
    $* \
     --latency-wait 60 \
    -j 50 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status_script.py \
    --default-resources "tmpdir='/ptmp/${USER}/tmp'" \
    --group-components tax_rank_group=500 \
    --rerun-incomplete \
    --restart-times 1 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/

    #--group-components kraken_bracken=25 \
    #--dry-run \
    #--group-components \
    # --dag \
    # | dot -Tsvg > dag.svg
