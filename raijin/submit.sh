#!/bin/bash
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=72:00:00
#PBS -l mem=5G
#PBS -l other=gdata1
#PBS -l wd
#PBS -M kevin@kdmurray.id.au
#PBS -m abe
#PBS -P xe2
set -xe

logdir=raijin/log
mkdir -p $logdir
mkdir -p data/log/
source raijin/modules.sh
export TMPDIR=${PBS_JOBFS:-$TMPDIR}

TARGET=${TARGET:-all}
SNAKEFILE=${SNAKEFILE:-Snakefile}

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name}"
QSUB="$QSUB -l wd -o $logdir -e $logdir -P {cluster.project}"

snakemake --unlock

snakemake                                                 \
    -j 100                                                \
    --cluster-config raijin/cluster.yaml                  \
    --local-cores ${PBS_NCPUS:-1}                         \
    --js raijin/jobscript.sh                              \
    --rerun-incomplete                                    \
    --keep-going                                          \
    --snakefile "$SNAKEFILE"                              \
    --cluster "$QSUB"                                     \
    "$TARGET"                                             \
    >data/log/submitter_${PBS_JOBID:-}_snakemake.log 2>&1 \

