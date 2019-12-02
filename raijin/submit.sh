#!/bin/bash
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=72:00:00
#PBS -l mem=5G
#PBS -l other=gdata1
#PBS -l wd
#PBS -j oe
#PBS -M kevin@kdmurray.id.au
#PBS -m abe
#PBS -P xe2

logdir=raijin/log
mkdir -p $logdir
mkdir -p data/log/
source raijin/modules.sh
export TMPDIR=${PBS_JOBFS:-$TMPDIR}

set -ueo pipefail
TARGET=${TARGET:-all}
SNAKEFILE=${SNAKEFILE:-Snakefile}

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name}"
QSUB="$QSUB -l wd -j oe -o $logdir -P {cluster.project}"

if [ "${RMTEMP:-yes}" == yes ]
then
	temp=''
else
	temp='--notemp'
fi

snakemake                                                          \
    -j 1000                                                        \
    --cluster-config raijin/cluster.yaml                           \
    --local-cores ${PBS_NCPUS:-1}                                  \
    --js raijin/jobscript.sh                                       \
    --nolock                                                       \
    --rerun-incomplete                                             \
    --keep-going                                                   \
    $temp                                                          \
    --snakefile "$SNAKEFILE"                                       \
    --cluster "$QSUB"                                              \
    "$TARGET"                                                      \
    |& tee data/log/submitter_${PBS_JOBID:-headnode}_snakemake.log \

