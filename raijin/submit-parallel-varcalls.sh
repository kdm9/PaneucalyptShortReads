set -u
ALIGNER=bwa
REF=grandisv2chl
REFPATH=/g/data1/xe2/references/eucalyptus/grandis_v2_chloro/Egrandis-v2-plus-chloro.fasta
SIZE=100000
for CALLER in mpileup freebayes
do
    for SAMPLESET in all_samples
    do
	# Don't process the whole way to filtering, causes race conditions. Filtering is quick and can happen later.
        qsub -N VC_${CALLER}_${SAMPLESET} -v SIZE=${SIZE},CALLER=${CALLER},ALIGNER=${ALIGNER},REF=${REF},REFPATH=${REFPATH},SAMPLESET=${SAMPLESET} raijin/parallel-varcall-one.pbs -l walltime=1:00:00 -q express
    done
done
