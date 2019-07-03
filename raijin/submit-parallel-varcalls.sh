set -u
ALIGNER=bwa
REF=grandisv2chl
REFPATH=/g/data1/xe2/references/eucalyptus/grandis_v2_chloro/Egrandis-v2-plus-chloro.fasta
REF=emel-scott
REFPATH=/g/data1/xe2/references/eucalyptus/emel_scott/Emelliodora_CSIROg1_SISH00000000.1.fasta
SIZE=100000
#for CALLER in mpileup freebayes
for CALLER in mpileup
do
    for SAMPLESET in LBEmel DP15-everything all_samples
    do
        if [ "$CALLER" == freebayes ]
        then
            req="-l walltime=24:00:00,ncpus=768,mem=3000G"
        else
            req="-l walltime=6:00:00,ncpus=512,mem=1000G"
        fi
        qsub $req -N VC_${CALLER}_${SAMPLESET} -v SIZE=${SIZE},CALLER=${CALLER},ALIGNER=${ALIGNER},REF=${REF},REFPATH=${REFPATH},SAMPLESET=${SAMPLESET} raijin/parallel-varcall-one.pbs
    done
done
