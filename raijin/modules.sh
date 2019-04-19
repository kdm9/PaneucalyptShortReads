module purge
source /short/xe2/kdm801/miniconda3/etc/profile.d/conda.sh
conda activate paneuc
module load pbs
module load adapterremoval nextgenmap/0.5.5 samtools bcftools htslib freebayes/v1.2.0
module load bwa plink19 mash khmer kwip pigz seqhax R/3.4.3 stampy sra-toolkit vt sourmash gatk4
module load angsd pcangsd
export TMPDIR=${PBS_JOBFS:-/tmp}
