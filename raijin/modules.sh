module purge
module load pbs
module load adapterremoval nextgenmap/0.5.5 samtools bcftools htslib freebayes/v1.2.0
module load plink19 mash pigz seqhax sra-toolkit vt gatk4
export TMPDIR=${PBS_JOBFS:-/tmp}
source /short/xe2/kdm801/miniconda3/etc/profile.d/conda.sh
conda activate paneuc
