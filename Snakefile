configfile: "config.yml"
import snkmk
RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")
SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv",
                                   setfile_glob="metadata/samplesets/*.txt")
VARCALL_REGIONS = snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"])
shell.prefix = "set -euo pipefail; "

wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",


#######################################################################
#                            Read-level QC                            #
#######################################################################

rule qc_runlib:
    input:
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib)
		for run, lib in RUNLIB2SAMP],

rule read_stats:
    input:
        "data/stats/reads/readnum_librun.tsv",
        "data/stats/reads/readnum_samples.tsv",

rule qc_samples:
    input:
        expand("data/reads/samples/{sample}.fastq.gz", sample=SAMP2RUNLIB),

rule reads:
    input:
        rules.qc_runlib.input,
        rules.read_stats.input,
        rules.qc_samples.input,

### Actual rules

#localrules: qcreads
ruleorder: qcreads_il > qcreads
rule qcreads:
    input:
        r1="rawdata/runs/{run}/{lib}_R1.fastq.gz",
        r2="rawdata/runs/{run}/{lib}_R2.fastq.gz",
    output:
        reads=temp("data/reads/runs/{run}/{lib}.fastq.gz"),
    log:
        log="data/log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        4
    params:
        adp1=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter1"],
        adp2=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter2"],
        minqual=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["minqual"],
    shell:
        "( AdapterRemoval"
        "   --file1 {input.r1}"
        "   --file2 {input.r2}"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --combined-output"
        "   --interleaved-output"
        "   --trimns"
        "   --trimqualities"
        "   --trimwindows 10"
        "   --minquality {params.minqual}"
        "   --threads {threads}"
        "   --settings {log.settings}"
        "   --output1 /dev/stdout"
        " | seqhax pairs"
        "   -l 20"
        "   -b >(pigz >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"

#localrules: qcreads
rule qcreads_il:
    input:
        il="rawdata/runs/{run}/{lib}.fastq.gz",
    output:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
    log:
        log="data/log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        1
    params:
        adp1=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter1"],
        adp2=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["adapter2"],
        minqual=lambda wc: config["qc"].get(wc.run, config["qc"]["_DEFAULT_"])["minqual"],
    shell:
        "( AdapterRemoval"
        "   --file1 {input.il}"
        "   --adapter1 {params.adp1}"
        "   --adapter2 {params.adp2}"
        "   --combined-output"
        "   --interleaved"
        "   --interleaved-output"
        "   --trimns"
        "   --trimqualities"
        "   --trimwindows 10"
        "   --minquality {params.minqual}"
        "   --threads {threads}"
        "   --settings {log.settings}"
        "   --output1 /dev/stdout"
        " | seqhax pairs"
        "   -l 20"
        "   -b >(gzip >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"

rule samplefastq:
    input:
        lambda wc: ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=r, lib=l) for r, l in SAMP2RUNLIB[wc.sample]],
    output: "data/reads/samples/{sample}.fastq.gz"
    log: "data/log/samplefastq/{sample}.log"
    threads: 1
    shell:
        "cat {input} > {output}"

rule read_count_librun:
    input:
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib)
		for run, lib in RUNLIB2SAMP],
    output:
        "data/stats/reads/readnum_librun.tsv",
    threads:
        32
    log:
        "data/log/readstats/seqhax-stats-librun.log",
    shell:
        "( seqhax stats"
        "    -t {threads}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"

rule read_count_sample:
    input:
    	expand("data/reads/samples/{sample}.fastq.gz", sample=SAMP2RUNLIB),
    output:
        "data/stats/reads/readnum_samples.tsv",
    threads:
        32
    log:
        "data/log/readstats/seqhax-stats-sample.log",
    shell:
        "( seqhax stats"
        "    -t {threads}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"


#######################################################################
#                      De-novo Distance analysis                      #
#######################################################################

rule mashsketch:
    input:
        lambda wc: expand("data/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        temp("data/mash/k{ksize}-s{sketchsize}/{set}.msh"),
    log:
        "data/log/mash/sketch/k{ksize}-s{sketchsize}-{set}.log"
    threads: 32
    shell:
        " mash sketch"
        "   -k {wildcards.ksize}"
        "   -s {wildcards.sketchsize}"
        "   -p {threads}"
        "   -o {output}"
        "   {input}"
        " >{log} 2>&1"


rule mashdist:
    input:
        "data/mash/k{ksize}-s{sketchsize}/{set}.msh"
    output:
        dist="data/mash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "data/log/mash/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads: 32
    shell:
        "mash dist"
        "   -p {threads}"
        "   -t" # tabular format
        "   {input} {input}" # needs input twice
        " >{output}"
        " 2>{log}"

rule countsketch:
    input:
        "data/reads/samples/{sample}.fastq.gz",
    output:
        ct=temp("data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz"),
        info="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info",
        tsv="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info.tsv",
    log:
        "data/log/kwip/sketch/k{ksize}-s{sketchsize}-{sample}.log"
    threads:
        4
    shell:
        "load-into-counting.py"
        "   -N 1"
        "   -x {wildcards.sketchsize}"
        "   -k {wildcards.ksize}"
        "   -b"
        "   -f"
        "   -s tsv"
        "   -T {threads}"
        "   {output.ct}"
        "   {input}"
        " >{log} 2>&1"

rule kwipdist:
    input:
        lambda wc: expand("data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        d="data/kwip/k{ksize}-s{sketchsize}/{set}.dist",
        k="data/kwip/k{ksize}-s{sketchsize}/{set}.kern",
    log:
        "data/log/kwip/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads:
        4
    shell:
        "kwip"
        " -d {output.d}"
        " -k {output.k}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"

rule unique_kmers:
    input:
        lambda wc: expand("data/reads/samples/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        "data/readstats/unique-kmers/{set}.tsv",
    threads:
        32
    params:
        kmersize=config["denovodist"]["ksize"],
    log:
        "data/log/readstats/unique-kmers/{set}.log",
    shell:
        "( kdm-unique-kmers.py"
        "    -t {threads}"
        "    -k {params.kmersize}"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"


rule sourmash_sketch:
    input:
        "data/reads/samples/{sample}.fastq.gz",
    output:
        temp("data/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh"),
    log:
        "data/log/sourmash/sketch/k{ksize}-s{sketchsize}-{sample}.log"
    shell:
        "( sourmash compute"
        "   --name '{wildcards.sample}'"
        "   -k {wildcards.ksize}"
        "   -n {wildcards.sketchsize}"
        "   -o {output}"
        "   {input}"
        ") >{log} 2>&1"

rule sourmash_dist:
    input:
        lambda wc: expand("data/sourmash/sketch/k{ksize}-s{sketchsize}/{sample}.smh",
                            ksize=wc.ksize, sketchsize=wc.sketchsize,
                            sample=SAMPLESETS[wc.set]),
    output:
        "data/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
    log:
        "data/log/sourmash/dist/k{ksize}-s{sketchsize}-{set}.log"
    threads: 1
    shell:
        "(sourmash compare -k {wildcards.ksize} -o {output} {input} ) >{log} 2>&1"


rule kwip:
    input:
        expand("data/kwip/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["kwip_sketchsize"],
               set=config["denovodist"]["kwip_sets"]),

rule sourmash:
    input:
        expand("data/sourmash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["sourmash_sketchsize"],
               set=config["denovodist"]["sourmash_sets"]),

rule mash:
    input:
        expand("data/mash/k{ksize}-s{sketchsize}/{set}.dist",
               ksize=config["denovodist"]["ksize"],
               sketchsize=config["denovodist"]["mash_sketchsize"],
               set=config["denovodist"]["mash_sets"]),

rule denovo:
    input:
        rules.kwip.input,
        rules.mash.input,
        rules.sourmash.input,

#######################################################################
#                       Alignment to Reference                        #
#######################################################################

rule ngmap:
    input:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/ngm/{ref}/byrun/{run}/{lib}.bam"),
    log:
        "data/log/ngm/{ref}/{run}/{lib}.log"
    threads:
        4
    params:
        sample=lambda wc: RUNLIB2SAMP.get((wc.run, wc.lib), "{}~{}".format(wc.run, wc.lib)),
    shell:
        "( ngm"
        "   -q {input.reads}"
        "   --paired --broken-pairs"
        "   -r {input.ref}"
        "   -t {threads}"
        "   --rg-id {wildcards.run}_{wildcards.lib}"
        "   --rg-sm {params.sample}"
        "   --sensitivity 0.5" # this is the mean from a bunch of different runs
        " | samtools view -Suh -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/ngm_{wildcards.run}_{wildcards.lib}"
        "   -@ {threads}"
        "   -m 1G"
        "   -o {output.bam}"
        "   -" # stdin
        " ) >{log} 2>&1"

rule bwamem:
    input:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/bwa/{ref}/byrun/{run}/{lib}.bam"),
    log:
        "data/log/bwa/{ref}/{run}/{lib}.log"
    threads:
        4
    params:
        sample=lambda wc: RUNLIB2SAMP.get((wc.run, wc.lib), "{}~{}".format(wc.run, wc.lib)),
    shell:
        "( bwa mem"
        "   -p" # paired input
        "   -t {threads}"
        "   -R '@RG\\tID:{wildcards.run}_{wildcards.lib}\\tSM:{params.sample}'"
        "   {input.ref}"
        "   {input.reads}"
        " | samtools view -Suh -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/bwa_{wildcards.run}_{wildcards.lib}"
        "   -@ {threads}"
        "   -m 1G"
        "   -o {output.bam}"
        "   -" # stdin
        " ) >{log} 2>&1"

rule stampy:
    input:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam="data/alignments/stampy/{ref}/byrun/{run}/{lib}.bam",
    log:
        "data/log/stampy/{ref}/{run}/{lib}.log"
    threads:
        16
    params:
        sample=lambda wc: RUNLIB2SAMP.get((wc.run, wc.lib), "{}~{}".format(wc.run, wc.lib)),
    shell:
        "( mkdir /dev/shm/stampyref_{wildcards.run}_{wildcards.lib}"
        " && cp {input.ref} {input.ref}.st* /dev/shm/stampyref_{wildcards.run}_{wildcards.lib}/"
        " && stampy.py"
        "   -t {threads}"
        "   --sensitive"
        "   --substitutionrate=0.05"
        "   -g /dev/shm/stampyref_{wildcards.run}_{wildcards.lib}/$(basename {input.ref})"
        "   -h /dev/shm/stampyref_{wildcards.run}_{wildcards.lib}/$(basename {input.ref})"
        "   -M {input.reads}"
        "   --readgroup='ID:{wildcards.run}_{wildcards.lib},SM:{params.sample}'"
        " | samtools view -Suh -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/stampy_{wildcards.run}_{wildcards.lib}"
        "   -@ {threads}"
        "   -m 1G"
        "   -o {output.bam}"
        "   -" # stdin
        " ; rm -rf /dev/shm/stampyref_{wildcards.run}_{wildcards.lib}"
        " ) >{log} 2>&1"


rule bamidx:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    shell:
        "samtools index {input}"

rule mergebam_samp:
    input:
        lambda wc: ["data/alignments/{aln}/{ref}/byrun/{run}/{lib}.bam".format(
                            run=r, lib=l, aln=wc.aligner, ref=wc.ref)
	                for r, l in SAMP2RUNLIB[wc.sample]]
    output:
        bam="data/alignments/{aligner}/{ref}/samples/{sample}.bam",
    log:
        "data/log/mergesamplebam/{aligner}/{ref}/{sample}.log"
    threads: 4
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   {output.bam}"
        "   {input}"
        " ) >{log} 2>&1"


localrules: bamlist
rule bamlist:
    input:
        lambda wc: expand("data/alignments/{aligner}/{ref}/samples/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        "data/bamlists/{aligner}/{ref}/{sampleset}.bamlist",
    run:
        with open(output[0], "w") as fh:
            for s in input:
                print(s, file=fh)


rule mergebam_set:
    input:
        lambda wc: expand("data/alignments/{aligner}/{ref}/samples/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
    log:
        "data/log/mergesetbam/{aligner}/{ref}/{sampleset}.log"
    threads: 4
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   {output.bam}"
        "   {input}"
        " ) >{log} 2>&1"

## Bam stats
localrules: mergestats
rule mergestats:
    input:
        lambda wc: [ "data/stats/{type}/{aligner}/{ref}/{run}~{lib}.tsv".format(
                            run=r, lib=l, aligner=wc.aligner, ref=wc.ref, type=wc.type)
                        for r, l in RUNLIB2SAMP],
    output:
        "data/stats/{type}-{aligner}~{ref}.tsv"
    shell:
        "cat {input} > {output}"

localrules: all_bamstats
rule all_bamstats:
    input:
        expand("data/stats/bamstats_sample/{aligner}~{ref}~{sample}.tsv",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMPLESETS["all_samples"]),
        expand("data/stats/summarynumbers_{aligner}~{ref}.tsv",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"]),

rule bamstat_summary_nums:
    input:
        lambda wc: expand("data/stats/bamstats_sample/{aligner}~{ref}~{sample}.tsv",
                          aligner=wc.aligner, ref=wc.ref,
                          sample=SAMPLESETS["all_samples"]),
    output:
        "data/stats/summarynumbers_{aligner}~{ref}.tsv"
    shell:
        "./scripts/extractsn.py {input} > {output}"


rule bamstat_samps:
    input:
        "data/alignments/{aligner}/{ref}/samples/{sample}.bam",
    output:
        "data/stats/bamstats_sample/{aligner}~{ref}~{sample}.tsv"
    log:
        "data/log/bamstats_sample/{aligner}~{ref}~{sample}.tsv"
    shell:
        "(samtools stats -i 5000 -x {input} >{output}) >{log}"

rule bamstats_insertionsize:
    input:
        "data/alignments/{aligner}/{ref}/byrun/{run}/{lib}.bam",
    output:
        "data/stats/insertsize/{aligner}/{ref}/{run}~{lib}.tsv"
    log:
        "data/log/insertionsize/{aligner}/{ref}/{run}~{lib}.log"
    shell:
        "(samtools stats -i 5000 -x  {input}"
        " | grep '^IS'"
        " | sed -e 's/^IS/{wildcards.run}\\~{wildcards.lib}/'"
        " > {output})"
        ">{log} 2>&1"

rule qualstat:
    input:
        bam="data/alignments/{aligner}/{ref}/byrun/{run}/{lib}.bam",
    output:
        "data/stats/qualstat/{aligner}/{ref}/{run}~{lib}.tsv"
    log:
        "data/log/qualstat/{aligner}/{ref}/{run}~{lib}.log"
    shell:
        "(samtools view {input} "
        "   | awk '{{print $5}}'"
        "   | seqhax clihist"
        "   | sed -e 's/^/{wildcards.run}~{wildcards.lib}	/'"
        "   > {output} ) >{log} 2>&1"

rule align_librun:
    input:
        lambda wc: ["data/alignments/{aln}/{ref}/byrun/{run}/{lib}.bam".
                        format(run=r, lib=l, aln=a, ref=ref)
                        for r, l in RUNLIB2SAMP
                        for a in config["mapping"]["aligners"]
                        for ref in config["mapping"]["refs"]],

rule align_samp:
    input:
        expand("data/alignments/{aligner}/{ref}/samples/{sample}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sample=SAMP2RUNLIB),

rule align_sampset:
    input:
        rules.align_samp.input,
        expand("data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=config["mapping"]["samplesets"]),
        expand("data/bamlists/{aligner}/{ref}/{sampleset}.bamlist",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=config["mapping"]["samplesets"]),

rule bamstats:
    input:
        expand("data/stats/{type}-{aligner}~{ref}.tsv",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               type=["insertsize", "qualstat"]),

rule align:
   input:
        rules.align_samp.input,
        rules.align_sampset.input,
        rules.bamstats.input,

#######################################################################
#                           Variant Calling                           #
#######################################################################

rule freebayes:
    input:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
        bai="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam.bai",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/freebayes~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.log"
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
    shell:
        "( samtools view"  # some versions of freebayes don't seek to region
        "   -u"            # also region format is zero-based in freebayes
        "    {input.bam}"  # so we extract the correct region from the BAM
        "   '{wildcards.region}'"
        " | freebayes"
        "   --theta {params.theta}"
        "   --use-best-n-alleles 4"
        "   --min-alternate-fraction 0"
        "   --min-alternate-count 1" # per sample
        "   --min-alternate-total 3" # across all samples
        "   --min-coverage 5" # across all samples
        "   --strict-vcf"
        "   --stdin"
        "   --fasta-reference {input.ref}"
        " | bcftools view"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule mpileup:
    input:
        bam="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam",
        bai="data/alignments/{aligner}/{ref}/sets/{sampleset}.bam.bai",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/mpileup~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.log"
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
    shell:
        "( samtools mpileup"
        "   --output-tags DP,AD,SP,INFO/AD" #output everything
        "   --region '{wildcards.region}'"
        "   --fasta-ref {input.ref}"
        "   --redo-BAQ"
        "   --BCF --uncompressed"
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{wildcards.region}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior {params.theta}"
        "   -O u"
        "   -o {output.bcf}"
        " ) >{log} 2>&1"


rule bcfnorm:
    input:
        bcf="data/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run multiple filters if a pipe
        bcf=temp("data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf"),
    log:
        "data/log/bcfnormalise/{caller}~{aligner}~{ref}~{sampleset}/{region}.log"
    shell:
        "( bcftools norm"
        "   --fasta-ref {input.ref}"
        "   -O u"
        "   {input.bcf}"
        " | vt decompose_blocksub + -o -" # decompose MNP to multipe SNPs
        " | bcftools norm" # Split multi-alleics
        "   --do-not-normalize"
        "   --multiallelics -snps"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule bcffilter:
    input:
        bcf="data/variants/norm_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        # Not a pipe! can't run all regions separately if this is a pipe into merge
        bcf=temp("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf"),
    log:
        "data/log/bcffilter/{caller}~{aligner}~{ref}~{sampleset}/{filter}/{region}.log"
    params:
        filtarg=lambda wc: config["varcall"]["filters"][wc.filter].replace('\n', ' ')
    shell:
        "( bcftools view"
        "   {params.filtarg}"
        "   -O u"
        "   {input.bcf}"
        " | bcftools norm" # We normalise here to re-join multi-allelic sites, after filtering with multi-allelics split
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics +snps" # Split multi-alleic sites
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"

localrules: bcfmerge_fofn
rule bcfmerge_fofn:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
    output:
        fofn=temp("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN"),
    run:
        with open(output[0], "w") as fh:
            for s in input:
                print(s, file=fh)

rule bcfmerge:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.ref])),
        fofn="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN",
    output:
        bcf="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf",
    log:
        "data/log/mergebcf/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}.log"
    threads: 4
    shell:
        "( bcftools concat"
        "   --threads {threads}"
        "   -O b"
        "   -o {output.bcf}"
        "   --file-list {input.fofn}"
        " ) >{log} 2>&1"


rule bcf2vcf:
    input:
        bcf="{path}.bcf",
    output:
        vcf="{path}.vcf.gz",
    log:
        "data/log/bcf2vcf/{path}.log"
    threads: 4
    shell:
        "( bcftools view"
        "   {input.bcf}"
        "   -O z"
        "   --threads {threads}"
        "   -o {output.vcf}"
        " ) >{log} 2>&1"

rule variantidx:
    input:
        "{path}"
    output:
        "{path}.csi"
    priority: 2
    shell:
        "bcftools index -f {input}"

rule stats:
    input:
        "data/variants/{path}"
    output:
        "data/stats/variants/{path}.varstats"
    shell:
        "bcftools stats -s - -d 0,1000,2 --threads {threads} {input} >{output}"

rule filtered_variants:
    input:
        expand("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.{ext}",
               ext=["bcf", "bcf.csi", "vcf.gz", "vcf.gz.csi"],
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),

rule varcall:
    input:
        rules.filtered_variants.input,


#######################################################################
#                              All rule                               #
#######################################################################


rule all:
    input:
        rules.denovo.input,
        rules.reads.input,
        rules.align.input,
        rules.varcall.input,
