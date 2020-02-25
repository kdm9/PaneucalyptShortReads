configfile: "config.yml"
import snkmk
RUNLIB2SAMP, SAMP2RUNLIB = snkmk.make_runlib2samp("metadata/sample2runlib.csv")
SAMPLESETS = snkmk.make_samplesets(s2rl_file="metadata/sample2runlib.csv",
                                   setfile_glob="metadata/samplesets/*.txt")

VARCALL_REGIONS = {
    vc: snkmk.make_regions(config["refs"], window=config["varcall"]["chunksize"][vc])
    for vc in config["varcall"]["chunksize"]
}

CHROMS = snkmk.make_chroms(config["refs"])
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
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib) for run, lib in RUNLIB2SAMP],
        expand("data/mash/everything/k{ksize}-s{sketchsize}_everything_librun.dist", 
               ksize=config["denovodist"]["ksize"], sketchsize=config["denovodist"]["mash_sketchsize"]),

rule read_stats:
    input:
        "data/stats/reads/readnum_librun.tsv",
        #"data/stats/reads/readnum_samples.tsv",

rule reads:
    input:
        rules.qc_runlib.input,
        rules.read_stats.input,

### Actual rules

ruleorder: qcreads_il > qcreads
rule qcreads:
    input:
        r1=ancient("rawdata/runs/{run}/{lib}_R1.fastq.gz"),
        r2=ancient("rawdata/runs/{run}/{lib}_R2.fastq.gz"),
    output:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
    log:
        log="data/log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        8
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
        "   -b >(pigz -p {threads} >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"

rule qcreads_il:
    input:
        il="rawdata/runs/{run}/{lib}.fastq.gz",
    output:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
    log:
        log="data/log/adapterremoval/{run}/{lib}.log",
        settings="data/stats/adapterremoval/{run}/{lib}.txt",
    threads:
        8
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
        "   -b >(pigz -p {threads} >{output.reads})"
        "   /dev/stdin"
        ") >{log.log} 2>&1"


rule samplefastqpipe:
    input:
        lambda wc: ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=r, lib=l) for r, l in SAMP2RUNLIB[wc.sample]],
    output: pipe("data/reads/samples_pipe/{sample}.fastq.gz")
    log: "data/log/samplefastqpipe/{sample}.log"
    threads: 1
    shell:
        "cat {input} > {output}"

rule samplefastqfile:
    input:
        lambda wc: ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=r, lib=l) for r, l in SAMP2RUNLIB[wc.sample]],
    output: "data/reads/samples/{sample}.fastq.gz"
    log: "data/log/samplefastqfile/{sample}.log"
    threads: 1
    shell:
        "cat {input} > {output}"


localrules: sample_fastqs
rule sample_fastqs:
    input:
        [expand("data/reads/samples/{sample}.fastq.gz", sample=SAMPLESETS[sset])
            for sset in config["persample_reads"]["samplesets"]]

#rule read_count_librun:
#    input:
#        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib)
#		for run, lib in RUNLIB2SAMP],
#    output:
#        "data/stats/reads/readnum_librun.tsv",
#    threads:
#        48
#    log:
#        "data/log/readstats/seqhax-stats-librun.log",
#    shell:
#        "( seqhax stats"
#        "    -t {threads}"
#        "    {input}"
#        "    >{output}"
#        " ) 2>{log}"

rule read_count_librun_indiv:
    input:
        "data/reads/runs/{run}/{lib}.fastq.gz"
    output:
        "data/stats/reads/readnum_librun/{run}~{lib}.tsv",
    log:
        "data/log/readstats/seqhax-stats-librun/{run}~{lib}.log",
    shell:
        "( seqhax stats"
        "    {input}"
        "    >{output}"
        " ) 2>{log}"

rule read_count_fromindiv:
    input:
        ["data/stats/reads/readnum_librun/{run}~{lib}.tsv".format(run=run, lib=lib)
          for run, lib in RUNLIB2SAMP],
    output:
        "data/stats/reads/readnum_librun.tsv",
    threads:
        1
    run:
        with open(output[0], "w") as fh:
            for i, tsv in enumerate(input):
                with open(tsv) as tsvfh:
                    if i > 0:
                        next(tsvfh)  # skip header on all but first file
                    for line in tsvfh:
                        fh.write(line)


rule read_count_sample:
    input:
    	expand("data/reads/samples_pipe/{sample}.fastq.gz", sample=SAMP2RUNLIB),
    output:
        "data/stats/reads/readnum_samples.tsv",
    threads:
        47
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

rule everything_mash_sketch:
    input:
        ["data/reads/runs/{run}/{lib}.fastq.gz".format(run=run, lib=lib) for run, lib in RUNLIB2SAMP],
    output:
        temp("data/mash/everything/k{ksize}-s{sketchsize}_everything_librun.msh"),
    log:
        "data/log/mash/everything/k{ksize}-s{sketchsize}.log"
    threads: 49
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
        "data/mash/everything/k{ksize}-s{sketchsize}_everything_librun.msh",
    output:
        dist="data/mash/everything/k{ksize}-s{sketchsize}_everything_librun.dist",
    log:
        "data/log/mash/everything/dist_k{ksize}-s{sketchsize}.log"
    threads: 48
    shell:
        "mash dist"
        "   -p {threads}"
        "   -t" # tabular format
        "   {input} {input}" # needs input twice
        " >{output}"
        " 2>{log}"

rule countsketch:
    input:
        "data/reads/samples_pipe/{sample}.fastq.gz",
    output:
        ct=temp("data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz"),
        info="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info",
        tsv="data/kwip/sketch/k{ksize}-s{sketchsize}/{sample}.ct.gz.info.tsv",
    log:
        "data/log/kwip/sketch/k{ksize}-s{sketchsize}-{sample}.log"
    threads:
        3
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
        lambda wc: expand("data/reads/samples_pipe/{sample}.fastq.gz",
                          sample=SAMPLESETS[wc.set]),
    output:
        "data/readstats/unique-kmers/{set}.tsv",
    threads:
        27
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
        "data/reads/samples_pipe/{sample}.fastq.gz",
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
        bam=temp("data/alignments/byrun.raw/ngm/{ref}/{run}/{lib}~{sample}.bam"),
    log:
        "data/log/ngm/{ref}/{run}/{lib}~{sample}.log"
    threads:
        8
    params:
        sensitivity=config["mapping"]["ngm"]["sensitivity"],
    shell:
        "( ngm"
        "   -q {input.reads}"
        "   --paired --broken-pairs"
        "   -r {input.ref}"
        "   -t {threads}"
        "   --rg-id {wildcards.run}_{wildcards.lib}_{wildcards.sample}"
        "   --rg-sm {wildcards.sample}"
        "   --sensitivity {params.sensitivity}" # this is the mean from a bunch of different runs
        "| samtools view -Suh - >{output.bam}"
        " ) >{log} 2>&1"

rule bwamem:
    input:
        reads="data/reads/runs/{run}/{lib}.fastq.gz",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/byrun.raw/bwa/{ref}/{run}/{lib}~{sample}.bam"),
    log: "data/log/bwa/{ref}/{run}/{lib}~{sample}.log"
    threads:
        12
    shell:
        "( bwa mem"
        "   -p" # paired input
        "   -t {threads}"
        "   -R '@RG\\tID:{wildcards.run}_{wildcards.lib}_{wildcards.sample}\\tSM:{wildcards.sample}'"
        "   {input.ref}"
        "   {input.reads}"
        "| samtools view -Suh - >{output.bam}"
        " ) >{log} 2>&1"

rule bam_markdups_sort:
    input:
        bam="data/alignments/byrun.raw/{aligner}/{ref}/{run}/{lib}~{sample}.bam",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bam=temp("data/alignments/byrun/{aligner}/{ref}/{run}/{lib}~{sample}.bam"),
    threads: 4
    log: "data/log/markdup/{aligner}/{ref}/{run}/{lib}~{sample}.log"
    shell:
        "( samtools fixmate "
        "   -m"
        "   -@ {threads}"
        "   --output-fmt bam,level=0"
        "   {input.bam}"
        "   -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.run}_{wildcards.lib}_sort_$RANDOM"
        "   --output-fmt bam,level=0"
        "   -@ {threads}"
        "   -m 1g"
        "   -"
        " | samtools markdup"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.run}_{wildcards.lib}_markdup_$RANDOM"
        "   -s" # report stats
        "   -@ {threads}"
        "   --output-fmt bam,level=3"
        "   -"
        "   {output.bam}"
        " ) >{log} 2>&1"



rule mergebam_samp:
    input:
        lambda wc: ["data/alignments/byrun/{aln}/{ref}/{run}/{lib}~{sample}.bam".format(
                            run=r, lib=l, aln=wc.aligner, ref=wc.ref, sample=wc.sample)
	                for r, l in SAMP2RUNLIB[wc.sample]]
    output:
        bam="data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    log:
        "data/log/mergesamplebam/{aligner}/{ref}/{sample}.log"
    threads: 8
    priority: 1 # so the temps get cleaned sooner
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   --output-fmt bam,level=4"
        "   {output.bam}"
        "   {input}"
        " ) >{log} 2>&1"


rule qualimap_samp:
    input:
        bam="data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    output:
        directory("data/alignments/qualimap/samples/{aligner}~{ref}~{sample}/"),
    log:
        "data/log/qualimap_sample/{aligner}~{ref}~{sample}.log"
    threads: 4
    shell:
        "( unset DISPLAY; qualimap bamqc"
        "   --java-mem-size=6G"
        "   -bam {input.bam}"
        "   -nr 10000"
        "   -nt {threads}"
        "   -outdir {output}"
        "   {input}"
        " ) >{log} 2>&1"


localrules: bamlist
rule bamlist:
    input:
        bam=lambda wc: expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),
        bai=lambda wc: expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam.bai",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        "data/alignments/bamlists/{aligner}~{ref}~{sampleset}.bamlist",
    run:
        with open(output[0], "w") as fh:
            for s in input.bam:
                print(s, file=fh)


rule mergebam_set:
    input:
        lambda wc: expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
                          aligner=wc.aligner, ref=wc.ref, sample=SAMPLESETS[wc.sampleset]),

    output:
        bam="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
        bai="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam.bai",
    log:
        "data/log/mergesetbam/{aligner}/{ref}/{sampleset}.log"
    threads: 16
    shell:
        "( samtools merge"
        "   --output-fmt bam,level=7"
        "   -@ {threads}"
        "   -"
        "   {input}"
        " | tee {output.bam}"
        " | samtools index - {output.bai}"  # indexing takes bloody ages, we may as well do this on the fly
        " ) >{log} 2>&1"


localrules: bamstat_samps
rule bamstat_samps:
    input:
        "data/alignments/samples/{aligner}/{ref}/{sample}.bam",
    output:
        "data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.tsv",
    log:
        "data/log/bamstats_sample/{aligner}~{ref}~{sample}.tsv"
    shell:
        "(samtools stats -i 5000 -x {input} >{output}) >{log}"



# Utils
ruleorder: mergebam_set > bamidx # as we index on the fly for sets
rule bamidx:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    log:
        "data/log/bamindex/{path}.log"
    shell:
        "samtools index {input}"

#ruleorder: mergebam_samp > sam2bam
#rule sam2bam:
#    input:
#        "{path}.sam"
#    output:
#        "{path}.bam"
#    shell:
#        "samtools view -Suh {input} >{output}"


### Align targets

rule align_librun:
    input:
        lambda wc: ["data/alignments/byrun/{aln}/{ref}/{run}/{lib}.bam".
                        format(run=r, lib=l, aln=a, ref=ref)
                        for r, l in RUNLIB2SAMP
                        for a in config["mapping"]["aligners"]
                        for ref in config["mapping"]["refs"]],

localrules: align_samples
rule align_samples:
    input:
        expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sample=SAMP2RUNLIB),

localrules: align_qualimap_samples
rule align_qualimap_samples:
    input:
        expand("data/alignments/qualimap/samples/{aligner}~{ref}~{sample}/",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMP2RUNLIB),

localrules: align_stats
rule align_stats:
    input:
        expand("data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.tsv",
               aligner=config["mapping"]["aligners"],
               ref=config["mapping"]["refs"],
               sample=SAMPLESETS["all_samples"]),
    output:
        expand("data/alnstats/everything_{type}.csv",
               type=["SN", "IS", "COV"])
    log: "data/log/bamstats/mergeallbamstats.log"
    shell:
        "python3 ./scripts/tidybamstat.py"
        "   -o data/alnstats/everything"  # prefix
        "   {input}"
        ">{log} 2>&1"


allsets = set(
    list(config["mapping"]["samplesets"]) +
    list(config["varcall"]["samplesets"]) +
    list(config["sample_sets"])
)
rule align_samplesets_all:
    # Currently there's no real need for this, as all variant calling uses the mega-bam
    # and ANGSD uses a list of filenames (which will be generated on the fly
    # for each ANGSD sample set)
    input:
        expand("data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=allsets),
        expand("data/alignments/bamlists/{aligner}~{ref}~all_samples.bamlist",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"],
               sampleset=allsets),

rule align:
   input:
        rules.align_samples.input,
        rules.align_stats.input,
        expand("data/alignments/sets/{aligner}~{ref}~all_samples.bam",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),
        expand("data/alignments/bamlists/{aligner}~{ref}~all_samples.bamlist",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),
        expand("data/alignments/sets/{aligner}~{ref}~all_samples.bam.bai",
               ref=config["mapping"]["refs"],
               aligner=config["mapping"]["aligners"]),

#######################################################################
#                           Variant Calling                           #
#######################################################################

# So I've rejigged the variant calling to use the megabam of all samples for
# each aligner/ref, as then we don't need 15 different massive subset bams.

rule freebayes:
    input:
        bam="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",  # use the megabam, see above
        bai="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam.bai",
        sset="data/samplelists/{sampleset}.txt",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/freebayes~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.log"
    benchmark:
        "data/log/freebayes/{aligner}~{ref}~{sampleset}/{region}.benchmark"
    priority: 1  # get them done earlier, normalisation is super quick
    params:
        theta=lambda wc: config["varcall"]["samplesets"][wc.sampleset].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    shell:
        "( freebayes"
        "   --theta {params.theta}"
        "   --samples {input.sset}"
        "   --ploidy 2"
        "   --use-best-n-alleles 4"
        "   --min-mapping-quality {params.minmq}"
        "   --min-base-quality {params.minbq}"
        "   --read-max-mismatch-fraction 0.15"
        "   --min-alternate-fraction 0.05"
        "   --min-alternate-count 3" # per sample
        "   --min-alternate-total 9" # across all samples
        "   --min-coverage 20" # across all samples
	"   --skip-coverage 100000"
        "   --prob-contamination 1e-3"
        "   --strict-vcf"
        "   --region '{wildcards.region}'"
        "   --fasta-reference {input.ref}"
        "   {input.bam}"
        " | bcftools view"
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"


rule mpileup:
    input:
        bam="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam",
        bai="data/alignments/sets/{aligner}~{ref}~{sampleset}.bam.bai",
        sset="data/samplelists/{sampleset}.txt",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="data/variants/raw_split/mpileup~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.log"
    benchmark:
        "data/log/mpileup/{aligner}~{ref}~{sampleset}/{region}.benchmark"
    params:
        theta=lambda wc: config["varcall"]["samplesets"][wc.sampleset].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    priority: 1  # get them done earlier, normalisation is super quick
    shell:
        "( bcftools mpileup"
        "   --redo-BAQ"
        "   --max-depth 100000" # the default per file max (250x) is insane, i.e. <1x for most sets. new limit of 20000x  equates to a max. of 20x across all samples.
        "   --min-MQ {params.minmq}"
        "   --min-BQ {params.minbq}"
        "   --fasta-ref {input.ref}"
        "   --samples-file {input.sset}"
        "   --annotate FORMAT/DP,FORMAT/AD,FORMAT/SP,INFO/AD" #output extra tags
        "   --region '{wildcards.region}'"
        "   --output-type u" # uncompressed bam
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{wildcards.region}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior {params.theta}"
        "   -O b" # compressed bam
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
        "   --multiallelics -snps"  # Split multi-alleics to filter each allele separately
        "   -O u  -o {output.bcf}"
        "   {input.bcf}"
        " ) >{log} 2>&1"

        #
        #  Old version with vt decompose_blocksub -- for some reason doesn't decompose block subsititions, and crashes on gadi
        #  "( bcftools norm"
        #  "   --fasta-ref {input.ref}"
        #  "   -O u"
        #  "   {input.bcf}"
        #  " | vt decompose_blocksub + -o -" # decompose MNP to multipe SNPs
        #  " | bcftools norm" # Split multi-alleics to filter each allele separately
        #  "   --fasta-ref {input.ref}"
        #  "   --do-not-normalize"
        #  "   --multiallelics -snps"
        #  "   -O u  -o {output.bcf}"
        #  " ) >{log} 2>&1"

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
        "   --multiallelics +snps" # re-join multi-alleic sites
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"

localrules: bcfmerge_fofn
rule bcfmerge_fofn:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.caller][wc.ref])),
    output:
        fofn=temp("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN"),
    run:
        with open(output[0], "w") as fh:
            for s in sorted(input):
                print(s, file=fh)

rule bcfmerge:
    input:
        bcf=lambda wc: expand("data/variants/filter_split/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}/{region}.bcf",
                              caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sampleset=wc.sampleset, filter=wc.filter,
                              region=sorted(VARCALL_REGIONS[wc.caller][wc.ref])),
        fofn="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf.INPUT_FOFN",
    output:
        bcf="data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.bcf",
    log:
        "data/log/mergebcf/{caller}~{aligner}~{ref}~{sampleset}_filtered~{filter}.log"
    threads: 12
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
    threads: 8
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
    shell:
        "bcftools index -f {input}"

rule bcfstats:
    input:
        "data/variants/{path}.bcf"
    output:
        "data/variants/{path}.bcf.stats"
    shell:
        "bcftools stats -s - -d 0,1000,1 --threads {threads} {input} >{output}"


def raw_variant_calls_input(wildcards):
    inputs = []
    for sampleset in config["varcall"]["samplesets"]:
        for caller in config["varcall"]["samplesets"][sampleset]["callers"]:
            for aligner in config["varcall"]["samplesets"][sampleset]["aligners"]:
                for ref in config["varcall"]["samplesets"][sampleset]["refs"]:
                    this_rawfiles = expand("data/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
                                           caller=caller, aligner=aligner, ref=ref, sampleset=sampleset, region=VARCALL_REGIONS[caller][ref])
                    inputs.extend(this_rawfiles)
    return inputs


rule raw_variant_calls:
    input: raw_variant_calls_input

localrules: filtered_variants
rule filtered_variants:
    input:
        [expand("data/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.{ext}",
               ext=["bcf", "bcf.csi", "vcf.gz", "vcf.gz.csi", "bcf.stats"],
               caller=config["varcall"]["samplesets"][sampleset]["callers"],
               aligner=config["varcall"]["samplesets"][sampleset]["aligners"],
               ref=config["varcall"]["samplesets"][sampleset]["refs"],
               filter=config["varcall"]["filters"],
               sampleset=sampleset)
        for sampleset in config["varcall"]["samplesets"]]



rule varcall:
    input:
        rules.filtered_variants.input,


## Haplotype Caller

#BSQR # persamp bam
#HC gvcf   #samp/region
#GenomicsDBImport  #samp/region
#GenotypeGVCFs  #region
#VariantRecalibrator, ApplyRecalibration  # per region
#Merge 




### ANGSD

# Angsd somewhat hacky verison of Rose's logic w/ hardcoded sites per step1 files

rule angsd_step2_maf_chrom:
    input:
        ref=lambda wc: config['refs'][wc.ref],
        sites=lambda wc: expand("rawdata/angsd-sites-rose/{chr}.sites_mm", chr=wc.chrom),
        bamlist="data/alignments/bamlists/{aligner}~{ref}~{sampleset}.bamlist",
    output:
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.arg",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.hwe.gz",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.mafs.gz",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.qs",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.saf.gz",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.saf.idx",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.saf.pos.gz",
        "data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.snpStat.gz",
    log:
        "data/log/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.log"
    shell:
        "angsd"
        "   -P 1"
        "   -bam {input.bamlist}"
        "   -out data/angsd/step2/maf_{wildcards.aligner}~{wildcards.ref}~{wildcards.sampleset}_{wildcards.chrom}"
        "   -ref {input.ref} -anc {input.ref}"
        "   -GL 2"
        "   -doMajorMinor 3"
        "   -skipTriallelic 1"
        "   -doHWE 1"
        "   -doMaf 1"
        "   -doSaf 1"
        "   -doCounts 1"
        "   -doQsDist 1"
        "   -doSnpStat 1"
        "   -C 50"
        "   -baq 1"
        "   -minMapQ 20"
        "   -minQ 20"
        "   -r {wildcards.chrom}"
        "   -sites {input.sites}"
        "   > {log} 2>&1"


rule all_angsd_step2_maf:
    input:
        expand("data/angsd/step2/maf_{aligner}~{ref}~{sampleset}_{chrom}.mafs.gz",
               aligner=config["angsd"]["aligners"],
               ref=config["angsd"]["refs"],
               sampleset=config["angsd"]["samplesets"],
               chrom=config["angsd"]["chroms"])


#######################################################################
#                              All rule                               #
#######################################################################


rule all:
    input:
        rules.denovo.input,
        rules.reads.input,
        rules.align.input,
        rules.varcall.input,
