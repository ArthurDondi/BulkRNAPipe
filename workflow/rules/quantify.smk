### Rules for gene-level read quantification with featureCounts (Subread)

rule FeatureCounts:
    input:
        bams = get_bam_files,
        gtf  = GTF,
    output:
        counts   = "quantify/counts.txt",
        summary  = "quantify/counts.txt.summary",
    params:
        strandedness = STRANDEDNESS,
        min_mq       = config['featureCounts']['min_mq'],
        paired_flag  = "-p --countReadPairs" if PAIRED else "",
        extra        = config['featureCounts'].get('extra', ''),
    threads: 8
    resources:
        mem_mb        = 16000,
        runtime       = 120,
        cpus_per_task = 8,
    conda:
        "../envs/subread.yaml"
    log:
        "logs/FeatureCounts/counts.log"
    benchmark:
        "benchmark/FeatureCounts/counts.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p quantify
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -s {params.strandedness} \
            -Q {params.min_mq} \
            {params.paired_flag} \
            {params.extra} \
            {input.bams}
        """
