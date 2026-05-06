### Rules for differential expression analysis with DESeq2

rule PCA:
    """PCA plot (PC1 vs PC2, % variance explained) across all samples."""
    input:
        counts = "quantify/counts.txt",
    output:
        pca_plot = "deseq2/pca.pdf",
    params:
        script  = f"{workflow.basedir}/scripts/pca.R",
        outdir  = "deseq2",
        # Use raw conditions (no contrast-specific remapping) so all samples
        # are represented with their original condition labels.
        samples = lambda wildcards: ",".join(
            f"{s}:{config['samples'][s]['condition']}" for s in SAMPLES
        ),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/PCA/pca.log"
    benchmark:
        "benchmark/PCA/pca.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --counts   {input.counts} \
            --outdir   {params.outdir} \
            --samples  {params.samples}
        """


rule DESeq2:
    input:
        counts = "quantify/counts.txt",
    output:
        results = "deseq2/{contrast}/results.csv",
        norm_counts = "deseq2/{contrast}/normalized_counts.csv",
        volcano     = "deseq2/{contrast}/volcano.pdf",
        ma_plot     = "deseq2/{contrast}/ma_plot.pdf",
    params:
        script       = f"{workflow.basedir}/scripts/deseq2.R",
        outdir       = "deseq2/{contrast}",
        contrast     = lambda wildcards: next(
            c for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
        padj_threshold = config['DESeq2']['padj_threshold'],
        lfc_threshold  = config['DESeq2']['lfc_threshold'],
        # Inline the sample → condition mapping as a compact string
        # Format: "sample1:condition1,sample2:condition2,..."
        # get_contrast_effective_condition() applies only the combine_conditions
        # groups that are referenced by this specific contrast, allowing the
        # same source condition to participate in different groups across contrasts.
        sample_conditions = lambda wildcards: (
            lambda ce: ",".join(
                f"{s}:{get_contrast_effective_condition(s, ce)}" for s in SAMPLES
            )
        )(next(c for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast)),
    threads: 2
    resources:
        mem_mb        = 8000,
        runtime       = 60,
        cpus_per_task = 2,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/DESeq2/{contrast}.log"
    benchmark:
        "benchmark/DESeq2/{contrast}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --counts       {input.counts} \
            --outdir       {params.outdir} \
            --contrast     {params.contrast[1]} {params.contrast[2]} \
            --samples      {params.sample_conditions} \
            --padj         {params.padj_threshold} \
            --lfc          {params.lfc_threshold}
        """
