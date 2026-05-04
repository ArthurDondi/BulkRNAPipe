### Rules for differential expression analysis with DESeq2

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
        sample_conditions = lambda _: ",".join(
            f"{s}:{config['samples'][s]['condition']}" for s in SAMPLES
        ),
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
