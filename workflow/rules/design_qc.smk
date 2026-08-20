### Library-level QC for confounder checking
###
### Answers "does something technical vary with group?" — library size,
### assignment rate, genes detected, reporter/transgene expression, and the
### sample-sample correlation and distance structure, all samples on one plot.

rule DesignQC:
    """Design-level QC across all samples together."""
    input:
        counts  = "quantify/counts.txt",
        summary = "quantify/counts.txt.summary",
    output:
        library_size = "design_qc/library_size.pdf",
        genes_det    = "design_qc/genes_detected.pdf",
        markers      = "design_qc/marker_expression.pdf",
        markers_csv  = "design_qc/marker_expression.csv",
        goi          = "design_qc/goi_expression.pdf",
        goi_csv      = "design_qc/goi_expression.csv",
        correlation  = "design_qc/sample_correlation.pdf",
        distance     = "design_qc/sample_distance.pdf",
        summary_csv  = "design_qc/design_qc_summary.csv",
    params:
        script  = f"{workflow.basedir}/scripts/design_qc.R",
        outdir  = "design_qc",
        samples = lambda wildcards: ",".join(
            f"{s}:{config['samples'][s]['condition']}" for s in SAMPLES
        ),
        marker_genes = ",".join(DESIGN_QC_GENES),
        goi_genes    = ",".join(DESIGN_QC_GOI),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/DesignQC/design_qc.log"
    benchmark:
        "benchmark/DesignQC/design_qc.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --counts       {input.counts} \
            --summary      {input.summary} \
            --outdir       {params.outdir} \
            --samples      {params.samples} \
            --marker_genes "{params.marker_genes}" \
            --goi_genes    "{params.goi_genes}"
        """
