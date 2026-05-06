### Rules for GO enrichment analysis (clusterProfiler) per DESeq2 contrast

rule GOEnrichment:
    """Run clusterProfiler enrichGO() for one contrast and one GO ontology.

    Gene symbols are mapped to Entrez IDs via org.Hs.eg.db, then enrichGO is
    run.  Outputs a CSV of enrichment results and a dotplot PDF.
    """
    input:
        results = "deseq2/{contrast}/results.csv",
    output:
        csv     = "go/{contrast}/go_{ont}_results.csv",
        dotplot = "go/{contrast}/go_{ont}_dotplot.pdf",
    params:
        script       = f"{workflow.basedir}/scripts/go.R",
        outdir       = "go/{contrast}",
        ontology     = lambda wildcards: wildcards.ont.upper(),
        padj_cutoff  = config.get('GO', {}).get('padj_cutoff', 0.05),
        min_gs_size  = config.get('GO', {}).get('min_gs_size', 10),
        max_gs_size  = config.get('GO', {}).get('max_gs_size', 500),
        padj_thr     = config['DESeq2']['padj_threshold'],
        lfc_thr      = config['DESeq2']['lfc_threshold'],
    threads: 2
    resources:
        mem_mb        = 8000,
        runtime       = 60,
        cpus_per_task = 2,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GOEnrichment/{contrast}/go_{ont}.log"
    benchmark:
        "benchmark/GOEnrichment/{contrast}/go_{ont}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results    {input.results} \
            --outdir     {params.outdir} \
            --ontology   {params.ontology} \
            --padj_cutoff {params.padj_cutoff} \
            --min_gs_size {params.min_gs_size} \
            --max_gs_size {params.max_gs_size} \
            --padj_thr   {params.padj_thr} \
            --lfc_thr    {params.lfc_thr}
        """
