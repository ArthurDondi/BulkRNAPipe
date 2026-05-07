### Rules for GO enrichment analysis (clusterProfiler) per DESeq2 contrast

rule GOEnrichment:
    """Run clusterProfiler enrichGO() for one contrast, ontology, and direction.

    Input IDs are mapped (or treated as ENTREZID), then enrichGO and simplify
    are run. Outputs raw/simplified CSVs and dotplots.
    """
    input:
        results = "deseq2/{contrast}/results.csv",
    output:
        csv_raw         = "go/{contrast}/go_{ont}_{dir}_results.csv",
        csv_simplified  = "go/{contrast}/go_{ont}_{dir}_results_simplified.csv",
        dotplot_raw     = "go/{contrast}/go_{ont}_{dir}_dotplot.pdf",
        dotplot_simplified = "go/{contrast}/go_{ont}_{dir}_dotplot_simplified.pdf",
        unmapped_universe = "go/{contrast}/go_{ont}_{dir}_unmapped_universe.csv",
        unmapped_sig      = "go/{contrast}/go_{ont}_{dir}_unmapped_sig.csv",
    params:
        script       = f"{workflow.basedir}/scripts/go.R",
        outdir       = "go/{contrast}",
        ontology     = lambda wildcards: wildcards.ont.upper(),
        direction    = lambda wildcards: wildcards.dir.lower(),
        padj_cutoff  = config.get('GO', {}).get('padj_cutoff', 0.05),
        min_gs_size  = config.get('GO', {}).get('min_gs_size', 10),
        max_gs_size  = config.get('GO', {}).get('max_gs_size', 500),
        padj_thr     = config['DESeq2']['padj_threshold'],
        lfc_thr      = config['DESeq2']['lfc_threshold'],
        gene_id_type = config.get('GO', {}).get('gene_id_type', 'SYMBOL'),
    threads: 2
    resources:
        mem_mb        = 8000,
        runtime       = 60,
        cpus_per_task = 2,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GOEnrichment/{contrast}/go_{ont}_{dir}.log"
    benchmark:
        "benchmark/GOEnrichment/{contrast}/go_{ont}_{dir}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results    {input.results} \
            --outdir     {params.outdir} \
            --ontology   {params.ontology} \
            --direction  {params.direction} \
            --padj_cutoff {params.padj_cutoff} \
            --min_gs_size {params.min_gs_size} \
            --max_gs_size {params.max_gs_size} \
            --padj_thr   {params.padj_thr} \
            --lfc_thr    {params.lfc_thr} \
            --gene_id_type {params.gene_id_type}
        """
