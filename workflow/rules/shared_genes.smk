### Genes shared by two DESeq2 contrasts with the same denominator
### Outputs: shared_genes/{name}/genes.csv and ora_{up,down}_{collection}.csv

rule SharedGenes:
    """Genes significant in the same direction in contrast_A and contrast_B,
    and an over-representation test on the genes shared in each direction.

    The report draws the Venn diagrams of these sets and the enrichment of
    their centres (section "Shared by ATRX_FL and ATRX_IFF").
    """
    input:
        results_a = lambda wildcards: "deseq2/{}/results.csv".format(
            get_shared_cfg(wildcards.name)['contrast_A']),
        results_b = lambda wildcards: "deseq2/{}/results.csv".format(
            get_shared_cfg(wildcards.name)['contrast_B']),
        kegg_gmt = ("resources/generated_gmts/kegg_current.gmt"
                    if "KEGG_CURRENT" in SHARED_COLLECTIONS else []),
    output:
        genes = "shared_genes/{name}/genes.csv",
        ora   = expand("shared_genes/{{name}}/ora_{dir}_{collection}.csv",
                       dir=["up", "down"],
                       collection=[_collection_to_slug(c) for c in SHARED_COLLECTIONS]),
    params:
        script      = f"{workflow.basedir}/scripts/shared_genes.R",
        outdir      = "shared_genes/{name}",
        padj_thr    = config['DESeq2']['padj_threshold'],
        lfc_thr     = config['DESeq2']['lfc_threshold'],
        exclude     = ",".join(_shared_cfg.get('exclude_genes') or []),
        collections = ",".join(SHARED_COLLECTIONS),
        min_size    = _shared_cfg.get('min_size', 10),
        max_size    = _shared_cfg.get('max_size', 500),
    threads: 1
    resources:
        mem_mb        = 8000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/SharedGenes/{name}.log"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results_a     {input.results_a} \
            --results_b     {input.results_b} \
            --padj_thr      {params.padj_thr} \
            --lfc_thr       {params.lfc_thr} \
            --exclude_genes "{params.exclude}" \
            --collections   "{params.collections}" \
            --kegg_gmt      "{input.kegg_gmt}" \
            --min_size      {params.min_size} \
            --max_size      {params.max_size} \
            --outdir        {params.outdir}
        """
