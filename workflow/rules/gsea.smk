### Rules for GSEA (fgsea) enrichment analysis per DESeq2 contrast

# ── Step 1: generate the HOX GMT once from the gene universe ─────────────────
rule GenerateHoxGmt:
    """Generate custom HOX gene sets (HOX_ALL and HOXB_ONLY) from the gene universe.

    Reads all gene symbols from the DESeq2 results of the first contrast and
    writes a GMT file to resources/generated_gmts/hox.gmt.
    """
    input:
        results = expand("deseq2/{contrast}/results.csv",
                         contrast=[CONTRASTS[0]]),
    output:
        gmt = "resources/generated_gmts/hox.gmt",
    params:
        script = f"{workflow.basedir}/scripts/generate_hox_gmt.R",
        outdir = "resources/generated_gmts",
    threads: 1
    resources:
        mem_mb        = 2000,
        runtime       = 10,
        cpus_per_task = 1,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GenerateHoxGmt/hox.log"
    benchmark:
        "benchmark/GenerateHoxGmt/hox.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results {input.results} \
            --outdir  {params.outdir}
        """


# ── Step 2: run fgsea per contrast and per collection ─────────────────────────
rule GSEA:
    """Run fgsea (fgseaMultilevel) for one contrast and one gene-set collection.

    Outputs a CSV of enrichment results and a dotplot PDF.
    """
    input:
        results = "deseq2/{contrast}/results.csv",
        hox_gmt = "resources/generated_gmts/hox.gmt",
    output:
        csv     = "gsea/{contrast}/{collection}_results.csv",
        dotplot = "gsea/{contrast}/{collection}_dotplot.pdf",
    params:
        script           = f"{workflow.basedir}/scripts/gsea.R",
        outdir           = "gsea/{contrast}",
        collection       = "{collection}",
        rank_metric      = config.get('GSEA', {}).get('rank_metric', 'stat'),
        min_size         = config.get('GSEA', {}).get('min_size', 15),
        max_size         = config.get('GSEA', {}).get('max_size', 500),
        nperm            = config.get('GSEA', {}).get('nperm', 1000),
        custom_gmt_files = lambda wildcards: ",".join(
            config.get('GSEA', {}).get('custom_gmt_files', [])
        ),
    threads: 2
    resources:
        mem_mb        = 8000,
        runtime       = 60,
        cpus_per_task = 2,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GSEA/{contrast}/{collection}.log"
    benchmark:
        "benchmark/GSEA/{contrast}/{collection}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results        {input.results} \
            --hox_gmt        {input.hox_gmt} \
            --outdir         {params.outdir} \
            --collection     {params.collection} \
            --rank_metric    {params.rank_metric} \
            --min_size       {params.min_size} \
            --max_size       {params.max_size} \
            --nperm          {params.nperm} \
            --custom_gmts    "{params.custom_gmt_files}"
        """
