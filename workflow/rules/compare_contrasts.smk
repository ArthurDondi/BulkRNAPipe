### Rules for comparative GSEA analysis between pairs of DESeq2 contrasts
### Outputs: ΔNES tables and residual-rank fgsea results

rule CompareContrasts:
    """Compute ΔNES and/or residual-rank GSEA for a configured contrast pair.

    For each ContrastComparisons entry this rule:
      - reads DESeq2 results for contrast_A and contrast_B
      - computes ΔNES (NES_A − NES_B) per pathway across all configured
        gene-set collections
      - optionally runs fgsea on residual ranks (rank_A − rank_B)
    Outputs are written to gsea_compare/{name}/.
    """
    input:
        results_a = lambda wildcards: (
            "deseq2/{contrast_a}/results.csv".format(
                contrast_a=get_comparison_cfg(wildcards.comparison)['contrast_A']
            )
        ),
        results_b = lambda wildcards: (
            "deseq2/{contrast_b}/results.csv".format(
                contrast_b=get_comparison_cfg(wildcards.comparison)['contrast_B']
            )
        ),
        hox_gmt = "resources/generated_gmts/hox.gmt",
        # Require that per-contrast GSEA is already complete for contrast_A
        gsea_a = lambda wildcards: expand(
            "gsea/{contrast}/{collection}_results.csv",
            contrast=get_comparison_cfg(wildcards.comparison)['contrast_A'],
            collection=GSEA_COLLECTIONS,
        ),
        # and for contrast_B
        gsea_b = lambda wildcards: expand(
            "gsea/{contrast}/{collection}_results.csv",
            contrast=get_comparison_cfg(wildcards.comparison)['contrast_B'],
            collection=GSEA_COLLECTIONS,
        ),
    output:
        summary        = "gsea_compare/{comparison}/delta_nes_summary.csv",
        residual_dir   = directory("gsea_compare/{comparison}/residual_rank"),
    params:
        script           = f"{workflow.basedir}/scripts/gsea_compare.R",
        outdir           = "gsea_compare/{comparison}",
        contrast_a       = lambda wildcards: get_comparison_cfg(wildcards.comparison)['contrast_A'],
        contrast_b       = lambda wildcards: get_comparison_cfg(wildcards.comparison)['contrast_B'],
        do_delta_nes     = lambda wildcards: str(get_comparison_cfg(wildcards.comparison).get('delta_nes', True)).upper(),
        do_residual_rank = lambda wildcards: str(get_comparison_cfg(wildcards.comparison).get('residual_rank', True)).upper(),
        gsea_dir_a       = lambda wildcards: "gsea/{c}".format(
            c=get_comparison_cfg(wildcards.comparison)['contrast_A']
        ),
        gsea_dir_b       = lambda wildcards: "gsea/{c}".format(
            c=get_comparison_cfg(wildcards.comparison)['contrast_B']
        ),
        rank_metric      = config.get('GSEA', {}).get('rank_metric', 'stat'),
        min_size         = config.get('GSEA', {}).get('min_size', 15),
        max_size         = config.get('GSEA', {}).get('max_size', 500),
        nperm            = config.get('GSEA', {}).get('nperm', 1000),
        collections      = lambda wildcards: ",".join(GSEA_COLLECTIONS),
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
        "logs/CompareContrasts/{comparison}.log"
    benchmark:
        "benchmark/CompareContrasts/{comparison}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results_a      {input.results_a} \
            --results_b      {input.results_b} \
            --hox_gmt        {input.hox_gmt} \
            --gsea_dir_a     {params.gsea_dir_a} \
            --gsea_dir_b     {params.gsea_dir_b} \
            --outdir         {params.outdir} \
            --contrast_a     "{params.contrast_a}" \
            --contrast_b     "{params.contrast_b}" \
            --collections    "{params.collections}" \
            --rank_metric    {params.rank_metric} \
            --min_size       {params.min_size} \
            --max_size       {params.max_size} \
            --nperm          {params.nperm} \
            --delta_nes      {params.do_delta_nes} \
            --residual_rank  {params.do_residual_rank} \
            --custom_gmts    "{params.custom_gmt_files}"
        """
