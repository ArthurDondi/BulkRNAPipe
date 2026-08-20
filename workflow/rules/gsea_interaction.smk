### GSEA on the interaction (difference-of-differences) contrasts
###
### Separate from the GSEA rule so that adding it cannot disturb the per-contrast
### enrichment already on disk: different rule, different output tree, and the
### existing rule's inputs and params are untouched.
###
### Why run GSEA here at all.  A plain contrast that crosses the clone or
### transduction boundary is dominated by whatever separates those groups - in
### practice the proliferation/inflammation axis, which is the largest source of
### variation in cultured cells and which the Hallmark collection reports through
### many overlapping sets.  The interaction cancels both nuisance effects by
### subtraction, so ranking on it asks the question the plain contrasts cannot:
### which pathways did the intervention fail to restore, or overshoot?

rule GSEAInteraction:
    """fgsea on one interaction contrast and one gene-set collection."""
    input:
        results = "deseq2_interaction/{interaction}/results.csv",
        hox_gmt = "resources/generated_gmts/hox.gmt",
    output:
        csv     = "gsea_interaction/{interaction}/{collection}_results.csv",
        dotplot = "gsea_interaction/{interaction}/{collection}_dotplot.pdf",
    params:
        script        = f"{workflow.basedir}/scripts/gsea.R",
        outdir        = "gsea_interaction/{interaction}",
        collection    = "{collection}",
        contrast_name = "{interaction}",
        # gsea.R labels the two halves of the dotplot with these strings.  For an
        # interaction there is no single numerator condition, so both sides are
        # spelled out as the comparison of the two effects:
        #   NES > 0  the A pair moved the genes further up than the B pair
        #   NES < 0  the A pair moved them less far (a shortfall)
        numerator = lambda wildcards: (
            lambda it: f"{it['group_A'][0]}/{it['group_A'][1]} > {it['group_B'][0]}/{it['group_B'][1]}"
        )(get_interaction_cfg(wildcards.interaction)),
        denominator = lambda wildcards: (
            lambda it: f"{it['group_A'][0]}/{it['group_A'][1]} < {it['group_B'][0]}/{it['group_B'][1]}"
        )(get_interaction_cfg(wildcards.interaction)),
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
        "logs/GSEAInteraction/{interaction}/{collection}.log"
    benchmark:
        "benchmark/GSEAInteraction/{interaction}/{collection}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results        {input.results} \
            --hox_gmt        {input.hox_gmt} \
            --outdir         {params.outdir} \
            --collection     {params.collection} \
            --contrast_name  "{params.contrast_name}" \
            --numerator      "{params.numerator}" \
            --denominator    "{params.denominator}" \
            --rank_metric    {params.rank_metric} \
            --min_size       {params.min_size} \
            --max_size       {params.max_size} \
            --nperm          {params.nperm} \
            --custom_gmts    "{params.custom_gmt_files}"
        """
