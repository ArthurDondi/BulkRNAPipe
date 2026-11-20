### Rules for differential expression analysis with DESeq2

rule DESeq2:
    input:
        counts = DOWNSTREAM_COUNTS,
    output:
        results       = "deseq2/{contrast}/results.csv",
        norm_counts   = "deseq2/{contrast}/normalized_counts.csv",
        volcano       = "deseq2/{contrast}/volcano.pdf",
        ma_plot       = "deseq2/{contrast}/ma_plot.pdf",
    params:
        script         = f"{workflow.basedir}/scripts/deseq2.R",
        outdir         = "deseq2/{contrast}",
        contrast       = lambda wildcards: next(
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
            --counts         {input.counts} \
            --outdir         {params.outdir} \
            --contrast       "{params.contrast[1]} {params.contrast[2]}" \
            --samples        {params.sample_conditions} \
            --padj           {params.padj_threshold} \
            --lfc            {params.lfc_threshold}
        """


rule DESeq2Interaction:
    """Difference of differences: (A1 - A2) - (B1 - B2).

    Used to compare across a boundary whose effect cannot be estimated as a
    covariate.  The nuisance effect cancels in the subtraction instead of being
    modelled, at the cost of assuming it is the same size in both pairs.
    """
    input:
        counts = DOWNSTREAM_COUNTS,
    output:
        results     = "deseq2_interaction/{interaction}/results.csv",
        components  = "deseq2_interaction/{interaction}/components.csv",
        norm_counts = "deseq2_interaction/{interaction}/normalized_counts.csv",
        volcano     = "deseq2_interaction/{interaction}/volcano.pdf",
        ma_plot     = "deseq2_interaction/{interaction}/ma_plot.pdf",
    params:
        script         = f"{workflow.basedir}/scripts/deseq2_interaction.R",
        outdir         = "deseq2_interaction/{interaction}",
        group_a        = lambda wildcards: " ".join(
            get_interaction_cfg(wildcards.interaction)['group_A']
        ),
        group_b        = lambda wildcards: " ".join(
            get_interaction_cfg(wildcards.interaction)['group_B']
        ),
        padj_threshold = lambda wildcards: float(
            get_interaction_cfg(wildcards.interaction).get(
                'padj_threshold', config['DESeq2']['padj_threshold'])
        ),
        lfc_threshold  = lambda wildcards: float(
            get_interaction_cfg(wildcards.interaction).get(
                'lfc_threshold', config['DESeq2']['lfc_threshold'])
        ),
        # Raw condition labels: an interaction names its four groups directly,
        # so combine_conditions remapping must not be applied here.
        sample_conditions = lambda wildcards: ",".join(
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
        "logs/DESeq2Interaction/{interaction}.log"
    benchmark:
        "benchmark/DESeq2Interaction/{interaction}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --counts  {input.counts} \
            --outdir  {params.outdir} \
            --samples {params.sample_conditions} \
            --name    "{wildcards.interaction}" \
            --group_a "{params.group_a}" \
            --group_b "{params.group_b}" \
            --padj    {params.padj_threshold} \
            --lfc     {params.lfc_threshold}
        """


rule DESeq2Proteomic:
    input:
        results = "deseq2/{contrast}/results.csv",
    output:
        volcano_proteomic = "deseq2/{contrast}/volcano_proteomic.pdf",
    params:
        script         = f"{workflow.basedir}/scripts/deseq2_proteomic.R",
        outdir         = "deseq2/{contrast}",
        contrast_name  = "{contrast}",
        contrast       = lambda wildcards: next(
            c for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
        padj_threshold = config['DESeq2']['padj_threshold'],
        lfc_threshold  = config['DESeq2']['lfc_threshold'],
        proteomics_xlsx = lambda wildcards: str(config.get('Proteomics', {}).get('limma_xlsx', "")),
        proteomics_sheet = lambda wildcards: str(config.get('Proteomics', {}).get('sheet', "limma result")),
        proteomics_gene_column = lambda wildcards: str(config.get('Proteomics', {}).get('gene_column', "")),
        proteomics_comparison_column = lambda wildcards: str(config.get('Proteomics', {}).get('comparison_column', "")),
        proteomics_fdr_column = lambda wildcards: str(config.get('Proteomics', {}).get('fdr_column', "")),
        proteomics_logfc_column = lambda wildcards: str(config.get('Proteomics', {}).get('logfc_column', "")),
        proteomics_fdr_threshold = lambda wildcards: float(config.get('Proteomics', {}).get('fdr_threshold', 0.05)),
        proteomics_lfc_threshold = lambda wildcards: float(config.get('Proteomics', {}).get('lfc_threshold', 1.0)),
        proteomics_comparison = lambda wildcards: str(get_proteomics_comparison(wildcards.contrast)),
    threads: 2
    resources:
        mem_mb        = 8000,
        runtime       = 60,
        cpus_per_task = 2,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/DESeq2_proteomic/{contrast}.log"
    benchmark:
        "benchmark/DESeq2_proteomic/{contrast}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results         {input.results} \
            --outdir          {params.outdir} \
            --contrast_name   "{params.contrast_name}" \
            --contrast        "{params.contrast[1]} {params.contrast[2]}" \
            --padj            {params.padj_threshold} \
            --lfc             {params.lfc_threshold} \
            --proteomics_xlsx "{params.proteomics_xlsx}" \
            --proteomics_sheet "{params.proteomics_sheet}" \
            --proteomics_gene_column "{params.proteomics_gene_column}" \
            --proteomics_comparison_column "{params.proteomics_comparison_column}" \
            --proteomics_fdr_column "{params.proteomics_fdr_column}" \
            --proteomics_logfc_column "{params.proteomics_logfc_column}" \
            --proteomics_fdr_threshold {params.proteomics_fdr_threshold} \
            --proteomics_lfc_threshold {params.proteomics_lfc_threshold} \
            --proteomics_comparison "{params.proteomics_comparison}"
        """
