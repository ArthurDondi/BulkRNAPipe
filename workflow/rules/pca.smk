### Sample-level PCA
###
### Runs a grid of (sample set x gene view) so that the effect you are trying to
### see is not permanently hidden behind a larger one:
###   views       all_genes | no_vector_genes | no_vector_no_transgene
###   sample sets all + every entry under PCA.subsets
### Outputs live in pca/{sampleset}/ and are named after the view.

import re


rule PCA:
    """PCA for one (sample set x gene view) combination.

    Writes PC1/PC2, PC3/PC4, a scree plot, the PC coordinates, the variance
    table, the list of excluded features, and a supplementary
    limma::removeBatchEffect panel that is for visual inspection only.
    """
    input:
        counts = "quantify/counts.txt",
    output:
        pc12          = "pca/{sampleset}/{view}_pc12.pdf",
        pc34          = "pca/{sampleset}/{view}_pc34.pdf",
        scree         = "pca/{sampleset}/{view}_scree.pdf",
        batch_removed = "pca/{sampleset}/{view}_pc12_batch_removed.pdf",
        coords        = "pca/{sampleset}/{view}_coords.csv",
        variance      = "pca/{sampleset}/{view}_variance_explained.csv",
        excluded      = "pca/{sampleset}/{view}_excluded_genes.csv",
    params:
        script                = f"{workflow.basedir}/scripts/pca.R",
        outdir                = "pca/{sampleset}",
        # Raw conditions (no contrast-specific remapping) so every sample keeps
        # its original label.
        samples               = lambda wildcards: get_pca_samples(wildcards.sampleset),
        batch                 = lambda wildcards: get_pca_batch(wildcards.sampleset),
        exclude_genes         = ",".join(PCA_EXCLUDE_GENES),
        exclude_gene_patterns = ",".join(PCA_EXCLUDE_GENE_PATTERNS),
        exclude_contigs       = ",".join(PCA_EXCLUDE_CONTIGS),
        transgene_genes       = ",".join(PCA_TRANSGENE_GENES),
        ntop                  = lambda wildcards: config.get('PCA', {}).get('ntop', 500),
    wildcard_constraints:
        sampleset = "|".join(re.escape(s) for s in PCA_SAMPLE_SETS),
        view      = "|".join(re.escape(v) for v in PCA_VIEWS),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/PCA/{sampleset}_{view}.log"
    benchmark:
        "benchmark/PCA/{sampleset}_{view}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --counts    {input.counts} \
            --outdir    {params.outdir} \
            --samples   {params.samples} \
            --sampleset "{wildcards.sampleset}" \
            --view      "{wildcards.view}" \
            --ntop      {params.ntop} \
            --batch                 "{params.batch}" \
            --exclude_genes         "{params.exclude_genes}" \
            --exclude_gene_patterns "{params.exclude_gene_patterns}" \
            --exclude_contigs       "{params.exclude_contigs}" \
            --transgene_genes       "{params.transgene_genes}"
        """
