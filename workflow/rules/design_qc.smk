### Library-level QC for confounder checking
###
### Answers "does something technical vary with group?" — library size,
### assignment rate, genes detected, reporter/transgene expression, and the
### sample-sample correlation and distance structure, all samples on one plot.

# ─── Pooled and tag-only transgene counts (no re-alignment) ──────────────────
# Both rules re-count the EXISTING BAMs with -Q 0 -M --fraction, which is the
# whole point: the main matrix drops every NH >= 2 fragment, and for a transgene
# that duplicates an endogenous gene that is most of them. Weighting each
# alignment by 1/NH and pooling the copies into one meta-feature puts each
# fragment back exactly once.
#
# Counts come out fractional. That is fine for a QC panel and NOT valid DESeq2
# input - these features never enter the main count matrix.

rule DesignQCTotalCounts:
    """Pool several gene_ids into one meta-feature and count ambiguous fragments once."""
    input:
        bams = get_bam_files,
        gtf  = GTF,
    output:
        gtf    = "design_qc/total_feature.gtf",
        counts = "design_qc/total_feature_counts.txt",
    params:
        genes        = ",".join(DESIGN_QC_TOTAL_GENES),
        name         = DESIGN_QC_TOTAL_NAME,
        strandedness = STRANDEDNESS,
        paired_flag  = "-p --countReadPairs" if PAIRED else "",
    threads: 8
    resources:
        mem_mb        = 16000,
        runtime       = 120,
        cpus_per_task = 8,
    conda:
        "../envs/subread.yaml"
    log:
        "logs/DesignQC/total_feature.log"
    benchmark:
        "benchmark/DesignQC/total_feature.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p design_qc
        # Keep only the listed genes' exon lines and give them all one gene_id.
        # Restricting to exon keeps the semantics identical to the main matrix
        # (intronic fragments stay unassigned); pooling the ids is what lets
        # --fraction reassemble a fragment split across copies.
        awk -v want="{params.genes}" -v newid="{params.name}" '
            BEGIN {{ n = split(want, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1 }}
            /^#/  {{ next }}
            $3 != "exon" {{ next }}
            {{
                if (match($0, /gene_id "[^"]*"/)) {{
                    g = substr($0, RSTART + 9, RLENGTH - 10)
                    if (g in keep) {{
                        sub(/gene_id "[^"]*"/, "gene_id \"" newid "\"")
                        print
                    }}
                }}
            }}' {input.gtf} > {output.gtf}

        if [ ! -s {output.gtf} ]; then
            echo "ERROR: no exon features matched DesignQC.total_feature.genes ({params.genes})." >&2
            echo "Check the gene_id values in {input.gtf}." >&2
            exit 1
        fi

        featureCounts \
            -T {threads} \
            -a {output.gtf} \
            -o {output.counts} \
            -s {params.strandedness} \
            -Q 0 -M --fraction \
            {params.paired_flag} \
            {input.bams}
        """

rule DesignQCTagCounts:
    """Count fragments over the vector tag window, pooled across constructs."""
    input:
        bams = get_bam_files,
    output:
        saf    = "design_qc/tag.saf",
        counts = "design_qc/tag_counts.txt",
    params:
        spec         = DESIGN_QC_TAG_SPEC,
        min_overlap  = DESIGN_QC_TAG_MIN_OV,
        strandedness = STRANDEDNESS,
        paired_flag  = "-p --countReadPairs" if PAIRED else "",
    threads: 8
    resources:
        mem_mb        = 16000,
        runtime       = 120,
        cpus_per_task = 8,
    conda:
        "../envs/subread.yaml"
    log:
        "logs/DesignQC/tag.log"
    benchmark:
        "benchmark/DesignQC/tag.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p design_qc
        awk -v spec='{params.spec}' 'BEGIN {{
            print "GeneID\tChr\tStart\tEnd\tStrand"
            n = split(spec, rows, ";")
            for (i = 1; i <= n; i++) {{
                split(rows[i], f, ",")
                print f[1] "\t" f[2] "\t" f[3] "\t" f[4] "\t" f[5]
            }}
        }}' > {output.saf}

        # --minOverlap keeps fragments that merely graze the window out of the
        # count; without it one shared base would be enough to call a fragment
        # tag-derived.
        featureCounts \
            -T {threads} \
            -a {output.saf} -F SAF \
            -o {output.counts} \
            -s {params.strandedness} \
            -Q 0 -M --fraction \
            --minOverlap {params.min_overlap} \
            {params.paired_flag} \
            {input.bams}
        """

rule DesignQCSubstituteTotal:
    """Adjust the count matrix downstream of design_qc: swap total_feature's
    pooled counts in for one gene and/or drop other genes outright.

    quantify/counts.txt is never touched; this only runs when
    DesignQC.total_feature.substitute_as and/or
    DesignQC.downstream_exclude_genes is set (DOWNSTREAM_COUNTS then points at
    the output here instead of the raw matrix).
    """
    input:
        counts = "quantify/counts.txt",
        total  = "design_qc/total_feature_counts.txt" if DESIGN_QC_SUBSTITUTE_AS else [],
    output:
        counts = "design_qc/counts_substituted.txt",
    params:
        script           = f"{workflow.basedir}/scripts/substitute_total_feature.R",
        total_name       = DESIGN_QC_TOTAL_NAME,
        target_gene      = DESIGN_QC_SUBSTITUTE_AS,
        drop_genes       = ",".join(DESIGN_QC_TOTAL_GENES),
        extra_drop_genes = ",".join(DESIGN_QC_DOWNSTREAM_EXCLUDE),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 15,
        cpus_per_task = 1,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/DesignQC/substitute_total.log"
    benchmark:
        "benchmark/DesignQC/substitute_total.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p design_qc
        Rscript {params.script} \
            --counts           {input.counts} \
            --total            "{input.total}" \
            --total_name       "{params.total_name}" \
            --target_gene      "{params.target_gene}" \
            --drop_genes       "{params.drop_genes}" \
            --extra_drop_genes "{params.extra_drop_genes}" \
            --output           {output.counts}
        """


rule DesignQC:
    """Design-level QC across all samples together."""
    input:
        counts  = "quantify/counts.txt",
        summary = "quantify/counts.txt.summary",
        total   = "design_qc/total_feature_counts.txt" if DESIGN_QC_TOTAL_GENES else [],
        tag     = "design_qc/tag_counts.txt" if DESIGN_QC_TAG_REGIONS else [],
    output:
        library_size = "design_qc/library_size.pdf",
        genes_det    = "design_qc/genes_detected.pdf",
        markers      = "design_qc/marker_expression.pdf",
        markers_lin  = "design_qc/marker_expression_linear.pdf",
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
        # Empty string when the feature is not configured; the script treats an
        # empty path as "not requested" and leaves the panel unchanged.
        total_counts = "design_qc/total_feature_counts.txt" if DESIGN_QC_TOTAL_GENES else "",
        total_name   = DESIGN_QC_TOTAL_NAME,
        tag_counts   = "design_qc/tag_counts.txt" if DESIGN_QC_TAG_REGIONS else "",
        tag_name     = DESIGN_QC_TAG_NAME,
        gene_labels  = DESIGN_QC_LABELS,
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
            --goi_genes    "{params.goi_genes}" \
            --total_counts "{params.total_counts}" \
            --total_name   "{params.total_name}" \
            --tag_counts   "{params.tag_counts}" \
            --tag_name     "{params.tag_name}" \
            --gene_labels  "{params.gene_labels}"
        """
