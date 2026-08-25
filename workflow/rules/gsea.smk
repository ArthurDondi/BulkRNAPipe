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


# ── Step 1b: fetch the CURRENT KEGG pathway definitions ──────────────────────
# Only built when 'KEGG_CURRENT' is among GSEA.collections. Needs outbound HTTPS
# to rest.kegg.jp; see workflow/scripts/generate_kegg_gmt.R for why the frozen
# MSigDB collection is not enough.

rule GenerateKeggGmt:
    """Build a GMT of current human KEGG pathways from the KEGG REST API."""
    output:
        gmt      = "resources/generated_gmts/kegg_current.gmt",
        pathways = "resources/generated_gmts/kegg_current_pathways.tsv",
    params:
        script = f"{workflow.basedir}/scripts/generate_kegg_gmt.R",
        outdir = "resources/generated_gmts",
    threads: 1
    retries: 2
    resources:
        mem_mb        = 4000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GenerateKeggGmt/kegg.log"
    benchmark:
        "benchmark/GenerateKeggGmt/kegg.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} --outdir {params.outdir}
        """


# ── Step 2: run fgsea per contrast and per collection ─────────────────────────
rule GSEA:
    """Run fgsea (fgseaMultilevel) for one contrast and one gene-set collection.

    Outputs a CSV of enrichment results and a dotplot PDF.
    """
    input:
        results = "deseq2/{contrast}/results.csv",
        hox_gmt = "resources/generated_gmts/hox.gmt",
        # Only the KEGG_CURRENT collection needs the fetched GMT, so the other
        # collections do not drag the network call into their dependency graph.
        kegg_gmt = lambda wildcards: (
            "resources/generated_gmts/kegg_current.gmt"
            if wildcards.collection == "KEGG_CURRENT" else []
        ),
    output:
        csv     = "gsea/{contrast}/{collection}_results.csv",
        dotplot = "gsea/{contrast}/{collection}_dotplot.pdf",
        audit   = "gsea/{contrast}/{collection}_pathway_audit.csv",
    params:
        script           = f"{workflow.basedir}/scripts/gsea.R",
        outdir           = "gsea/{contrast}",
        collection       = "{collection}",
        contrast_name    = "{contrast}",
        numerator        = lambda wildcards: next(
            c[1] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
        denominator      = lambda wildcards: next(
            c[2] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
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
            --kegg_gmt       "{input.kegg_gmt}" \
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


# ── Step 3: per-KEGG-category panels (re-plots step 2, no fgsea re-run) ───────

rule GSEAKeggCategories:
    """Slice the KEGG GSEA results by BRITE category and plot every pathway."""
    input:
        results    = f"gsea/{{contrast}}/{KEGG_SOURCE_SLUG}_results.csv",
        audit      = f"gsea/{{contrast}}/{KEGG_SOURCE_SLUG}_pathway_audit.csv",
        categories = KEGG_CATEGORIES_FILE,
    output:
        pdfs = expand("gsea/{{contrast}}/kegg/{slug}_dotplot.pdf",
                      slug=KEGG_CATEGORIES),
        csvs = expand("gsea/{{contrast}}/kegg/{slug}_results.csv",
                      slug=KEGG_CATEGORIES),
    params:
        script      = f"{workflow.basedir}/scripts/gsea_kegg.R",
        outdir      = "gsea/{contrast}/kegg",
        select      = ",".join(KEGG_CATEGORIES),
        padj_cutoff = lambda wildcards: config.get('GSEA', {}).get('kegg_padj_cutoff', 0.05),
        numerator   = lambda wildcards: next(
            c[1] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
        denominator = lambda wildcards: next(
            c[2] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 20,
        cpus_per_task = 1,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GSEAKeggCategories/{contrast}.log"
    benchmark:
        "benchmark/GSEAKeggCategories/{contrast}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results     {input.results} \
            --audit       {input.audit} \
            --categories  {input.categories} \
            --select      "{params.select}" \
            --outdir      {params.outdir} \
            --numerator   "{params.numerator}" \
            --denominator "{params.denominator}" \
            --padj_cutoff {params.padj_cutoff}
        """


# ── Step 4: GO keyword panels (same categories, matched on term names) ───────

rule GSEAGoCategories:
    """Group GO terms by keyword into the KEGG categories and plot each."""
    input:
        results    = f"gsea/{{contrast}}/{GO_SOURCE_SLUG}_results.csv",
        categories = GO_CATEGORIES_FILE,
    output:
        pdfs = expand("gsea/{{contrast}}/go/{slug}_dotplot.pdf",
                      slug=GO_CATEGORIES),
        csvs = expand("gsea/{{contrast}}/go/{slug}_results.csv",
                      slug=GO_CATEGORIES),
    params:
        script      = f"{workflow.basedir}/scripts/gsea_go_categories.R",
        outdir      = "gsea/{contrast}/go",
        select      = ",".join(GO_CATEGORIES),
        top_n       = GO_CAT_TOP_N,
        padj_cutoff = lambda wildcards: config.get('GSEA', {}).get('kegg_padj_cutoff', 0.05),
        numerator   = lambda wildcards: next(
            c[1] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
        denominator = lambda wildcards: next(
            c[2] for c in config['DESeq2']['contrasts'] if c[0] == wildcards.contrast
        ),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 20,
        cpus_per_task = 1,
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GSEAGoCategories/{contrast}.log"
    benchmark:
        "benchmark/GSEAGoCategories/{contrast}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results     {input.results} \
            --categories  {input.categories} \
            --select      "{params.select}" \
            --outdir      {params.outdir} \
            --top_n       {params.top_n} \
            --numerator   "{params.numerator}" \
            --denominator "{params.denominator}" \
            --padj_cutoff {params.padj_cutoff}
        """
