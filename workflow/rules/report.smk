### Rendered analysis report
###
### Deliberately NOT part of the default target: the report is a communication
### artefact, written and rewritten by hand, and rebuilding it on every pipeline
### run would overwrite edits and add a job to every invocation. Build it when
### you want it:
###
###   snakemake -s workflow/Snakefile --configfile config/<your>.yaml \
###       --use-conda --conda-frontend conda --cores 1 report/analysis_report.html
###
### The Rmd reads the pipeline's CSV outputs and rebuilds every figure, so it
### renders whatever exists and prints a note for whatever does not.

rule Report:
    """Render report/analysis_report.Rmd against this run's outputs."""
    input:
        # Only the sample sheet is required; every other section degrades to a
        # note if its input is missing, so the report can be built early.
        qc = "design_qc/design_qc_summary.csv",
    output:
        html = "report/analysis_report.html",
    params:
        rmd    = f"{workflow.basedir}/../report/analysis_report.Rmd",
        outdir = lambda wildcards: config['User']['output_dir'],
        padj   = config['DESeq2']['padj_threshold'],
        lfc    = config['DESeq2']['lfc_threshold'],
    threads: 1
    resources:
        mem_mb        = 8000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/report.yaml"
    log:
        "logs/Report/analysis_report.log"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p report
        Rscript -e 'rmarkdown::render(
            "{params.rmd}",
            params = list(outdir = "{params.outdir}",
                          padj   = {params.padj},
                          lfc    = {params.lfc}),
            output_file = normalizePath("{output.html}", mustWork = FALSE),
            intermediates_dir = tempdir(), knit_root_dir = tempdir())'
        """
