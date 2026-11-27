### Rendered analysis reports
###
### Two reports, from report/<name>_report.Rmd:
###   cell_lines - the SK-N-SH ATRX-rescue experiment (this pipeline's outputs)
###   patients   - the neuroblastoma patient cohort: cohort and groups, ATRX
###                interactors (with the SK-N-SH lines next to the patients),
###                differential expression between the ATRX/MYCN groups
### Both source report/report_common.R for their setup and shared helpers.
###
### Deliberately NOT part of the default target: the reports are communication
### artefacts, written and rewritten by hand, and rebuilding them on every
### pipeline run would overwrite edits and add jobs to every invocation. Build
### them when you want them:
###
###   snakemake -s workflow/Snakefile --configfile config/<your>.yaml \
###       --use-conda --conda-frontend conda --cores 1 \
###       report/cell_lines_report.html report/patients_report.html
###
### Each Rmd reads the pipeline's CSV outputs and rebuilds every figure, so it
### renders whatever exists and prints a note for whatever does not.

REPORT_NAMES = ["cell_lines", "patients"]
# Resolved at parse time: inside a params lambda, workflow.basedir is not
# guaranteed to still point at workflow/.
REPORT_SRC = f"{workflow.basedir}/../report"

rule Report:
    """Render report/{report}_report.Rmd against this run's outputs."""
    input:
        # Only the sample sheet is required; every other section degrades to a
        # note if its input is missing, so the report can be built early.
        qc     = "design_qc/design_qc_summary.csv",
        # The shared setup, so editing it re-renders both reports.
        common = f"{REPORT_SRC}/report_common.R",
    output:
        html = "report/{report}_report.html",
    wildcard_constraints:
        report = "|".join(REPORT_NAMES),
    params:
        rmd    = lambda wildcards: f"{REPORT_SRC}/{wildcards.report}_report.Rmd",
        outdir = lambda wildcards: config['User']['output_dir'],
        padj   = config['DESeq2']['padj_threshold'],
        lfc    = config['DESeq2']['lfc_threshold'],
        # rmarkdown::render resolves a RELATIVE output_file against the input
        # .Rmd's own directory, not the working directory - which put the report
        # inside the repo's report/ folder and failed with
        # "The directory 'report' does not exist". Passing output_dir and a bare
        # filename is the documented way to control the destination.
        # The per-plot PDFs and per-table CSVs go next to the HTML:
        # report/figures/<report>/<section>/ and report/tables/<report>/<section>/.
        report_dir = "report",
        html_name  = lambda wildcards: f"{wildcards.report}_report.html",
    threads: 1
    resources:
        # patients: 64 GB, 8 GB ran out in the patient section (MSigDB through
        # msigdbr, GSVA on GO:BP per analysis, 24 comparisons kept in memory) on
        # the ~250-sample cohort. The cell-line report rendered within 8 GB.
        mem_mb        = lambda wildcards: 64000 if wildcards.report == "patients" else 8000,
        # patients: 4 limma models and 72 fgsea analyses (24 comparisons x 3
        # collections; GO:BP is the slow one). 120 min is also the tinyq limit
        # (profile/slurm/config.yaml).
        runtime       = lambda wildcards: 120 if wildcards.report == "patients" else 60,
        cpus_per_task = 1,
    conda:
        "../envs/report.yaml"
    log:
        "logs/Report/{report}_report.log"
    shell:
        r"""
        exec > {log} 2>&1
        # The report uses non-ASCII characters (en dashes, +/-, superscripts). Under
        # the C locale R cannot translate them and floods the log with encoding
        # warnings, so force UTF-8 for the render.
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        mkdir -p {params.report_dir}
        # report_dir is resolved before render(): render() evaluates params
        # lazily, after it has moved to knit_root_dir, where a relative path
        # would land the exported figures/ and tables/ in a temporary folder.
        Rscript -e 'report_dir <- normalizePath("{params.report_dir}")
          rmarkdown::render(
            "{params.rmd}",
            params = list(outdir = "{params.outdir}",
                          padj   = {params.padj},
                          lfc    = {params.lfc},
                          export_dir = report_dir),
            output_dir  = report_dir,
            output_file = "{params.html_name}",
            intermediates_dir = tempdir(), knit_root_dir = tempdir())'
        """
