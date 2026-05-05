### Rules for read quality control (FastQC + MultiQC)

# ── Per-sample FastQC on raw reads ──────────────────────────────────────────

if PAIRED:
    rule FastQC_raw:
        input:
            r1 = get_raw_fastq_r1,
            r2 = get_raw_fastq_r2,
        output:
            html_r1 = "QC/raw/fastqc/{sample}/{sample}_R1_fastqc.html",
            html_r2 = "QC/raw/fastqc/{sample}/{sample}_R2_fastqc.html",
            zip_r1  = "QC/raw/fastqc/{sample}/{sample}_R1_fastqc.zip",
            zip_r2  = "QC/raw/fastqc/{sample}/{sample}_R2_fastqc.zip",
        params:
            outdir = "QC/raw/fastqc/{sample}",
        threads: 2
        resources:
            mem_mb  = 4000,
            runtime = 60,
            cpus_per_task = 2,
        conda:
            "../envs/fastqc.yaml"
        log:
            "logs/FastQC_raw/{sample}.log"
        benchmark:
            "benchmark/FastQC_raw/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            fastqc --threads {threads} --outdir {params.outdir} {input.r1} {input.r2}
            b1=$(basename {input.r1}); b1="${{b1%.fastq.gz}}"; b1="${{b1%.fq.gz}}"
            b2=$(basename {input.r2}); b2="${{b2%.fastq.gz}}"; b2="${{b2%.fq.gz}}"
            mv {params.outdir}/${{b1}}_fastqc.html {output.html_r1}
            mv {params.outdir}/${{b1}}_fastqc.zip  {output.zip_r1}
            mv {params.outdir}/${{b2}}_fastqc.html {output.html_r2}
            mv {params.outdir}/${{b2}}_fastqc.zip  {output.zip_r2}
            """
else:
    rule FastQC_raw:
        input:
            r1 = get_raw_fastq_r1,
        output:
            html_r1 = "QC/raw/fastqc/{sample}/{sample}_R1_fastqc.html",
            zip_r1  = "QC/raw/fastqc/{sample}/{sample}_R1_fastqc.zip",
        params:
            outdir = "QC/raw/fastqc/{sample}",
        threads: 1
        resources:
            mem_mb  = 4000,
            runtime = 60,
            cpus_per_task = 1,
        conda:
            "../envs/fastqc.yaml"
        log:
            "logs/FastQC_raw/{sample}.log"
        benchmark:
            "benchmark/FastQC_raw/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            fastqc --threads {threads} --outdir {params.outdir} {input.r1}
            b1=$(basename {input.r1}); b1="${{b1%.fastq.gz}}"; b1="${{b1%.fq.gz}}"
            mv {params.outdir}/${{b1}}_fastqc.html {output.html_r1}
            mv {params.outdir}/${{b1}}_fastqc.zip  {output.zip_r1}
            """

# ── Aggregate raw FastQC with MultiQC ───────────────────────────────────────

rule MultiQC_raw:
    input:
        expand("QC/raw/fastqc/{sample}/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("QC/raw/fastqc/{sample}/{sample}_R2_fastqc.zip", sample=SAMPLES) if PAIRED else [],
    output:
        report = "QC/raw/multiqc/multiqc_report.html",
    params:
        indir  = "QC/raw/fastqc",
        outdir = "QC/raw/multiqc",
    resources:
        mem_mb  = 4000,
        runtime = 30,
        cpus_per_task = 1,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/MultiQC_raw/multiqc.log"
    benchmark:
        "benchmark/MultiQC_raw/multiqc.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        multiqc --force --outdir {params.outdir} {params.indir}
        """

# ── Per-sample FastQC on trimmed reads ──────────────────────────────────────

if PAIRED:
    rule FastQC_trimmed:
        input:
            r1 = "trim/{sample}/{sample}_val_1.fq.gz",
            r2 = "trim/{sample}/{sample}_val_2.fq.gz",
        output:
            html_r1 = "QC/trimmed/fastqc/{sample}/{sample}_val_1_fastqc.html",
            html_r2 = "QC/trimmed/fastqc/{sample}/{sample}_val_2_fastqc.html",
            zip_r1  = "QC/trimmed/fastqc/{sample}/{sample}_val_1_fastqc.zip",
            zip_r2  = "QC/trimmed/fastqc/{sample}/{sample}_val_2_fastqc.zip",
        params:
            outdir = "QC/trimmed/fastqc/{sample}",
        threads: 2
        resources:
            mem_mb  = 4000,
            runtime = 60,
            cpus_per_task = 2,
        conda:
            "../envs/fastqc.yaml"
        log:
            "logs/FastQC_trimmed/{sample}.log"
        benchmark:
            "benchmark/FastQC_trimmed/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            fastqc --threads {threads} --outdir {params.outdir} {input.r1} {input.r2}
            """
else:
    rule FastQC_trimmed:
        input:
            r1 = "trim/{sample}/{sample}_trimmed.fq.gz",
        output:
            html_r1 = "QC/trimmed/fastqc/{sample}/{sample}_trimmed_fastqc.html",
            zip_r1  = "QC/trimmed/fastqc/{sample}/{sample}_trimmed_fastqc.zip",
        params:
            outdir = "QC/trimmed/fastqc/{sample}",
        threads: 1
        resources:
            mem_mb  = 4000,
            runtime = 60,
            cpus_per_task = 1,
        conda:
            "../envs/fastqc.yaml"
        log:
            "logs/FastQC_trimmed/{sample}.log"
        benchmark:
            "benchmark/FastQC_trimmed/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            fastqc --threads {threads} --outdir {params.outdir} {input.r1}
            """

# ── Aggregate trimmed FastQC with MultiQC ───────────────────────────────────

rule MultiQC_trimmed:
    input:
        expand("QC/trimmed/fastqc/{sample}/{sample}_val_1_fastqc.zip", sample=SAMPLES) if PAIRED else
        expand("QC/trimmed/fastqc/{sample}/{sample}_trimmed_fastqc.zip", sample=SAMPLES),
        expand("trim/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES) if PAIRED else
        expand("trim/{sample}/{sample}_trimmed.fq.gz", sample=SAMPLES),
    output:
        report = "QC/trimmed/multiqc/multiqc_report.html",
    params:
        indirs = ["QC/trimmed/fastqc", "trim"],
        outdir = "QC/trimmed/multiqc",
    resources:
        mem_mb  = 4000,
        runtime = 30,
        cpus_per_task = 1,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/MultiQC_trimmed/multiqc.log"
    benchmark:
        "benchmark/MultiQC_trimmed/multiqc.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        multiqc --force --outdir {params.outdir} {params.indirs}
        """

# ── Aggregate alignment logs with MultiQC ────────────────────────────────────

rule MultiQC_align:
    input:
        expand("align/{sample}/{sample}.Log.final.out", sample=SAMPLES),
    output:
        report = "QC/align/multiqc/multiqc_report.html",
    params:
        indir  = "align",
        outdir = "QC/align/multiqc",
    resources:
        mem_mb  = 4000,
        runtime = 30,
        cpus_per_task = 1,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/MultiQC_align/multiqc.log"
    benchmark:
        "benchmark/MultiQC_align/multiqc.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        multiqc --force --outdir {params.outdir} {params.indir}
        """
