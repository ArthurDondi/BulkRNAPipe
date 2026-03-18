### Rules for adapter and quality trimming with Trim Galore

if PAIRED:
    rule TrimGalore:
        input:
            r1 = get_raw_fastq_r1,
            r2 = get_raw_fastq_r2,
        output:
            r1 = "trim/{sample}/{sample}_R1_val_1.fq.gz",
            r2 = "trim/{sample}/{sample}_R2_val_2.fq.gz",
        params:
            outdir     = "trim/{sample}",
            min_length = config['TrimGalore']['min_length'],
            quality    = config['TrimGalore']['quality'],
            extra      = config['TrimGalore'].get('extra', ''),
        threads: 4
        resources:
            mem_mb        = 8000,
            runtime       = 120,
            cpus_per_task = 4,
        conda:
            "../envs/trimgalore.yaml"
        log:
            "logs/TrimGalore/{sample}.log"
        benchmark:
            "benchmark/TrimGalore/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            trim_galore \
                --paired \
                --quality {params.quality} \
                --length {params.min_length} \
                --cores {threads} \
                --gzip \
                --basename {wildcards.sample} \
                {params.extra} \
                --output_dir {params.outdir} \
                {input.r1} {input.r2}
            """
else:
    rule TrimGalore:
        input:
            r1 = get_raw_fastq_r1,
        output:
            r1 = "trim/{sample}/{sample}_trimmed.fq.gz",
        params:
            outdir     = "trim/{sample}",
            min_length = config['TrimGalore']['min_length'],
            quality    = config['TrimGalore']['quality'],
            extra      = config['TrimGalore'].get('extra', ''),
        threads: 4
        resources:
            mem_mb        = 8000,
            runtime       = 120,
            cpus_per_task = 4,
        conda:
            "../envs/trimgalore.yaml"
        log:
            "logs/TrimGalore/{sample}.log"
        benchmark:
            "benchmark/TrimGalore/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            trim_galore \
                --quality {params.quality} \
                --length {params.min_length} \
                --cores {threads} \
                --gzip \
                --basename {wildcards.sample} \
                {params.extra} \
                --output_dir {params.outdir} \
                {input.r1}
            """
