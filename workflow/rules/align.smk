### Rules for genome indexing and read alignment with STAR

# sjdbOverhang should be max(read_length) - 1 (STAR documentation recommendation)
SJDB_OVERHANG = READ_LENGTH - 1

# ── STAR genome index ────────────────────────────────────────────────────────

rule STARindex:
    input:
        fasta = GENOME,
        gtf   = GTF,
    output:
        directory(STAR_INDEX),
    params:
        extra = config['STAR'].get('index_extra', ''),
    threads: 16
    resources:
        mem_mb        = 64000,
        runtime       = 240,
        cpus_per_task = 16,
    conda:
        "../envs/star.yaml"
    log:
        "logs/STARindex/genome.log"
    benchmark:
        "benchmark/STARindex/genome.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {SJDB_OVERHANG} \
            {params.extra}
        """

# ── STAR alignment (2-pass mode via STARsolo default) ────────────────────────

if PAIRED:
    rule STARalign:
        input:
            index = STAR_INDEX,
            reads = get_star_input,
        output:
            bam       = "align/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            log_final = "align/{sample}/{sample}.Log.final.out",
        params:
            prefix  = "align/{sample}/{sample}.",
            extra   = config['STAR'].get('align_extra', ''),
        threads: 8
        resources:
            mem_mb        = 32000,
            runtime       = 180,
            cpus_per_task = 8,
        conda:
            "../envs/star.yaml"
        log:
            "logs/STARalign/{sample}.log"
        benchmark:
            "benchmark/STARalign/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p align/{wildcards.sample}
            STAR \
                --runThreadN {threads} \
                --genomeDir {input.index} \
                --readFilesIn {input.reads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI AS NM MD \
                --outFileNamePrefix {params.prefix} \
                --outReadsUnmapped Fastx \
                --twopassMode Basic \
                {params.extra}
            """
else:
    rule STARalign:
        input:
            index = STAR_INDEX,
            reads = get_star_input,
        output:
            bam       = "align/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            log_final = "align/{sample}/{sample}.Log.final.out",
        params:
            prefix = "align/{sample}/{sample}.",
            extra  = config['STAR'].get('align_extra', ''),
        threads: 8
        resources:
            mem_mb        = 32000,
            runtime       = 180,
            cpus_per_task = 8,
        conda:
            "../envs/star.yaml"
        log:
            "logs/STARalign/{sample}.log"
        benchmark:
            "benchmark/STARalign/{sample}.benchmark.txt"
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p align/{wildcards.sample}
            STAR \
                --runThreadN {threads} \
                --genomeDir {input.index} \
                --readFilesIn {input.reads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI AS NM MD \
                --outFileNamePrefix {params.prefix} \
                --outReadsUnmapped Fastx \
                --twopassMode Basic \
                {params.extra}
            """

# ── Index BAM files for downstream tools ─────────────────────────────────────

rule SamtoolsIndex:
    input:
        bam = "align/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai = "align/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    threads: 4
    resources:
        mem_mb        = 4000,
        runtime       = 60,
        cpus_per_task = 4,
    conda:
        "../envs/star.yaml"
    log:
        "logs/SamtoolsIndex/{sample}.log"
    benchmark:
        "benchmark/SamtoolsIndex/{sample}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        samtools index -@ {threads} {input.bam}
        """
