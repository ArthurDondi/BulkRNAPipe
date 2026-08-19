### Signature reversal — scoring a rescue without crossing a confounded boundary
###
### Compares the log2 fold-changes of two contrasts that each sit entirely
### inside one genetic background: a perturbation signature and the response to
### an intervention.  A working intervention reverses the signature, which shows
### up as a negative slope.  Outputs are written to signature_reversal/{name}/.

rule SignatureReversal:
    input:
        results_signature = lambda wildcards: "deseq2/{c}/results.csv".format(
            c=get_signature_reversal_cfg(wildcards.comparison)['signature_contrast']
        ),
        results_response = lambda wildcards: "deseq2/{c}/results.csv".format(
            c=get_signature_reversal_cfg(wildcards.comparison)['response_contrast']
        ),
    output:
        scatter  = "signature_reversal/{comparison}/reversal_scatter.pdf",
        recovery = "signature_reversal/{comparison}/recovery_distribution.pdf",
        stats    = "signature_reversal/{comparison}/reversal_stats.csv",
        genes    = "signature_reversal/{comparison}/reversal_genes.csv",
    params:
        script = f"{workflow.basedir}/scripts/signature_reversal.R",
        outdir = "signature_reversal/{comparison}",
        # Human-readable "numerator / denominator" labels for the plot axes, so
        # the sign convention is visible on the figure itself.
        signature_label = lambda wildcards: (
            lambda sides: f"{sides[0]} / {sides[1]}"
        )(get_contrast_sides(
            get_signature_reversal_cfg(wildcards.comparison)['signature_contrast'])),
        response_label = lambda wildcards: (
            lambda sides: f"{sides[0]} / {sides[1]}"
        )(get_contrast_sides(
            get_signature_reversal_cfg(wildcards.comparison)['response_contrast'])),
        padj = lambda wildcards: float(
            get_signature_reversal_cfg(wildcards.comparison).get(
                'padj_threshold', config['DESeq2']['padj_threshold'])),
        lfc = lambda wildcards: float(
            get_signature_reversal_cfg(wildcards.comparison).get('lfc_threshold', 0.0)),
        restored_fraction = lambda wildcards: float(
            get_signature_reversal_cfg(wildcards.comparison).get(
                'restored_fraction', 0.5)),
        label_n = lambda wildcards: int(
            get_signature_reversal_cfg(wildcards.comparison).get('label_n', 20)),
        n_perm = lambda wildcards: int(
            get_signature_reversal_cfg(wildcards.comparison).get('n_perm', 1000)),
        # The conditions each contrast is built from, so the script can detect a
        # shared group — which would make the score uninterpretable.
        signature_groups = lambda wildcards: ",".join(get_contrast_sides(
            get_signature_reversal_cfg(wildcards.comparison)['signature_contrast'])),
        response_groups = lambda wildcards: ",".join(get_contrast_sides(
            get_signature_reversal_cfg(wildcards.comparison)['response_contrast'])),
    threads: 1
    resources:
        mem_mb        = 4000,
        runtime       = 30,
        cpus_per_task = 1,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/SignatureReversal/{comparison}.log"
    benchmark:
        "benchmark/SignatureReversal/{comparison}.benchmark.txt"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --results_signature {input.results_signature} \
            --results_response  {input.results_response} \
            --outdir            {params.outdir} \
            --name              "{wildcards.comparison}" \
            --signature_label   "{params.signature_label}" \
            --response_label    "{params.response_label}" \
            --padj              {params.padj} \
            --lfc               {params.lfc} \
            --restored_fraction {params.restored_fraction} \
            --label_n           {params.label_n} \
            --n_perm            {params.n_perm} \
            --signature_groups  "{params.signature_groups}" \
            --response_groups   "{params.response_groups}"
        """
