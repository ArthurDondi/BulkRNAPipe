# ATRX rescue experiment: what the PCA is telling us, and what we do about it

*A short note for everyone working on the Epicode bulk RNA-seq dataset.*

---

## 1. The setup

Five groups, three replicates each:

| Group | What it is | Parental clone | Vector |
| --- | --- | --- | --- |
| `TP53` | TP53 KO | TP53 clone | none |
| `E6` | TP53 KO + ATRX KO | **E6 clone** | none |
| `EmptyVector` | E6 + empty EGFP/mCherry vector | E6 clone | empty |
| `ATRX_FL` | E6 + ATRX full length (rescue) | E6 clone | ATRX-FL |
| `ATRX_IFF` | E6 + ATRX in-frame fusion (semi-rescue) | E6 clone | ATRX-IFF |

The hope was that `ATRX_FL` would look like `TP53` (rescue worked) and that `E6`
would look like `EmptyVector` (the vector does nothing).

## 2. What we actually see

Replicates cluster tightly — the data quality is good. But:

- **PC1 (77%)** separates `TP53` from everything else. That is the E6 clone.
- **PC2 (10%)** separates `E6` from the three transduced groups, with `ATRX_FL`
  pushed furthest the other way.

So the two things we wanted to read off the PCA are both sitting behind larger
effects that have nothing to do with the question.

## 3. Why "they should cluster together" is not a question the data can answer

Every group is **one clone**. `TP53` is one clone; `E6` and everything derived
from it is a second clone. So the difference between them is:

> ATRX genotype **+** clonal background **+** everything that happened during
> the derivation of that clone

and there is no sample anywhere in the experiment where one of those varies
while the others are held fixed. They are **confounded**: not "hard to
separate", but *the same column of numbers*. Three replicates make each group
mean more precise; they do not create a group that separates the terms.

The same applies, more mildly, to PC2: `E6` is untransduced and the other three
are transduced, so "vector effect" and "which construct" are tangled in the same
way at the `E6` → `EmptyVector` boundary.

## 4. Why we cannot fix this with a covariate — the bit people ask about

The natural instinct is to add a term to the model: `~ vector + condition`, or a
`clone` covariate, or run RUV / ComBat / `removeBatchEffect` and carry on. Here
is why that does not work, using one gene with made-up group means to keep it
concrete (log2 scale):

```
TP53 = 10    E6 = 6    EmptyVector = 7    ATRX_FL = 9    ATRX_IFF = 7.5
```

Adding a `vector` term means asking the model to satisfy:

```
TP53 :  Int                    = 10
E6   :  Int + b_E6             =  6
EV   :  Int + v + b_EV         =  7
FL   :  Int + v + b_FL         =  9
IFF  :  Int + v + b_IFF        =  7.5
```

Five equations, six unknowns. Pick **any** value for the vector effect `v` and
the rest absorbs it:

| `v` | `b_EV` | `b_FL` | `b_IFF` | fits the data? |
| --- | --- | --- | --- | --- |
| 0.0 | −3.0 | −1.0 | −2.5 | **exactly** |
| +2.0 | −5.0 | −3.0 | −4.5 | **exactly** |
| −5.0 | +2.0 | +4.0 | +2.5 | **exactly** |

All three reproduce every observed value perfectly. One says the vector does
nothing, one says it drives expression up 4-fold, one says it drops it 32-fold.
**The data cannot distinguish them.** DESeq2 detects this and refuses to fit:
*"the model matrix is not full rank."*

The counting rule behind it: our experiment measures **five numbers** — one mean
per group. A model with six parameters describing five numbers has no unique
answer. Replicates improve the *precision* of those five numbers; they never add
a sixth.

**So there is nothing to tune, and no better software to reach for.** The
quantity we would want to subtract is not small, or noisy — in this design it is
**not a measurable quantity at all**.

### But can't we just measure the vector effect and subtract it?

Yes — and we already do. This is worth being precise about, because "confounded"
is easy to over-apply.

Within the E6 lineage the vector effect **is** estimable: `EmptyVector` minus
`E6` measures it directly. And subtracting it from the rescue is not a new
analysis, it is the one already in the config. The algebra collapses:

```
(ATRX_FL − E6) − (EmptyVector − E6)  =  ATRX_FL − EmptyVector
 ^ rescue vs untransduced KO   ^ vector effect       ^ what we already run
```

Subtracting the vector effect **is** using `EmptyVector` as the denominator.
They are the same operation, and `EmptyVector_vs_ATRX_FL` already does it. There
is nothing extra to gain, and the `ATRX_FL` / `ATRX_IFF` comparisons were never
affected by the confounding in the first place — both sides carry the vector, so
it cancels whether or not you think of it as a subtraction.

So the confounding is narrower than "the whole experiment". It bites in exactly
two places:

- **Anything involving `TP53`.** The clone difference has no control group
  anywhere in the design, so it cannot be measured, let alone subtracted. This
  is the one that is genuinely unfixable.
- **Decomposing `E6` → `EmptyVector` into sub-causes** — vector insertion vs
  reporter-protein burden vs the selection/expansion the transduced cells went
  through. We can measure the total, not the parts. In practice we don't need
  the parts.

What we cannot do is take the vector effect measured in the E6 background and
apply it to `TP53`. That step assumes the vector would do the same thing in a
different clone, and nothing in this experiment tests that assumption — it is
exactly the additivity assumption behind the interaction contrast in §5c.

### What about RUV?

RUV was suggested and it is a reasonable thing to ask. Our conclusion: **use it
for QC if you like, but not in the DE models.**

RUV estimates a hidden factor `W` and adds it to the design. That only works if
`W` is at least partly independent of the thing you care about. Here it cannot
be — every candidate technical factor is a deterministic function of the group
label. The difference from the `~ vector + condition` case above is that `W` is
*continuous*, so it is not an exact duplicate of a group column: the model
appears full rank, DESeq2 runs happily, and `W` quietly absorbs part of the
biology. You get a shorter, cleaner-looking gene list and **no error message**.
That is worse than the clean failure, not better.

If you do run it, the diagnostic is: regress `W_1` on the group labels. If R² is
near 1, `W` is the condition variable wearing a hat.

The one honest use: run RUV or `removeBatchEffect` as a *picture*, clearly
labelled, to show that a vector effect exists. The pipeline now produces exactly
that panel, and stamps a warning onto the figure itself so it cannot be
mistaken for evidence when it turns up in a slide deck.

## 5. What we do instead

The good news: the questions we actually care about **are** answerable from
these 15 samples. They just have to be asked as comparisons that never cross a
confounded boundary. Four things are now in the pipeline.

### 5a. PCAs that are not dominated by the design — `pca/`

A grid of sample sets × gene views, so no single effect can hide the others:

| | `all_genes` | `no_vector_genes` | `no_vector_no_transgene` |
| --- | --- | --- | --- |
| `all` | the original plot | minus EGFP/mCherry | minus those + ATRX |
| `E6_derived` | drops TP53 → the clone axis is gone | | |
| `transduced_only` | also drops E6 → construct is the only difference | | |

Two things this fixes. **EGFP and mCherry are rows in our count matrix** (we
align to `hg38_plus_vectors`) — they are zero in `TP53`/`E6` and very high in the
transduced groups, which puts them straight into the top-variance gene set *and*
shifts the library-size normalisation for every other gene. **ATRX itself** is
the largest single difference between the groups by construction, so a PCA
including it partly re-plots the design rather than its consequences.

Each combination also gets PC3/PC4 and a scree plot — when PC1 is spoken for,
the interesting structure often sits lower down.

### 5b. Rescue scored as signature reversal — `signature_reversal/` ⭐

**This is the primary rescue readout, replacing "do they cluster".**

Take the KO signature from `TP53_vs_E6` and ask what the rescue does to those
same genes, measured as `ATRX_FL` vs `EmptyVector`. Both contrasts sit entirely
inside one background, so the clonal difference never enters. Plot one log2FC
against the other:

- a working rescue pushes KO-responsive genes back the other way → **negative slope**
- reported as a **reversal score**: 1 = complete reversal, 0 = no effect

`ATRX_IFF` gets the same treatment, so "semi-rescue" becomes a number rather
than an impression. The baseline comes from a **permutation null** computed inside each comparison
(the response log2FCs are shuffled across genes), not from an assumed 50%.
Read `excess_over_null` and `perm_p`.

One trap worth knowing: the obvious negative control — `EmptyVector` vs `E6` —
is invalid here. It shares the `E6` group with the signature, so E6's sampling
noise enters one axis positively and the other negatively and manufactures a
negative correlation from nothing (~66% "reversed" on pure noise). The script
now detects any shared condition, marks the comparison `interpretable = FALSE`,
and stamps the warning onto the figure.

On simulated data with a known answer the scores recovered the truth in the
right order and roughly the right size (0.59 for a true 0.70, 0.28 for a true
0.35, 0.13 for a true 0.00). The mild shrinkage toward zero is expected, which
is exactly why the negative control matters.

### 5c. Interaction contrasts — `deseq2_interaction/`

`(ATRX_FL / EmptyVector) − (TP53 / E6)`: how far the rescue actually moved each
gene, minus how far it had to move to undo the knockout. The first term is
measured inside the transduced background, the second inside the untransduced
one, so the clone and vector effects **cancel in the subtraction** instead of
being modelled — no covariate, no rank problem. This is the honest replacement
for the plain `TP53_vs_ATRX_FL` contrast.

Read it as a **shortfall**: `log2FC = 0` is a complete rescue, negative means the
rescue fell short, positive means it overshot the TP53 level. Note the direction
of the second pair — it is `TP53 / E6` (the restoration target), not `E6 / TP53`
(the knockout effect). Flipping it reports the sum of the two effects instead of
the residual, and a perfect rescue would score twice the effect size rather than
zero.

It assumes the nuisance effects are the same size in both pairs (additivity).
That assumption is untestable in this design, so it goes in the methods section —
but it is far weaker and far more visible than the assumption buried inside a
batch-correction step.

### 5d. Design QC — `design_qc/`

Library size, featureCounts assignment rate, genes detected, sample-sample
correlation and distance, and the reporter/transgene expression levels — all
samples on one plot, because a group-linked technical shift is invisible one
sample at a time. Read alongside the MultiQC reports in `QC/`.

## 6. What the ATRX numbers already tell us

From the normalised counts (endogenous ATRX locus):

| Group | mean | % of TP53 | replicate spread |
| --- | --- | --- | --- |
| `TP53` | 2705 | 100% | 1.0× |
| `E6` | 636 | 24% | 1.2× |
| `EmptyVector` | 771 | 29% | 1.2× |
| `ATRX_FL` | 2351 | **87%** | 1.2× |
| `ATRX_IFF` | 4068 | 150% | **2.8×** |

Three things worth knowing:

1. **`ATRX_FL` restores ATRX to near-physiological level** (87% of TP53), not a
   gross overexpression. An earlier worry — that FL sits far from TP53 on PC2
   simply because the transgene is massively overexpressed — is **not supported**.
   Whatever separates them on PC2, it is not ATRX dosage.
2. **`ATRX_IFF` is over-expressed and inconsistent** — 150% on average but
   varying 2.8-fold across replicates (2316 to 6504). Worth flagging before any
   IFF result is interpreted as a property of the fusion rather than of dosage.
3. **The KO is not transcript-null** — `E6` retains ~24% of TP53's ATRX signal.
   Expected for a frameshift/exon-disrupting KO, but worth stating explicitly
   rather than describing the line as ATRX-negative.

One caveat to check: the vector-borne `ATRX_FL` / `ATRX_IFF` features carry far
fewer counts than the endogenous `ATRX` row, which suggests most transgene reads
are being assigned to the endogenous locus. That would make the `ATRX` row
"endogenous + transgene" in the transduced samples. It does not change the
conclusions above, but someone should confirm how the custom GTF defines those
features.

## 7. What would actually fix the design

For the next round, the single most valuable addition is a **`TP53` + empty
vector** group. It is the one missing cell: it lets the vector effect be
measured in a background where ATRX is intact, which both pins the effect down
and makes the additivity assumption in §5c *testable* instead of assumed.

Second: **two or three independent ATRX-KO clones**. That is what breaks the
clone/genotype confounding — and it is the only thing that does.

## 8. What to say when asked

> "TP53 and ATRX_FL don't cluster on the PCA because they are different clonal
> backgrounds with different transduction status, and in this design those are
> inseparable from the ATRX genotype. So we assess the rescue by signature
> reversal within the E6 background instead, where the confounding does not
> apply."

Not: "we corrected for it."

---

*Implementation: `workflow/rules/pca.smk`, `signature_reversal.smk`,
`design_qc.smk`, and the `DESeq2Interaction` rule. Configured under `PCA:`,
`SignatureReversal:`, `DESeq2.interactions:` and `DesignQC:` in
`config/config_epicode.yaml`.*
