# Decisions and the machine rules behind them

The model proposes a decision; the server checks it and may override it (`machine_override` in the result). Report the
final decision as the server returns it.

## Existing gene models

| decision | meaning | machine check |
|---|---|---|
| `retain` | the current model is complete and supported | must be a complete ORF, otherwise it becomes `unresolved` |
| `recommend_candidate` | adopt a structure you built | the candidate must come from this session and be materializable; it is re-validated before it is written |
| `split` | the model joins two or more genes; replace it by the parts | >= 2 materializable parts from `split_part`, no overlap between parts, each part mostly inside the model's span; refused when the whole model is supported as one gene (full-chain long read, or an intact homolog of its full length) |
| `genome_frameshift` | a protein clearly longer than the model needs a frameshift or stop that the reads confirm in this genome | needs a read-supported site inside the model, needed by a protein longer than the model; otherwise `unresolved` |
| `assembly_conflict` | assembly-level problem (gap, misjoin, or reads disagree with the assembly) | - |
| `unresolved` | conflicting or insufficient evidence | - |

## Loci without a model

| decision | meaning | machine check |
|---|---|---|
| `add_gene` | add the gene you built | materializable candidate from this session, re-validated |
| `no_gene` | not a protein-coding gene (transposon, non-coding transcript, insufficient evidence) | - |
| `genome_frameshift` | the gene is disabled in this genome | needs a read-supported site from `readcheck` |
| `unresolved` | cannot decide | - |

## What makes a candidate materializable

- Complete ORF: ATG start, terminal stop, length a multiple of 3, no internal stop, no N; canonical splice sites
  (GT-AG, GC-AG, AT-AC).
- Built only from an observed transcript chain, or (new genes and split parts) from an intact protein projection.
- No same-strand CDS overlap >= 30 bp with another gene unless that gene is merged.
- At least one evidence gate:
  - a full-chain long read covering >= 95 % of the CDS;
  - every intron with >= 3 junction reads or reproduced by a projection, and start, stop and length (0.9-1.1x)
    anchored by an independent or related-species projection;
  - the dominant ORF across samples.
- Fix only what is broken: a complete model that already passes a gate is replaced only with a truncation proof.
  A truncation proof means an intact protein of about the candidate's length (0.9-1.1x) exists, the original is
  < 90 % of that protein, and every candidate intron is supported; the candidate's start and stop are anchored either
  by a projection or, since skill 3.8, by being the dominant ORF in five or more samples (such an ORF may also be
  measured against a protein that projects here with a frameshift or stop, since the genome evidently encodes the ORF).
- Anti-truncation: a candidate shorter than the dominant ORF and supported by fewer samples is refused.
