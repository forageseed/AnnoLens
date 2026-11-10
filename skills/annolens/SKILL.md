---
name: annolens
description: Curate protein-coding gene models from evidence with the AnnoLens MCP tools (AnnoLens 1.2.0). Use it to review a suspicious gene model (truncated, split, fused, missing start or stop, frameshift), to check whether a gene is missing at a locus, or to work through an AnnoLens review queue.
---

# AnnoLens 1.2.0

AnnoLens reviews ONE gene model or ONE locus at a time through the `annolens` MCP server. You reason over evidence
images and tables; the server builds candidates only from observed transcripts or intact protein projections, enforces
the budgets and re-validates every change before it is written. It never edits the genome or the input annotation.

## Before you start
If no `annolens_*` tools are available, AnnoLens is not connected yet: follow `references/setup.md`.

## Workflow
1. Existing gene model: open with `annolens_open_gene(gene_id)` and follow `references/gene-policy.md` (skill v3.8).
2. Locus without a model (possible missing gene): open with `annolens_open_locus(locus_id)` and follow
   `references/locus-policy.md` (skill v3.6).
3. Several genes or a whole queue: `annolens_list_queue(name)`, then one session per id; for hundreds of ids suggest the
   batch command `annolens review --manifest <annolens.toml> --queue <name>` (it uses the model configured in the manifest).
4. Always finish a session with the decide tool and report the server's final decision, including any machine
   override. What each decision means: `references/decisions.md`. All tools: `references/tools.md`.

The policies describe actions as JSON objects; in this MCP setting call the matching tools instead
(`references/tools.md` maps them).

## Without the MCP server
`scripts/check_gff.py genes.gff3 genome.fa` runs the deterministic model checks (ATG, terminal stop, length % 3,
internal stops, phases, splice motifs) on any GFF3 using only the Python standard library.
