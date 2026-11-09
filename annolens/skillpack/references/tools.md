# AnnoLens MCP tools

| policy action | existing gene | locus without a model |
|---|---|---|
| list a queue | `annolens_list_queue(name, offset, limit)` | `annolens_list_queue("missing_loci")` |
| open | `annolens_open_gene(gene_id)` | `annolens_open_locus(locus_id)` |
| `view` | `annolens_gene_view(session, start, end)` | `annolens_locus_view(session, start, end[, candidate])` |
| `orfs` | `annolens_gene_orfs(session, transcript)` | `annolens_locus_orfs(session, transcript)` |
| `build` | `annolens_gene_build(session, transcript[, start][, merge])` | `annolens_locus_build(session, transcript=T# or projection=H#[, start])` |
| `split_part` | `annolens_gene_split_part(session, transcript=T# or projection=H#[, start])` | - |
| `readcheck` | `annolens_readcheck(session, position[, kind])` | `annolens_readcheck(session, position[, kind])` |
| `decide` | `annolens_gene_decide(session, status, candidate_id?, parts?, confidence, rationale, ...)` | `annolens_locus_decide(session, status, candidate_id?, confidence, rationale, ...)` |

Every open/view/build returns text (FACTS, TRANSCRIPTS, PROJECTIONS tables) and an evidence image. Coordinates are
1-based in everything you see and send. The server limits views, builds and read checks per session (see FACTS.budgets).

Prompts: `annolens_review_gene(gene_id)`, `annolens_add_gene(locus_id)`.
Resources: `annolens://skill/gene`, `annolens://skill/locus`, `annolens://server/info`.
