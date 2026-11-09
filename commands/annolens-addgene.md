---
description: Check a locus for a missing gene with AnnoLens and add it if the evidence supports it
argument-hint: <locus_id or chrom:start-end+/->
---

Review the locus `$ARGUMENTS` with the AnnoLens MCP tools, following the `annolens` skill
(references/locus-policy.md). Open one session with `annolens_open_locus` (if the id is not in the missing_loci queue,
list the queue first with `annolens_list_queue("missing_loci")` and use the matching id), build from a transcript ORF or
an intact projection when justified, and finish with `annolens_locus_decide`. Report the server's final decision and
the new gene's protein length if one was added.
