<p align="center"><img src="docs/logo.svg" alt="AnnoLens" width="420"></p>


<p align="center">Review and repair gene models from the evidence, the way a human curator does: an agent plus a skill, served over MCP, for any species</p>

<p align="center"><a href="README.md">中文</a> | <b>English</b></p>

<p align="center">
<a href="https://www.npmjs.com/package/annolens-mcp"><img alt="npm" src="https://img.shields.io/npm/v/annolens-mcp?style=flat-square&color=cb3837&label=npm"></a>
<a href="server.json"><img alt="MCP Registry" src="https://img.shields.io/badge/MCP-registry-5c4ee5?style=flat-square"></a>
<a href="https://skills.sh/forageseed/AnnoLens"><img alt="skills.sh" src="https://skills.sh/b/forageseed/AnnoLens?style=flat-square"></a>
<a href="npm/package.json"><img alt="Node.js" src="https://img.shields.io/badge/Node.js-%E2%89%A518-339933?style=flat-square&logo=nodedotjs&logoColor=white"></a>
<a href="pyproject.toml"><img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.11-3776ab?style=flat-square&logo=python&logoColor=white"></a>
<a href="LICENSE"><img alt="License" src="https://img.shields.io/badge/license-Apache--2.0-0e7c86?style=flat-square"></a>
</p>

<p align="center"><code>npx -y annolens-mcp install --host codex --ssh &lt;server&gt; --manifest &lt;path/annolens.toml&gt;</code></p>


---

## What it is

Automatic annotation pipelines such as EVM, BRAKER and MAKER always leave some models wrong. Some are truncated, some are split in two, some fuse two genes, some lack a start or stop codon, some have no UTRs, and some genes are missing entirely. Human curators fix these one by one in a genome browser, checking transcripts, homologous proteins and splice reads. AnnoLens gives this job to a vision-capable language model and keeps it from making arbitrary changes:

- **The agent judges evidence images.** Each gene or locus gets an image and a text table. The image shows the current model, observed transcripts per sample, homologous protein projections onto the genome, and splice-junction reads.
- **A versioned skill defines how to review.** It covers how to settle gene extent before structure, which evidence counts, and when to retain, change, add a gene or leave it unresolved.
- **An engine enforces the hard rules.** A candidate can only come from an observed transcript or an intact protein projection; the model never types coordinates. The engine checks every candidate, and the model has no say over these rules:
  - complete ORF and canonical splice sites
  - no conflict with neighbouring genes
  - every intron supported by reads or by a projection
  - length consistent with homologs
  - no shortening below the ORF that most samples agree on
- **It proposes and never edits.** Outputs are review records and candidate GFF3 files. The input annotation and the genome are never modified.

## Features

| Feature | Details |
|---|---|
| Genome-wide rule scan | `annolens scan` checks every gene and adds the problem genes to prioritised queues. It covers:<br>• coding integrity and UTRs<br>• overlaps between genes<br>• evidence support<br>• length against an independent annotation, splits and fusions<br>• read-supported frameshifts<br>• repeat overlap<br>It also finds loci that have no gene model. |
| Repair existing models | Extend truncated models, merge split ones, complete missing start/stop codons, take UTRs from observed transcripts |
| Split fused models | When one model joins two or more genes, the agent builds each gene and each is validated like a new gene; a model that is supported as one gene is never split |
| Add missing genes | Build a new gene from transcripts or an intact projection at loci that have no model |
| Read checks | At any position, count genomic long reads (e.g. HiFi) and RNA reads that agree or disagree with the assembly, to separate assembly errors from real frameshifts and stops |
| Rule-based UTRs | Add UTRs in bulk from transcripts whose intron chain matches the model, with conflict checks against neighbouring genes |
| Batch review | Review thousands of genes in parallel with any vision model behind an OpenAI- or Anthropic-compatible API |
| Conversational use | Ask in plain words inside Codex, Claude Code, WorkBuddy, Cursor or any MCP host |

## Install

On the machine that holds the data:

```bash
pip install git+https://github.com/forageseed/AnnoLens.git
# external tools: miniprot (projections); StringTie when you supply BAMs without GTFs
```

## Quick start

1. Copy [`examples/annolens.toml`](examples/annolens.toml) and fill in the genome, annotation, RNA alignments and protein sets. The long reads, HiFi, repeat annotation and model endpoint fields are optional.
2. Build the work directory and run the genome-wide scan:

```bash
annolens prepare --manifest annolens.toml --threads 16
annolens scan    --manifest annolens.toml
```

3. Run the batch reviews. The API key is read from an environment variable only:

```bash
export ANNOLENS_API_KEY=...        # the variable named by api_key_env in the manifest
annolens review  --manifest annolens.toml --queue P1 --workers 4
annolens addgene --manifest annolens.toml --queue missing_high --workers 4
```

Results go to `<work_dir>/reviews/`: one review record (JSON) per gene. Validated changes and new genes are also written as GFF3.

### Protein roles

| role | meaning | used for |
|---|---|---|
| `independent` | another annotation of the same species, made independently | anchoring start, stop and length; main reference for truncation |
| `related` | proteome of a related species | anchoring start, stop and length |
| `paralog` | proteins from the same pipeline as the reviewed annotation (e.g. homeologs) | supporting introns only |

Do not mark a proteome that was used as evidence for the reviewed annotation as `independent`: it repeats the same truncations.

## Use it in conversation (Codex, Claude Code, WorkBuddy, ...)

In any MCP host you **just talk** and ask in plain words. The host model calls AnnoLens on its own: it looks at the evidence, builds candidates and checks reads, then answers with a decision and the evidence images. You do not need to know any commands or tool names.

> **You**: Is this gene model right? gene_123
>
> **Assistant**: The current model is 199 aa and only covers the 3′ part of the gene. Transcripts from 7 samples and the independent protein support a 646 aa structure, but the end of exon 1 carries a 1-bp frameshift. 26 HiFi reads and 170 RNA reads confirm that this is the real sequence of this genome. Decision: genome_frameshift; the model is kept and flagged. (evidence images attached)

> **You**: Is a gene missing at chr3:1,204,000–1,212,000?
>
> **Assistant**: Yes. IsoSeq transcript T3 encodes a complete 606 aa ORF whose length matches both protein sets, and every intron has junction reads. The candidate was built and validated. Decision: add_gene. (evidence images and candidate GFF3 attached)

> **You**: Review the first 20 genes of the review queue and summarise.

Other typical questions:

- "Were these two neighbouring genes originally one?"
- "Is this gene really two genes joined together?"
- "Why is this gene so much shorter than its homologs?"
- "Is this frameshift an assembly error?"

Every decision is saved as a review record (JSON). Changed and new genes are also written as GFF3. The input annotation and the genome are never modified.

### One-time setup

Run this once on the computer you work from (Node.js required). It **registers the MCP server and installs the skill** together. The skill text is read from the annolens installed on the server, so it always matches the server's rules:

```bash
npx -y annolens-mcp install --host codex  --ssh <server> --manifest /path/on/server/annolens.toml
npx -y annolens-mcp install --host claude --ssh <server> --manifest /path/on/server/annolens.toml
npx -y annolens-mcp install --host cursor --ssh <server> --manifest /path/on/server/annolens.toml
npx -y annolens-mcp install --host generic --ssh <server> --manifest /path/on/server/annolens.toml   # WorkBuddy and others: prints mcpServers JSON
```

You can also **install from the chat**: tell your agent "install AnnoLens; my data are on <server> at /path/annolens.toml", and it runs the command above.

- Leave out `--ssh` when the data are on this machine; add `--dry-run` to see what would change.
- Host configuration files are backed up as `*.annolens.bak` before they are changed.
- To check the connection: `npx -y annolens-mcp doctor --ssh <server> --manifest ...`
- Without Node.js, run `annolens install` on the server to do the same.

### Install channels

| channel | command / usage |
|---|---|
| npm (recommended) | `npx -y annolens-mcp install --host <codex\|claude\|cursor\|generic> --ssh <server> --manifest <path>` |
| from the chat | tell your agent "install AnnoLens; my data are on <server> at <path>" |
| skills.sh (skill only) | `npx skills add forageseed/AnnoLens`; on first use the skill guides you through connecting the MCP server |
| Claude Code plugin | `/plugin marketplace add forageseed/AnnoLens`, then `/plugin install annolens@annolens`; set the environment variable `ANNOLENS_MANIFEST` (optionally `ANNOLENS_SSH`) |
| Codex plugin | `.codex-plugin/` in this repository; same environment variables |
| MCP Registry | `io.github.forageseed/annolens` (see [`server.json`](server.json)) |

### What the skill package contains

An install gives the host a complete skill folder, plus sub-agents and slash commands (Claude Code plugin):

```
skills/annolens/
  SKILL.md                    when to use it, how to connect, the workflow
  references/gene-policy.md   policy for existing gene models (incl. split)
  references/locus-policy.md  policy for loci without a model (add gene)
  references/decisions.md     decisions and the machine rules behind them
  references/tools.md         MCP tool reference
  references/setup.md         install and connection
  scripts/check_gff.py        deterministic checks without MCP (any GFF3 + FASTA)
agents/annolens-curator.md    sub-agent: review one gene or locus
agents/annolens-batch.md      sub-agent: review a queue and summarise
commands/annolens-review.md   /annolens-review <gene>
commands/annolens-addgene.md  /annolens-addgene <locus>
commands/annolens-queue.md    /annolens-queue <queue> [N]
```

The skill folder is generated by `annolens skill-files` from the same versioned texts the server uses; the repository copy is checked with `python scripts/build_skill.py --check`.

<details>
<summary>Tool reference (developers)</summary>

| step | existing gene | locus without a model |
|---|---|---|
| list a queue | `annolens_list_queue` | `annolens_list_queue` |
| open | `annolens_open_gene` | `annolens_open_locus` |
| other window | `annolens_gene_view` | `annolens_locus_view` |
| ORFs | `annolens_gene_orfs` | `annolens_locus_orfs` |
| build | `annolens_gene_build` | `annolens_locus_build` (transcript or projection) |
| split (one part at a time) | `annolens_gene_split_part` | - |
| read check | `annolens_readcheck` | `annolens_readcheck` |
| decide | `annolens_gene_decide` | `annolens_locus_decide` |

The server also provides:

- the prompts `annolens_review_gene` and `annolens_add_gene`;
- the resources `annolens://skill/gene` and `annolens://skill/locus`.

It sets how many views, builds and read checks each session may use, and it re-validates every change or new gene before accepting it.

</details>

## Decisions

| existing gene | meaning |
|---|---|
| `retain` | the current model is complete and supported |
| `recommend_candidate` | adopt a structure the engine has validated |
| `split` | the model is a fusion; it is replaced by two or more validated genes |
| `genome_frameshift` | reads confirm that this genome really has the frameshift or stop the independent protein needs; the model is kept and flagged |
| `assembly_conflict` | assembly-level problem (gap, misjoin, or reads disagree with the assembly) |
| `unresolved` | conflicting evidence |

| locus without a model | meaning |
|---|---|
| `add_gene` | add the gene the engine has validated |
| `no_gene` | not a protein-coding gene (transposon, non-coding transcript, insufficient evidence) |
| `genome_frameshift` | the gene is disabled in this genome |
| `unresolved` | cannot decide |

## Design notes

- **Structure only from observed evidence.** The model chooses which transcript, which ATG or which projection to use; it never writes coordinates.
- **Evidence gates.** A change is accepted only if at least one of these holds:
  - a full-chain long read supports it;
  - its introns have read support and its start, stop and length are anchored by a projection;
  - it is the dominant ORF across samples.
- **Fix only what is broken.** A complete model that already has supporting evidence is never replaced. The one exception is when the dominant ORF and an intact independent protein both show that the model is truncated.
- **Anti-truncation.** A candidate is never accepted if it is shorter than the dominant ORF and has fewer samples supporting it.
- **Provenance.** Every candidate records a content hash, the engine hash and an input fingerprint. If the evidence or the code changes, the old candidate becomes invalid.

## Known limitations

- The host model must accept images; with the text tables alone, its judgement is weaker.
- Rule-based UTRs require the transcript's intron chain to match the model exactly; unclear cases go to the agent.

## Layout

```
annolens/            Python package: config, engine, tools, addgene, readcheck, scan, prepare, session, agent, mcp_server, install, render
annolens/skills/     versioned skill texts (existing genes gene_v3.7 with split; previous gene_v3.5 kept; missing genes locus_v3.6)
examples/            example project manifest
scripts/             release gate check_sensitive.py; skill file generator/check build_skill.py
skills/annolens/     skill file (used by skills.sh and plugins; generated by build_skill.py)
.claude-plugin/ .codex-plugin/ server.json   Claude Code and Codex plugin manifests, MCP Registry entry
npm/                 npm package annolens-mcp: MCP launcher and one-step installer
```

## Citation and license

A manuscript is in preparation. License: [Apache-2.0](LICENSE).
