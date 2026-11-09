# AGENTS.md

AnnoLens is an MCP server plus an agent skill for evidence-driven gene-model curation (review, repair, split and add
protein-coding gene models). This file tells a coding agent how to install it for a user and how the repository is
organised. Everything runs on the user's own machine.

## Install or update for the user

When the user asks to install or update this plugin:

1. Check Python: `python --version` must be 3.11 or newer. The engine needs Linux, macOS or WSL (it reads BAM files
   with `pysam`, which has no native Windows build). On Windows, run the next steps inside WSL.
2. Install or update the package:
   `python -m pip install --upgrade git+https://github.com/forageseed/AnnoLens.git`
3. Register the MCP server and install the skill for the host you are running in:
   `annolens install --host <claude|codex|cursor|generic> [--manifest /path/to/annolens.toml]`
   - `--manifest` is the user's project file. Ask for its path if the user has one; never guess. Without it the server
     uses `./annolens.toml` of the directory the host starts it in.
   - Existing host configuration is backed up as `*.annolens.bak`. Add `--dry-run` to preview.
4. Run `annolens check` and report every line marked FAIL with its fix.
5. Tell the user to restart the agent session so the new MCP server is loaded, then to ask in plain words, for example
   "这个基因模型对吗？gene_123".

A new project: copy `examples/annolens.toml`, fill in genome, annotation, RNA alignments and protein sets, then
`annolens prepare` and `annolens scan`. Model API keys (only for batch review) are read from environment variables; never
write them into files or chat.

## What runs

- Python package `annolens` (`annolens/`), Python 3.11+. CLI: `annolens` (`annolens mcp` is the stdio MCP server).
- MCP tools named `annolens_*`; standard client configuration in `.mcp.json`.
- The engine validates every change; the host model chooses among observed evidence and never writes coordinates.
  The genome and the input annotation are never edited.

## Host adapters

| Host | Files |
|---|---|
| Any MCP client | `.mcp.json` (`annolens mcp`, stdio) |
| Claude Code | `.claude-plugin/` (plugin + marketplace), `skills/`, `agents/`, `commands/` |
| Codex | `.codex-plugin/plugin.json`, `.mcp.json`, `skills/` |
| npm | `npm/` (`npx annolens-mcp` runs `annolens mcp`; `install` / `check` are forwarded) |
| MCP Registry | `server.json` |

## Development rules

- Skill texts are versioned in `annolens/skills/`; change a skill by adding a new version and keep the old ones.
- `skills/annolens/` is generated: `python scripts/build_skill.py`, verified with `--check`.
- Before every push run `python scripts/check_sensitive.py`; it must report 0.
- Keep version numbers in `pyproject.toml`, `annolens/__init__.py`, `npm/package.json`, `server.json` and the plugin
  manifests equal.
