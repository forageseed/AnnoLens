"""One-step install into an agent host: register the local MCP server (`annolens mcp`) AND install the skill (built from
the same versioned texts the server serves).

  annolens install --host claude  [--manifest ./annolens.toml]
  annolens install --host codex   [--manifest ./annolens.toml]
  annolens install --host cursor  [--manifest ./annolens.toml]
  annolens install --host generic [--manifest ./annolens.toml]    (prints mcpServers JSON, writes the skill folder)
With --manifest the project path is stored as ANNOLENS_MANIFEST in the host config; without it the server uses
./annolens.toml of the directory the host starts it in. Config files are backed up to *.annolens.bak first.
"""
import json, re, shutil, subprocess, sys
from pathlib import Path

from . import __version__
from .session import SKILL, skill_text

PACK = Path(__file__).with_name('skillpack')


def skill_md():
    gv, lv = SKILL['gene'][1], SKILL['locus'][1]
    return f"""---
name: annolens
description: Curate protein-coding gene models from evidence with the AnnoLens MCP tools (AnnoLens {__version__}). Use it to review a suspicious gene model (truncated, split, fused, missing start or stop, frameshift), to check whether a gene is missing at a locus, or to work through an AnnoLens review queue.
---

# AnnoLens {__version__}

AnnoLens reviews ONE gene model or ONE locus at a time through the `annolens` MCP server. You reason over evidence
images and tables; the server builds candidates only from observed transcripts or intact protein projections, enforces
the budgets and re-validates every change before it is written. It never edits the genome or the input annotation.

## Before you start
If no `annolens_*` tools are available, AnnoLens is not connected yet: follow `references/setup.md`.

## Workflow
1. Existing gene model: open with `annolens_open_gene(gene_id)` and follow `references/gene-policy.md` (skill v{gv}).
2. Locus without a model (possible missing gene): open with `annolens_open_locus(locus_id)` and follow
   `references/locus-policy.md` (skill v{lv}).
3. Several genes or a whole queue: `annolens_list_queue(name)`, then one session per id; for hundreds of ids suggest the
   batch command `annolens review --manifest <annolens.toml> --queue <name>` (it uses the model configured in the manifest).
4. Always finish a session with the decide tool and report the server's final decision, including any machine
   override. What each decision means: `references/decisions.md`. All tools: `references/tools.md`.

The policies describe actions as JSON objects; in this MCP setting call the matching tools instead
(`references/tools.md` maps them).

## Without the MCP server
`scripts/check_gff.py genes.gff3 genome.fa` runs the deterministic model checks (ATG, terminal stop, length % 3,
internal stops, phases, splice motifs) on any GFF3 using only the Python standard library.
"""


def skill_files():
    """The whole skill folder as {relative path: text}; one source for every installer and for the repository copy."""
    files = {'SKILL.md': skill_md(),
             'references/gene-policy.md': f"# Policy for existing gene models (skill v{SKILL['gene'][1]})\n\n" + skill_text('gene') + '\n',
             'references/locus-policy.md': f"# Policy for loci without a model (skill v{SKILL['locus'][1]})\n\n" + skill_text('locus') + '\n'}
    for p in sorted(PACK.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            files[p.relative_to(PACK).as_posix()] = p.read_text(encoding='utf-8')
    return files


def write_skill_dir(root, write):
    for rel, text in skill_files().items():
        write(Path(root) / rel, text)


def install(a):
    cmd = [shutil.which('annolens') or 'annolens', 'mcp']
    env = {'ANNOLENS_MANIFEST': str(Path(a.manifest).expanduser().resolve())} if a.manifest else {}

    def write(path, text):
        path = Path(path).expanduser()
        if a.dry_run: print(f'[dry-run] would write {path} ({len(text)} chars)'); return
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists(): shutil.copy2(path, str(path) + '.annolens.bak')
        path.write_text(text, encoding='utf-8'); print('wrote', path)

    entry = {'command': cmd[0], 'args': cmd[1:], **({'env': env} if env else {})}
    if a.host == 'claude':
        reg = ['claude', 'mcp', 'add', '--scope', 'user'] + [x for k, v in env.items() for x in ('-e', f'{k}={v}')] + ['annolens', '--'] + cmd
        print(('[dry-run] ' if a.dry_run else '') + 'register: ' + ' '.join(reg))
        if not a.dry_run:
            if not shutil.which('claude'): sys.exit('claude CLI not found; register the server with the command above')
            subprocess.run(['claude', 'mcp', 'remove', '--scope', 'user', 'annolens'], capture_output=True)
            subprocess.run(reg, check=True)
        write_skill_dir(Path(a.skills_dir or '~/.claude/skills') / 'annolens', write)
    elif a.host == 'codex':
        cfg = Path('~/.codex/config.toml').expanduser(); old = cfg.read_text(encoding='utf-8') if cfg.exists() else ''
        old = re.sub(r'(?m)^\[mcp_servers\.annolens(?:\.env)?\]\n(?:(?!\[).*(?:\n|$))*', '', old).rstrip() + '\n'
        block = (f'\n[mcp_servers.annolens]\ncommand = {json.dumps(cmd[0])}\nargs = {json.dumps(cmd[1:])}\n'
                 'startup_timeout_sec = 300\ntool_timeout_sec = 900\n')
        if env: block += '\n[mcp_servers.annolens.env]\n' + ''.join(f'{k} = {json.dumps(v)}\n' for k, v in env.items())
        write(cfg, old + block)
        write_skill_dir(Path(a.skills_dir or '~/.codex/skills') / 'annolens', write)
    elif a.host == 'cursor':
        cfg = Path('~/.cursor/mcp.json').expanduser(); j = json.loads(cfg.read_text(encoding='utf-8') or '{}') if cfg.exists() else {}
        j.setdefault('mcpServers', {})['annolens'] = entry
        write(cfg, json.dumps(j, indent=2) + '\n')
        write_skill_dir(Path(a.skills_dir or '~/.cursor/skills') / 'annolens', write)
    else:
        print(json.dumps({'mcpServers': {'annolens': entry}}, indent=2))
        write_skill_dir(Path(a.skills_dir or 'annolens-skill'), write)
    print('done: MCP server "annolens" registered for', a.host, 'and skill installed; next: annolens check')
    return 0
