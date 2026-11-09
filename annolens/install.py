"""One-step install into an agent host: register the MCP server AND install the skill (same text the server serves).

  annolens install --host codex  --ssh myserver --manifest /data/proj/annolens.toml
  annolens install --host claude --manifest ./annolens.toml          (local server, no ssh)
  annolens install --host generic --ssh myserver --manifest ...      (prints mcpServers JSON, writes SKILL.md)
Config files are backed up to *.annolens.bak before they are changed.
"""
import json, re, shlex, shutil, subprocess, sys
from pathlib import Path

from . import __version__
from .session import skill_text


def skill_md():
    return f"""---
name: annolens
description: Curate protein-coding gene models with evidence images and read checks (AnnoLens {__version__}). Use it to review a suspicious gene model (truncated, split, fused, missing start/stop, frameshift) or a locus where a gene may be missing, through the annolens MCP tools.
---

# AnnoLens {__version__}

AnnoLens reviews ONE gene model or ONE locus at a time through the `annolens` MCP server. You reason over the
evidence; the server builds candidates only from observed transcripts or intact protein projections, enforces the
budgets and re-validates every change before it is written. It never edits the genome or the input annotation.

## Tools
| step | existing gene | locus without a model |
|---|---|---|
| list a queue | `annolens_list_queue(name)` | `annolens_list_queue("missing_loci")` |
| open | `annolens_open_gene(gene_id)` | `annolens_open_locus(locus_id)` |
| look further | `annolens_gene_view(session,start,end)` | `annolens_locus_view(session,start,end)` |
| ORFs | `annolens_gene_orfs(session,transcript)` | `annolens_locus_orfs(session,transcript)` |
| build | `annolens_gene_build(session,transcript[,start][,merge])` | `annolens_locus_build(session,transcript=T# or projection=H#)` |
| read check | `annolens_readcheck(session,position)` | `annolens_readcheck(session,position)` |
| finish | `annolens_gene_decide(session,status,...)` | `annolens_locus_decide(session,status,...)` |

Call the tools instead of writing the JSON actions described below. Always finish a session with decide.

## Existing gene models
{skill_text('gene')}

## Loci without a model
{skill_text('locus')}
"""


def install(a):
    cmd = ['ssh', '-T', '-o', 'BatchMode=yes', '-o', 'ServerAliveInterval=30', a.ssh, f'{a.remote_cmd} mcp --manifest {shlex.quote(a.manifest)}'] if a.ssh \
        else [a.remote_cmd, 'mcp', '--manifest', a.manifest]
    md = skill_md()

    def write(path, text):
        path = Path(path).expanduser()
        if a.dry_run: print(f'[dry-run] would write {path} ({len(text)} chars)'); return
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists(): shutil.copy2(path, str(path) + '.annolens.bak')
        path.write_text(text, encoding='utf-8'); print('wrote', path)

    if a.host == 'claude':
        reg = ['claude', 'mcp', 'add', '--scope', 'user', 'annolens', '--'] + cmd
        print(('[dry-run] ' if a.dry_run else '') + 'register: ' + ' '.join(reg))
        if not a.dry_run:
            if not shutil.which('claude'): sys.exit('claude CLI not found; register the server with the command above')
            subprocess.run(reg, check=True)
        write(Path(a.skills_dir or '~/.claude/skills') / 'annolens' / 'SKILL.md', md)
    elif a.host == 'codex':
        cfg = Path('~/.codex/config.toml').expanduser(); old = cfg.read_text(encoding='utf-8') if cfg.exists() else ''
        old = re.sub(r'(?m)^\[mcp_servers\.annolens\]\n(?:(?!\[).*(?:\n|$))*', '', old).rstrip() + '\n'
        write(cfg, old + f'\n[mcp_servers.annolens]\ncommand = {json.dumps(cmd[0])}\nargs = {json.dumps(cmd[1:])}\nstartup_timeout_sec = 300\ntool_timeout_sec = 900\n')
        write(Path(a.skills_dir or '~/.codex/skills') / 'annolens' / 'SKILL.md', md)
    else:
        print(json.dumps({'mcpServers': {'annolens': {'command': cmd[0], 'args': cmd[1:]}}}, indent=2))
        write(Path(a.skills_dir or 'annolens-skill') / 'SKILL.md', md)
    print('done: MCP server "annolens" registered for', a.host, 'and skill installed')
