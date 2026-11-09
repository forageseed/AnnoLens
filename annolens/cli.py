"""annolens command line. Everything runs on this machine.

  annolens check    [--manifest annolens.toml]                 is this machine ready? (Python, modules, programs, inputs, key)
  annolens prepare  [--manifest annolens.toml] [--threads 8]   build the work dir (junctions, projections, homology, read sites)
  annolens scan     [--manifest annolens.toml]                 rule-based check of every gene + loci without a model
  annolens review   [--manifest annolens.toml] --queue P1      batch review of existing models with the configured model
  annolens addgene  [--manifest annolens.toml] [--queue missing_high]   batch review of loci without a model
  annolens mcp      [--manifest annolens.toml]                 MCP server (stdio)
  annolens install  --host claude|codex|cursor|generic [--manifest annolens.toml]   register the MCP server + install the skill
  annolens skill-md | skill-files                              print SKILL.md / the whole skill folder as JSON
The manifest is taken from --manifest, else the environment variable ANNOLENS_MANIFEST, else ./annolens.toml.
"""
import argparse, json, sys
from pathlib import Path

from .check import find_manifest


def main(argv=None):
    ap = argparse.ArgumentParser(prog='annolens'); sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('check', 'prepare', 'scan', 'review', 'addgene', 'mcp'):
        p = sub.add_parser(name); p.add_argument('--manifest'); p.add_argument('--out')
        if name == 'prepare': p.add_argument('--threads', type=int, default=8)
        if name in ('review', 'addgene'):
            p.add_argument('--queue', default='P1' if name == 'review' else 'missing_high'); p.add_argument('--ids', help='file with ids (overrides --queue)')
            p.add_argument('--workers', type=int, default=4)
    sub.add_parser('skill-md'); sub.add_parser('skill-files')
    p = sub.add_parser('install')
    p.add_argument('--host', required=True, choices=['claude', 'codex', 'cursor', 'generic'])
    p.add_argument('--manifest', help='project manifest; stored as ANNOLENS_MANIFEST in the host config')
    p.add_argument('--skills-dir'); p.add_argument('--dry-run', action='store_true')
    a = ap.parse_args(argv)
    if a.cmd == 'skill-files':
        from .install import skill_files
        sys.stdout.write(json.dumps(skill_files())); return 0
    if a.cmd == 'skill-md':
        from .install import skill_md
        sys.stdout.write(skill_md()); return 0
    if a.cmd == 'install':
        from .install import install
        return install(a)
    if a.cmd == 'check':
        from .check import run
        return run(a.manifest)
    mp = find_manifest(a.manifest)
    if a.cmd == 'mcp':
        from .mcp_server import serve
        if mp is None: return serve(None, None)            # starts anyway; tools explain how to point it at a project
        from .config import Manifest
        M = Manifest(mp); return serve(M, Path(a.out) if a.out else M.w('reviews'))
    if mp is None:
        sys.exit('annolens: no project manifest; pass --manifest, set ANNOLENS_MANIFEST or put annolens.toml here '
                 '(template: examples/annolens.toml). Run `annolens check` for details.')
    from .config import Manifest
    M = Manifest(mp); out = Path(a.out) if a.out else M.w('reviews')
    if a.cmd == 'prepare':
        from . import prepare; prepare.run(M, a.threads)
    elif a.cmd == 'scan':
        from . import scan; print(json.dumps(scan.run(M), indent=1))
    elif a.cmd in ('review', 'addgene'):
        from . import agent
        ids = Path(a.ids).read_text().split() if a.ids else M.w('queues', a.queue + '.txt').read_text().split()
        agent.run(M.path, 'gene' if a.cmd == 'review' else 'locus', ids, out, a.workers)
    return 0


if __name__ == '__main__':
    sys.exit(main())
