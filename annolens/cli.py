"""annolens command line.

  annolens prepare  --manifest annolens.toml [--threads 8]     build the work dir (junctions, projections, homology, read sites)
  annolens scan     --manifest annolens.toml                   rule-based check of every gene + loci without a model
  annolens review   --manifest annolens.toml --queue P1        batch review of existing models with the configured model
  annolens addgene  --manifest annolens.toml [--queue missing_high]   batch review of loci without a model
  annolens mcp      --manifest annolens.toml                   MCP server (stdio)
  annolens install  --host claude|codex|generic --ssh HOST --manifest /remote/annolens.toml [--remote-cmd annolens]
  annolens skill-md                                            print SKILL.md
  annolens skill-files                                         print the whole skill folder as JSON (used by the npm installer)
"""
import argparse, json, sys
from pathlib import Path

from .config import Manifest


def main(argv=None):
    ap = argparse.ArgumentParser(prog='annolens'); sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('prepare', 'scan', 'review', 'addgene', 'mcp'):
        p = sub.add_parser(name); p.add_argument('--manifest', required=True); p.add_argument('--out')
        if name == 'prepare': p.add_argument('--threads', type=int, default=8)
        if name in ('review', 'addgene'):
            p.add_argument('--queue', default='P1' if name == 'review' else 'missing_high'); p.add_argument('--ids', help='file with ids (overrides --queue)')
            p.add_argument('--workers', type=int, default=4)
    sub.add_parser('skill-md'); sub.add_parser('skill-files')
    p = sub.add_parser('install')
    p.add_argument('--host', required=True, choices=['claude', 'codex', 'generic']); p.add_argument('--ssh', help='ssh host of the data server (omit for local)')
    p.add_argument('--manifest', required=True, help='manifest path on the server'); p.add_argument('--remote-cmd', default='annolens')
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
    M = Manifest(a.manifest); out = Path(a.out) if a.out else M.w('reviews')
    if a.cmd == 'prepare':
        from . import prepare; prepare.run(M, a.threads)
    elif a.cmd == 'scan':
        from . import scan; print(json.dumps(scan.run(M), indent=1))
    elif a.cmd in ('review', 'addgene'):
        from . import agent
        ids = Path(a.ids).read_text().split() if a.ids else M.w('queues', a.queue + '.txt').read_text().split()
        agent.run(M.path, 'gene' if a.cmd == 'review' else 'locus', ids, out, a.workers)
    elif a.cmd == 'mcp':
        from .mcp_server import serve; serve(M, out)


if __name__ == '__main__':
    sys.exit(main())
