#!/usr/bin/env python3
"""Write the skill folder skills/annolens/ from the package (single source), or check that it is up to date.
  python scripts/build_skill.py          # write
  python scripts/build_skill.py --check  # exit 1 if skills/annolens/ differs from the package"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from annolens.install import skill_files

D = ROOT / 'skills' / 'annolens'; files = skill_files()
have = {p.relative_to(D).as_posix() for p in D.rglob('*') if p.is_file()} if D.exists() else set()
if '--check' in sys.argv:
    bad = [r for r, t in files.items() if not (D / r).exists() or (D / r).read_text(encoding='utf-8') != t] + sorted(have - set(files))
    print('skills/annolens/', 'up to date' if not bad else 'OUT OF DATE: ' + ', '.join(bad)); sys.exit(1 if bad else 0)
for r in sorted(have - set(files)): (D / r).unlink(); print('removed', r)
for r, t in files.items():
    (D / r).parent.mkdir(parents=True, exist_ok=True); (D / r).write_text(t, encoding='utf-8', newline='\n'); print('wrote', (D / r).relative_to(ROOT))
