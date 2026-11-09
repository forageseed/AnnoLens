#!/usr/bin/env python3
"""Write skills/annolens/SKILL.md from the package's skill texts (single source), or check it is up to date.
  python scripts/build_skill.py          # write
  python scripts/build_skill.py --check  # exit 1 if skills/annolens/SKILL.md differs from the package"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from annolens.install import skill_md

out = ROOT / 'skills' / 'annolens' / 'SKILL.md'; text = skill_md()
if '--check' in sys.argv:
    ok = out.exists() and out.read_text(encoding='utf-8') == text
    print('skills/annolens/SKILL.md', 'up to date' if ok else 'OUT OF DATE: run python scripts/build_skill.py'); sys.exit(0 if ok else 1)
out.parent.mkdir(parents=True, exist_ok=True); out.write_text(text, encoding='utf-8', newline='\n'); print('wrote', out.relative_to(ROOT))
