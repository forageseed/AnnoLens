#!/usr/bin/env python3
"""Release gate: fail if the repository contains site-, person- or project-specific strings or anything that looks
like a credential. Run before every push:  python scripts/check_sensitive.py
Extra patterns (one regex per line) can be kept locally in .sensitive-patterns (git-ignored)."""
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = [
    r'/home/[a-z]', r'[A-Z]:[\\/](work|Users)', r'\b10\.\d+\.\d+\.\d+\b', r'\b192\.168\.\d+\.\d+\b',   # paths, private IPs
    r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[a-z]{2,}',                                               # e-mail addresses
    r'\b(sk|pk|rk)-[A-Za-z0-9]{16,}', r'\b[0-9a-f]{32}\.[A-Za-z0-9]{16}\b', r'AKIA[0-9A-Z]{16}',     # API keys / tokens
    r'(?i)(api[_-]?key|secret|token|password)\s*[=:]\s*["\'][^"\']{8,}',
]
local = ROOT / '.sensitive-patterns'
if local.exists(): PATTERNS += [l.strip() for l in local.read_text().splitlines() if l.strip() and not l.startswith('#')]
ALLOW = {r'noreply@anthropic\.com'}
bad = 0
for p in ROOT.rglob('*'):
    if p.is_dir() or '.git' in p.parts or p.suffix in ('.png', '.jpg', '.ico') or p.name == '.sensitive-patterns': continue
    try: text = p.read_text(encoding='utf-8')
    except UnicodeDecodeError: continue
    for n, line in enumerate(text.splitlines(), 1):
        for pat in PATTERNS:
            m = re.search(pat, line)
            if m and not any(re.search(a, m.group(0)) for a in ALLOW):
                print(f'{p.relative_to(ROOT)}:{n}: {pat}: {line.strip()[:120]}'); bad += 1
print('sensitive strings:', bad)
sys.exit(1 if bad else 0)
