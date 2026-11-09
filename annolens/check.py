"""`annolens check`: is this machine ready to run AnnoLens locally?

Checks, in order: Python version, operating system, Python dependencies, external programs, the project manifest and
its input files, the prepared work dir, and the model API key (only whether the variable is set; the value is never
printed). Exit code 0 when nothing blocks the next step.
"""
import importlib, os, platform, shutil, sys
from pathlib import Path

from . import __version__


def find_manifest(arg=None):
    """Project manifest: explicit argument, then ANNOLENS_MANIFEST, then ./annolens.toml."""
    for c in (arg, os.environ.get('ANNOLENS_MANIFEST'), 'annolens.toml'):
        if c and Path(c).expanduser().is_file(): return Path(c).expanduser().resolve()
    return None


def run(manifest_arg=None):
    rows, blocking = [], 0

    def row(ok, what, detail='', block=True):
        nonlocal blocking
        mark = 'ok  ' if ok else ('FAIL' if block else 'warn')
        if not ok and block: blocking += 1
        rows.append(f'[{mark}] {what}' + (f': {detail}' if detail else ''))

    row(sys.version_info >= (3, 11), 'Python >= 3.11', platform.python_version())
    unix = os.name == 'posix'
    row(unix, 'operating system', platform.system() + ('' if unix else ' (use Linux, macOS or WSL on Windows)'))
    for mod, why in (('pysam', 'reads BAM files'), ('matplotlib', 'draws evidence images')):
        try: importlib.import_module(mod); row(True, f'python module {mod}')
        except Exception as ex: row(False, f'python module {mod}', f'{why}; pip install {mod} ({type(ex).__name__})')

    mp = find_manifest(manifest_arg)
    if mp is None:
        row(False, 'project manifest', 'none found: pass --manifest, set ANNOLENS_MANIFEST, or put annolens.toml in this '
                                       'directory (template: examples/annolens.toml)')
    else:
        try:
            from .config import Manifest
            M = Manifest(mp); row(True, 'project manifest', str(mp))
            try: M.check(); row(True, 'input files')
            except Exception as ex: row(False, 'input files', str(ex))
            for tool, why in (('miniprot', 'protein projections (prepare)'), ('samtools', 'indexing (prepare)'),
                              ('stringtie', 'only when a sample has a BAM but no GTF')):
                if tool == 'stringtie' and not any(r['bam'] and not r['gtf'] for r in M.rna): continue
                exe = M.tools[tool]; found = shutil.which(exe) or Path(exe).is_file()
                row(bool(found), f'program {tool}', why if not found else exe, block=False)
            prepared = M.w('homology.jsonl').exists()
            row(prepared, 'work dir prepared', str(M.work) if prepared else f'run: annolens prepare --manifest {mp}', block=False)
            key = M.model.get('api_key_env')
            if key: row(bool(os.environ.get(key)), f'model key {key}', 'set' if os.environ.get(key) else
                        'not set (only needed for batch review; in a chat host the host model does the review)', block=False)
        except Exception as ex:
            row(False, 'project manifest', f'{mp}: {ex}')

    print(f'AnnoLens {__version__}')
    print('\n'.join(rows))
    print('ready' if not blocking else f'{blocking} blocking problem(s)')
    return 0 if not blocking else 1
