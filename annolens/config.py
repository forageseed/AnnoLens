"""Project manifest (annolens.toml). Everything species- or site-specific lives here; the code has no paths.

[project]   name, work_dir
[genome]    fasta (indexed .fai), annotation (GFF3 with gene/mRNA/exon/CDS)
[[rna]]     any number of short-read samples: name, bam, gtf (optional; built with StringTie when missing)
[longread]  optional: name, bam, gtf, index (full-length transcripts, e.g. IsoSeq / ONT cDNA; index = .bai/.csi stored elsewhere)
[[proteins]] any number of protein sets: name, fasta, role = independent | related | paralog
             independent  an annotation of the same species made independently of the one under review: anchors
                          start/stop/length and is the reference for truncation
             related      proteome of a related species: anchors start/stop/length
             paralog      proteins from the same pipeline (homeologs/paralogs): support introns only
[hifi]      optional: bam (genomic long reads, for single-base read checks)
[repeats]   optional: gff3 (TE annotation), seqid_map (two-column map when the GFF uses renamed sequences)
[tools]     miniprot, stringtie, samtools (executables; default: found on PATH)
[model]     provider = anthropic | openai, base_url, model, api_key_env (name of the environment variable)
"""
import tomllib
from pathlib import Path

ROLES = {'independent': 'ind', 'related': 'pot', 'paralog': 'hom'}     # internal kind codes used by the engine
KIND_LABEL = {'ind': 'independent_annotation', 'pot': 'related_species', 'hom': 'homeolog_or_paralog'}


class Manifest:
    def __init__(self, path):
        self.path = Path(path).resolve(); base = self.path.parent
        d = tomllib.loads(self.path.read_text(encoding='utf-8'))
        P = lambda p: (base / p).resolve() if p and not Path(p).is_absolute() else (Path(p) if p else None)
        self.name = d.get('project', {}).get('name', 'project')
        self.work = P(d.get('project', {}).get('work_dir', 'annolens_work'))
        g = d.get('genome', {})
        if not g.get('fasta') or not g.get('annotation'): raise ValueError('[genome] needs fasta and annotation')
        self.genome, self.annotation = P(g['fasta']), P(g['annotation'])
        self.rna = [{'name': r['name'], 'bam': P(r.get('bam')), 'gtf': P(r.get('gtf'))} for r in d.get('rna', [])]
        names = [r['name'] for r in self.rna]
        if len(set(names)) != len(names): raise ValueError('[[rna]] names must be unique')
        lr = d.get('longread')
        self.longread = {'name': lr.get('name', 'longread'), 'bam': P(lr.get('bam')), 'gtf': P(lr.get('gtf')), 'index': P(lr.get('index'))} if lr else None
        self.proteins = []
        for p in d.get('proteins', []):
            if p.get('role') not in ROLES: raise ValueError(f"protein set {p.get('name')}: role must be one of {sorted(ROLES)}")
            self.proteins.append({'name': p['name'], 'fasta': P(p['fasta']), 'role': p['role'], 'kind': ROLES[p['role']]})
        self.hifi_bam = P(d.get('hifi', {}).get('bam'))
        rep = d.get('repeats', {})
        self.repeats, self.repeats_map = P(rep.get('gff3')), P(rep.get('seqid_map'))
        t = d.get('tools', {})
        self.tools = {k: t.get(k, k) for k in ('miniprot', 'stringtie', 'samtools')}
        self.model = d.get('model', {})

    # derived paths inside the work dir
    def w(self, *p): return self.work.joinpath(*p)

    @property
    def sources(self):
        """Transcript assemblies by source name (each sample is one source; the long-read assembly is one source)."""
        out = {r['name']: self.w('transcripts', f"{r['name']}.gtf") for r in self.rna}
        if self.longread: out[self.longread['name']] = self.w('transcripts', f"{self.longread['name']}.gtf")
        return out

    def check(self):
        missing = [str(p) for p in [self.genome, self.annotation, Path(str(self.genome) + '.fai')] if not p.exists()]
        for r in self.rna:
            if not (r['bam'] and r['bam'].exists()) and not (r['gtf'] and r['gtf'].exists()): missing.append(f"rna {r['name']}: bam or gtf")
        for p in self.proteins:
            if not p['fasta'].exists(): missing.append(str(p['fasta']))
        if missing: raise FileNotFoundError('missing inputs: ' + '; '.join(missing))
        if not self.rna and not self.longread: raise ValueError('at least one [[rna]] sample or [longread] is required')
        return True
