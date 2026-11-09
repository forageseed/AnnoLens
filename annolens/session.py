"""Review sessions shared by the batch agent and the MCP server: one place for actions, budgets and machine rules.

GeneSession   an existing gene model (skill gene_v3.5): open, view, orfs, build, readcheck, decide
LocusSession  a locus without a model (skill locus_v3.6): open, view, orfs, build (transcript or projection), readcheck,
              decide
Every action returns {'text': <dict for the model>, 'image': <png path or None>}. Decisions are machine-checked and
written to <out>/<mode>/reviews/*.json; the genome and the input annotation are never edited.
"""
import csv, hashlib, json, time
from pathlib import Path

from . import genome as G
from .render import render

SKILLS = Path(__file__).with_name('skills')
SKILL = {'gene': (SKILLS / 'gene_v3.5.txt', '3.5'), 'locus': (SKILLS / 'locus_v3.6.txt', '3.6')}
BUDGET = {'gene': {'views': 4, 'builds': 2, 'readchecks': 3}, 'locus': {'views': 4, 'builds': 3, 'readchecks': 3}}
GENE_STATUSES = {'retain', 'recommend_candidate', 'unresolved', 'assembly_conflict', 'genome_frameshift'}
LOCUS_STATUSES = {'add_gene', 'no_gene', 'genome_frameshift', 'unresolved'}


def skill_text(mode): return SKILL[mode][0].read_text(encoding='utf-8')


class Workspace:
    """Lazily loaded engines for one project manifest."""
    def __init__(self, manifest, out):
        self.M, self.out = manifest, Path(out); self.out.mkdir(parents=True, exist_ok=True); self._g = self._l = self._rc = None

    def gene(self):
        if self._g is None:
            from .engine import GeneEngine
            from .tools import GeneTools
            e = GeneEngine(self.M, self.out / 'gene' / 'engine_out'); sp = self.M.w('readcheck', 'sites.jsonl')
            sites = {json.loads(l)['gene_id']: json.loads(l)['sites'] for l in open(sp)} if sp.exists() else {}
            self._g = (e, GeneTools(e), sites)
        return self._g

    def locus(self):
        if self._l is None:
            from . import addgene
            e = addgene.LocusEngine(addgene.config_from_manifest(self.M, self.out / 'locus' / 'engine_out'))
            lt = self.M.w('loci.tsv')
            if lt.exists():
                for r in csv.DictReader(open(lt), delimiter='\t'):
                    e.add_locus(r['locus_id'], r['chrom'], int(r['start_1based']) - 1, int(r['end_1based']), r['strand'])
            self._l = (e, addgene.LocusTools(e))
        return self._l

    def readcheck(self):
        if self._rc is None:
            from .readcheck import ReadCheck
            self._rc = ReadCheck(self.M.w('readcheck', 'config.json'))
        return self._rc


class _Session:
    mode = None

    def __init__(self, ws, key):
        self.ws, self.key = ws, key; self.used = {k: 0 for k in BUDGET[self.mode]}; self.n_img = 0
        self.builds, self.readchecks, self.trace, self.closed = {}, [], [], False

    def _budget(self, k):
        if self.closed: raise ValueError('session already decided')
        if self.used[k] >= BUDGET[self.mode][k]: raise ValueError(f'{k} budget used ({BUDGET[self.mode][k]})')
        self.used[k] += 1

    def _img(self, data):
        p = self.ws.out / self.mode / 'images' / f"{self.key.replace(':', '_')}.{self.n_img}.png"; p.parent.mkdir(parents=True, exist_ok=True)
        render(data, p); self.n_img += 1
        self.trace.append({'image': str(p), 'image_sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'view': [data['start0'] + 1, data['end0']]})
        return str(p)

    def readcheck(self, position, kind='any'):
        self._budget('readchecks'); r = self.ws.readcheck().check(self.chrom, int(position), kind)
        self.readchecks.append(r); self.trace.append({'readcheck': r}); return {'text': {'READCHECK': r}, 'image': None}

    def _write(self, rec):
        p = self.ws.out / self.mode / 'reviews' / f"{self.key.replace(':', '_')}.json"; p.parent.mkdir(parents=True, exist_ok=True)
        rec.update(trace=self.trace, budgets_used=self.used, skill_version=SKILL[self.mode][1],
                   policy_sha256=hashlib.sha256(skill_text(self.mode).encode()).hexdigest(), written=time.strftime('%Y-%m-%d %H:%M:%S'),
                   genome_edits=[], input_annotation_modified=False)
        p.write_text(json.dumps(rec, indent=1)); self.closed = True; return str(p)


class GeneSession(_Session):
    mode = 'gene'

    def open(self):
        e, T, sites = self.ws.gene(); m = e.model(self.key); self.chrom = m['seqid']
        self.checks = e.describe(m, m['exons'], m['cds']); self.sites = [dict(s) for s in sites.get(self.key, [])][:6]
        lo, hi = T.overview_region(self.key); d = T.view_data(self.key, lo, hi)
        facts = {'gene_id': self.key, 'seqid': m['seqid'], 'strand': m['strand'],
                 'current_model': {'exons': len(m['exons']), 'cds_span_1based': [m['cds'][0][0] + 1, m['cds'][-1][1]], 'aa': d['current']['aa'],
                                   'start_codon': self.checks['start_codon'], 'terminal_stop': self.checks['terminal_stop'],
                                   'internal_stops': self.checks['internal_stops'], 'coding_complete': e.complete(self.checks)},
                 'dominant_orf': None if not d['dominant_orf'] else {k: d['dominant_orf'][k] for k in ('aa', 'sources', 'adoptable')},
                 'overview_window_1based': [lo + 1, hi], 'budgets': BUDGET['gene'],
                 'read_checked_sites': [{k: s.get(k) for k in ('kind', 'position_1based', 'protein_pos', 'protein_length', 'context', 'hifi_reads',
                                                               'hifi_reads_disagreeing', 'rna_reads', 'rna_reads_disagreeing', 'call')} for s in self.sites]}
        return {'text': {'FACTS': facts, 'TRANSCRIPTS': T.transcript_table(d)}, 'image': self._img(d)}

    def view(self, start, end):
        self._budget('views'); e, T, _ = self.ws.gene(); d = T.view_data(self.key, int(start) - 1, int(end))
        return {'text': {'VIEW': [int(start), int(end)], 'TRANSCRIPTS': T.transcript_table(d)}, 'image': self._img(d)}

    def orfs(self, transcript):
        return {'text': {'ORFS': self.ws.gene()[1].orf_options(self.key, transcript)}, 'image': None}

    def build(self, transcript, start=None, merge=False):
        self._budget('builds'); e, T, _ = self.ws.gene()
        res = T.build_from_transcript(self.key, transcript, start, bool(merge)); self.trace.append({'build': res})
        if not res.get('built'): return {'text': {'BUILD': res}, 'image': None}
        self.builds[res['candidate_id']] = res
        d = T.view_data(self.key, *T.overview_region(self.key), candidate_id=res['candidate_id'])
        return {'text': {'BUILD': res}, 'image': self._img(d)}

    def decide(self, a):
        e, T, _ = self.ws.gene(); st = a.get('status')
        if self.closed: raise ValueError('session already decided')
        if st not in GENE_STATUSES: raise ValueError(f'status must be one of {sorted(GENE_STATUSES)}')
        if st == 'recommend_candidate' and not (self.builds.get(a.get('candidate_id')) or {}).get('materializable'):
            ok = [k for k, b in self.builds.items() if b.get('materializable')]
            raise ValueError(f'recommend_candidate needs a materializable candidate_id built in this session (available: {ok or "none - build first"})')
        final = {k: a.get(k) for k in ('status', 'candidate_id', 'confidence', 'rationale', 'evidence_gaps', 'visual_findings')}
        proposal = cds = None
        if st == 'recommend_candidate':
            v = e.validate_annotation_edit(self.key, final['candidate_id'])
            final['machine_validation'] = {k: v.get(k) for k in ('materializable_proposal', 'errors', 'evidence_gate', 'cross_tissue_support_sources')}
            if v['materializable_proposal']: proposal = e.materialize_proposal(self.key, final['candidate_id']); cds = e.candidate(self.key, final['candidate_id'])['new_cds']
            else: final.update(status='unresolved', machine_override='validator_refused:' + ','.join(v['errors']))
        if final['status'] == 'genome_frameshift' and not any(x.get('call') == 'genome_sequence_supported' for x in self.sites + self.readchecks):
            final.update(status='unresolved', machine_override='genome_frameshift_without_read_supported_site')
        if final['status'] == 'retain' and not e.complete(self.checks):
            final.update(status='unresolved', machine_override='retain_to_unresolved_for_incomplete_CDS')
        path = self._write({'gene_id': self.key, 'decision': final, 'emitted_cds': cds, 'proposal': proposal,
                            'read_checked_sites': self.sites, 'readchecks': self.readchecks})
        return {'text': {'FINAL': final, 'written': path}, 'image': None}


class LocusSession(_Session):
    mode = 'locus'

    def open(self):
        e, T = self.ws.locus()
        if self.key not in e.loci: raise ValueError('unknown locus_id')
        L = e.loci[self.key]; self.chrom = L['chrom']; lo, hi = T.region(self.key); d = T.view_data(self.key, lo, hi); t, h = T.tables(d)
        facts = {'locus_id': self.key, 'seqid': L['chrom'], 'strand': L['strand'], 'locus_1based': [L['s'] + 1, L['e']], 'current_model': None,
                 'overview_window_1based': [lo + 1, hi], 'budgets': BUDGET['locus']}
        return {'text': {'FACTS': facts, 'TRANSCRIPTS': t, 'PROJECTIONS': h}, 'image': self._img(d)}

    def view(self, start, end, candidate=None):
        self._budget('views'); e, T = self.ws.locus(); d = T.view_data(self.key, int(start) - 1, int(end), candidate=candidate); t, h = T.tables(d)
        return {'text': {'VIEW': [int(start), int(end)], 'TRANSCRIPTS': t, 'PROJECTIONS': h}, 'image': self._img(d)}

    def orfs(self, transcript):
        return {'text': {'ORFS': self.ws.locus()[1].orf_options(self.key, transcript)}, 'image': None}

    def build(self, transcript=None, projection=None, start=None):
        self._budget('builds'); e, T = self.ws.locus()
        res = T.build_from_projection(self.key, projection) if projection else T.build_from_transcript(self.key, transcript, start)
        short = {k: res[k] for k in ('candidate_id', 'materializable', 'errors', 'source')}
        short.update(protein_aa=res['detail']['protein_aa'], homology=res['detail']['homology'],
                     unsupported_introns=res['detail']['unsupported_introns'][:5], overlaps=res['detail']['overlaps'])
        self.builds[res['candidate_id']] = res; self.trace.append({'build': short})
        lo, hi = T.region(self.key)
        d = T.view_data(self.key, max(0, min(lo, res['exons'][0][0] - 500)), max(hi, res['exons'][-1][1] + 500), candidate=res['candidate_id'])
        return {'text': {'BUILD': short}, 'image': self._img(d)}

    def decide(self, a):
        e, T = self.ws.locus(); L = e.loci[self.key]; st = a.get('status')
        if self.closed: raise ValueError('session already decided')
        if st not in LOCUS_STATUSES: raise ValueError(f'status must be one of {sorted(LOCUS_STATUSES)}')
        final = {k: a.get(k) for k in ('status', 'candidate_id', 'confidence', 'rationale', 'evidence_gaps', 'visual_findings')}; new = None
        if st == 'add_gene':
            b = self.builds.get(final['candidate_id'])
            if not b or not b['materializable']:
                ok = [k for k, x in self.builds.items() if x['materializable']]
                raise ValueError(f'add_gene needs a materializable candidate_id built in this session (available: {ok or "none - build first"})')
            errs, det, prot = e.validate(self.key, b['exons'], b['cds'], b['source'])
            if errs: final.update(status='unresolved', machine_override='revalidation_failed:' + ','.join(errs))
            else: new = {'chrom': L['chrom'], 'strand': L['strand'], 'exons': b['exons'], 'cds': b['cds'], 'protein': prot, 'source': b['source'], 'homology': det['homology']}
        if final['status'] == 'genome_frameshift' and not any(x.get('call') == 'genome_sequence_supported' for x in self.readchecks):
            final.update(status='unresolved', machine_override='genome_frameshift_without_read_supported_site')
        path = self._write({'locus_id': self.key, 'decision': final, 'new_gene': new, 'readchecks': self.readchecks})
        return {'text': {'FINAL': final, 'written': path}, 'image': None}
