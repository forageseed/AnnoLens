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
SKILL = {'gene': (SKILLS / 'gene_v3.7.txt', '3.7'), 'locus': (SKILLS / 'locus_v3.6.txt', '3.6')}
BUDGET = {'gene': {'views': 4, 'builds': 3, 'readchecks': 3}, 'locus': {'views': 4, 'builds': 3, 'readchecks': 3}}
GENE_STATUSES = {'retain', 'recommend_candidate', 'split', 'unresolved', 'assembly_conflict', 'genome_frameshift'}
SPLIT_IN_SPAN, SPLIT_LEN_RATIO, FS_LONGER = 0.5, (0.8, 1.2), 1.1


def _in_span(cds, m):
    """Fraction of a CDS inside the current model's span (first to last CDS base)."""
    return G.overlap(cds, [[m['cds'][0][0], m['cds'][-1][1]]]) / max(1, sum(b - a for a, b in cds))
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
        self.t0 = time.time()   # wall time from opening the session to the written decision

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
        rec.update(trace=self.trace, budgets_used=self.used, skill_version=SKILL[self.mode][1], elapsed_s=round(time.time() - self.t0, 1),
                   policy_sha256=hashlib.sha256(skill_text(self.mode).encode()).hexdigest(), written=time.strftime('%Y-%m-%d %H:%M:%S'),
                   genome_edits=[], input_annotation_modified=False)
        p.write_text(json.dumps(rec, indent=1)); self.closed = True; return str(p)


class GeneSession(_Session):
    mode = 'gene'

    def open(self):
        e, T, sites = self.ws.gene(); m = e.model(self.key); self.chrom = m['seqid']
        self.checks = e.describe(m, m['exons'], m['cds']); self.sites = [dict(s) for s in sites.get(self.key, [])][:6]
        lo, hi = T.overview_region(self.key); d = T.view_data(self.key, lo, hi)
        self.split_locus = 'split:' + self.key
        le, lt = self.ws.locus(); le.add_locus(self.split_locus, m['seqid'], lo, hi, m['strand'], ignore={self.key})
        self.parts = {}; d['homology'] = lt.projection_rows(self.split_locus, lo, hi)
        facts = {'gene_id': self.key, 'seqid': m['seqid'], 'strand': m['strand'],
                 'current_model': {'exons': len(m['exons']), 'cds_span_1based': [m['cds'][0][0] + 1, m['cds'][-1][1]], 'aa': d['current']['aa'],
                                   'start_codon': self.checks['start_codon'], 'terminal_stop': self.checks['terminal_stop'],
                                   'internal_stops': self.checks['internal_stops'], 'coding_complete': e.complete(self.checks)},
                 'dominant_orf': None if not d['dominant_orf'] else {k: d['dominant_orf'][k] for k in ('aa', 'sources', 'adoptable')},
                 'overview_window_1based': [lo + 1, hi], 'budgets': BUDGET['gene'], 'scan_flags': self._scan_flags(),
                 'read_checked_sites': [{k: s.get(k) for k in ('kind', 'position_1based', 'protein_pos', 'protein_length', 'context', 'hifi_reads',
                                                               'hifi_reads_disagreeing', 'rna_reads', 'rna_reads_disagreeing', 'call')} for s in self.sites]}
        return {'text': {'FACTS': facts, 'TRANSCRIPTS': T.transcript_table(d), 'PROJECTIONS': self.ws.locus()[1].projection_table(d['homology'])}, 'image': self._img(d)}

    def _scan_flags(self):
        """Findings of `annolens scan` for this gene (rule-based hints, not decisions)."""
        if not hasattr(self.ws, '_scan'):
            p = self.ws.M.w('scan', 'gene_scan.tsv'); self.ws._scan = {}
            if p.exists():
                for r in csv.DictReader(open(p), delimiter='\t'): self.ws._scan[r['gene_id']] = [x for x in r['problems'].split(';') if x]
        return self.ws._scan.get(self.key, [])

    def view(self, start, end):
        self._budget('views'); e, T, _ = self.ws.gene(); d = T.view_data(self.key, int(start) - 1, int(end))
        lt = self.ws.locus()[1]; d['homology'] = lt.projection_rows(self.split_locus, int(start) - 1, int(end))
        return {'text': {'VIEW': [int(start), int(end)], 'TRANSCRIPTS': T.transcript_table(d), 'PROJECTIONS': lt.projection_table(d['homology'])}, 'image': self._img(d)}

    def split_part(self, transcript=None, projection=None, start=None):
        """Build ONE part of a split (a separate gene inside the current model) from a transcript ORF or a projection."""
        self._budget('builds'); e, T, _ = self.ws.gene(); le, lt = self.ws.locus(); m = e.model(self.key)
        if projection: res = lt.build_from_projection(self.split_locus, projection)
        else:
            t = T.labels.get(self.key, {}).get(transcript)
            if not t: raise ValueError('unknown transcript label; view first')
            res = lt.build_from_tx(self.split_locus, f'{t[0]}:{t[1]}', transcript, start)
        if _in_span(res['cds'], m) < SPLIT_IN_SPAN:
            res = {**res, 'errors': res['errors'] + ['part_outside_current_model'], 'materializable': False}
        self.parts[res['candidate_id']] = res
        short = {k: res[k] for k in ('candidate_id', 'materializable', 'errors', 'source')}
        short.update(protein_aa=res['detail']['protein_aa'], homology=res['detail']['homology'], unsupported_introns=res['detail']['unsupported_introns'][:5])
        self.trace.append({'split_part': short})
        lo, hi = T.overview_region(self.key); d = T.view_data(self.key, lo, hi); d['homology'] = lt.projection_rows(self.split_locus, lo, hi)
        d['candidate'] = {'candidate_id': res['candidate_id'], 'exons': res['exons'], 'cds': res['cds']}
        return {'text': {'SPLIT_PART': short}, 'image': self._img(d)}

    def _supported_as_one_gene(self, m):
        """A model supported as one gene must not be split: a full-chain long read, or an intact protein projection
        of matching length (0.8-1.2 x) covering >= 50 % of its CDS."""
        e, _, _ = self.ws.gene(); le, _ = self.ws.locus(); why = []
        if e.evidence(m, m['exons'], m['cds'], {'kind': 'original'})['full_length_chain_observed']: why.append('full_chain_long_read')
        aa = G.aa_len(m['cds']); tot = max(1, sum(b - a for a, b in m['cds']))
        for x in le.projections(m['seqid'], m['cds'][0][0], m['cds'][-1][1], m['strand']):
            if x['frameshift'] or x['stop'] or x['cov'] < .9 or not x['protein_aa']: continue
            if G.overlap(m['cds'], x['cds']) / tot >= .5 and SPLIT_LEN_RATIO[0] <= aa / x['protein_aa'] <= SPLIT_LEN_RATIO[1]:
                why.append('intact_homolog_of_full_length:' + x['protein']); break
        return why

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
        proposal = cds = None; split = None
        if st == 'split':
            ids = list(dict.fromkeys(a.get('parts') or []))
            missing = [i for i in ids if i not in self.parts or not self.parts[i]['materializable']]
            if len(ids) < 2 or missing:
                ok = [k for k, x in self.parts.items() if x['materializable']]
                raise ValueError(f'split needs >= 2 materializable parts built in this session (available: {ok or "none - use split_part"})')
            final['parts'] = ids; split = self._decide_split(final, ids)
        if st == 'recommend_candidate':
            v = e.validate_annotation_edit(self.key, final['candidate_id'])
            final['machine_validation'] = {k: v.get(k) for k in ('materializable_proposal', 'errors', 'evidence_gate', 'cross_tissue_support_sources')}
            if v['materializable_proposal']: proposal = e.materialize_proposal(self.key, final['candidate_id']); cds = e.candidate(self.key, final['candidate_id'])['new_cds']
            else: final.update(status='unresolved', machine_override='validator_refused:' + ','.join(v['errors']))
        if final['status'] == 'genome_frameshift':
            m = e.model(self.key); aa = G.aa_len(m['cds']); lo, hi = m['cds'][0][0] - 3000, m['cds'][-1][1] + 3000
            ok = [x for x in self.sites if x.get('call') == 'genome_sequence_supported' and lo <= x['position_1based'] <= hi
                  and (x.get('protein_length') or 0) > FS_LONGER * aa]
            ok += [x for x in self.readchecks if x.get('call') == 'genome_sequence_supported' and lo <= x['position_1based'] <= hi]
            if not ok: final.update(status='unresolved', machine_override='genome_frameshift_without_read_supported_site_of_a_longer_protein_in_this_model')
        if final['status'] == 'retain' and not e.complete(self.checks):
            final.update(status='unresolved', machine_override='retain_to_unresolved_for_incomplete_CDS')
        path = self._write({'gene_id': self.key, 'decision': final, 'emitted_cds': cds, 'proposal': proposal, 'split': split,
                            'read_checked_sites': self.sites, 'readchecks': self.readchecks})
        return {'text': {'FINAL': final, 'written': path}, 'image': None}


    def _decide_split(self, final, ids):
        """Machine checks of a split; writes one GFF3 with the parts as separate genes. Returns the split record."""
        e, _, _ = self.ws.gene(); le, _ = self.ws.locus(); m = e.model(self.key); errors = []
        parts = [self.parts[i] for i in ids]
        for i in ids:   # re-validate every part exactly as it will be written
            r = self.parts[i]; errs, det, prot = le.validate(self.split_locus, r['exons'], r['cds'], r['source'])
            if _in_span(r['cds'], m) < SPLIT_IN_SPAN: errs = errs + ['part_outside_current_model']
            if errs: errors.append(f'{i}:' + ','.join(errs))
        for x in range(len(parts)):
            for y in range(x + 1, len(parts)):
                if G.overlap(parts[x]['exons'], parts[y]['exons']) or G.overlap(parts[x]['cds'], parts[y]['cds']): errors.append('parts_overlap')
        one = self._supported_as_one_gene(m)
        if one: errors.append('current_model_supported_as_one_gene:' + ';'.join(one))
        if errors:
            final.update(status='unresolved', machine_override='split_refused:' + ' | '.join(sorted(set(errors))))
            return {'parts': ids, 'errors': errors}
        parts = sorted(parts, key=lambda r: r['cds'][0][0] if m['strand'] == '+' else -r['cds'][-1][1])
        d = self.ws.out / 'gene' / 'proposals'; d.mkdir(parents=True, exist_ok=True)
        gff = d / f"{G.sha(self.key)[:16]}.split.gff3"
        rows = ['##gff-version 3', f'# AnnoLens split of {self.key}; genome unchanged; terminal stop included in CDS']
        R = lambda typ, a, b, ph, at: rows.append('\t'.join(map(str, [m['seqid'], 'AnnoLens', typ, a + 1, b, '.', m['strand'], ph, at])))
        for k, r in enumerate(parts, 1):
            gid, tid = f'{self.key}.s{k}', f'{m["transcript_id"]}.s{k}'
            R('gene', r['exons'][0][0], r['exons'][-1][1], '.', f'ID={gid};split_from={self.key}')
            R('mRNA', r['exons'][0][0], r['exons'][-1][1], '.', f'ID={tid};Parent={gid}')
            for i, (a, b) in enumerate(r['exons'], 1): R('exon', a, b, '.', f'ID={tid}.exon{i};Parent={tid}')
            for i, ((a, b), ph) in enumerate(zip(r['cds'], G.phases(r['cds'], m['strand'])), 1): R('CDS', a, b, ph, f'ID={tid}.cds{i};Parent={tid}')
        gff.write_text('\n'.join(rows) + '\n')
        return {'parts': ids, 'genes': [{'gene_id': f'{self.key}.s{k}', 'exons': r['exons'], 'cds': r['cds'], 'protein': r['protein'],
                                         'source': r['source'], 'homology': r['detail']['homology']} for k, r in enumerate(parts, 1)],
                'proposed_gff': str(gff)}


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
