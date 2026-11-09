#!/usr/bin/env python3
"""AnnoLens add-gene: locus engine + tools for loci where the annotation has no gene model.

A locus (from screen_v6: unannotated_locus_*) has no current model. Candidates come only from evidence:
  transcript   the complete ORF of an observed transcript (StringTie of 7 tissues or IsoSeq); exons = transcript exons
  projection   an intact protein projection (miniprot, independent annotation or related species); exons = CDS
Validation of a new gene (all must hold, otherwise not materializable):
  complete ORF (ATG .. terminal stop, length % 3 == 0, no internal stop, no N); canonical splice sites (GT-AG,
  GC-AG, AT-AC) for every intron; no overlap >= 30 bp with any existing CDS on either strand; protein >= MIN_AA;
  homology: an intact alignment (frameshift 0, stop 0, coverage >= 0.9) overlapping >= 50 % of the candidate CDS
  with candidate length 0.8-1.2 x that protein; introns: every CDS intron has >= 3 RNA junction reads or equals an
  intron of an intact projection (independent annotation or related species; for intron existence, unlike model
  length, the related-species proteome is valid evidence); TE: < 50 % of the CDS inside intact LTR retrotransposons.
All coordinates are genomic 0-based half-open. Nothing here edits the genome or the release annotation.
"""
import bisect, collections, hashlib, json, re
from pathlib import Path

from .config import KIND_LABEL

STOPS = {'TAA', 'TAG', 'TGA'}
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
rc = lambda s: s.translate(COMP)[::-1]
MIN_AA, OV_BP, MIN_JUNC, COV, LEN_RATIO, LTR_MAX = 50, 30, 3, 0.9, (0.8, 1.2), 0.5


class Idx:
    """Sorted interval index: items (a, b, payload)."""
    def __init__(self, items):
        self.v = sorted(items, key=lambda x: (x[0], x[1])); self.s = [x[0] for x in self.v]; self.m = max((x[1] - x[0] for x in self.v), default=0)
    def q(self, a, b):
        i = bisect.bisect_left(self.s, a - self.m)
        while i < len(self.v) and self.v[i][0] < b:
            if self.v[i][1] > a: yield self.v[i]
            i += 1


def intr(bl): return [(bl[i][1], bl[i + 1][0]) for i in range(len(bl) - 1)]
def ov(x, y): return sum(max(0, min(b, d) - max(a, c)) for a, b in x for c, d in y)


class LocusEngine:
    def __init__(self, cfg):
        cfg = json.loads(Path(cfg).read_text()) if not isinstance(cfg, dict) else cfg
        self.cfg = cfg
        import sys
        if cfg.get('pysam_dir'): sys.path.insert(0, cfg['pysam_dir'])
        import pysam
        self.fa = pysam.FastaFile(cfg['genome'])
        # existing models (CDS) on both strands
        cds = collections.defaultdict(list); tx2g = {}
        for line in open(cfg['gff']):
            if line.startswith('#') or not line.strip(): continue
            f = line.rstrip('\n').split('\t'); at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
            if f[2] == 'mRNA': tx2g[at['ID']] = at['Parent']
            elif f[2] == 'CDS': cds[f[0]].append((int(f[3]) - 1, int(f[4]), (tx2g[at['Parent'].split(',')[0]], f[6])))
        self.cds = {c: Idx(v) for c, v in cds.items()}
        # transcripts
        tx = {}
        def gtf(path, tag):
            for line in open(path):
                if line.startswith('#') or not line.strip(): continue
                f = line.rstrip('\n').split('\t'); at = dict(re.findall(r'(\S+) "([^"]*)"', f[8])); t = at.get('transcript_id')
                if not t: continue
                x = tx.setdefault(tag + ':' + t, {'src': tag, 'chrom': f[0], 'strand': f[6], 'ex': [], 'cov': 0.0})
                if f[2] == 'exon': x['ex'].append((int(f[3]) - 1, int(f[4])))
                elif f[2] == 'transcript': x['cov'] = float(at.get('cov', 0) or 0)
        for t, p in cfg['transcripts'].items():
            if Path(p).exists(): gtf(p, t)
        self.longread_name = cfg.get('longread_name')
        by = collections.defaultdict(list)
        for k, x in tx.items():
            x['ex'].sort(); x['id'] = k
            if x['ex']: by[x['chrom']].append((x['ex'][0][0], x['ex'][-1][1], k))
        self.tx = tx; self.txi = {c: Idx(v) for c, v in by.items()}
        # junctions (pooled over tissues)
        j = collections.Counter()
        for t, p in cfg['junction_beds'].items():
            for line in open(p):
                c = line.split(); j[(c[0], int(c[1]), int(c[2]), c[5])] += int(c[4])
        self.junc = j
        jb = collections.defaultdict(list)
        for (c, a, b, s), n in j.items(): jb[c].append((a, b, (s, n)))
        self.jidx = {c: Idx(v) for c, v in jb.items()}
        # protein projections
        self.aln = {}
        for kind, (gff, faa) in cfg['projections'].items():
            L = {}; cur = None
            for line in open(faa):
                if line[0] == '>': cur = line[1:].split()[0].split('|')[-1]; L[cur] = 0
                else: L[cur] += len(line.strip())
            per = collections.defaultdict(list); cur = None
            for line in open(gff):
                if line.startswith('#') or not line.strip(): continue
                f = line.rstrip('\n').split('\t')
                if f[2] == 'mRNA':
                    at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x); t = at['Target'].split(); q = t[0].split('|')[-1]
                    cur = {'kind': kind, 'protein': q, 'protein_aa': L.get(q, 0), 'chrom': f[0], 'strand': f[6],
                           'identity': float(at.get('Identity', 0)), 'positive': float(at.get('Positive', 0)),
                           'frameshift': int(at.get('Frameshift', 0)), 'stop': int(at.get('StopCodon', 0)),
                           'cov': (int(t[2]) - int(t[1]) + 1) / max(1, L.get(q, 0)), 'cds': [], 'id': kind + ':' + at['ID']}
                    per[f[0]].append(cur)
                elif f[2] == 'CDS' and cur is not None: cur['cds'].append((int(f[3]) - 1, int(f[4])))
            for c, v in per.items():
                for x in v: x['cds'].sort()
                self.aln.setdefault(c, []).extend(v)
        self.alni = {c: Idx([(x['cds'][0][0], x['cds'][-1][1], x) for x in v if x['cds']]) for c, v in self.aln.items()}
        ltr = collections.defaultdict(list)
        smap = dict(l.rstrip('\n').split('\t')[:2] for l in open(cfg['ltr_seqid_map'])) if cfg.get('ltr_seqid_map') else {}
        for line in (open(cfg['ltr_intact_gff3']) if cfg.get('ltr_intact_gff3') else []):
            if line.startswith('#') or not line.strip(): continue
            f = line.split('\t')
            # any repeat feature except structural sub-parts counts as repeat sequence
            if len(f) > 4 and f[2] not in ('long_terminal_repeat', 'target_site_duplication', 'repeat_region', 'primer_binding_site'):
                ltr[smap.get(f[0], f[0])].append((int(f[3]) - 1, int(f[4]), None))
        self.ltr = {c: Idx(v) for c, v in ltr.items()}
        self.loci = {}
        self.out = Path(cfg.get('out', '.'))

    # ---------- loci ----------
    def add_locus(self, locus_id, chrom, s, e, strand):
        self.loci[locus_id] = {'id': locus_id, 'chrom': chrom, 's': s, 'e': e, 'strand': strand}

    def seq(self, chrom, blocks, strand):
        s = ''.join(self.fa.fetch(chrom, a, b).upper() for a, b in blocks)
        return s if strand == '+' else rc(s)

    def transcripts(self, chrom, a, b, strand):
        return [self.tx[k] for x, y, k in self.txi.get(chrom, Idx([])).q(a, b) if self.tx[k]['strand'] in (strand, '.')]

    def projections(self, chrom, a, b, strand):
        return [x for p, q, x in self.alni.get(chrom, Idx([])).q(a, b) if x['strand'] == strand]

    def models(self, chrom, a, b):
        out = collections.defaultdict(list)
        for x, y, (g, st) in self.cds.get(chrom, Idx([])).q(a, b): out[(g, st)].append([x, y])
        return out

    def junctions(self, chrom, a, b, strand):
        return [[x, y, n] for x, y, (s, n) in self.jidx.get(chrom, Idx([])).q(a, b) if s in (strand, '.') and x >= a and y <= b]

    # ---------- ORFs ----------
    def orfs(self, chrom, exons, strand, min_aa=MIN_AA):
        """Complete ORFs (ATG..stop) of a spliced transcript; returns [{'aa', 'cds'}] longest first."""
        mrna = self.seq(chrom, exons, strand); pos = []   # mRNA offset -> genomic base (0-based)
        for a, b in (exons if strand == '+' else exons[::-1]):
            pos += list(range(a, b)) if strand == '+' else list(range(b - 1, a - 1, -1))
        res = []
        for fr in range(3):
            start = None
            for i in range(fr, len(mrna) - 2, 3):
                cod = mrna[i:i + 3]
                if start is None and cod == 'ATG': start = i
                if cod in STOPS:
                    if start is not None and (i - start) // 3 >= min_aa:
                        g = sorted(pos[start:i + 3]); blocks = []
                        for x in g:
                            if blocks and x == blocks[-1][1]: blocks[-1][1] += 1
                            else: blocks.append([x, x + 1])
                        res.append({'aa': (i - start) // 3, 'cds': blocks})
                    start = None
        return sorted(res, key=lambda r: -r['aa'])

    # ---------- validation ----------
    def validate(self, locus, exons, cds, source):
        L = self.loci[locus]; c, st = L['chrom'], L['strand']; errors = []; detail = {}
        s = self.seq(c, cds, st); aa = len(s) // 3 - 1
        prot = ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) // 3 * 3, 3))
        if s[:3] != 'ATG': errors.append('no_ATG_start')
        if s[-3:] not in STOPS: errors.append('no_terminal_stop')
        if len(s) % 3: errors.append('cds_not_mod3')
        if '*' in prot[:-1]: errors.append('internal_stop')
        if 'N' in s: errors.append('ambiguous_bases')
        if aa < MIN_AA: errors.append('protein_too_short')
        bad = []
        for a, b in intr(exons):
            m = self.fa.fetch(c, a, a + 2).upper() + '-' + self.fa.fetch(c, b - 2, b).upper()
            if st == '-': m = rc(self.fa.fetch(c, b - 2, b).upper()) + '-' + rc(self.fa.fetch(c, a, a + 2).upper())
            if m not in ('GT-AG', 'GC-AG', 'AT-AC'): bad.append([a, b, m])
        if bad: errors.append('noncanonical_splice')
        hits = collections.Counter()
        for x, y, (g, gst) in self.cds.get(c, Idx([])).q(cds[0][0], cds[-1][1]):
            hits[(g, gst)] += ov(cds, [(x, y)])
        clash = {f'{g}({gst})': n for (g, gst), n in hits.items() if n >= OV_BP}
        if clash: errors.append('overlaps_existing_CDS')
        # homology support
        hom = []
        for x in self.projections(c, cds[0][0], cds[-1][1], st):
            if x['frameshift'] or x['stop'] or x['cov'] < COV or not x['protein_aa']: continue
            o = ov(cds, x['cds']) / max(1, sum(b - a for a, b in cds)); r = aa / x['protein_aa']
            if o >= .5 and LEN_RATIO[0] <= r <= LEN_RATIO[1]: hom.append({'protein': x['kind'] + '|' + x['protein'], 'ratio': round(r, 3), 'overlap': round(o, 3)})
        if not hom: errors.append('no_intact_homolog_of_matching_length')
        proj_introns = {tuple(i) for x in self.projections(c, cds[0][0], cds[-1][1], st)
                        if not x['frameshift'] and not x['stop'] and x['cov'] >= COV for i in intr(x['cds'])}
        unsup = []
        for a, b in intr(cds):
            n = self.junc.get((c, a, b, st), 0) + self.junc.get((c, a, b, '.'), 0)
            if n < MIN_JUNC and (a, b) not in proj_introns: unsup.append([a, b, n])
        if unsup: errors.append('intron_without_RNA_or_intact_projection')
        tot = sum(b - a for a, b in cds)
        lt = sum(ov(cds, [(x, y)]) for x, y, _ in self.ltr.get(c, Idx([])).q(cds[0][0], cds[-1][1])) / max(1, tot)
        if lt >= LTR_MAX: errors.append('cds_inside_repeat')
        if not (cds[-1][1] > L['s'] and cds[0][0] < L['e']): errors.append('outside_locus')
        detail.update(protein_aa=aa, homology=hom[:3], unsupported_introns=unsup, noncanonical=bad, overlaps=clash,
                      ltr_fraction=round(lt, 2), source=source)
        return errors, detail, prot.rstrip('*')

    def build(self, locus, exons, cds, source):
        exons = [list(x) for x in sorted(exons)]; cds = [list(x) for x in sorted(cds)]
        errors, detail, prot = self.validate(locus, exons, cds, source)
        cid = 'new_' + hashlib.sha256(json.dumps([locus, exons, cds]).encode()).hexdigest()[:20]
        rec = {'candidate_id': cid, 'locus': locus, 'exons': exons, 'cds': cds, 'source': source, 'errors': errors,
               'materializable': not errors, 'detail': detail, 'protein': prot, 'genome_edits': []}
        p = self.out / 'candidates'; p.mkdir(parents=True, exist_ok=True)
        (p / f'{cid}.json').write_text(json.dumps(rec))
        return rec


class LocusTools:
    """Agent-facing tools; transcript labels T1..Tn and projection labels H1..Hn per locus."""
    def __init__(self, e): self.e = e; self.labels = {}; self.hlabels = {}

    def region(self, locus, pad=3000):
        L = self.e.loci[locus]; return max(0, L['s'] - pad), L['e'] + pad

    def view_data(self, locus, a, b, candidate=None):
        L = self.e.loci[locus]; c, st = L['chrom'], L['strand']
        if b - a > 200000 or b <= a: raise ValueError('view must be 1..200000 bp')
        txs = sorted(self.e.transcripts(c, a, b, st), key=lambda x: (x['ex'][0][0], x['id']))
        lab = self.labels.setdefault(locus, {})
        for x in txs:
            if x['id'] not in lab.values(): lab[f'T{len(lab) + 1}'] = x['id']
        inv = {v: k for k, v in lab.items()}
        trs = []
        for x in txs:
            o = self.e.orfs(c, x['ex'], st)
            trs.append({'label': inv[x['id']], 'source': ('longread:' if x['src'] == self.e.longread_name else 'rna:') + x['src'],
                        'cov': round(x['cov'], 1), 'exons': [list(e) for e in x['ex']],
                        'orf_aa': o[0]['aa'] if o else 0, 'orf_cds': o[0]['cds'] if o else []})
        hl = self.hlabels.setdefault(locus, {}); hom = []
        for x in sorted(self.e.projections(c, a, b, st), key=lambda x: -x['cov'] * x['identity'])[:8]:
            if x['id'] not in hl.values(): hl[f'H{len(hl) + 1}'] = x['id']
            k = {v: k for k, v in hl.items()}[x['id']]
            hom.append({'label': k, 'kind': KIND_LABEL.get(x['kind'].split(':')[0], x['kind']),
                        'protein': x['protein'], 'protein_aa': x['protein_aa'], 'identity': x['identity'], 'coverage': round(x['cov'], 3),
                        'frameshift': x['frameshift'], 'stop_codons': x['stop'], 'blocks': [list(y) for y in x['cds']]})
        others = {f'{g}({gst})': sorted(v) for (g, gst), v in self.e.models(c, a, b).items()}
        cand = None
        if candidate:
            r = json.loads((self.e.out / 'candidates' / f'{candidate}.json').read_text())
            cand = {'candidate_id': candidate, 'exons': r['exons'], 'cds': r['cds']}
        return {'gene_id': 'NEW.' + c + '.' + str(L['s'] + 1), 'seqid': c, 'strand': st, 'start0': a, 'end0': b,
                'current': {'exons': [], 'cds': [], 'aa': 0}, 'other_models': others, 'transcripts': trs,
                'junctions': [j for j in self.e.junctions(c, a, b, st) if j[2] >= MIN_JUNC], 'homology': hom,
                'dominant_orf': None, 'candidate': cand}

    def tables(self, data):
        t = [{'id': x['label'], 'source': x['source'], 'cov': x['cov'], 'span_1based': [x['exons'][0][0] + 1, x['exons'][-1][1]],
              'exons': len(x['exons']), 'longest_ORF_aa': x['orf_aa']} for x in data['transcripts']]
        h = [{'id': x['label'], 'kind': x['kind'], 'protein_aa': x['protein_aa'], 'identity': x['identity'], 'coverage': x['coverage'],
              'frameshift': x['frameshift'], 'stop_codons': x['stop_codons'],
              'span_1based': [x['blocks'][0][0] + 1, x['blocks'][-1][1]], 'cds_exons': len(x['blocks'])} for x in data['homology']]
        return t, h

    def orf_options(self, locus, label, top=6):
        L = self.e.loci[locus]; x = self.e.tx[self.labels[locus][label]]
        return [{'aa': o['aa'], 'start_1based': (o['cds'][0][0] + 1) if L['strand'] == '+' else o['cds'][-1][1],
                 'span_1based': [o['cds'][0][0] + 1, o['cds'][-1][1]]} for o in self.e.orfs(L['chrom'], x['ex'], L['strand'])[:top]]

    def build_from_transcript(self, locus, label, start_1based=None):
        L = self.e.loci[locus]
        if label not in self.labels.get(locus, {}): raise ValueError(f'unknown transcript label {label}; view first')
        x = self.e.tx[self.labels[locus][label]]
        os_ = self.e.orfs(L['chrom'], x['ex'], L['strand'])
        if start_1based:
            os_ = [o for o in os_ if ((o['cds'][0][0] + 1) if L['strand'] == '+' else o['cds'][-1][1]) == int(start_1based)]
            if not os_: raise ValueError('no complete ORF starts at that position; use orfs to list them')
        if not os_: raise ValueError('transcript has no complete ORF >= 50 aa')
        # exons: transcript exons (UTRs included) clipped so that they contain the chosen CDS
        return self.e.build(locus, x['ex'], os_[0]['cds'], {'kind': 'transcript', 'label': label, 'transcript': x['id']})

    def build_from_projection(self, locus, label):
        L = self.e.loci[locus]
        if label not in self.hlabels.get(locus, {}): raise ValueError(f'unknown homology label {label}; view first')
        x = next(y for y in self.e.aln[L['chrom']] if y['id'] == self.hlabels[locus][label])
        cds = [list(b) for b in x['cds']]
        # miniprot CDS may end before the stop codon: include it when the next codon is a stop
        if L['strand'] == '+':
            if self.e.fa.fetch(L['chrom'], cds[-1][1], cds[-1][1] + 3).upper() in STOPS: cds[-1][1] += 3
        elif cds[0][0] >= 3 and rc(self.e.fa.fetch(L['chrom'], cds[0][0] - 3, cds[0][0]).upper()) in STOPS: cds[0][0] -= 3
        return self.e.build(locus, cds, cds, {'kind': 'projection', 'label': label, 'protein': x['kind'] + '|' + x['protein']})


def config_from_manifest(M, out):
    """Locus-engine configuration from the project manifest (no site-specific paths in code)."""
    return {'genome': str(M.genome), 'gff': str(M.annotation), 'transcripts': {k: str(v) for k, v in M.sources.items()},
            'longread_name': M.longread['name'] if M.longread else None,
            'junction_beds': {p.stem: str(p) for p in sorted(M.w('junctions').glob('*.bed'))},
            'projections': {p['kind'] + ':' + p['name']: [str(M.w('proteins', p['name'] + '.miniprot.gff')), str(p['fasta'])] for p in M.proteins},
            'ltr_intact_gff3': str(M.repeats) if M.repeats else None, 'ltr_seqid_map': str(M.repeats_map) if M.repeats_map else None,
            'pysam_dir': str(M.w('deps')), 'out': str(out)}
