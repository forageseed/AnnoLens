"""Agent-facing tools for EXISTING gene models: what the agent can look at and which edits it can request.
Read-only on the genome and the input annotation. A structure can only be built from an observed exon chain
(a transcript of an RNA sample or of the long-read assembly) and one of its ORFs; the engine validates it.
Intervals are 0-based half-open internally; the agent sees 1-based closed coordinates."""
import bisect
from collections import defaultdict

from . import genome as G
from .config import KIND_LABEL

MAX_VIEW_BP, OVERVIEW_FLANK = 200000, 5000


class GeneTools:
    def __init__(self, engine):
        self.e = engine; self.labels = {}       # gene -> {label: (source, tid, exons, strand, cov)}

    def _transcripts(self, chrom, a, b):
        rows, starts = self.e._tx.get(chrom, []), self.e._starts.get(chrom, [])
        return [r for r in rows[bisect.bisect_left(starts, a - 300000):bisect.bisect_left(starts, b)] if r[1] > a and r[0] < b]

    def _label(self, gene, chrom, a, b, strand):
        """Stable per-gene labels T1..Tn for same-strand transcripts, most covered first."""
        known = self.labels.setdefault(gene, {}); have = {(v[0], v[1]) for v in known.values()}
        for t in sorted((t for t in self._transcripts(chrom, a, b) if t[4] in (strand, '.')), key=lambda t: -t[5]):
            if (t[2], t[3]) not in have:
                known[f'T{len(known) + 1}'] = (t[2], t[3], t[6], t[4], t[5]); have.add((t[2], t[3]))
        return known

    def _orf(self, m, exons, start_genomic=None, gene=None):
        """Default: the cross-sample dominant ORF if this chain encodes it, else the longest ORF overlapping the CDS."""
        seq = self.e.fa.spliced(m['seqid'], exons, m['strand']); olen = max(1, sum(y - x for x, y in m['cds'])); options = []
        for ts, te in G.find_orfs(seq, min_aa=1):
            cds = G.project(exons, m['strand'], ts, te); start = cds[0][0] if m['strand'] == '+' else cds[-1][1] - 1
            options.append({'cds': cds, 'aa': (te - ts) // 3 - 1, 'start_1based': start + 1,
                            'overlap_original_fraction': round(G.overlap(cds, m['cds']) / olen, 3)})
        if start_genomic is not None:
            hit = [o for o in options if o['start_1based'] == int(start_genomic)]
            return (max(hit, key=lambda o: o['aa']) if hit else None), options
        if gene is not None:
            dom = self.e.cross_tissue_orfs(gene)['dominant']
            if dom is not None:
                hit = [o for o in options if o['cds'] == dom['cds']]
                if hit: return hit[0], options
        over = [o for o in options if o['overlap_original_fraction'] > 0]
        return (over[0] if over else None), options

    def overview_region(self, gene):
        m = self.e.model(gene); a, b = m['exons'][0][0], m['exons'][-1][1]; lo, hi = a, b
        for s, e, src, tid, strand, cov, ex in self._transcripts(m['seqid'], a, b):
            if strand in (m['strand'], '.') and G.overlap(ex, m['cds']) > 0: lo, hi = min(lo, s), max(hi, e)
        for al in self.e.homology.get(gene, []):
            if al['strand'] == m['strand'] and al['cds_blocks']:
                lo = min(lo, min(x for x, y in al['cds_blocks'])); hi = max(hi, max(y for x, y in al['cds_blocks']))
        lo, hi = max(0, lo - OVERVIEW_FLANK), hi + OVERVIEW_FLANK
        if hi - lo > MAX_VIEW_BP: mid = (a + b) // 2; lo, hi = max(0, mid - MAX_VIEW_BP // 2), mid + MAX_VIEW_BP // 2
        return lo, min(hi, self.e.fa.length(m['seqid']))

    def view_data(self, gene, start0, end0, candidate_id=None):
        m = self.e.model(gene); L = self.e.fa.length(m['seqid'])
        start0, end0 = max(0, int(start0)), min(L, int(end0))
        if end0 - start0 > MAX_VIEW_BP or end0 <= start0: raise ValueError(f'view must be 1..{MAX_VIEW_BP} bp')
        tx = []
        for lab, (src, tid, ex, strand, cov) in self._label(gene, m['seqid'], start0, end0, m['strand']).items():
            if ex[-1][1] > start0 and ex[0][0] < end0:
                best, _ = self._orf(m, ex, gene=gene)
                tx.append({'label': lab, 'source': src, 'cov': round(cov, 1), 'exons': ex, 'orf_aa': best['aa'] if best else 0, 'orf_cds': best['cds'] if best else []})
        models = defaultdict(list)
        for g2, x, y in self.e.conn.execute('select gene,a,b from coding where chrom=? and strand=? and a<? and b>?', (m['seqid'], m['strand'], end0, start0)):
            models[g2].append([x, y])
        juncs = defaultdict(int)
        for a, b, n in self.e.conn.execute('select a,b,n from junctions where chrom=? and a>=? and b<=? and strand in (?,?)', (m['seqid'], start0, end0, m['strand'], '.')):
            juncs[(a, b)] += n
        hom = [{'kind': KIND_LABEL.get(al['kind'], al['kind']), 'protein': al['protein'], 'protein_aa': al['protein_length'], 'identity': al['identity'],
                'coverage': round((al['target_end_aa'] - al['target_start_aa'] + 1) / max(1, al['protein_length']), 3), 'blocks': [list(x) for x in al['cds_blocks']]}
               for al in self.e.homology.get(gene, []) if al['strand'] == m['strand'] and al['cds_blocks']]
        dom = self.e.cross_tissue_orfs(gene)['dominant']; cand = None
        if candidate_id:
            c = self.e.candidate(gene, candidate_id); cand = {'candidate_id': candidate_id, 'exons': c['new_exons'], 'cds': c['new_cds']}
        return {'gene_id': gene, 'seqid': m['seqid'], 'strand': m['strand'], 'start0': start0, 'end0': end0,
                'current': {'exons': m['exons'], 'cds': m['cds'], 'aa': G.aa_len(m['cds'])},
                'other_models': {g: sorted(v) for g, v in models.items() if g != gene},
                'transcripts': sorted(tx, key=lambda t: int(t['label'][1:])),
                'junctions': [[a, b, n] for (a, b), n in sorted(juncs.items()) if n >= 3], 'homology': hom,
                'dominant_orf': None if not dom else {'cds': dom['cds'], 'aa': G.aa_len(dom['cds']), 'sources': sorted(dom['sources']), 'adoptable': dom.get('adoptable', True)},
                'candidate': cand}

    @staticmethod
    def transcript_table(data):
        return [{'id': t['label'], 'source': t['source'], 'cov': t['cov'], 'span_1based': [t['exons'][0][0] + 1, t['exons'][-1][1]],
                 'exons': len(t['exons']), 'longest_ORF_aa': t['orf_aa']} for t in data['transcripts']]

    def orf_options(self, gene, label, top=6):
        m = self.e.model(gene); t = self.labels.get(gene, {}).get(label)
        if not t: raise ValueError('unknown transcript label; request a view that shows it first')
        _, options = self._orf(m, t[2], gene=gene)
        return [{'start_1based': o['start_1based'], 'aa': o['aa'], 'overlap_original_fraction': o['overlap_original_fraction']} for o in options[:top]]

    def build_from_transcript(self, gene, label, start_1based=None, merge=False):
        m = self.e.model(gene); t = self.labels.get(gene, {}).get(label)
        if not t: raise ValueError('unknown transcript label; request a view that shows it first')
        exons = [list(x) for x in t[2]]
        orf, _ = self._orf(m, exons, start_1based, gene)
        if orf is None: return {'built': False, 'reason': 'no complete ORF (with that start) overlapping the current CDS on this transcript'}
        merge_with = sorted(self.e.conflicting_genes(m, orf['cds'], gene)) if merge else []
        try: cid = self.e.add_candidate(gene, exons, orf['cds'], t[0], t[1], merge_with)
        except (ValueError, AssertionError) as ex: return {'built': False, 'reason': str(ex), 'orf_aa': orf['aa']}
        v = self.e.validate_annotation_edit(gene, cid)
        return {'built': True, 'candidate_id': cid, 'orf_aa': orf['aa'], 'source_transcript': label, 'merge_with': merge_with,
                'materializable': v['materializable_proposal'], 'errors': v['errors'], 'evidence_gate': v['evidence_gate'],
                'cross_tissue_support_sources': v['cross_tissue_support_sources'], 'overlapping_models': sorted(self.e.conflicting_genes(m, orf['cds'], gene))}
