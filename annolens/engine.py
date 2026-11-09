"""AnnoLens gene engine: deterministic evidence, candidate records and validation for EXISTING gene models.

Species- and site-independent: every input comes from the manifest (annolens.toml) and the work dir built by
`annolens prepare`. The rules are those of the engine chain used by skill v3.5 (labels in brackets name the
version that introduced each rule):
  candidate      only an exon chain that was observed (a transcript of some RNA sample or of the long-read assembly)
                 with one of its complete ORFs overlapping the current CDS [3.1]; optional merge of neighbouring
                 models that the new CDS overlaps [3.1]
  validation     complete ORF (ATG .. terminal stop, %3, no internal stop), no N, canonical splice sites, inside the
                 widened locus [v5], CDS inside exons, no same-strand CDS overlap >= 30 bp unless merged [v3/3.1],
                 single-exon structures need coding support [v3]
  evidence gate  one of: full-chain long read covering >= 95 % of the CDS [v3]; RNA/homology gate: every intron
                 >= 3 junction reads or reproduced by a projection, start/stop/length anchored by an independent or
                 related-species projection (0.9-1.1 x protein) [v4, 3.4]; cross-sample dominant ORF [v5]
  fix only what is broken: a complete original that already passes a gate is not replaced [v4], except when the
                 dominant ORF and an intact independent protein prove the original truncated [3.4]
  anti-truncation: no candidate shorter than a strong dominant ORF with fewer sources [v5]; a dominant ORF never
                 shortens a complete original [v5]; ambiguous coverage vote -> no dominant ORF [3.3]
The genome and the input annotation are never edited; outputs are proposals.
"""
import bisect, fcntl, hashlib, json, sqlite3, sys
from collections import defaultdict
from pathlib import Path

from . import genome as G

ENGINE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
MIN_RNA_JUNCTION, MIN_HOMOLOGY_IDENTITY, END_TOLERANCE_AA, LENGTH_RATIO = 3, 0.7, 3, (0.9, 1.1)
LOCUS_FLANK, MAX_LOCUS_EXTENSION, CONFLICT_BP = 10000, 100000, 30
DOM_MIN_SOURCES, TIE_COVERAGE_RATIO, COVERAGE_OVERRIDE_RATIO = 3, 1.1, 2.0
ANCHOR_KINDS = ('pot', 'ind')
IND_COV, TRUNC = 0.9, 0.9
GATE_ERRORS = {'no_exact_full_chain_longread_support', 'intron_without_RNA_or_homology_support', 'no_anchoring_projection',
               'start_not_supported_by_projection', 'stop_not_supported_by_projection', 'length_outside_projection_range',
               'single_exon_coding_support_insufficient'}
ANCHOR_ERRORS = {'no_anchoring_projection', 'start_not_supported_by_projection', 'stop_not_supported_by_projection',
                 'length_outside_projection_range'}


def _stamp(paths):
    return G.sha(json.dumps([[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in paths if p and Path(p).exists()]))


class GeneEngine:
    def __init__(self, manifest, out=None):
        self.M = manifest; self.fa = G.Fasta(manifest.genome)
        self.out = Path(out or manifest.w('engine_out')); (self.out / 'candidates').mkdir(parents=True, exist_ok=True)
        if str(manifest.w('deps')) not in sys.path: sys.path.insert(0, str(manifest.w('deps')))
        import pysam
        self._pysam = pysam
        self.models = G.read_gff_models(manifest.annotation)
        self._index()
        hp = manifest.w('homology.jsonl')
        self.homology = {}
        if hp.exists():
            for line in hp.open():
                r = json.loads(line); self.homology[r['gene_id']] = r['alignments']
        self.homology_sha256 = hashlib.sha256(hp.read_bytes()).hexdigest() if hp.exists() else None
        self._tx = defaultdict(list)
        for src, path in manifest.sources.items():
            if not path.exists(): continue
            for chrom, tid, strand, cov, ex in G.read_gtf(path, src):
                try: ex = G.normalize(ex)
                except AssertionError: continue
                self._tx[chrom].append((ex[0][0], ex[-1][1], src, tid, strand, cov, ex))
        for v in self._tx.values(): v.sort(key=lambda r: r[0])
        self._starts = {c: [r[0] for r in v] for c, v in self._tx.items()}
        self._dom_cache, self._locus_cache = {}, {}
        self.evidence_sha256 = G.sha(json.dumps([self.inputs_stamp, self.homology_sha256]))

    # ---------- index (coding blocks of all genes, junction counts per sample) ----------
    def _index(self):
        db = self.M.w('index.sqlite'); jb = sorted(self.M.w('junctions').glob('*.bed')) if self.M.w('junctions').exists() else []
        self.inputs_stamp = _stamp([self.M.annotation] + jb)
        lock = self.M.w('index.lock'); lock.parent.mkdir(parents=True, exist_ok=True)
        with lock.open('a') as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            ok = False
            if db.exists():
                try:
                    with sqlite3.connect(db) as c: ok = dict(c.execute('select k,v from metadata')).get('stamp') == self.inputs_stamp
                except sqlite3.DatabaseError: ok = False
            if not ok:
                tmp = Path(str(db) + '.tmp'); tmp.unlink(missing_ok=True); c = sqlite3.connect(tmp)
                c.executescript('create table metadata(k text primary key,v text);'
                                'create table coding(gene text,chrom text,a integer,b integer,strand text);create index coding_lookup on coding(chrom,a,b);'
                                'create table junctions(chrom text,a integer,b integer,strand text,tissue text,n integer);create index junction_lookup on junctions(chrom,a,b,strand);')
                c.execute('insert into metadata values (?,?)', ('stamp', self.inputs_stamp))
                c.executemany('insert into coding values (?,?,?,?,?)', ((g, m['seqid'], a, b, m['strand']) for g, m in self.models.items() for a, b in m['cds']))
                for p in jb:
                    batch = []
                    for line in p.open():
                        x = line.split(); batch.append((x[0], int(x[1]), int(x[2]), x[5], p.stem, int(x[4])))
                        if len(batch) >= 20000: c.executemany('insert into junctions values (?,?,?,?,?,?)', batch); batch = []
                    if batch: c.executemany('insert into junctions values (?,?,?,?,?,?)', batch)
                c.commit(); c.close(); tmp.replace(db)
        self.conn = sqlite3.connect(db, check_same_thread=False)

    # ---------- basics ----------
    def model(self, gene):
        m = self.models.get(gene); assert m, 'unknown_gene'
        return json.loads(json.dumps(m))

    def reference_hash(self, m):
        a, b = self.locus(m); return G.sha(self.fa.get(m['seqid'], a, b))

    def splice_context(self, m, exons):
        out = []
        for a, b in G.introns(exons):
            left, right = self.fa.get(m['seqid'], a, a + 2), self.fa.get(m['seqid'], b - 2, b)
            motif = left + '-' + right if m['strand'] == '+' else G.rc(right) + '-' + G.rc(left)
            out.append({'intron_start0': a, 'intron_end0': b, 'motif': motif, 'recognized': motif in ('GT-AG', 'GC-AG', 'AT-AC')})
        return out

    def describe(self, m, exons, cds):
        seq = self.fa.spliced(m['seqid'], cds, m['strand']); aa = G.translate(seq)
        return {'exons': exons, 'cds': cds, 'cds_phases_genomic_order': G.phases(cds, m['strand']),
                'cds_sequence_transcript_orientation': seq, 'cds_sha256': G.sha(seq), 'translation': aa,
                'start_codon': seq[:3], 'terminal_stop': seq[-3:], 'terminal_stop_included_in_CDS': seq[-3:] in G.STOPS,
                'cds_length_mod3': len(seq) % 3, 'internal_stops': aa[:-1].count('*') if aa.endswith('*') else aa.count('*'),
                'splice_sites': self.splice_context(m, exons)}

    @staticmethod
    def complete(d): return d['start_codon'] == 'ATG' and d['terminal_stop_included_in_CDS'] and not d['cds_length_mod3'] and not d['internal_stops']

    def locus(self, m):
        """Model +-10 kb, widened to same-strand observed transcripts overlapping the CDS (<= 100 kb further) [v5]."""
        key = (m['seqid'], m['transcript_id'])
        if key in self._locus_cache: return self._locus_cache[key]
        L = self.fa.length(m['seqid'])
        a, b = max(0, m['exons'][0][0] - LOCUS_FLANK), min(L, m['exons'][-1][1] + LOCUS_FLANK)
        lo, hi = a, b
        rows, starts = self._tx.get(m['seqid'], []), self._starts.get(m['seqid'], [])
        for s, e, src, tid, strand, cov, ex in rows[bisect.bisect_left(starts, a - 300000):bisect.bisect_left(starts, b)]:
            if e > a and strand in (m['strand'], '.') and G.overlap(ex, m['cds']) > 0: lo, hi = min(lo, s), max(hi, e)
        r = (max(0, lo, a - MAX_LOCUS_EXTENSION), min(L, hi, b + MAX_LOCUS_EXTENSION))
        self._locus_cache[key] = r; return r

    def junction_reads(self, chrom, a, b, strand):
        return sum(r[0] for r in self.conn.execute('select n from junctions where chrom=? and a=? and b=? and strand in (?,?)', (chrom, a, b, strand, '.')))

    # ---------- evidence ----------
    def evidence(self, m, exons, cds, source):
        js = G.introns(exons); counts = [0] * len(js); samples = [set() for _ in js]
        for i, (a, b) in enumerate(js):
            for t, n in self.conn.execute('select tissue,n from junctions where chrom=? and a=? and b=? and strand in (?,?)', (m['seqid'], a, b, m['strand'], '.')):
                counts[i] += n
                if n: samples[i].add(t)
        full = []
        lr = self.M.longread
        if lr and lr.get('bam') and Path(lr['bam']).exists():
            with self._pysam.AlignmentFile(str(lr['bam']), 'rb', index_filename=str(lr['index']) if lr.get('index') else None) as f:
                for r in f.fetch(m['seqid'], exons[0][0], exons[-1][1]):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
                    pos = r.reference_start; rj = []
                    for op, n in r.cigartuples or []:
                        if op == 3: rj.append([pos, pos + n])
                        if op in (0, 2, 3, 7, 8): pos += n
                    if [j for j in rj if exons[0][0] <= j[0] and j[1] <= exons[-1][1]] != js: continue
                    cover = G.overlap(r.get_blocks(), cds) / sum(b - a for a, b in cds)
                    if cover >= .95: full.append({'read_id': r.query_name, 'mapq': r.mapping_quality, 'coding_fraction': round(cover, 5)})
                    if len(full) >= 20: break
        return {'source_track': source, 'junction_read_counts': counts, 'supporting_samples': [sorted(x) for x in samples],
                'exact_chain_full_length_longreads': full, 'full_length_chain_observed': bool(full),
                'RNA_chain_all_junctions_ge3': bool(js) and all(x >= 3 for x in counts)}

    def protein_cds_fraction(self, m):
        """Best fraction of the model CDS covered by one protein projection (coding support of single-exon models)."""
        tot = max(1, sum(b - a for a, b in m['cds'])); best = 0.0
        for al in self.homology.get(m['gene_id'], []):
            if al['strand'] == m['strand'] and al['identity'] >= MIN_HOMOLOGY_IDENTITY and al['cds_blocks']:
                best = max(best, G.overlap(m['cds'], al['cds_blocks']) / tot)
        return best

    # ---------- cross-sample dominant ORF [v5, 3.3] ----------
    def cross_tissue_orfs(self, gene):
        if gene in self._dom_cache: return self._dom_cache[gene]
        m = self.model(gene); a, b = self.locus(m); ocds = m['cds']; olen = max(1, sum(y - x for x, y in ocds))
        rows, starts = self._tx.get(m['seqid'], []), self._starts.get(m['seqid'], [])
        groups = {}
        for s, e, src, tid, strand, cov, ex in rows[bisect.bisect_left(starts, a - 300000):bisect.bisect_left(starts, b)]:
            if e <= a or strand not in (m['strand'], '.') or ex[0][0] < a or ex[-1][1] > b or G.overlap(ex, ocds) == 0: continue
            seq = self.fa.spliced(m['seqid'], ex, m['strand']); best = None
            for ts, te in G.find_orfs(seq, min_aa=1):
                cds = G.project(ex, m['strand'], ts, te)
                if G.overlap(cds, ocds) / olen >= .5: best = cds; break
            if best is None: continue
            g = groups.setdefault(json.dumps(best), {'cds': best, 'sources': set(), 'cov': 0.0, 'chains': {}})
            g['sources'].add(src); g['cov'] += cov
            ck = json.dumps(ex)
            if ck not in g['chains'] or g['chains'][ck]['cov'] < cov: g['chains'][ck] = {'exons': ex, 'cov': cov, 'src': src, 'tid': tid}
        ranked = sorted(groups.values(), key=lambda g: (-len(g['sources']), -g['cov']))
        dom = None
        if ranked and len(ranked[0]['sources']) >= DOM_MIN_SOURCES:
            d = ranked[0]
            tied = len(ranked) > 1 and len(ranked[1]['sources']) == len(d['sources']) and d['cov'] < TIE_COVERAGE_RATIO * ranked[1]['cov']
            juncs = all(self.junction_reads(m['seqid'], x, y, m['strand']) >= MIN_RNA_JUNCTION for x, y in G.introns(d['cds']))
            oc = self.describe(m, m['exons'], m['cds'])
            if not tied and juncs:
                dom = d; dom['adoptable'] = not (self.complete(oc) and sum(y - x for x, y in d['cds']) < olen)
        info = {'groups': ranked, 'dominant': dom}
        if dom is not None:     # [3.3] ambiguous vote: a rival with >= 2 sources and >= 2x coverage cancels the dominant ORF
            rival = [g for g in ranked if g is not dom and len(g['sources']) >= 2 and g['cov'] >= COVERAGE_OVERRIDE_RATIO * dom['cov']]
            if rival: info = {**info, 'dominant': None, 'coverage_conflict': {'rival_cov': round(rival[0]['cov'], 1), 'dominant_cov': round(dom['cov'], 1)}}
        self._dom_cache[gene] = info
        return info

    # ---------- RNA / homology gate [v4, 3.4] ----------
    @staticmethod
    def _end(cds, strand, which):
        if which == 'start': return cds[0][0] if strand == '+' else cds[-1][1] - 1
        return cds[-1][1] - 1 if strand == '+' else cds[0][0]

    def rna_homology_gate(self, gene, m, exons, cds):
        errors, detail = [], {}
        als = self.homology.get(gene, [])
        anchors = [al for al in als if al['kind'] in ANCHOR_KINDS and al['strand'] == m['strand'] and al['identity'] >= MIN_HOMOLOGY_IDENTITY and al['cds_blocks']]
        hom_introns = {tuple(x) for al in als if al['strand'] == m['strand'] and al['identity'] >= MIN_HOMOLOGY_IDENTITY and al['cds_blocks']
                       for x in G.introns(G.normalize(sorted(al['cds_blocks'])))}
        unsupported = [[a, b, n] for a, b in G.introns(exons) for n in [self.junction_reads(m['seqid'], a, b, m['strand'])]
                       if n < MIN_RNA_JUNCTION and (a, b) not in hom_introns]
        if unsupported: errors.append('intron_without_RNA_or_homology_support')
        detail['unsupported_introns'] = unsupported
        if not anchors:
            errors.append('no_anchoring_projection'); return errors, detail
        cstart, cstop = self._end(cds, m['strand'], 'start'), self._end(cds, m['strand'], 'stop')
        od = self.describe(m, m['exons'], m['cds'])
        start_ok = cstart == self._end(m['cds'], m['strand'], 'start') and od['start_codon'] == 'ATG'
        stop_ok = cstop == self._end(m['cds'], m['strand'], 'stop') and od['terminal_stop_included_in_CDS']
        aa = sum(y - x for x, y in cds) // 3 - 1; ratios = []
        for al in anchors:
            bl = G.normalize(sorted(al['cds_blocks'])); k = al['target_start_aa']
            if k <= END_TOLERANCE_AA and abs(self._end(bl, m['strand'], 'start') - cstart) <= 3 * (k - 1): start_ok = True
            if al['protein_length'] and al['target_end_aa'] >= al['protein_length'] - END_TOLERANCE_AA and abs(self._end(bl, m['strand'], 'stop') - cstop) <= 3: stop_ok = True
            if al['protein_length']: ratios.append(aa / al['protein_length'])
        if not start_ok: errors.append('start_not_supported_by_projection')
        if not stop_ok: errors.append('stop_not_supported_by_projection')
        if not any(LENGTH_RATIO[0] <= r <= LENGTH_RATIO[1] for r in ratios): errors.append('length_outside_projection_range')
        detail.update(start_supported=start_ok, stop_supported=stop_ok, length_ratios=[round(r, 3) for r in ratios],
                      anchor_kinds=sorted({al['kind'] for al in anchors}))
        return errors, detail

    # ---------- candidates [3.1] ----------
    def conflicting_genes(self, m, cds, gene):
        hit = defaultdict(int)
        for x0, y0 in cds:
            for g2, x, y in self.conn.execute('select gene,a,b from coding where chrom=? and strand=? and a<? and b>? and gene!=?', (m['seqid'], m['strand'], y0, x0, gene)):
                hit[g2] += max(0, min(y0, y) - max(x0, x))
        return {g: n for g, n in hit.items() if n >= CONFLICT_BP}

    def observed(self, m, src, tid, exons):
        return any(s2 == src and t2 == tid and [list(x) for x in ex] == [list(x) for x in exons] for s, e, s2, t2, st, cov, ex in self._tx.get(m['seqid'], []))

    def _cand_path(self, gene): return self.out / 'candidates' / (G.sha(gene) + '.json')

    def add_candidate(self, gene, exons, cds, src, tid, merge_with=()):
        m = self.model(gene); rh = self.reference_hash(m)
        exons = G.normalize([list(x) for x in exons]); cds = G.normalize([list(x) for x in cds])
        if G.overlap(cds, m['cds']) == 0: raise ValueError('ORF does not overlap the current CDS; that is a different gene')
        if not self.observed(m, src, tid, exons): raise ValueError('exon chain is not an observed transcript')
        source = {'kind': 'observed_transcript', 'model_id': f'{src}:{tid}', 'exact_transcript_exons': exons, 'source_strand': None}
        cand = {'engine_sha256': ENGINE_SHA256, 'evidence_sha256': self.evidence_sha256, 'gene_id': gene, 'transcript_id': m['transcript_id'],
                'seqid': m['seqid'], 'strand': m['strand'], 'reference_locus_sha256': rh, 'old_exons': m['exons'], 'old_cds': m['cds'],
                'new_exons': exons, 'new_cds': cds, 'checks': self.describe(m, exons, cds), 'evidence': self.evidence(m, exons, cds, source),
                'genome_edits': [], 'merge_with': sorted(merge_with), 'UTR': G.utr_intervals(exons, cds, m['strand'])}
        cand['candidate_id'] = 'cand_' + G.sha(json.dumps(cand, sort_keys=True))[:24]
        p = self._cand_path(gene)
        data = json.loads(p.read_text()) if p.exists() else {'reference_locus_sha256': rh, 'candidates': []}
        if data['reference_locus_sha256'] != rh: data = {'reference_locus_sha256': rh, 'candidates': []}
        if not any(c['candidate_id'] == cand['candidate_id'] for c in data['candidates']): data['candidates'].append(cand)
        p.write_text(json.dumps(data))
        return cand['candidate_id']

    def candidate(self, gene, ident):
        p = self._cand_path(gene); assert p.exists(), 'unknown_candidate_id'
        data = json.loads(p.read_text()); m = self.model(gene)
        assert data['reference_locus_sha256'] == self.reference_hash(m), 'reference_hash_mismatch'
        hits = [x for x in data['candidates'] if x['candidate_id'] == ident]; assert len(hits) == 1, 'unknown_candidate_id'
        c = hits[0]; payload = {k: v for k, v in c.items() if k != 'candidate_id'}
        assert ident == 'cand_' + G.sha(json.dumps(payload, sort_keys=True))[:24], 'candidate_content_hash_mismatch'
        assert c['engine_sha256'] == ENGINE_SHA256, 'engine_code_hash_mismatch'
        assert c['evidence_sha256'] == self.evidence_sha256, 'evidence_changed_since_candidate'
        return c

    # ---------- validation ----------
    def validate_annotation_edit(self, gene, candidate_id):
        m = self.model(gene); c = self.candidate(gene, candidate_id); errors = []; gate, gate_detail = None, {}
        if c.get('genome_edits'): errors.append('genome_edits_forbidden')
        if c['seqid'] != m['seqid'] or c['strand'] != m['strand'] or c['transcript_id'] != m['transcript_id']: errors.append('identity_or_strand_change')
        if c['old_cds'] != m['cds'] or c['old_exons'] != m['exons']: errors.append('original_annotation_mismatch')
        a, b = self.locus(m); merge = set(c.get('merge_with') or [])
        try:
            exons = G.normalize(c['new_exons']); cds = G.normalize(c['new_cds'])
            if not a <= exons[0][0] < exons[-1][1] <= b: errors.append('outside_locus')
            if G.overlap(cds, exons) != sum(y - x for x, y in cds): errors.append('CDS_outside_exons')
            desc = self.describe(m, exons, cds)
            if not self.complete(desc): errors.append('invalid_complete_ORF')
            if 'N' in desc['cds_sequence_transcript_orientation']: errors.append('ambiguous_CDS')
            if any(not x['recognized'] for x in desc['splice_sites']): errors.append('noncanonical_requires_special_review')
            src, tid = c['evidence']['source_track']['model_id'].split(':', 1)
            if not self.observed(m, src, tid, exons): errors.append('unsupported_or_invented_exon_chain')
            conflicts = self.conflicting_genes(m, cds, gene)
            if conflicts and not merge: errors.append('overlaps_other_same_strand_CDS_requires_merge_review')
            orig_complete = self.complete(self.describe(m, m['exons'], m['cds']))
            if c['evidence'].get('full_length_chain_observed'):
                gate = 'longread_full_chain'
                if orig_complete and self.evidence(m, m['exons'], m['cds'], {'kind': 'original'})['full_length_chain_observed']:
                    errors.append('original_model_already_has_full_chain_longread'); gate = None
            else:
                gerr, gate_detail = self.rna_homology_gate(gene, m, exons, cds)
                if orig_complete and not self.rna_homology_gate(gene, m, m['exons'], m['cds'])[0]:
                    gerr.append('original_model_already_passes_rna_homology_gate')
                if gerr: errors.extend(gerr)
                else: gate = 'rna_homology'
        except (AssertionError, KeyError, ValueError):
            errors.append('invalid_intervals'); cds = c['new_cds']; exons = c['new_exons']
        if len(c['new_exons']) == 1 and gate != 'rna_homology' and self.protein_cds_fraction(m) < .8:
            errors.append('single_exon_coding_support_insufficient')
        if gate is None and not c['evidence'].get('full_length_chain_observed'): errors.append('no_exact_full_chain_longread_support')
        # [v5] dominant ORF: gate bypass, never shorten a complete original, anti-truncation
        info = self.cross_tissue_orfs(gene); dom = info['dominant']
        support = next((len(g['sources']) for g in info['groups'] if g['cds'] == cds), 0)
        if dom is not None and cds == dom['cds']:
            if not dom.get('adoptable', True): errors.append('dominant_ORF_shorter_than_complete_original')
            elif [e for e in errors if e in GATE_ERRORS]:
                errors = [e for e in errors if e not in GATE_ERRORS]; gate = 'cross_sample_dominant'
        elif dom is not None:
            if sum(y - x for x, y in cds) < sum(y - x for x, y in dom['cds']) and support < len(dom['sources']):
                errors.append('shorter_than_cross_sample_dominant_ORF')
        # [3.1] merge of neighbouring models
        if merge:
            problems = []
            conflicts = self.conflicting_genes(m, cds, gene)
            if not set(conflicts) <= merge: problems.append('merge_list_misses_overlapping_models')
            for g2 in merge:
                rows = list(self.conn.execute('select a,b from coding where gene=?', (g2,)))
                if not rows or min(x for x, y in rows) < exons[0][0] or max(y for x, y in rows) > exons[-1][1]: problems.append('merged_model_not_inside_new_transcript')
            strong = (dom is not None and dom.get('adoptable', True) and cds == dom['cds']) or bool(c['evidence'].get('full_length_chain_observed'))
            if not strong: problems.append('merge_requires_dominant_ORF_or_full_chain_longread')
            errors += sorted(set(problems))
        # [3.4] narrow lift of the "already supported" rule when an independent protein proves truncation
        if errors == ['original_model_already_has_full_chain_longread'] and dom is not None and dom.get('adoptable', True) and cds == dom['cds']:
            gerr, det = self.rna_homology_gate(gene, m, G.normalize(c['new_exons']), cds)
            orig_aa = G.aa_len(m['cds'])
            ind = [al for al in self.homology.get(gene, []) if al['kind'] == 'ind' and al['strand'] == m['strand'] and al['frameshift'] == 0
                   and al['stop_codon_found'] == 0 and al['protein_length'] and al['target_end_aa'] - al['target_start_aa'] + 1 >= IND_COV * al['protein_length']
                   and orig_aa < TRUNC * al['protein_length']]
            if not set(gerr) & ANCHOR_ERRORS and ind:
                errors = []; gate = 'dominant_ORF_plus_independent_annotation'
                gate_detail = {**det, 'lifted': 'original_model_already_has_full_chain_longread', 'independent_protein': ind[0]['protein'],
                               'independent_length': ind[0]['protein_length'], 'original_aa': orig_aa}
        noop = c['new_cds'] == m['cds'] and c['new_exons'] == m['exons']
        return {'candidate_id': candidate_id, 'is_noop': noop, 'materializable_proposal': not errors and not noop, 'errors': errors,
                'evidence_gate': gate, 'gate_detail': gate_detail, 'cross_tissue_support_sources': support, 'merge_with': sorted(merge),
                'genome_unchanged': True}

    def materialize_proposal(self, gene, candidate_id):
        c = self.candidate(gene, candidate_id); v = self.validate_annotation_edit(gene, candidate_id)
        assert v['materializable_proposal'], 'candidate_not_materializable:' + ','.join(v['errors'])
        m = self.model(gene); d = self.out / 'proposals'; d.mkdir(parents=True, exist_ok=True)
        gff = d / f'{G.sha(gene)[:16]}.{candidate_id}.gff3'
        rows = ['##gff-version 3', '# AnnoLens proposal; genome unchanged; terminal stop included in CDS']
        R = lambda typ, a, b, ph, at: rows.append('\t'.join(map(str, [m['seqid'], 'AnnoLens', typ, a + 1, b, '.', m['strand'], ph, at])))
        ex, cds, t = c['new_exons'], c['new_cds'], m['transcript_id']
        R('gene', ex[0][0], ex[-1][1], '.', f'ID={gene}'); R('mRNA', ex[0][0], ex[-1][1], '.', f'ID={t};Parent={gene}')
        for i, (a, b) in enumerate(ex, 1): R('exon', a, b, '.', f'ID={t}.exon{i};Parent={t}')
        for i, ((a, b), ph) in enumerate(zip(cds, G.phases(cds, m['strand'])), 1): R('CDS', a, b, ph, f'ID={t}.cds{i};Parent={t}')
        for typ in ('five_prime_UTR', 'three_prime_UTR'):
            for i, (a, b) in enumerate(c['UTR'][typ], 1): R(typ, a, b, '.', f'ID={t}.{typ}{i};Parent={t}')
        gff.write_text('\n'.join(rows) + '\n')
        rec = d / f'{G.sha(gene)[:16]}.{candidate_id}.json'; rec.write_text(json.dumps({**c, 'validation': v, 'proposed_gff': str(gff)}))
        return {'candidate_id': candidate_id, 'proposed_gff': str(gff), 'record': str(rec), 'genome_unchanged': True}
