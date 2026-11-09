"""`annolens scan`: rule-based check of every gene model and of loci without a model (step 1, before the agent).

Per gene (first mRNA):
  coding    missing_start, missing_stop, internal_stop, cds_not_mod3, phase_error, noncanonical_splice, cds_has_N,
            micro_intron (< 50 bp), tiny_internal_exon (< 10 bp), very_short_protein (< 50 aa)
  utr       no UTR in the model; a rule-based UTR from an observed transcript whose introns inside the CDS equal the
            model's (long-read transcripts first) -> utr_rule_* ; utr_conflict when that transcript reaches another
            same-strand gene's CDS
  overlap   same_strand_cds_overlap / antisense_cds_overlap (>= 30 bp)
  evidence  introns_without_rna (an intron with < 3 junction reads), little_external_evidence (< 20 % of the CDS
            covered by transcripts and by protein projections)
  homology  vs intact, full-length independent-annotation proteins: short_vs_independent (model <= 0.75 x protein),
            split_vs_independent (>= 2 models under one protein), fused_vs_independent (one model >= 1.5 x the
            longest of >= 2 proteins), frameshift_vs_independent (>= 99 % identity but the alignment needs a
            frameshift/stop and the model is short)
  reads     read-supported frameshift/stop sites (from `annolens prepare`)
  repeats   >= 50 % of the CDS inside annotated repeats
Routes: agent | genome_frameshift_list | possible_false_gene | utr_rule | utr_unsupported | ok.
Loci without a model: intact projections (independent or related) with no model on that strand, clustered;
  high = an intact alignment, >= 100 aa, and agreement: both protein roles align there when both are configured,
  otherwise >= 2 alignments of the one configured role.
Outputs (work dir): scan/gene_scan.tsv, scan/utr_rule.gff3, scan/summary.json, loci.tsv, queues/{P1,P2}.txt
"""
import bisect, collections, csv, json, re
from pathlib import Path

from . import genome as G

MICRO_INTRON, TINY_EXON, OV, COV, SHORT, LONG, FS_ID, MIN_AA = 50, 10, 30, 0.9, 0.75, 1.5, 0.99, 100
HARD = {'missing_start', 'missing_stop', 'internal_stop', 'cds_not_mod3', 'phase_error', 'noncanonical_splice', 'cds_has_N',
        'micro_intron', 'tiny_internal_exon', 'very_short_protein', 'same_strand_cds_overlap',
        'short_vs_independent', 'split_vs_independent', 'fused_vs_independent', 'frameshift_vs_independent'}


class Idx:
    def __init__(self, items):
        self.v = sorted(items, key=lambda x: (x[0], x[1])); self.s = [x[0] for x in self.v]; self.m = max((x[1] - x[0] for x in self.v), default=0)

    def q(self, a, b):
        i = bisect.bisect_left(self.s, a - self.m)
        while i < len(self.v) and self.v[i][0] < b:
            if self.v[i][1] > a: yield self.v[i]
            i += 1


def read_alignments(gff, faa, kind):
    L, cur = {}, None
    for line in open(faa):
        if line[0] == '>': cur = line[1:].split()[0]; L[cur] = 0
        elif cur: L[cur] += len(line.strip())
    out, cur = [], None
    for line in open(gff):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if f[2] == 'mRNA':
            at = G.gff_attributes(f[8]); t = at['Target'].split(); q = t[0]
            cur = {'kind': kind, 'q': q, 'qlen': L.get(q, 0), 'c': f[0], 'st': f[6], 'cov': (int(t[2]) - int(t[1]) + 1) / max(1, L.get(q, 0)),
                   'id': float(at.get('Identity', 0)), 'pos': float(at.get('Positive', 0)), 'fs': int(at.get('Frameshift', 0)),
                   'stop': int(at.get('StopCodon', 0)), 'score': int(f[5]) if f[5] != '.' else 0, 'blocks': []}
            out.append(cur)
        elif f[2] == 'CDS' and cur is not None: cur['blocks'].append((int(f[3]) - 1, int(f[4])))
    for x in out: x['blocks'].sort()
    return [x for x in out if x['blocks']]


def run(M):
    fa = G.Fasta(M.genome); models = G.read_gff_models(M.annotation)
    cds_idx = collections.defaultdict(list)
    for g, m in models.items():
        for a, b in m['cds']: cds_idx[m['seqid']].append((a, b, (g, m['strand'])))
    cds_idx = {c: Idx(v) for c, v in cds_idx.items()}
    # transcripts
    lr = M.longread['name'] if M.longread else None; tx = []
    for src, p in M.sources.items():
        if p.exists():
            for chrom, tid, strand, cov, ex in G.read_gtf(p, src): tx.append((chrom, src, tid, strand, cov, ex))
    tx_idx = collections.defaultdict(list)
    for t in tx: tx_idx[t[0]].append((t[5][0][0], t[5][-1][1], t))
    tx_idx = {c: Idx(v) for c, v in tx_idx.items()}
    junc = collections.Counter()
    for p in sorted(M.w('junctions').glob('*.bed')) if M.w('junctions').exists() else []:
        for line in open(p):
            x = line.split(); junc[(x[0], int(x[1]), int(x[2]), x[5])] += int(x[4])
    alns = []
    for p in M.proteins:
        gp = M.w('proteins', p['name'] + '.miniprot.gff')
        if gp.exists(): alns += read_alignments(gp, p['fasta'], p['kind'])
    al_idx = collections.defaultdict(list)
    for x in alns: al_idx[(x['c'], x['st'])].append((x['blocks'][0][0], x['blocks'][-1][1], x))
    al_idx = {k: Idx(v) for k, v in al_idx.items()}
    rep = collections.defaultdict(list)
    if M.repeats and M.repeats.exists():
        smap = dict(l.rstrip('\n').split('\t')[:2] for l in open(M.repeats_map)) if M.repeats_map else {}
        for line in open(M.repeats):
            if line.startswith('#') or not line.strip(): continue
            f = line.split('\t')
            if len(f) > 4 and f[2] not in ('long_terminal_repeat', 'target_site_duplication', 'repeat_region', 'primer_binding_site'):
                rep[smap.get(f[0], f[0])].append((int(f[3]) - 1, int(f[4]), None))
    rep = {c: Idx(v) for c, v in rep.items()}
    sites = {}
    sp = M.w('readcheck', 'sites.jsonl')
    if sp.exists():
        for line in open(sp): r = json.loads(line); sites[r['gene_id']] = collections.Counter(s['call'] for s in r['sites'])

    def under(c, st, blocks):
        return sorted({g for x, y, (g, gst) in (cds_idx[c].q(blocks[0][0], blocks[-1][1]) if c in cds_idx else []) if gst == st and G.overlap(blocks, [(x, y)]) > 0
                       and G.overlap(blocks, models[g]['cds']) >= OV})

    # ---- homology routes vs intact independent proteins ----
    route = collections.defaultdict(set); intact = lambda x: x['fs'] == 0 and x['stop'] == 0 and x['cov'] >= COV and x['pos'] >= .9
    per_gene = collections.defaultdict(list)
    for x in alns:
        if x['kind'] != 'ind': continue
        gs = under(x['c'], x['st'], x['blocks'])
        if intact(x):
            for g in gs: per_gene[g].append(x)
            if len(gs) == 1 and G.aa_len(models[gs[0]]['cds']) <= SHORT * x['qlen']: route[gs[0]].add('short_vs_independent')
            elif len(gs) >= 2 and all(G.aa_len(models[g]['cds']) <= SHORT * x['qlen'] for g in gs):
                for g in gs: route[g].add('split_vs_independent')
        elif x['fs'] + x['stop'] > 0 and x['cov'] >= COV and x['id'] >= FS_ID:
            for g in gs:
                if G.aa_len(models[g]['cds']) <= SHORT * x['qlen']: route[g].add('frameshift_vs_independent')
    for g, xs in per_gene.items():
        keep = []
        for x in sorted(xs, key=lambda x: -x['score']):
            if all(x['blocks'][-1][1] <= y['blocks'][0][0] or x['blocks'][0][0] >= y['blocks'][-1][1] for y in keep): keep.append(x)
        if len(keep) >= 2 and G.aa_len(models[g]['cds']) >= LONG * max(x['qlen'] for x in keep): route[g].add('fused_vs_independent')
    for g in route:
        if 'short_vs_independent' in route[g]: route[g].discard('frameshift_vs_independent')

    # ---- per-gene scan ----
    rows, utr = [], []
    for g, m in models.items():
        c, st, cds = m['seqid'], m['strand'], m['cds']; p = []
        seq = fa.spliced(c, cds, st); ph0 = (m['original_phases'] if st == '+' else m['original_phases'][::-1])[0]
        coding = seq[ph0:]; aa = G.translate(coding[:len(coding) // 3 * 3])
        if coding[:3] != 'ATG': p.append('missing_start')
        if coding[-3:] not in G.STOPS: p.append('missing_stop')
        if aa[:-1].count('*'): p.append('internal_stop')
        if len(coding) % 3: p.append('cds_not_mod3')
        order = cds if st == '+' else cds[::-1]; phs = m['original_phases'] if st == '+' else m['original_phases'][::-1]; cum = 0; pe = 0
        for k, ((a, b), ph) in enumerate(zip(order, phs)):
            if k and ph != (3 - ((cum - ph0) % 3)) % 3: pe += 1
            cum += b - a
        if pe: p.append('phase_error')
        if 'N' in seq: p.append('cds_has_N')
        ints = G.introns(cds); bad = 0
        for a, b in ints:
            mot = fa.get(c, a, a + 2) + '-' + fa.get(c, b - 2, b) if st == '+' else G.rc(fa.get(c, b - 2, b)) + '-' + G.rc(fa.get(c, a, a + 2))
            bad += mot not in ('GT-AG', 'GC-AG', 'AT-AC')
        if bad: p.append('noncanonical_splice')
        if any(b - a < MICRO_INTRON for a, b in ints): p.append('micro_intron')
        if len(cds) > 2 and any(b - a < TINY_EXON for a, b in cds[1:-1]): p.append('tiny_internal_exon')
        prot = len(coding) // 3 - (1 if coding[-3:] in G.STOPS else 0)
        if prot < 50: p.append('very_short_protein')
        coding_problem = bool(p)
        same = anti = 0
        for a, b in cds:
            for x, y, (g2, st2) in (cds_idx[c].q(a, b) if c in cds_idx else []):
                if g2 == g: continue
                if st2 == st: same += min(b, y) - max(a, x)
                else: anti += min(b, y) - max(a, x)
        if same >= OV: p.append('same_strand_cds_overlap')
        if anti >= OV: p.append('antisense_cds_overlap')
        if ints and any(junc.get((c, a, b, st), 0) + junc.get((c, a, b, '.'), 0) < 3 for a, b in ints): p.append('introns_without_rna')
        tot = max(1, sum(b - a for a, b in cds))
        txs = [t for _, _, t in (tx_idx[c].q(cds[0][0], cds[-1][1]) if c in tx_idx else []) if t[3] in (st, '.')]
        txcov = max((G.overlap(cds, t[5]) / tot for t in txs), default=0)
        prcov = max((G.overlap(cds, x['blocks']) / tot for _, _, x in (al_idx[(c, st)].q(cds[0][0], cds[-1][1]) if (c, st) in al_idx else [])), default=0)
        little = txcov < .2 and prcov < .2
        if little: p.append('little_external_evidence')
        r6 = sorted(route.get(g, ())); p += r6
        sc = sites.get(g, collections.Counter())
        if sc['genome_sequence_supported']: p.append('read_supported_frameshift_or_stop')
        rf = min(1.0, sum(G.overlap(cds, [(x, y)]) for x, y, _ in (rep[c].q(cds[0][0], cds[-1][1]) if c in rep else [])) / tot)
        if rf >= .5: p.append('cds_in_repeat')
        # UTR rule
        ci = G.introns(cds); s, e = cds[0][0], cds[-1][1]; best, conflict = None, False
        for t in txs:
            ex = t[5]
            if ex[0][0] > s or ex[-1][1] < e: continue
            ti = G.introns(ex)
            if [y for y in ti if y[0] >= s and y[1] <= e] != ci or any(y[0] < s < y[1] or y[0] < e < y[1] for y in ti): continue
            if not any(a <= s < b for a, b in ex) or not any(a < e <= b for a, b in ex) or (ex[0][0] == s and ex[-1][1] == e): continue
            other = {g2 for a, b in ex for x, y, (g2, st2) in cds_idx[c].q(a, b) if st2 == st and g2 != g}
            if other: conflict = True; continue
            key = (t[1] == lr, t[4])
            if best is None or key > best[0]: best = (key, t)
        t = best[1] if best else None
        has_utr = m['exons'][0][0] < s or m['exons'][-1][1] > e
        if has_utr: u = 'has_utr'
        elif t:
            lo, hi = t[5][0][0] < s, t[5][-1][1] > e; u5, u3 = (lo, hi) if st == '+' else (hi, lo)
            u = 'utr_rule_' + ('both' if u5 and u3 else '5only' if u5 else '3only')
            mid = 'utr.' + m['transcript_id']
            utr.append(f"{c}\tAnnoLens_UTR_rule\tmRNA\t{t[5][0][0] + 1}\t{t[5][-1][1]}\t.\t{st}\t.\tID={mid};Parent={g};utr_source={t[1]}:{t[2]}")
            utr += [f"{c}\tAnnoLens_UTR_rule\texon\t{a + 1}\t{b}\t.\t{st}\t.\tParent={mid}" for a, b in t[5]]
            utr += [f"{c}\tAnnoLens_UTR_rule\tCDS\t{a + 1}\t{b}\t.\t{st}\t{ph}\tParent={mid}" for (a, b), ph in zip(cds, G.phases(cds, st))]
        else: u = 'utr_conflict' if conflict else 'utr_unsupported'
        if coding_problem or r6 or 'same_strand_cds_overlap' in p or u == 'utr_conflict' or ('introns_without_rna' in p and not t and not little): rt = 'agent'
        elif 'read_supported_frameshift_or_stop' in p: rt = 'genome_frameshift_list'
        elif little and not r6 and (rf >= .5 or prot < 100): rt = 'possible_false_gene'
        elif u.startswith('utr_rule'): rt = 'utr_rule'
        elif u == 'has_utr': rt = 'ok'
        else: rt = 'utr_unsupported'
        rows.append({'gene_id': g, 'chrom': c, 'start_1based': cds[0][0] + 1, 'end_1based': cds[-1][1], 'strand': st, 'protein_aa': prot,
                     'route': rt, 'utr': u, 'repeat_cds_fraction': round(rf, 2), 'problems': ';'.join(p)})
    out = M.w('scan'); out.mkdir(parents=True, exist_ok=True)
    with open(out / 'gene_scan.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
    (out / 'utr_rule.gff3').write_text('##gff-version 3\n' + '\n'.join(utr) + '\n')
    q = M.w('queues'); q.mkdir(parents=True, exist_ok=True)
    p1 = [r['gene_id'] for r in rows if r['route'] == 'agent' and (set(r['problems'].split(';')) & HARD or r['utr'] == 'utr_conflict')]
    p1s = set(p1); p2 = [r['gene_id'] for r in rows if r['route'] == 'agent' and r['gene_id'] not in p1s]
    (q / 'P1.txt').write_text('\n'.join(p1) + '\n'); (q / 'P2.txt').write_text('\n'.join(p2) + '\n')

    # ---- loci without a model ----
    cand = [x for x in alns if x['kind'] in ('ind', 'pot') and x['cov'] >= COV and x['pos'] >= .7 and not under(x['c'], x['st'], x['blocks'])]
    cand.sort(key=lambda x: (x['c'], x['st'], x['blocks'][0][0])); loci = []
    for x in cand:
        Lc = loci[-1] if loci else None
        if Lc and Lc['c'] == x['c'] and Lc['st'] == x['st'] and x['blocks'][0][0] <= Lc['e']: Lc['e'] = max(Lc['e'], x['blocks'][-1][1]); Lc['aln'].append(x)
        else: loci.append({'c': x['c'], 'st': x['st'], 's': x['blocks'][0][0], 'e': x['blocks'][-1][1], 'aln': [x]})
    both_roles = {'ind', 'pot'} <= {p['kind'] for p in M.proteins}   # both roles configured: both must agree
    lrows = []
    for Lc in loci:
        good = [x for x in Lc['aln'] if x['fs'] == 0 and x['stop'] == 0]; b = max(good or Lc['aln'], key=lambda x: x['score'])
        kinds = {x['kind'] for x in Lc['aln']}
        agree = len(kinds) == 2 if both_roles else len(Lc['aln']) >= 2
        tier = ('high' if agree and b['qlen'] >= MIN_AA else 'low') if good else 'disrupted'
        lrows.append({'locus_id': f"{Lc['c']}:{Lc['s'] + 1}-{Lc['e']}{Lc['st']}", 'chrom': Lc['c'], 'start_1based': Lc['s'] + 1, 'end_1based': Lc['e'],
                      'strand': Lc['st'], 'tier': tier, 'best_protein': b['q'], 'expected_aa': b['qlen'], 'alignments': len(Lc['aln'])})
    with open(M.w('loci.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['locus_id', 'chrom', 'start_1based', 'end_1based', 'strand', 'tier', 'best_protein', 'expected_aa', 'alignments'], delimiter='\t')
        w.writeheader(); w.writerows(lrows)
    (q / 'missing_high.txt').write_text('\n'.join(r['locus_id'] for r in lrows if r['tier'] == 'high') + '\n')
    summ = {'genes': len(rows), 'routes': dict(collections.Counter(r['route'] for r in rows)), 'utr': dict(collections.Counter(r['utr'] for r in rows)),
            'problems': dict(collections.Counter(x for r in rows for x in r['problems'].split(';') if x).most_common()),
            'queues': {'P1': len(p1), 'P2': len(p2)}, 'loci': dict(collections.Counter(r['tier'] for r in lrows))}
    (out / 'summary.json').write_text(json.dumps(summ, indent=1))
    return summ
