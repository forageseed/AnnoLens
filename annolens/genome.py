"""Sequence helpers. All intervals are genomic 0-based half-open; strands '+'/'-'."""
import hashlib, itertools, re
from pathlib import Path

STOPS = ('TAA', 'TAG', 'TGA')
CODONS = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)),
                  'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def sha(x): return hashlib.sha256(x.encode()).hexdigest()
def rc(s): return s.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
def translate(s): return ''.join(CODONS.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def normalize(iv):
    iv = [list(map(int, x)) for x in iv]
    assert iv and all(0 <= a < b for a, b in iv), 'invalid_intervals'
    assert iv == sorted(iv), 'intervals_must_be_genomically_sorted'
    assert all(iv[i][1] <= iv[i + 1][0] for i in range(len(iv) - 1)), 'overlapping_intervals'
    return iv


def introns(exons): return [[exons[i][1], exons[i + 1][0]] for i in range(len(exons) - 1)]
def overlap(a, b): return sum(max(0, min(y, v) - max(x, u)) for x, y in a for u, v in b)
def aa_len(cds): return sum(y - x for x, y in cds) // 3 - 1


def phases(cds, strand):
    order = cds if strand == '+' else cds[::-1]; cum = 0; out = {}
    for a, b in order: out[(a, b)] = (3 - cum % 3) % 3; cum += b - a
    return [out[tuple(x)] for x in cds]


def project(exons, strand, start, end):
    """Project transcript-oriented [start, end) to genomic CDS blocks."""
    exons = normalize(exons); out = []; pos = 0
    for a, b in (exons if strand == '+' else exons[::-1]):
        left = max(start, pos); right = min(end, pos + b - a)
        if left < right:
            x, y = left - pos, right - pos
            out.append([a + x, a + y] if strand == '+' else [b - y, b - x])
        pos += b - a
    assert 0 <= start < end <= pos, 'transcript_interval_out_of_bounds'
    return sorted(out)


def utr_intervals(exons, cds, strand):
    """Exonic sequence before/after a continuous spliced ORF; the stop codon belongs to the CDS.
    Observed transcript ends do not prove complete UTRs or TSS/TTS."""
    exons = normalize(exons); cds = normalize(cds)
    first = cds[0][0] if strand == '+' else cds[-1][1] - 1
    start = 0; found = False
    for a, b in (exons if strand == '+' else exons[::-1]):
        if a <= first < b:
            start += (first - a if strand == '+' else b - 1 - first); found = True; break
        start += b - a
    assert found, 'CDS_outside_exons'
    end = start + sum(b - a for a, b in cds); total = sum(b - a for a, b in exons)
    assert project(exons, strand, start, end) == cds, 'CDS_not_continuous_in_spliced_transcript'
    return {'five_prime_UTR': project(exons, strand, 0, start) if start else [],
            'three_prime_UTR': project(exons, strand, end, total) if end < total else [],
            'endpoint_status': 'observed_transcript_ends_not_proven_complete_UTRs_TSS_TTS'}


def find_orfs(seq, min_aa=30):
    """Complete ATG-to-first-in-frame-stop ORFs (every ATG), terminal stop included; longest first."""
    result = []
    for frame in range(3):
        starts = []
        for i in range(frame, len(seq) - 2, 3):
            cod = seq[i:i + 3]
            if cod == 'ATG': starts.append(i)
            if cod in STOPS:
                for start in starts:
                    if (i - start) // 3 >= min_aa: result.append((start, i + 3))
                starts = []
    return sorted(result, key=lambda z: -(z[1] - z[0]))


class Fasta:
    """Random access through the .fai index (no full genome in memory)."""
    def __init__(self, path):
        self.path = Path(path); self.index = {}
        for l in Path(str(path) + '.fai').read_text().splitlines():
            c = l.split(); self.index[c[0]] = list(map(int, c[1:5]))

    def length(self, chrom): return self.index[chrom][0]

    def get(self, chrom, a, b):
        L, off, lb, lw = self.index[chrom]
        assert 0 <= a < b <= L, 'reference_bounds'
        with self.path.open('rb') as f:
            f.seek(off + a // lb * lw + a % lb)
            s = f.read((b - a) // lb * lw + lw + 5).replace(b'\n', b'').replace(b'\r', b'')[:b - a].decode().upper()
        s = re.sub('[^ACGTN]', 'N', s)
        assert len(s) == b - a, 'reference_read_failed'
        return s

    def spliced(self, chrom, iv, strand):
        return ''.join(self.get(chrom, a, b) if strand == '+' else rc(self.get(chrom, a, b))
                       for a, b in (iv if strand == '+' else iv[::-1]))


def gff_attributes(s): return dict(x.split('=', 1) for x in s.strip().split(';') if '=' in x)


def read_gff_models(path):
    """First mRNA of every gene: {gene: {'gene_id','transcript_id','seqid','strand','exons','cds','original_phases'}}.
    Exons default to the CDS when the GFF has none."""
    tx2g, first, models = {}, {}, {}
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9: continue
        at = gff_attributes(c[8]); typ = c[2]
        if typ in ('mRNA', 'transcript'):
            g = at.get('Parent', at['ID']).split(',')[0]; tx2g[at['ID']] = g
            if g not in first:
                first[g] = at['ID']
                models[g] = {'gene_id': g, 'transcript_id': at['ID'], 'seqid': c[0], 'strand': c[6], 'exons': [], 'cds': [], 'original_phases': []}
        elif typ in ('exon', 'CDS'):
            for p in at.get('Parent', '').split(','):
                g = tx2g.get(p)
                if g is None or first[g] != p: continue
                m = models[g]; a, b = int(c[3]) - 1, int(c[4])
                if typ == 'exon': m['exons'].append([a, b])
                else: m['cds'].append([a, b]); m['original_phases'].append(int(c[7]) if c[7] in '012' else 0)
    out = {}
    for g, m in models.items():
        if not m['cds']: continue
        paired = sorted(zip(m['cds'], m['original_phases'])); m['cds'] = [x for x, _ in paired]; m['original_phases'] = [p for _, p in paired]
        ex = sorted(m['exons']) or [list(x) for x in m['cds']]
        merged = []
        for a, b in ex:
            if merged and a < merged[-1][1]: merged[-1][1] = max(merged[-1][1], b)   # overlapping only; adjacency kept
            else: merged.append([a, b])
        m['exons'] = merged; out[g] = m
    return out


def read_gtf(path, source):
    """Transcripts of a GTF: [(chrom, tid, strand, cov, exons)]."""
    tx = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        t = tx.setdefault((f[0], m.group(1)), {'ex': [], 'strand': f[6], 'cov': 0.0})
        t['ex'].append([int(f[3]) - 1, int(f[4])])
        c = re.search(r'cov "([0-9.eE+-]+)"', f[8])
        if c: t['cov'] = max(t['cov'], float(c.group(1)))
    return [(chrom, tid, t['strand'], t['cov'], sorted(t['ex'])) for (chrom, tid), t in tx.items()]
