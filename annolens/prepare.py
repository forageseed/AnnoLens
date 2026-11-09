"""`annolens prepare`: build the work directory from the manifest. Every step is skipped when its output exists.

  transcripts/<sample>.gtf     copied/linked from the manifest, or assembled with StringTie from the BAM
  junctions/<sample>.bed       spliced-alignment counts per intron (primary, MAPQ >= 10): chrom start0 end0 . count strand
  proteins/<set>.miniprot.gff  miniprot alignments of every protein set to the genome
  homology.jsonl               per gene, the best projection of every protein set overlapping its CDS (>= 30 bp)
  readcheck/config.json        configuration of the single-position read checker
  readcheck/sites.jsonl        frameshift / in-frame-stop sites that independent proteins need, with read counts
"""
import collections, json, os, shutil, subprocess, sys
from pathlib import Path

from . import genome as G


def log(*a): print('[prepare]', *a, file=sys.stderr, flush=True)


def _pysam(M):
    try:
        import pysam; return pysam
    except ImportError:
        sys.path.insert(0, str(M.w('deps'))); import pysam; return pysam


def transcripts(M, threads):
    d = M.w('transcripts'); d.mkdir(parents=True, exist_ok=True)
    items = [(r['name'], r['bam'], r['gtf'], False) for r in M.rna] + ([(M.longread['name'], M.longread['bam'], M.longread['gtf'], True)] if M.longread else [])
    for name, bam, gtf, lr in items:
        out = d / f'{name}.gtf'
        if out.exists(): continue
        if gtf and gtf.exists():
            try: os.symlink(gtf, out)
            except OSError: shutil.copy(gtf, out)
        elif bam and bam.exists():
            log('StringTie', name)
            cmd = [M.tools['stringtie'], str(bam), '-p', str(threads), '-o', str(out)] + (['-L'] if lr else [])
            subprocess.run(cmd, check=True)
        else:
            raise FileNotFoundError(f'{name}: needs gtf or bam')


def junctions(M):
    pysam = _pysam(M); d = M.w('junctions'); d.mkdir(parents=True, exist_ok=True)
    for r in M.rna:
        out = d / f"{r['name']}.bed"
        if out.exists() or not (r['bam'] and r['bam'].exists()): continue
        log('junctions', r['name']); cnt = collections.Counter()
        with pysam.AlignmentFile(str(r['bam']), 'rb') as f:
            for x in f.fetch(until_eof=True):
                if x.is_unmapped or x.is_secondary or x.is_supplementary or x.mapping_quality < 10 or not x.cigartuples: continue
                st = '.'
                if x.has_tag('XS'): st = x.get_tag('XS')
                elif x.has_tag('ts'): st = {'+': '-', '-': '+'}[x.get_tag('ts')] if x.is_reverse else x.get_tag('ts')
                pos = x.reference_start
                for op, n in x.cigartuples:
                    if op == 3: cnt[(x.reference_name, pos, pos + n, st)] += 1
                    if op in (0, 2, 3, 7, 8): pos += n
        tmp = Path(str(out) + '.tmp')
        with tmp.open('w') as fh:
            for (c, a, b, s), n in sorted(cnt.items()): fh.write(f'{c}\t{a}\t{b}\t.\t{n}\t{s}\n')
        tmp.replace(out)


def proteins(M, threads):
    d = M.w('proteins'); d.mkdir(parents=True, exist_ok=True); idx = M.w('proteins', 'genome.mpi')
    need = [p for p in M.proteins if not (d / f"{p['name']}.miniprot.gff").exists()]
    if need and not idx.exists():
        log('miniprot index'); subprocess.run([M.tools['miniprot'], '-t', str(threads), '-d', str(idx), str(M.genome)], check=True)
    for p in M.proteins:
        out = d / f"{p['name']}.miniprot.gff"
        if out.exists(): continue
        log('miniprot', p['name']); tmp = Path(str(out) + '.tmp')
        with tmp.open('w') as fh:
            subprocess.run([M.tools['miniprot'], '-t', str(threads), '--gff', '-N', '3', str(idx), str(p['fasta'])], stdout=fh, check=True)
        tmp.replace(out)


def homology(M):
    out = M.w('homology.jsonl')
    if out.exists(): return
    models = G.read_gff_models(M.annotation); by = collections.defaultdict(list)
    for g, m in models.items(): by[(m['seqid'], m['strand'])].append((m['cds'][0][0], m['cds'][-1][1], g))
    for v in by.values(): v.sort()
    import bisect
    starts = {k: [x[0] for x in v] for k, v in by.items()}; maxlen = max((b - a for v in by.values() for a, b, _ in v), default=0)
    best = collections.defaultdict(dict)
    for p in M.proteins:
        L, cur = {}, None
        for line in open(p['fasta']):
            if line[0] == '>': cur = line[1:].split()[0]; L[cur] = 0
            elif cur: L[cur] += len(line.strip())
        aln, cur = [], None
        for line in open(M.w('proteins', f"{p['name']}.miniprot.gff")):
            if line.startswith('#') or not line.strip(): continue
            f = line.rstrip('\n').split('\t')
            if f[2] == 'mRNA':
                at = G.gff_attributes(f[8]); t = at['Target'].split()
                cur = {'alignment_id': p['name'] + ':' + at['ID'], 'kind': p['kind'], 'protein': p['name'] + '|' + t[0], 'protein_length': L.get(t[0], 0),
                       'target_start_aa': int(t[1]), 'target_end_aa': int(t[2]), 'identity': float(at.get('Identity', 0)),
                       'positive': float(at.get('Positive', 0)), 'frameshift': int(at.get('Frameshift', 0)), 'stop_codon_found': int(at.get('StopCodon', 0)),
                       'seqid': f[0], 'strand': f[6], 'cds_blocks': [], '_score': int(f[5]) if f[5] != '.' else 0}
                aln.append(cur)
            elif f[2] == 'CDS' and cur is not None: cur['cds_blocks'].append([int(f[3]) - 1, int(f[4])])
        for al in aln:
            if not al['cds_blocks']: continue
            al['cds_blocks'].sort(); s, e = al['cds_blocks'][0][0], al['cds_blocks'][-1][1]; key = (al['seqid'], al['strand'])
            v = by.get(key, []); i = bisect.bisect_left(starts.get(key, []), s - maxlen)
            while i < len(v) and v[i][0] < e:
                gs, ge, g = v[i]; i += 1
                if ge <= s or G.overlap(al['cds_blocks'], models[g]['cds']) < 30: continue
                if p['name'] not in best[g] or al['_score'] > best[g][p['name']]['_score']: best[g][p['name']] = al
    tmp = Path(str(out) + '.tmp')
    with tmp.open('w') as fh:
        for g in sorted(best):
            fh.write(json.dumps({'gene_id': g, 'alignments': [{k: v for k, v in al.items() if k != '_score'} for al in best[g].values()]}) + '\n')
    tmp.replace(out); log('homology genes', len(best))


def readcheck_setup(M):
    d = M.w('readcheck'); d.mkdir(parents=True, exist_ok=True)
    cfg = {'genome': str(M.genome), 'hifi_bam': str(M.hifi_bam) if M.hifi_bam else None,
           'rna_bams': {r['name']: str(r['bam']) for r in M.rna if r['bam'] and r['bam'].exists()},
           'miniprot': M.tools['miniprot'], 'pysam_dir': str(M.w('deps'))}
    (d / 'config.json').write_text(json.dumps(cfg, indent=1))
    if M.hifi_bam and not (d / 'sites.jsonl').exists():
        from . import readcheck
        faa = {p['name']: p['fasta'] for p in M.proteins if p['role'] == 'independent'}
        if faa:
            log('read-check sites'); readcheck.precompute(d / 'config.json', M.w('homology.jsonl'), faa, d / 'sites.jsonl')


def run(M, threads=8):
    M.check(); M.work.mkdir(parents=True, exist_ok=True)
    transcripts(M, threads); junctions(M); proteins(M, threads); homology(M); readcheck_setup(M)
    log('work dir ready:', M.work)
