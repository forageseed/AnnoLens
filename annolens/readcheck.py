#!/usr/bin/env python3
"""Read-level check of single genome positions. Shared by `annolens prepare` (precompute) and the agent tool.

check(chrom, pos, kind) counts, in a +-WIN bp window around a 1-based position:
  HiFi  primary MAPQ>=20 reads that span the window; 'alt' = read has an insertion/deletion in the window
        (kind frameshift) or, for kind in_frame_stop/any, also a base different from the reference
  RNA   the same for every tissue BAM; reads whose alignment skips (intron) over the window are not counted
Calls (HiFi decides; RNA is reported as support):
  assembly_error_candidate   HiFi depth >= MIN_DEPTH and alt >= 50 %   (the assembly disagrees with its own reads)
  genome_sequence_supported  HiFi depth >= MIN_DEPTH and alt <= 10 %   (the reads confirm the assembly)
  ambiguous                  otherwise
Sites of an independent-annotation protein are located with miniprot --aln ('!' or '$' = frameshift,
internal '*' = in-frame stop).

CLI (precompute for every packet gene with an 'ind' alignment that has a frameshift or stop):
  readcheck.py precompute CONFIG.json HOMOLOGY.jsonl IND.faa OUT.jsonl
CONFIG.json: {"genome": fa, "hifi_bam": bam, "rna_bams": {tissue: bam}, "miniprot": path, "pysam_dir": dir}
"""
import json, re, subprocess, sys, tempfile
from pathlib import Path

WIN, MIN_DEPTH, FLANK = 6, 8, 2000


class ReadCheck:
    def __init__(self, config):
        cfg = json.loads(Path(config).read_text()) if not isinstance(config, dict) else config
        if cfg.get('pysam_dir') and cfg['pysam_dir'] not in sys.path: sys.path.insert(0, cfg['pysam_dir'])
        import pysam
        self.cfg = cfg; self.pysam = pysam
        self.ref = pysam.FastaFile(cfg['genome'])
        self.hifi = pysam.AlignmentFile(cfg['hifi_bam'], 'rb') if cfg.get('hifi_bam') else None
        self.rna = {t: pysam.AlignmentFile(b, 'rb') for t, b in cfg.get('rna_bams', {}).items()}

    def _count(self, bf, chrom, a, b, mismatch, mapq, pos0):
        """a, b: 0-based half-open window around pos0 (0-based). Returns (depth, alt). A spliced read is excluded
        only when its intron covers pos0 itself; sites at exon edges keep their spliced RNA reads."""
        refseq = self.ref.fetch(chrom, a, b).upper(); tot = alt = 0
        for r in bf.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < mapq: continue
            if r.reference_start > a or r.reference_end < b: continue
            rp, hit, skipped = r.reference_start, False, False
            for op, n in r.cigartuples:
                if op == 3 and rp <= pos0 < rp + n: skipped = True
                if op in (1, 2) and a <= rp <= b: hit = True
                if op in (0, 2, 3, 7, 8): rp += n
            if skipped: continue
            if mismatch and not hit:
                q = r.query_sequence
                for qp, rpos in r.get_aligned_pairs(matches_only=True):
                    if a <= rpos < b and q[qp].upper() != refseq[rpos - a]: hit = True; break
            tot += 1; alt += hit
        return tot, alt

    def check(self, chrom, pos, kind='frameshift'):
        a, b = max(0, pos - 1 - WIN), pos + WIN
        mism = kind != 'frameshift'
        hd, ha = self._count(self.hifi, chrom, a, b, mism, 20, pos - 1) if self.hifi else (0, 0)
        rna = {t: self._count(bf, chrom, a, b, mism, 10, pos - 1) for t, bf in self.rna.items()}
        rd, ra = sum(x[0] for x in rna.values()), sum(x[1] for x in rna.values())
        frac = ha / hd if hd else 0.0
        call = ('ambiguous' if hd < MIN_DEPTH else 'assembly_error_candidate' if frac >= .5
                else 'genome_sequence_supported' if frac <= .1 else 'ambiguous')
        return {'position_1based': pos, 'kind': kind, 'window_bp': [a + 1, b],
                'context': self.ref.fetch(chrom, max(0, pos - 9), pos + 8).upper(),
                'hifi_reads': hd, 'hifi_reads_disagreeing': ha, 'rna_reads': rd, 'rna_reads_disagreeing': ra,
                'call': call}


def locate_sites(miniprot, region_seq, chrom, rstart, protein_id, protein_seq):
    """miniprot --aln of one protein to one region; returns ([(kind, pos_1based, protein_pos)], attrs)."""
    with tempfile.TemporaryDirectory() as d:
        Path(d, 'r.fa').write_text(f'>{chrom}\n{region_seq}\n'); Path(d, 'p.fa').write_text(f'>{protein_id}\n{protein_seq}\n')
        txt = subprocess.run([miniprot, '--gff', '--aln', '-t', '1', f'{d}/r.fa', f'{d}/p.fa'], capture_output=True, text=True).stdout
    lines = txt.splitlines(); mrna = [l for l in lines if '\tmRNA\t' in l]
    if not mrna: return [], None
    f = mrna[0].split('\t'); ms, me, st = int(f[3]), int(f[4]), f[6]
    atn = next((l.split('\t', 1)[1] for l in lines if l.startswith('##ATN')), None)
    ata = next((l.split('\t', 1)[1] for l in lines if l.startswith('##ATA')), None)
    if atn is None or ata is None: return [], f[8]
    ext_n, ext_a = [], []
    for m in re.finditer(r'~(\d+)~|.', atn):
        tok = m.group(0); j = m.start()
        if tok.startswith('~'): ext_n.append('n' * int(m.group(1))); ext_a.append(' ' * int(m.group(1)))
        else: ext_n.append(tok); ext_a.append(ata[j:j + 1] or ' ')
    atn2, ata2 = ''.join(ext_n), ''.join(ext_a)
    last = len(ata2.rstrip()) - 3
    res, nbase, aa = [], 0, 0
    for col, (cn, ca) in enumerate(zip(atn2, ata2)):
        if ca.isalpha(): aa += 1
        g = (ms + nbase) if st == '+' else (me - nbase)
        if ca in '!$': res.append(('frameshift', rstart + g - 1, aa))
        elif ca == '*' and col < last: res.append(('in_frame_stop', rstart + g - 1, aa))
        if cn != '-': nbase += 1
    return res, f[8]


def precompute(config, homology, faa, out):
    """faa: {protein_set_name: fasta}; homology protein ids are 'set|id'."""
    rc = ReadCheck(config)
    prot, cur = {}, None
    for name, path in faa.items():
        for line in open(path):
            if line[0] == '>': cur = name + '|' + line[1:].split()[0]; prot[cur] = []
            elif cur: prot[cur].append(line.strip())
    n_genes = n_sites = 0
    with open(out, 'w') as fh:
        for line in open(homology):
            rec = json.loads(line); sites = []
            for al in rec['alignments']:
                if al['kind'] != 'ind' or al['frameshift'] + al['stop_codon_found'] == 0 or not al['cds_blocks']: continue
                q = al['protein']
                if q not in prot: continue
                s, e = al['cds_blocks'][0][0], al['cds_blocks'][-1][1]
                rs = max(1, s + 1 - FLANK); seq = rc.ref.fetch(al['seqid'], rs - 1, e + FLANK)
                found, _ = locate_sites(rc.cfg['miniprot'], seq, al['seqid'], rs, q, ''.join(prot[q]))
                for kind, pos, ppos in found:
                    sites.append({'protein': al['protein'], 'protein_length': al['protein_length'], 'protein_pos': ppos,
                                  **rc.check(al['seqid'], pos, kind)})
            if sites:
                fh.write(json.dumps({'gene_id': rec['gene_id'], 'sites': sites}) + '\n'); n_genes += 1; n_sites += len(sites)
    print(json.dumps({'genes_with_sites': n_genes, 'sites': n_sites}))

