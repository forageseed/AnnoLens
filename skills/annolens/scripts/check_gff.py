#!/usr/bin/env python3
"""Deterministic checks of gene models, no MCP needed: for every mRNA in a GFF3, against the genome FASTA:
ATG start, terminal stop codon, CDS length % 3, internal stops, CDS phases, splice motifs (GT-AG, GC-AG, AT-AC).
Only the Python standard library is used. Usage: check_gff.py genes.gff3 genome.fa [--tsv out.tsv]
"""
import sys
from collections import defaultdict

STOPS = {'TAA', 'TAG', 'TGA'}
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
rc = lambda s: s.translate(comp)[::-1]


def read_fasta(path, need):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line[0] == '>':
            if name in need: seqs[name] = ''.join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif name in need: buf.append(line.strip())
    if name in need: seqs[name] = ''.join(buf).upper()
    return seqs


def main():
    gff, fa = sys.argv[1], sys.argv[2]; out = sys.argv[sys.argv.index('--tsv') + 1] if '--tsv' in sys.argv else None
    tx, cds = {}, defaultdict(list)
    for line in open(gff):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] in ('mRNA', 'transcript'): tx[at['ID']] = (f[0], f[6])
        elif f[2] == 'CDS':
            for p in at.get('Parent', '').split(','): cds[p].append((int(f[3]) - 1, int(f[4]), f[7]))
    seqs = read_fasta(fa, {c for c, _ in tx.values()})
    rows = []
    for t, (c, st) in tx.items():
        bl = sorted(cds.get(t, []))
        if not bl: rows.append((t, 'no_CDS')); continue
        order = bl if st == '+' else bl[::-1]; g = seqs[c]
        s = ''.join(g[a:b] if st == '+' else rc(g[a:b]) for a, b, _ in order)
        prob = []
        if s[:3] != 'ATG': prob.append('no_ATG')
        if s[-3:] not in STOPS: prob.append('no_terminal_stop')
        if len(s) % 3: prob.append('length_not_mod3')
        if any(s[i:i + 3] in STOPS for i in range(0, len(s) - 3, 3)): prob.append('internal_stop')
        cum = 0
        for k, (a, b, ph) in enumerate(order):
            if ph in '012' and int(ph) != (3 - cum % 3) % 3: prob.append(f'phase_cds{k + 1}'); break
            cum += b - a
        for (a1, b1, _), (a2, b2, _) in zip(bl, bl[1:]):
            m = g[b1:b1 + 2] + '-' + g[a2 - 2:a2] if st == '+' else rc(g[a2 - 2:a2]) + '-' + rc(g[b1:b1 + 2])
            if m not in ('GT-AG', 'GC-AG', 'AT-AC'): prob.append(f'splice_{m}')
        rows.append((t, ';'.join(prob) or 'ok', len(s) // 3 - (1 if s[-3:] in STOPS else 0)))
    ok = sum(r[1] == 'ok' for r in rows)
    print(f'{ok}/{len(rows)} transcripts pass')
    for r in rows:
        if r[1] != 'ok': print('\t'.join(map(str, r)))
    if out:
        with open(out, 'w') as fh:
            fh.write('transcript\tresult\tprotein_aa\n')
            for r in rows: fh.write('\t'.join(map(str, r)) + ('\t\n' if len(r) == 2 else '\n'))


if __name__ == '__main__':
    main()
