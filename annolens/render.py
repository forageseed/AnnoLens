"""Evidence image for one view (dict from GeneTools.view_data or LocusTools.view_data) -> PNG bytes on disk."""


def short(name, n=28):
    """Compact label for long gene ids: keep the last two dot-separated fields when present."""
    p = name.split('.')
    s = '.'.join(p[-2:]) if len(p) >= 2 else name
    return s if len(s) <= n else s[-n:]


def render(d, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Arc
    a, b = d['start0'], d['end0']; span = b - a
    rows = []   # (label, blocks_cds, blocks_exon, color, note)
    if d['current']['cds']:
        rows.append(('CURRENT ' + short(d['gene_id']), d['current']['cds'], d['current']['exons'], '#d62728', f"{d['current']['aa']} aa"))
    if d.get('candidate'):
        c = d['candidate']; rows.append(('CANDIDATE', c['cds'], c['exons'], '#ff7f0e', c['candidate_id'][:12]))
    if d.get('dominant_orf'):
        o = d['dominant_orf']; rows.append(('dominant ORF', o['cds'], o['cds'], '#8c564b', f"{o['aa']} aa, {len(o['sources'])} src"))
    for h in d['homology']:
        lab = (h.get('label', '') + ' ' if h.get('label') else '') + {'related_species': 'related-species protein', 'independent_annotation': 'independent-annotation protein'}.get(h['kind'], 'homeolog/paralog protein')
        rows.append((lab, h['blocks'], h['blocks'], '#17becf' if h['kind'] == 'independent_annotation' else '#1f77b4',
                     f"{h.get('protein_aa', '?')} aa, id {h['identity']:.2f} cov {h['coverage']:.2f}"))
    for g, blocks in sorted(d['other_models'].items()):
        rows.append(('model ' + short(g), blocks, blocks, '#7f7f7f', 'other gene'))
    for t in d['transcripts']:
        rows.append((f"{t['label']} {t['source']}", t['orf_cds'], t['exons'], '#2ca02c' if t['source'] == 'IsoSeq' else '#9467bd',
                     f"cov {t['cov']:.0f}, ORF {t['orf_aa']} aa"))
    n = len(rows)
    fig, ax = plt.subplots(figsize=(16, 0.36 * n + 2.6))
    y = n + 2
    for lab, cds, ex, col, note in rows:
        y -= 1
        if ex:
            ax.plot([max(a, ex[0][0]), min(b, ex[-1][1])], [y, y], color=col, lw=0.7)
            for s, e in ex:
                if e > a and s < b: ax.add_patch(Rectangle((s, y - 0.12), e - s, 0.24, color=col, alpha=0.5, lw=0))
            for s, e in cds:
                if e > a and s < b: ax.add_patch(Rectangle((s, y - 0.3), e - s, 0.6, color=col, lw=0))
        ax.text(a - span * 0.005, y, lab, ha='right', va='center', fontsize=7)
        ax.text(b + span * 0.005, y, note, ha='left', va='center', fontsize=6.5)
    y -= 1.8
    mx = max([n for *_, n in d['junctions']] or [1])
    for s, e, c in d['junctions']:
        ax.add_patch(Arc(((s + e) / 2, y), e - s, 1.6, theta1=0, theta2=180, lw=0.3 + 2.5 * c / mx, color='#555555'))
        if c >= 0.05 * mx: ax.text((s + e) / 2, y + 0.85, str(c), ha='center', fontsize=5)
    ax.text(a - span * 0.005, y + 0.3, 'RNA junctions (all samples)', ha='right', va='center', fontsize=7)
    ax.set_xlim(a, b); ax.set_ylim(y - 0.6, n + 1.6); ax.set_yticks([])
    ticks = [a + i * span / 8 for i in range(9)]
    ax.set_xticks(ticks); ax.set_xticklabels([f'{int(t) + 1:,}' for t in ticks], fontsize=7)
    ax.set_title(f"{d['gene_id']}  {d['seqid']}:{a + 1:,}-{b:,} ({d['strand']} strand)   thick = CDS / ORF, thin = exon", fontsize=9)
    for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.13, right=0.88, top=0.95, bottom=0.05)
    fig.savefig(out, dpi=110); plt.close(fig)
    return out
