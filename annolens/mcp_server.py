"""AnnoLens MCP server: stdio, JSON-RPC 2.0, one message per line, no SDK dependency.

tools      gene sessions (existing models) and locus sessions (missing genes); budgets and every machine rule are
           enforced in annolens.session, never by the host model
prompts    annolens_review_gene, annolens_add_gene (versioned skill text + tool mapping)
resources  annolens://skill/gene, annolens://skill/locus, annolens://server/info
Run: annolens mcp [--manifest annolens.toml] [--out DIR]   (else ANNOLENS_MANIFEST, else ./annolens.toml; without a
           project the server still starts and every tool says how to set one)
"""
import base64, csv, json, sys, traceback, uuid
from pathlib import Path

from . import __version__
from .session import BUDGET, GENE_STATUSES, LOCUS_STATUSES, GeneSession, LocusSession, Workspace, skill_text

PROTOCOL = '2025-06-18'
STATE = {'sessions': {}}
NO_PROJECT = ('No AnnoLens project is set. Point the server at a manifest: set ANNOLENS_MANIFEST=/path/annolens.toml in the '
              'MCP config (annolens install --host <host> --manifest /path/annolens.toml) or start it in the directory that '
              'holds annolens.toml, then reconnect. `annolens check` shows what is missing.')

MAP = {'gene': ('In this MCP session the actions are tools: open = annolens_open_gene(gene_id); view = annolens_gene_view(session,start,end); '
                'orfs = annolens_gene_orfs(session,transcript); build = annolens_gene_build(session,transcript[,start][,merge]); '
                'split_part = annolens_gene_split_part(session, transcript=T# | projection=H#[, start]); '
                'readcheck = annolens_readcheck(session,position); decide = annolens_gene_decide(session,status,...). '
                'Call the tools instead of writing JSON actions. Budgets and validation are enforced by the server.'),
       'locus': ('In this MCP session the actions are tools: open = annolens_open_locus(locus_id); view = annolens_locus_view; '
                 'orfs = annolens_locus_orfs; build = annolens_locus_build(session, transcript=T# | projection=H#[, start]); '
                 'readcheck = annolens_readcheck; decide = annolens_locus_decide. Call the tools instead of writing JSON actions.')}


def content(res):
    out = [{'type': 'text', 'text': json.dumps(res['text'], separators=(',', ':'))}]
    if res.get('image'): out.append({'type': 'image', 'data': base64.b64encode(Path(res['image']).read_bytes()).decode(), 'mimeType': 'image/png'})
    return out


def sess(sid, mode=None):
    s = STATE['sessions'].get(sid)
    if not s or (mode and s.mode != mode): raise ValueError(f'unknown {mode or ""} session {sid}; open one first')
    return s


def open_(mode, key):
    s = (GeneSession if mode == 'gene' else LocusSession)(STATE['ws'], key); res = s.open()
    sid = uuid.uuid4().hex[:12]; STATE['sessions'][sid] = s; res['text'] = {'session': sid, **res['text']}
    return content(res)


def list_queue(name, offset=0, limit=50):
    q = {**{p.stem: p for p in STATE['M'].w('queues').glob('*.txt')}} if STATE['M'].w('queues').exists() else {}
    if STATE['M'].w('loci.tsv').exists(): q['missing_loci'] = STATE['M'].w('loci.tsv')
    if name not in q: raise ValueError(f'queue must be one of {sorted(q)}')
    p = q[name]
    items = [r['locus_id'] for r in csv.DictReader(open(p), delimiter='\t')] if p.suffix == '.tsv' else p.read_text().split()
    sel = items[int(offset):int(offset) + min(200, int(limit))]
    return [{'type': 'text', 'text': json.dumps({'queue': name, 'total': len(items), 'offset': int(offset), 'items': sel})}]


def schema(**props):
    req = [k for k, v in props.items() if not v.pop('optional', False)]
    return {'type': 'object', 'properties': props, 'required': req}


S = lambda d, opt=False: {'type': 'string', 'description': d, **({'optional': True} if opt else {})}
I = lambda d, opt=False: {'type': 'integer', 'description': d, **({'optional': True} if opt else {})}


def decide_schema(statuses, parts=False):
    extra = {'parts': {'type': 'array', 'items': {'type': 'string'}, 'description': 'split only: >= 2 materializable part ids from split_part', 'optional': True}} if parts else {}
    return schema(session=S('session id'), status={'type': 'string', 'enum': sorted(statuses)},
                  candidate_id=S('candidate from a materializable build of this session', True), confidence={'type': 'number', 'description': '0..1'}, **extra,
                  rationale={'type': 'array', 'items': {'type': 'string'}},
                  evidence_gaps={'type': 'array', 'items': {'type': 'string'}, 'optional': True},
                  visual_findings={'type': 'array', 'items': {'type': 'string'}, 'optional': True})


TOOLS = [
    ('annolens_list_queue', 'List IDs of a review queue (from `annolens scan`), paged.', schema(name=S('queue name'), offset=I('start', True), limit=I('max 200', True)),
     lambda a: list_queue(a['name'], a.get('offset', 0), a.get('limit', 50))),
    ('annolens_open_gene', 'Open a review session for an EXISTING gene model: facts, read-checked sites, transcript table, overview evidence image.',
     schema(gene_id=S('gene id from the annotation')), lambda a: open_('gene', a['gene_id'])),
    ('annolens_gene_view', 'View another window (1-based, <= 200 kb).', schema(session=S('session id'), start=I('start'), end=I('end')),
     lambda a: content(sess(a['session'], 'gene').view(a['start'], a['end']))),
    ('annolens_gene_orfs', 'Complete ORFs of an observed transcript (T#).', schema(session=S('session id'), transcript=S('T#')),
     lambda a: content(sess(a['session'], 'gene').orfs(a['transcript']))),
    ('annolens_gene_build', 'Build a candidate from an observed transcript (optional ATG start; merge=true absorbs split neighbours); the engine validates it.',
     schema(session=S('session id'), transcript=S('T#'), start=I('1-based ATG', True), merge={'type': 'boolean', 'description': 'merge neighbours', 'optional': True}),
     lambda a: content(sess(a['session'], 'gene').build(a['transcript'], a.get('start'), a.get('merge', False)))),
    ('annolens_gene_split_part', 'Build ONE gene of a split (the current model joins several genes) from a transcript ORF (transcript=T#) or an intact projection (projection=H#); validated like a new gene.',
     schema(session=S('session id'), transcript=S('T#', True), projection=S('H#', True), start=I('1-based ATG', True)),
     lambda a: content(sess(a['session'], 'gene').split_part(a.get('transcript'), a.get('projection'), a.get('start')))),
    ('annolens_gene_decide', 'Finish: retain | recommend_candidate | split (parts=[ids]) | unresolved | assembly_conflict | genome_frameshift. The server re-validates and writes the review.',
     decide_schema(GENE_STATUSES, parts=True), lambda a: content(sess(a['session'], 'gene').decide(a))),
    ('annolens_open_locus', 'Open a review session for a LOCUS WITHOUT a gene model (possible missing gene).', schema(locus_id=S('locus id from the missing_loci queue')),
     lambda a: open_('locus', a['locus_id'])),
    ('annolens_locus_view', 'View another window of the locus (1-based, <= 200 kb).', schema(session=S('session id'), start=I('start'), end=I('end'), candidate=S('candidate to draw', True)),
     lambda a: content(sess(a['session'], 'locus').view(a['start'], a['end'], a.get('candidate')))),
    ('annolens_locus_orfs', 'Complete ORFs of an observed transcript (T#) at the locus.', schema(session=S('session id'), transcript=S('T#')),
     lambda a: content(sess(a['session'], 'locus').orfs(a['transcript']))),
    ('annolens_locus_build', 'Build a new-gene candidate from a transcript ORF (transcript=T#) or an intact protein projection (projection=H#).',
     schema(session=S('session id'), transcript=S('T#', True), projection=S('H#', True), start=I('1-based ATG', True)),
     lambda a: content(sess(a['session'], 'locus').build(a.get('transcript'), a.get('projection'), a.get('start')))),
    ('annolens_locus_decide', 'Finish: add_gene | no_gene | genome_frameshift | unresolved. The server re-validates and writes the review.',
     decide_schema(LOCUS_STATUSES), lambda a: content(sess(a['session'], 'locus').decide(a))),
    ('annolens_readcheck', 'Genomic long reads and RNA reads that agree/disagree with the assembly at one 1-based position (gene or locus session).',
     schema(session=S('session id'), position=I('1-based position'), kind={'type': 'string', 'enum': ['any', 'frameshift'], 'optional': True}),
     lambda a: content(sess(a['session']).readcheck(a['position'], a.get('kind', 'any')))),
]
FN = {n: f for n, _, _, f in TOOLS}
PROMPTS = {'annolens_review_gene': ('Review one existing gene model with AnnoLens.', 'gene_id', 'gene'),
           'annolens_add_gene': ('Review one locus without a gene model and add the gene if the evidence supports it.', 'locus_id', 'locus')}


def resources():
    return {'annolens://skill/gene': ('AnnoLens skill for existing gene models', lambda: skill_text('gene')),
            'annolens://skill/locus': ('AnnoLens skill for loci without a model (add gene)', lambda: skill_text('locus')),
            'annolens://server/info': ('Server summary', lambda: json.dumps({'version': __version__, 'project': STATE['M'].name if STATE.get('M') else None, 'budgets': BUDGET}, indent=1))}


def handle(msg):
    m, p = msg.get('method'), msg.get('params') or {}
    if m == 'initialize':
        return {'protocolVersion': p.get('protocolVersion', PROTOCOL), 'serverInfo': {'name': 'annolens', 'version': __version__},
                'capabilities': {'tools': {}, 'prompts': {}, 'resources': {}},
                'instructions': 'AnnoLens gene-model curation. Use the prompts annolens_review_gene / annolens_add_gene or read annolens://skill/*. '
                                'The server enforces budgets and validation and never edits the genome or the input annotation.'}
    if m == 'ping': return {}
    if m == 'tools/list': return {'tools': [{'name': n, 'description': d, 'inputSchema': json.loads(json.dumps(s))} for n, d, s, _ in TOOLS]}
    if m == 'tools/call':
        name, args = p.get('name'), p.get('arguments') or {}
        if name not in FN: raise KeyError(f'unknown tool {name}')
        if STATE.get('M') is None:
            return {'content': [{'type': 'text', 'text': NO_PROJECT}], 'isError': True}
        try: return {'content': FN[name](args), 'isError': False}
        except Exception as ex:
            traceback.print_exc(file=sys.stderr); return {'content': [{'type': 'text', 'text': f'TOOL ERROR {type(ex).__name__}: {ex}'}], 'isError': True}
    if m == 'prompts/list':
        return {'prompts': [{'name': n, 'description': d, 'arguments': [{'name': arg, 'description': arg, 'required': True}]} for n, (d, arg, _) in PROMPTS.items()]}
    if m == 'prompts/get':
        d, arg, mode = PROMPTS[p['name']]; v = (p.get('arguments') or {}).get(arg, '')
        text = MAP[mode] + '\n\n' + skill_text(mode) + f'\n\nStart with annolens_open_{mode}({arg}="{v}").'
        return {'description': d, 'messages': [{'role': 'user', 'content': {'type': 'text', 'text': text}}]}
    if m == 'resources/list': return {'resources': [{'uri': u, 'name': n, 'mimeType': 'text/plain'} for u, (n, _) in resources().items()]}
    if m == 'resources/read': u = p['uri']; return {'contents': [{'uri': u, 'mimeType': 'text/plain', 'text': resources()[u][1]()}]}
    raise KeyError(f'method not found: {m}')


def serve(manifest, out):
    out_stream = sys.stdout; sys.stdout = sys.stderr      # library output never reaches the JSON-RPC channel
    STATE.update(M=manifest, ws=Workspace(manifest, out) if manifest else None)
    for line in sys.stdin:
        line = line.strip()
        if not line: continue
        try: msg = json.loads(line)
        except json.JSONDecodeError:
            out_stream.write(json.dumps({'jsonrpc': '2.0', 'id': None, 'error': {'code': -32700, 'message': 'parse error'}}) + '\n'); out_stream.flush(); continue
        if 'id' not in msg: continue
        try: resp = {'jsonrpc': '2.0', 'id': msg['id'], 'result': handle(msg)}
        except KeyError as ex: resp = {'jsonrpc': '2.0', 'id': msg['id'], 'error': {'code': -32601, 'message': str(ex)}}
        except Exception as ex:
            traceback.print_exc(file=sys.stderr); resp = {'jsonrpc': '2.0', 'id': msg['id'], 'error': {'code': -32603, 'message': f'{type(ex).__name__}: {ex}'}}
        out_stream.write(json.dumps(resp) + '\n'); out_stream.flush()
