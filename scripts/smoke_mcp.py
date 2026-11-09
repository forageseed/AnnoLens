#!/usr/bin/env python3
"""Minimal MCP stdio client: scripted tool calls against an AnnoLens server, no model involved.
Usage: smoke_mcp.py "<server command>" GENE_ID LOCUS_ID [READCHECK_POSITION]"""
import json, shlex, subprocess, sys

cmd, gene, locus = sys.argv[1], sys.argv[2], sys.argv[3]; pos = int(sys.argv[4]) if len(sys.argv) > 4 else None
p = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
n = 0


def call(method, params=None, notify=False):
    global n
    msg = {'jsonrpc': '2.0', 'method': method, **({'params': params} if params is not None else {})}
    if not notify: n += 1; msg['id'] = n
    p.stdin.write(json.dumps(msg) + '\n'); p.stdin.flush()
    if notify: return None
    r = json.loads(p.stdout.readline())
    if 'error' in r: print('ERROR', method, r['error']); return None
    return r['result']


def show(res, label):
    parts = []
    for c in res['content']:
        parts.append(c['text'][:260] if c['type'] == 'text' else f"<image {len(c['data'])} b64 chars>")
    print(label, 'isError=' + str(res['isError']), '|', ' || '.join(parts))


r = call('initialize', {'protocolVersion': '2025-06-18', 'capabilities': {}, 'clientInfo': {'name': 'smoke', 'version': '0'}})
print('init', r['serverInfo'], r['protocolVersion']); call('notifications/initialized', notify=True)
print('tools', [t['name'] for t in call('tools/list')['tools']])
print('prompts', [x['name'] for x in call('prompts/list')['prompts']])
print('resources', [x['uri'] for x in call('resources/list')['resources']])
pr = call('prompts/get', {'name': 'annolens_review_gene', 'arguments': {'gene_id': gene}})
print('prompt chars', len(pr['messages'][0]['content']['text']))
show(call('tools/call', {'name': 'annolens_list_queue', 'arguments': {'name': 'P1', 'limit': 3}}), 'queue')
res = call('tools/call', {'name': 'annolens_open_gene', 'arguments': {'gene_id': gene}}); show(res, 'open_gene')
sid = json.loads(res['content'][0]['text'])['session']
show(call('tools/call', {'name': 'annolens_gene_decide', 'arguments': {'session': sid, 'status': 'recommend_candidate', 'candidate_id': 'nope', 'confidence': 0.9, 'rationale': ['test']}}), 'decide_without_build(expect error)')
res = call('tools/call', {'name': 'annolens_gene_build', 'arguments': {'session': sid, 'transcript': 'T1'}}); show(res, 'build T1')
if pos: show(call('tools/call', {'name': 'annolens_readcheck', 'arguments': {'session': sid, 'position': pos}}), 'readcheck')
show(call('tools/call', {'name': 'annolens_gene_decide', 'arguments': {'session': sid, 'status': 'unresolved', 'confidence': 0.1, 'rationale': ['smoke test only']}}), 'decide')
res = call('tools/call', {'name': 'annolens_open_locus', 'arguments': {'locus_id': locus}}); show(res, 'open_locus')
lsid = json.loads(res['content'][0]['text'])['session']
show(call('tools/call', {'name': 'annolens_locus_build', 'arguments': {'session': lsid, 'transcript': 'T3'}}), 'locus build T3')
show(call('tools/call', {'name': 'annolens_locus_decide', 'arguments': {'session': lsid, 'status': 'unresolved', 'confidence': 0.1, 'rationale': ['smoke test only']}}), 'locus decide')
p.stdin.close(); p.wait(timeout=60)
