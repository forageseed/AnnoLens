"""Batch agent: a vision-language model reviews genes or loci through GeneSession / LocusSession.

Model access is configured in the manifest [model] table and works with any OpenAI-compatible (chat/completions)
or Anthropic-compatible (messages) endpoint:
  provider = "openai" | "anthropic", base_url, model, api_key_env (environment variable holding the key),
  max_tokens (default 1500), extra = {...} (merged into the request body, e.g. to disable reasoning)
The key is read from the environment only; it is never written to disk or logs.
"""
import base64, json, os, sys, time, traceback, urllib.error, urllib.request
from multiprocessing import Pool
from pathlib import Path

from .session import GeneSession, LocusSession, Workspace, skill_text

MAX_REPLIES = 7
W = {}


class Model:
    def __init__(self, cfg):
        self.cfg = cfg; self.provider = cfg.get('provider', 'openai'); self.model = cfg['model']
        self.base = cfg['base_url'].rstrip('/'); self.key = os.environ.get(cfg.get('api_key_env', 'ANNOLENS_API_KEY'))
        if not self.key: raise SystemExit(f"environment variable {cfg.get('api_key_env', 'ANNOLENS_API_KEY')} is not set")
        self.max_tokens = int(cfg.get('max_tokens', 1500)); self.extra = cfg.get('extra', {}); self.last = 0.0

    @staticmethod
    def _anthropic_messages(messages):
        system = '\n'.join(m['content'] for m in messages if m['role'] == 'system'); out = []
        for m in messages:
            if m['role'] == 'system': continue
            parts = m['content'] if isinstance(m['content'], list) else [{'type': 'text', 'text': m['content']}]
            conv = []
            for p in parts:
                if p['type'] == 'text': conv.append({'type': 'text', 'text': p['text'] or ' '})
                else:
                    head, data = p['image_url']['url'].split(',', 1)
                    conv.append({'type': 'image', 'source': {'type': 'base64', 'media_type': head[5:].split(';')[0], 'data': data}})
            if out and out[-1]['role'] == m['role']: out[-1]['content'] += conv
            else: out.append({'role': m['role'], 'content': conv})
        return system, out

    def ask(self, messages):
        gap = 1.0 - (time.time() - self.last)
        if gap > 0: time.sleep(gap)
        if self.provider == 'anthropic':
            system, msgs = self._anthropic_messages(messages)
            url, body = self.base + '/v1/messages', {'model': self.model, 'system': system, 'messages': msgs, 'temperature': 0, 'max_tokens': self.max_tokens}
            headers = {'x-api-key': self.key, 'anthropic-version': '2023-06-01'}
        else:
            url, body = self.base + '/chat/completions', {'model': self.model, 'messages': messages, 'temperature': 0, 'max_tokens': self.max_tokens}
            headers = {'Authorization': 'Bearer ' + self.key}
        body.update(self.extra)
        tries = 6
        for attempt in range(tries):
            wait = min(60, 3 * 2 ** attempt)
            try:
                req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={**headers, 'Content-Type': 'application/json'})
                with urllib.request.urlopen(req, timeout=180) as r:
                    self.last = time.time(); d = json.load(r)
                if self.provider == 'anthropic':
                    u = d.get('usage', {}); return ''.join(c.get('text', '') for c in d.get('content', []) if c.get('type') == 'text'), u.get('input_tokens', 0) + u.get('output_tokens', 0)
                return d['choices'][0]['message'].get('content') or '', d.get('usage', {}).get('total_tokens', 0)
            except urllib.error.HTTPError as ex:
                if ex.code in (401, 402, 403): raise RuntimeError(f'provider auth/balance HTTP {ex.code}') from None
                if attempt == tries - 1: raise RuntimeError(f'provider HTTP {ex.code}') from None
                ra = ex.headers.get('Retry-After') if ex.headers else None   # rate limits (429) / overload (5xx)
                if ra and ra.isdigit(): wait = min(120, int(ra))
            except (TimeoutError, urllib.error.URLError):
                if attempt == tries - 1: raise RuntimeError('provider timeout') from None
            time.sleep(wait)


def _user(res, prefix=''):
    text = prefix + '\n'.join(f'{k} ' + json.dumps(v, separators=(',', ':')) for k, v in res['text'].items())
    if not res.get('image'): return {'role': 'user', 'content': text}
    b64 = base64.b64encode(Path(res['image']).read_bytes()).decode()
    return {'role': 'user', 'content': [{'type': 'text', 'text': text}, {'type': 'image_url', 'image_url': {'url': 'data:image/png;base64,' + b64}}]}


def _parse(raw):
    body = raw.strip().strip('`').removeprefix('json').strip(); i = body.find('{')
    act, end = json.JSONDecoder().raw_decode(body[i:]) if i >= 0 else (json.loads(body), 0)
    if not isinstance(act, dict): raise ValueError('not an object')
    return act, body[i + end:].strip()


def review(key):
    ws, model, mode = W['ws'], W['model'], W['mode']
    s = (GeneSession if mode == 'gene' else LocusSession)(ws, key); tokens = 0
    messages = [{'role': 'system', 'content': skill_text(mode)}, _user(s.open())]
    for turn in range(MAX_REPLIES):
        if turn == MAX_REPLIES - 1: messages.append({'role': 'user', 'content': 'Budget reached: reply now with {"action":"decide",...}.'})
        raw, t = model.ask(messages); tokens += t
        messages.append({'role': 'assistant', 'content': raw})
        for m in messages[1:-1]:          # keep only the newest image in context
            if isinstance(m.get('content'), list): m['content'] = [c for c in m['content'] if c['type'] == 'text'] + [{'type': 'text', 'text': '[earlier image omitted]'}]
        try:
            act, rest = _parse(raw)
            if rest: s.trace.append({'turn': turn, 'trailing_text_ignored': rest[:300]}); messages[-1] = {'role': 'assistant', 'content': json.dumps(act)}
        except Exception:
            s.trace.append({'turn': turn, 'invalid_json': raw[:300]}); messages.append({'role': 'user', 'content': 'Reply with ONE valid JSON object.'}); continue
        s.trace.append({'turn': turn, 'action': act}); a = act.get('action')
        try:
            if a == 'view': res = s.view(act['start'], act['end'])
            elif a == 'orfs': res = s.orfs(act['transcript'])
            elif a == 'readcheck': res = s.readcheck(act['position'], act.get('kind', 'any'))
            elif a == 'build':
                res = s.build(act['transcript'], act.get('start'), act.get('merge', False)) if mode == 'gene' else \
                    s.build(act.get('transcript'), act.get('projection'), act.get('start'))
            elif a == 'decide': res = s.decide(act); return key, res['text']['FINAL']['status'], tokens
            else: raise ValueError('unknown action')
            messages.append(_user(res))
        except Exception as ex:
            messages.append({'role': 'user', 'content': f'TOOL ERROR {type(ex).__name__}: {ex}'})
    res = s.decide({'status': 'unresolved', 'confidence': 0, 'rationale': ['No valid decision within budget.']})
    return key, res['text']['FINAL']['status'], tokens


def _safe(key):
    try: return review(key)
    except Exception:
        p = W['ws'].out / W['mode'] / 'failures' / f"{key.replace(':', '_')}.txt"; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(traceback.format_exc())
        return key, 'failed', 0


def _init(manifest_path, out, mode):
    from .config import Manifest
    M = Manifest(manifest_path); W.update(ws=Workspace(M, out), model=Model(M.model), mode=mode)


def run(manifest_path, mode, items, out, workers=4):
    done_dir = Path(out) / mode / 'reviews'
    done = {p.stem for p in done_dir.glob('*.json')} if done_dir.exists() else set()
    todo = [k for k in items if k.replace(':', '_') not in done]
    print(json.dumps({'mode': mode, 'items': len(items), 'todo': len(todo)}), flush=True); tot = 0
    with Pool(workers, initializer=_init, initargs=(manifest_path, out, mode)) as pool:
        for k, st, t in pool.imap_unordered(_safe, todo):
            tot += t; print(json.dumps({'id': k, 'status': st, 'tokens': t}), flush=True)
    print(json.dumps({'done': len(todo), 'total_tokens': tot}), flush=True)
