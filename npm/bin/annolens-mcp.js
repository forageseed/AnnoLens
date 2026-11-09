#!/usr/bin/env node
'use strict';
/*
 * annolens-mcp: npm launcher and installer for the AnnoLens MCP server.
 * The engine is the Python package `annolens` on the machine that holds the data; this launcher only connects to it.
 *
 *   npx -y annolens-mcp --manifest /data/proj/annolens.toml [--ssh myserver] [--remote-cmd annolens]
 *       run as an MCP stdio server (hosts start it themselves after `install`)
 *   npx -y annolens-mcp install --host codex|claude|cursor|generic --manifest ... [--ssh myserver] [--remote-cmd annolens] [--dry-run]
 *       register the MCP server in the host AND install the skill (SKILL.md is read from the installed annolens,
 *       so the skill always matches the server's rules)
 *   npx -y annolens-mcp doctor --manifest ... [--ssh myserver]
 *       check that the server side is reachable and annolens is installed
 * Settings can also come from ANNOLENS_MANIFEST, ANNOLENS_SSH, ANNOLENS_REMOTE_CMD.
 */
const { spawn, spawnSync } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');

function parse(argv) {
  const a = { _: [] };
  for (let i = 0; i < argv.length; i++) {
    const t = argv[i];
    if (t === '--dry-run') a.dryRun = true;
    else if (t.startsWith('--')) { a[t.slice(2).replace(/-([a-z])/g, (_, c) => c.toUpperCase())] = argv[++i]; }
    else a._.push(t);
  }
  a.manifest = a.manifest || process.env.ANNOLENS_MANIFEST;
  a.ssh = a.ssh || process.env.ANNOLENS_SSH;
  a.remoteCmd = a.remoteCmd || process.env.ANNOLENS_REMOTE_CMD || 'annolens';
  return a;
}

const q = (s) => (/^[A-Za-z0-9_./:=@%+-]+$/.test(s) ? s : "'" + String(s).replace(/'/g, "'\\''") + "'");

// command that runs `annolens <args>` where the data are (over ssh, or locally)
function target(a, args) {
  const remote = [a.remoteCmd, ...args.map(q)].join(' ');
  if (a.ssh) return ['ssh', ['-T', '-o', 'BatchMode=yes', '-o', 'ServerAliveInterval=30', a.ssh, remote]];
  const parts = a.remoteCmd.split(/\s+/);
  return [parts[0], [...parts.slice(1), ...args]];
}

function need(a) {
  if (!a.manifest) { console.error('annolens-mcp: --manifest (or ANNOLENS_MANIFEST) is required: path of annolens.toml on the data machine'); process.exit(2); }
}

function serve(a) {
  need(a);
  const [cmd, args] = target(a, ['mcp', '--manifest', a.manifest]);
  const child = spawn(cmd, args, { stdio: ['pipe', 'pipe', 'inherit'] });
  process.stdin.pipe(child.stdin); child.stdout.pipe(process.stdout);
  child.on('error', (e) => { console.error('annolens-mcp: cannot start', cmd, '-', e.message); process.exit(1); });
  child.on('exit', (code) => process.exit(code ?? 1));
}

function launcherArgs(a) {
  const args = ['-y', 'annolens-mcp', '--manifest', a.manifest];
  if (a.ssh) args.push('--ssh', a.ssh);
  if (a.remoteCmd !== 'annolens') args.push('--remote-cmd', a.remoteCmd);
  return args;
}

function write(file, text, dry) {
  if (dry) { console.log(`[dry-run] would write ${file} (${text.length} chars)`); return; }
  fs.mkdirSync(path.dirname(file), { recursive: true });
  if (fs.existsSync(file)) fs.copyFileSync(file, file + '.annolens.bak');
  fs.writeFileSync(file, text); console.log('wrote', file);
}

function skillMd(a) {
  const [cmd, args] = target(a, ['skill-md']);
  const r = spawnSync(cmd, args, { encoding: 'utf8' });
  if (r.status !== 0 || !r.stdout.includes('name: annolens')) {
    console.error('annolens-mcp: could not read the skill from the server side (is `annolens` installed there? try `doctor`)');
    console.error(r.stderr || r.error || ''); process.exit(1);
  }
  return r.stdout;
}

function install(a) {
  need(a);
  const host = a.host;
  if (!['codex', 'claude', 'cursor', 'generic'].includes(host)) { console.error('annolens-mcp: --host must be codex, claude, cursor or generic'); process.exit(2); }
  const md = skillMd(a); const home = os.homedir();
  // Windows hosts cannot spawn npx directly (it is a .cmd shim): go through cmd /c
  const entry = process.platform === 'win32' ? { command: 'cmd', args: ['/c', 'npx', ...launcherArgs(a)] } : { command: 'npx', args: launcherArgs(a) };
  if (host === 'codex') {
    const cfg = path.join(home, '.codex', 'config.toml');
    let old = fs.existsSync(cfg) ? fs.readFileSync(cfg, 'utf8') : '';
    old = old.replace(/^\[mcp_servers\.annolens\]\n(?:(?!\[).*(?:\n|$))*/m, '').trimEnd() + '\n';
    const block = `\n[mcp_servers.annolens]\ncommand = ${JSON.stringify(entry.command)}\nargs = ${JSON.stringify(entry.args)}\nstartup_timeout_sec = 300\ntool_timeout_sec = 900\n`;
    write(cfg, old + block, a.dryRun);
    write(path.join(a.skillsDir || path.join(home, '.codex', 'skills'), 'annolens', 'SKILL.md'), md, a.dryRun);
  } else if (host === 'claude') {
    const reg = ['mcp', 'add', '--scope', 'user', 'annolens', '--', entry.command, ...entry.args];
    console.log((a.dryRun ? '[dry-run] ' : '') + 'claude ' + reg.join(' '));
    if (!a.dryRun) {
      const r = spawnSync('claude', reg, { stdio: 'inherit', shell: process.platform === 'win32' });
      if (r.status !== 0) { console.error('annolens-mcp: `claude mcp add` failed; is Claude Code installed?'); process.exit(1); }
    }
    write(path.join(a.skillsDir || path.join(home, '.claude', 'skills'), 'annolens', 'SKILL.md'), md, a.dryRun);
  } else if (host === 'cursor') {
    const cfg = path.join(home, '.cursor', 'mcp.json');
    const j = fs.existsSync(cfg) ? JSON.parse(fs.readFileSync(cfg, 'utf8') || '{}') : {};
    j.mcpServers = j.mcpServers || {}; j.mcpServers.annolens = entry;
    write(cfg, JSON.stringify(j, null, 2) + '\n', a.dryRun);
    write(path.join(a.skillsDir || path.join(process.cwd(), 'annolens-skill'), 'SKILL.md'), md, a.dryRun);
  } else {
    console.log(JSON.stringify({ mcpServers: { annolens: entry } }, null, 2));
    write(path.join(a.skillsDir || path.join(process.cwd(), 'annolens-skill'), 'SKILL.md'), md, a.dryRun);
  }
  console.log(`done: MCP server "annolens" registered for ${host} and skill installed`);
}

function doctor(a) {
  need(a);
  const [cmd, args] = target(a, ['skill-md']);
  const r = spawnSync(cmd, args, { encoding: 'utf8' });
  const ok = r.status === 0 && r.stdout.includes('name: annolens');
  console.log(ok ? `ok: annolens reachable${a.ssh ? ' on ' + a.ssh : ' locally'}; skill ${r.stdout.match(/AnnoLens ([0-9.]+)/)?.[1] || ''}`
                 : `not ok: ${(r.stderr || String(r.error || '')).trim().slice(0, 400)}\nOn the data machine: pip install git+https://github.com/forageseed/AnnoLens.git`);
  process.exit(ok ? 0 : 1);
}

const a = parse(process.argv.slice(2));
const sub = a._[0];
if (sub === 'install') install(a);
else if (sub === 'doctor') doctor(a);
else if (sub === 'help' || a.help !== undefined) console.log(fs.readFileSync(__filename, 'utf8').split('*/')[0]);
else serve(a);
