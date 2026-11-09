#!/usr/bin/env node
'use strict';
/*
 * annolens-mcp: npm launcher for the AnnoLens MCP server. Everything runs on this machine; the engine is the Python
 * package `annolens` (pip install git+https://github.com/forageseed/AnnoLens.git).
 *
 *   npx -y annolens-mcp [--manifest annolens.toml]            run the MCP server (stdio) = `annolens mcp`
 *   npx -y annolens-mcp install --host codex|claude|cursor|generic [--manifest annolens.toml] [--dry-run]
 *                                                             register the server in the host and install the skill
 *   npx -y annolens-mcp check [--manifest annolens.toml]      is this machine ready?
 * The project manifest can also come from ANNOLENS_MANIFEST, or ./annolens.toml in the working directory.
 */
const { spawn } = require('child_process');

const argv = process.argv.slice(2);
const sub = ['install', 'check'].includes(argv[0]) ? argv[0] : 'mcp';
const args = sub === 'mcp' ? ['mcp', ...argv] : argv;
const exe = process.env.ANNOLENS_CMD || 'annolens';

const child = spawn(exe, args, { stdio: 'inherit', shell: process.platform === 'win32' });
child.on('error', (e) => {
  console.error(`annolens-mcp: cannot run \`${exe}\` (${e.message}).`);
  console.error('Install the engine first:  python -m pip install git+https://github.com/forageseed/AnnoLens.git');
  process.exit(1);
});
child.on('exit', (code) => process.exit(code ?? 1));
