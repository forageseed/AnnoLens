# annolens-mcp

Launcher and one-step installer for the [AnnoLens](https://github.com/forageseed/AnnoLens) MCP server. AnnoLens reviews and repairs gene models from evidence (transcripts, protein projections, read checks) the way a human curator does.

The engine runs where your data are, as the Python package `annolens`. This npm package only connects your agent host to it, locally or over ssh.

## Install into your agent (one command)

```bash
npx -y annolens-mcp install --host codex  --ssh myserver --manifest /data/proj/annolens.toml
npx -y annolens-mcp install --host claude --ssh myserver --manifest /data/proj/annolens.toml
npx -y annolens-mcp install --host cursor --ssh myserver --manifest /data/proj/annolens.toml
npx -y annolens-mcp install --host generic --ssh myserver --manifest /data/proj/annolens.toml   # prints mcpServers JSON
```

This registers the MCP server and installs the AnnoLens skill. The skill text is read from the `annolens` installed with your data, so it always matches the server's rules. Leave out `--ssh` when the data are on this machine. Existing host configuration is backed up as `*.annolens.bak`. Add `--dry-run` to see what would change.

You can also ask your agent in plain words, for example "install AnnoLens; my data are on myserver at /data/proj/annolens.toml". It will run the command above.

## Check the connection

```bash
npx -y annolens-mcp doctor --ssh myserver --manifest /data/proj/annolens.toml
```

## On the data machine

```bash
pip install git+https://github.com/forageseed/AnnoLens.git
annolens prepare --manifest annolens.toml && annolens scan --manifest annolens.toml
```

Settings can also come from the environment variables `ANNOLENS_MANIFEST`, `ANNOLENS_SSH` and `ANNOLENS_REMOTE_CMD`.

License: Apache-2.0.
