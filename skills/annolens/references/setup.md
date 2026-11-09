# Connecting AnnoLens

AnnoLens runs as the Python package `annolens` on the machine that holds the data; the agent host reaches it through
the MCP server, locally or over ssh.

1. Ask the user where the data are: the ssh host of the data machine (or "this machine") and the path of the project
   manifest `annolens.toml` there. Do not guess paths.
2. Register the server and install this skill in one step:
   `npx -y annolens-mcp install --host <codex|claude|cursor|generic> --manifest <path> [--ssh <host>]`
   (without Node.js: `annolens install ...` on the data machine). Configuration files are backed up first.
3. Check the connection: `npx -y annolens-mcp doctor --manifest <path> [--ssh <host>]`.
4. Ask the user to restart the agent session so the new MCP server is loaded.

Plugin installs (Claude Code, Codex) read the same settings from the environment variables `ANNOLENS_MANIFEST`,
`ANNOLENS_SSH` and `ANNOLENS_REMOTE_CMD`.

On a new data machine:
`pip install git+https://github.com/forageseed/AnnoLens.git`, write `annolens.toml` (see `examples/annolens.toml` in
the repository), then `annolens prepare --manifest annolens.toml` and `annolens scan --manifest annolens.toml`.
Model API keys are read from environment variables only; never write them into files or chat.
