# Connecting AnnoLens

AnnoLens runs on this machine as the Python package `annolens`; the agent host talks to it through the MCP server
`annolens mcp`.

1. Install the engine if `annolens` is not on PATH: `python -m pip install git+https://github.com/forageseed/AnnoLens.git`
   (Python 3.11+, Linux, macOS or WSL).
2. Find the project manifest `annolens.toml`; ask the user for its path if it is not in the working directory. Do not
   guess paths. A new project starts from `examples/annolens.toml` in the repository, then `annolens prepare` and
   `annolens scan`.
3. Register the server and install this skill: `annolens install --host <claude|codex|cursor|generic> --manifest <path>`.
   Configuration files are backed up first.
4. Run `annolens check` and report anything marked FAIL.
5. Ask the user to restart the agent session so the new MCP server is loaded.

Model API keys are read from environment variables only; never write them into files or chat.
