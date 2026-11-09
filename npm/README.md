# annolens-mcp

npm launcher for the [AnnoLens](https://github.com/forageseed/AnnoLens) MCP server. AnnoLens reviews and repairs gene
models from evidence (transcripts, protein projections, read checks) the way a human curator does. Everything runs on
your own machine; the engine is the Python package `annolens`.

```bash
python -m pip install git+https://github.com/forageseed/AnnoLens.git   # engine (Linux, macOS or WSL)
npx -y annolens-mcp install --host codex|claude|cursor|generic --manifest ./annolens.toml
npx -y annolens-mcp check
```

`npx -y annolens-mcp` alone runs the MCP server (the same as `annolens mcp`). The project manifest comes from
`--manifest`, the environment variable `ANNOLENS_MANIFEST`, or `./annolens.toml`.

License: Apache-2.0.
