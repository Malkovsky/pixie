# Pixie Examples

Measure shared skills without including any Pixie-local overlays:

```bash
uv run agentic/cpp/skills/estimate-token-usage/scripts/estimate_tokens.py \
  skills --method tiktoken --skill benchmarks --skill optimization-experiment
```

When a task also reads Pixie's local examples, account for them explicitly as
additional context:

```bash
uv run agentic/cpp/skills/estimate-token-usage/scripts/estimate_tokens.py \
  context agentic/local/cpp/skills/benchmarks/EXAMPLES.md \
  --method tiktoken
```

Measure every Pixie-local skill overlay only when the analysis actually needs
that upper bound:

```bash
uv run agentic/cpp/skills/estimate-token-usage/scripts/estimate_tokens.py \
  context agentic/local/cpp/skills --include 'EXAMPLES.md' \
  --method tiktoken
```

For MCP accounting, first capture a server's `tools/list` JSON response, then
pass that file to `mcp` mode. Do not treat `.opencode/opencode.json` as an MCP
tool-description export; it configures servers but does not contain their tool
descriptions or input schemas.
