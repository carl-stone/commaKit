# commaBot — Agent Bootstrap

This document is the setup guide for initializing or re-initializing commaBot's memory and identity. It's for whoever (Engels or Carl) needs to bootstrap commaBot from scratch or fix a corrupted memory state.

## Who commaBot is

commaBot is a dedicated Letta agent that processes GitHub webhook events for the commaKit R package. It acts as the engineering team — making decisions, running code, filing issues, creating PRs — not just acknowledging notifications.

## commaBot's environment

commaBot runs inside a Docker container on Carl's home PC (WSL2). Key facts:

| Item | Value |
|------|-------|
| Container name | `commabot` |
| R version | 4.5.2 (pinned — S4Vectors C API break on 4.6.0) |
| R packages | 192 via renv.lock + devtools/roxygen2/testthat/pkgdown |
| Python | 3.12 + Flask/PyJWT/cryptography |
| Node.js | 22 + @letta-ai/letta-code-sdk 0.25.9 |
| gh CLI | 2.92.0 (authenticate before use) |
| letta CLI | 0.25.9 (at /usr/local/bin/letta) |
| commaKit repo | `/home/commabot/commaKit/` (mounted from host, read-write) |
| memfs volume | `/home/commabot/.letta/` (Docker named volume, persists across restarts) |
| App scripts | `/app/` (webhook-listener.py, recovery.py, commabot-github-token.py, entrypoint.sh) |
| SDK relay | `/app/sdk/` (index.mjs, node_modules) |

## Container constraints

- Runtime installs (R packages, Python packages, apt packages) are **ephemeral** — lost on `docker compose up --build`. If commaBot needs a package baked in, add it to `dev/DOCKER_PACKAGES.md` and Carl will rebuild.
- commaBot **cannot** edit the Dockerfile, compose.yaml, or rebuild the image. Those live on the host at `/home/carls/commabot-infrastructure/`.
- The commaKit repo mount is the only persistent writable path (besides the memfs volume).
- No access to the host's `~/.letta/` or other host directories.

## GitHub authentication

commaBot authenticates as the `commabot[bot]` GitHub App:

```bash
GH_TOKEN=$(python3 /app/commabot-github-token.py --token) && gh auth login --with-token <<< "$GH_TOKEN"
```

Token expires after 1 hour. Re-run to refresh. The App ID is 3718770, installation ID is 132464445.

## Memory initialization

commaBot's memory blocks should contain:

1. **persona** — commaBot's identity: who it is, what it does, its relationship to Carl and the commaKit package
2. **human** — Carl's identity and preferences (can be a subset of Engels' knowledge)
3. **comma-context** — commaKit package architecture, conventions, and gotchas (from `dev/AGENTS.md`)
4. **environment** — the container environment, tools, constraints, and paths (from this doc)
5. **webhook-playbook** — per-event decision framework (the `handling-webhooks` skill content)

## Memory contamination warning

commaBot's memory was contaminated during the LETTA_MEMORY_SCOPE incident (2026-05-14). Engels' memory blocks leaked into commaBot's agent. If you see Engels' persona, Carl's full identity, or project references (seascape, epi-csc, megan-sc) in commaBot's memory, those are contamination — remove them.

## Host-side maintenance (for Engels/Carl)

The container is managed from the host. Key commands:

```bash
# Rebuild and restart (after Dockerfile changes)
cd /home/carls/commabot-infrastructure && docker compose up -d --build

# Restart without rebuild (after env file changes)
cd /home/carls/commabot-infrastructure && docker compose up -d

# View logs
docker logs -f commabot

# Shell into container
docker exec -it commabot bash

# Check health
curl http://localhost:8080/health

# Inspect memfs volume
docker run --rm -v commabot-infrastructure_commabot-letta:/data alpine ls -la /data
```

## Syncing renv.lock

When commaBot updates `renv.lock` in the mounted commaKit repo, the build context copy must be synced before the next image rebuild:

```bash
cp /home/carls/commaKit/renv.lock /home/carls/commabot-infrastructure/commaKit/
cp /home/carls/commaKit/renv/settings.json /home/carls/commabot-infrastructure/commaKit/renv/
cp /home/carls/commaKit/renv/activate.R /home/carls/commabot-infrastructure/commaKit/renv/
```

## Files reference

| Host path | Container path | Purpose |
|-----------|---------------|---------|
| `~/commabot-infrastructure/Dockerfile` | — | Image definition |
| `~/commabot-infrastructure/compose.yaml` | — | Container orchestration |
| `~/commabot-infrastructure/commabot.env` | — | Secrets and config (chmod 600) |
| `~/commabot-infrastructure/entrypoint.sh` | `/app/entrypoint.sh` | Container startup |
| `~/commabot-infrastructure/webhook-listener.py` | `/app/webhook-listener.py` | GitHub webhook receiver |
| `~/commabot-infrastructure/webhook-sdk/index.mjs` | `/app/sdk/index.mjs` | SDK relay |
| `~/commabot-infrastructure/recovery.py` | `/app/recovery.py` | Missed delivery recovery |
| `~/commabot-infrastructure/commabot-github-token.py` | `/app/commabot-github-token.py` | GitHub App token generator |
| `~/commabot-infrastructure/commabot.private-key.pem` | `/app/commabot.private-key.pem` | GitHub App private key |
| `~/webhook-conversation-prompt.md` | — | commaBot's webhook conversation prompt |
| `~/commaKit/` | `/home/commabot/commaKit/` | commaKit R package (mounted) |
| `~/commaKit/dev/DOCKER_PACKAGES.md` | `/home/commabot/commaKit/dev/DOCKER_PACKAGES.md` | Package request manifest |
