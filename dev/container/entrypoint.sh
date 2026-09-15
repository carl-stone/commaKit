#!/usr/bin/env bash
set -euo pipefail

export HOME="${HOME:-/tmp/commakit-home}"
export TMPDIR="${TMPDIR:-/tmp}"

mkdir -p "$HOME" "$TMPDIR"

exec "$@"
