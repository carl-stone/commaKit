#!/usr/bin/env bash
set -euo pipefail

export HOME="${HOME:-/tmp/commakit-home}"
export TMPDIR="${TMPDIR:-/tmp}"

mkdir -p "$HOME" "$TMPDIR"

# The image is keyed by renv.lock, not DESCRIPTION. Check the current checkout
# at runtime so adding a package dependency fails clearly until the locked
# development image is updated. An interactive diagnostic shell remains usable
# even if dependencies are missing; scripted shell commands still get checked.
if [[ "${COMMAKIT_CONTAINER:-}" == "1" ]] &&
  [[ "${1:-}" != "bash" || $# -ne 1 ]]; then
  Rscript --vanilla /workspace/commaKit/dev/container/check-dependencies.R
fi

exec "$@"
