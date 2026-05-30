#!/bin/sh
set -eu

# Run R 4.5 from a temporary framework home whose wrappers do not point at
# /Library/Frameworks/R.framework/Current. This keeps Bioconductor 3.22 source
# builds from accidentally compiling or linking against R 4.6 on local macOS.

real_r_home=${COMMA_R45_HOME:-}
if [ -z "$real_r_home" ]; then
  for candidate in \
    /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources \
    /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources \
    /Library/Frameworks/R.framework/Versions/4.5/Resources
  do
    if [ -x "$candidate/bin/exec/R" ]; then
      real_r_home=$candidate
      break
    fi
  done
fi

if [ -z "$real_r_home" ] || [ ! -x "$real_r_home/bin/exec/R" ]; then
  cat >&2 <<'EOF'
Could not find a local R 4.5 framework.
Install R 4.5 or set COMMA_R45_HOME to its Resources directory.
EOF
  exit 1
fi

tmp_root=${TMPDIR:-/tmp}
if [ -n "${COMMA_R45_OVERLAY_HOME:-}" ]; then
  overlay_r_home=$COMMA_R45_OVERLAY_HOME
  cleanup_overlay=false
  case $(basename "$overlay_r_home") in
    comma-r-4.5-home.*) ;;
    *)
      cat >&2 <<'EOF'
COMMA_R45_OVERLAY_HOME must point to a scratch directory named comma-r-4.5-home.*.
Refusing to populate or remove an arbitrary path.
EOF
      exit 1
      ;;
  esac
  if [ -e "$overlay_r_home" ] && [ -n "$(find "$overlay_r_home" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    cat >&2 <<'EOF'
COMMA_R45_OVERLAY_HOME already exists and is not empty.
Choose an empty scratch directory or unset COMMA_R45_OVERLAY_HOME.
EOF
    exit 1
  fi
else
  overlay_r_home=$(mktemp -d "$tmp_root/comma-r-4.5-home.XXXXXX")
  cleanup_overlay=true
fi

cleanup() {
  if [ "$cleanup_overlay" = true ]; then
    rm -rf "$overlay_r_home"
  fi
}
trap cleanup EXIT HUP INT TERM

mkdir -p "$overlay_r_home/bin" "$overlay_r_home/etc"

for item in "$real_r_home"/*; do
  name=$(basename "$item")
  case "$name" in
    bin|etc) ;;
    *) ln -s "$item" "$overlay_r_home/$name" ;;
  esac
done

for item in "$real_r_home/bin"/*; do
  name=$(basename "$item")
  ln -s "$item" "$overlay_r_home/bin/$name"
done

rm "$overlay_r_home/bin/R"
sed "s#/Library/Frameworks/R.framework/Resources#$overlay_r_home#g" \
  "$real_r_home/bin/R" > "$overlay_r_home/bin/R"
chmod +x "$overlay_r_home/bin/R"

for item in "$real_r_home/etc"/*; do
  name=$(basename "$item")
  ln -s "$item" "$overlay_r_home/etc/$name"
done

rm "$overlay_r_home/etc/Makeconf"
sed \
  -e "s#/Library/Frameworks/R.framework/Resources#$overlay_r_home#g" \
  -e 's#^LIBR = .*#LIBR = -L"$(R_HOME)/lib$(R_ARCH)" -lR#' \
  "$real_r_home/etc/Makeconf" > "$overlay_r_home/etc/Makeconf"

export R_HOME=$overlay_r_home
export R_SHARE_DIR=$overlay_r_home/share
export R_INCLUDE_DIR=$overlay_r_home/include
export R_DOC_DIR=$overlay_r_home/doc

set +e
if [ "${1:-}" = "CMD" ]; then
  "$overlay_r_home/bin/R" "$@"
else
  "$overlay_r_home/bin/exec/R" --no-echo --no-restore "$@"
fi
status=$?
set -e

exit "$status"
