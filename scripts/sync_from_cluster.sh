#!/usr/bin/env bash
# Pull runs from the cluster into data/, bringing only the file types in sync_filter.txt (next to this script).
#
# Edit the settings, then:   bash scripts/sync_from_cluster.sh
# It first shows what would be copied (sizes, the total, and every file over BIG_MB), then asks before copying.
set -euo pipefail

# ============================================================================== settings
CLUSTER="bc889@rupc-01.rutgers.edu"                                   # ssh host or alias, e.g. "bc889@cluster.example.edu"; required
REMOTE="/home/bc889/mnt/w/axion/Y2Ir2O7/phonon/base/u_3.0/192344bc889/"   # job folder on the cluster; it is created, by name, inside LOCAL
LOCAL="/Users/treycole/Repos/axion-phonon/data/Y2Ir2O7/base/soc/u_3.0/output/"  # the folder to put it in: a full path, or relative to the repository root
BIG_MB=100                                   # list every file larger than this in the preview
# ======================================================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"
FILTER="$HERE/sync_filter.txt"
if [[ -z "$CLUSTER" ]]; then
    echo "Set CLUSTER at the top of $0 (the ssh host you rsync from)." >&2
    exit 1
fi
SOURCE="$CLUSTER:${REMOTE%/}"   # no trailing /: rsync then copies the folder itself (e.g. 192346/) into LOCAL
if [[ "$LOCAL" == /* ]]; then TARGET="$LOCAL"; else TARGET="$REPO/$LOCAL"; fi

# The rules are read here and passed one by one: with --filter="merge FILE", macOS's rsync (openrsync) sends
# the rule to the cluster, whose rsync then looks for FILE on the cluster.
OPTIONS=(-a --prune-empty-dirs)
while IFS= read -r rule || [[ -n "$rule" ]]; do
    [[ -z "${rule//[[:space:]]/}" || "$rule" == \#* ]] && continue
    OPTIONS+=(--filter="$rule")
done < "$FILTER"

echo "from: $SOURCE"
echo "to:   $TARGET"
echo "files that would be copied (bytes, path):"
if ! LISTING="$(rsync "${OPTIONS[@]}" --dry-run --out-format='%l %n' "$SOURCE" "$TARGET")"; then
    echo "rsync failed (message above); nothing was copied." >&2
    exit 1
fi
PREVIEW="$(printf '%s\n' "$LISTING" | grep -v '/$' | grep -v '^[[:space:]]*$' || true)"
if [[ -z "$PREVIEW" ]]; then
    echo "  nothing new"
    exit 0
fi
echo "$PREVIEW" | awk '{printf "  %10.1f MB  %s\n", $1 / 1e6, substr($0, index($0, $2))}'
echo "$PREVIEW" | awk -v big="$BIG_MB" '
    { total += $1; n += 1; if ($1 > big * 1e6) { large = large sprintf("    %8.1f MB  %s\n", $1 / 1e6, $2) } }
    END {
        printf "%d files, %.2f GB in total\n", n, total / 1e9
        if (large != "") { printf "files over %d MB:\n%s", big, large }
    }'

read -r -p "copy these? [y/N] " answer
if [[ "$answer" == [yY] ]]; then
    mkdir -p "$TARGET"
    rsync "${OPTIONS[@]}" --progress "$SOURCE" "$TARGET"   # --progress: macOS openrsync has no --info=progress2
fi
