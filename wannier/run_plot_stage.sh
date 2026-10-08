#!/usr/bin/env bash
# Stage 1, final step: write H(R) and A(R) with Jae-Mo's Wannier90 fork.
#
#   wannier/run_plot_stage.sh RUN_DIR [SEED] [WANNIER90_X]
#
#   RUN_DIR      a finished Wannier90 run holding SEED.win, SEED.chk and SEED.eig
#   SEED         seedname; default = the only *.win in RUN_DIR
#   WANNIER90_X  default ~/Repos/wannier90-JaeMo/build/wannier90.x
#                (branch plan5-write-ndegen-applied)
#
# Everything is written into RUN_DIR/jaemo/, so the production files are never
# touched:  SEED_hr.dat  SEED_r.dat  SEED_wsvec.dat  SEED_centres.xyz  SEED.wout
# `restart = plot` never rewrites the .chk, so the symlinked one is safe.
# Keyword names are matched in lower case, as Wannier90 files are normally written.
set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

RUN_DIR=$(cd "$1" && pwd)
W90X=${3:-$HOME/Repos/wannier90-JaeMo/build/wannier90.x}

if [[ -n "${2:-}" ]]; then
    SEED=$2
else
    shopt -s nullglob
    wins=("$RUN_DIR"/*.win)
    [[ ${#wins[@]} -eq 1 ]] || { echo "need exactly one *.win in $RUN_DIR (or pass SEED)"; exit 1; }
    SEED=$(basename "${wins[0]}" .win)
fi

for ext in win chk eig; do
    [[ -f "$RUN_DIR/$SEED.$ext" ]] || { echo "missing $RUN_DIR/$SEED.$ext"; exit 1; }
done
[[ -x "$W90X" ]] || { echo "wannier90.x not found or not executable: $W90X"; exit 1; }

OUT="$RUN_DIR/jaemo"
mkdir -p "$OUT"
ln -sfn "../$SEED.chk" "$OUT/$SEED.chk"
ln -sfn "../$SEED.eig" "$OUT/$SEED.eig"
WIN="$OUT/$SEED.win"
cp "$RUN_DIR/$SEED.win" "$WIN"

# set KEY = VALUE, replacing an existing line or appending one (never duplicating:
# Wannier90 rejects a keyword that appears twice)
setkey() {
    local key=$1 value=$2
    if grep -qE "^[[:space:]]*${key}[[:space:]]*[=:]" "$WIN"; then
        sed -i.bak -E "s|^[[:space:]]*${key}[[:space:]]*[=:].*|${key} = ${value}|" "$WIN"
        rm -f "$WIN.bak"
    else
        printf '%s = %s\n' "$key" "$value" >> "$WIN"
    fi
}

setkey restart               plot
setkey use_ws_distance       true
setkey transl_inv_full       true
setkey write_ndegen_applied  true
setkey write_hr              true
setkey write_rmn             true
setkey bands_plot            false
setkey write_tb              false

echo "running $W90X $SEED in $OUT"
( cd "$OUT" && "$W90X" "$SEED" > "$SEED.stdout" 2>&1 ) || {
    echo "wannier90.x failed; see $OUT/$SEED.wout and $SEED.werr"; exit 1; }
[[ ! -s "$OUT/$SEED.werr" ]] || { echo "wannier90 wrote an error file:"; cat "$OUT/$SEED.werr"; exit 1; }

head -1 "$OUT/${SEED}_wsvec.dat" | grep -q 'write_ndegen_applied=.true.' \
    || { echo "WARNING: ${SEED}_wsvec.dat does not carry write_ndegen_applied=.true."; exit 1; }

ls -lh "$OUT/${SEED}_hr.dat" "$OUT/${SEED}_r.dat" "$OUT/${SEED}_wsvec.dat"
echo "done.  next:  python wannier/check_outputs.py $OUT"
