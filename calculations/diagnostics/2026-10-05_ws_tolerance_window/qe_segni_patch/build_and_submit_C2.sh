#!/bin/bash
# Build pw.x with the PAW segni_rad lsign guard (paw_segni_lsign.patch) in a COPY of ~/q-e,
# then rerun NSCF test C (193145 fixed SCF save, 60 k) with it as test C2.
# Does not modify ~/q-e or any existing run folder. Stops on the first error.
set -euo pipefail
SRC=$HOME/q-e
DST=$HOME/q-e-segnifix
PATCH=$HOME/paw_segni_lsign.patch
T=$HOME/mnt/w/axion/Y2Ir2O7/phonon/base/u_3.0/test_nscf_sym_2026-10-05
EXPECTED_MD5=c41c409e8cb00733e864085801b5f6fc

echo "== source tree"
[ -e "$DST" ] && { echo "$DST already exists; not overwriting"; exit 1; }
[ -e "$T/C2" ] && { echo "$T/C2 already exists; not overwriting"; exit 1; }
md5=$(md5sum "$SRC/bin/pw.x" | cut -d' ' -f1)
[ "$md5" = "$EXPECTED_MD5" ] || { echo "~/q-e/bin/pw.x md5 $md5 != $EXPECTED_MD5"; exit 1; }
git -C "$SRC" log -1 --format='%h %ad %s' --date=short
du -sh "$SRC"
df -h "$HOME" | tail -1

echo "== copy"
cp -a "$SRC" "$DST"
cd "$DST"
echo "absolute paths to the source tree in make.inc before:"; grep -n "$SRC\b" make.inc || true
sed -i "s#$SRC\([/ ]\|\$\)#$DST\1#g" make.inc
echo "after:"; grep -n "$DST\|$SRC\b" make.inc || true

echo "== patch"
patch -p1 --dry-run < "$PATCH"
patch -p1 < "$PATCH"
git checkout -q -b fix/paw-segni-lsign
git -c user.name="Trey Cole" -c user.email="treco713@gmail.com" commit -q -am \
  "PAW noncollinear GGA: use ux for segni_rad only if lsign, as compute_rho does"
git log -1 --format='%h %s'

echo "== build"
module load intel/2024 intel/ompi
nice make -j4 pw > build_segnifix.log 2>&1 || { tail -40 build_segnifix.log; exit 1; }
tail -3 build_segnifix.log
ls -la bin/pw.x PW/src/pw.x
NEW_MD5=$(md5sum PW/src/pw.x | cut -d' ' -f1)
echo "new pw.x md5 $NEW_MD5"
[ "$NEW_MD5" != "$EXPECTED_MD5" ] || { echo "binary unchanged?!"; exit 1; }

echo "== test C2"
mkdir "$T/C2"
cp -p "$T/C/tmp.tar.gz" "$T"/C/*.UPF "$T/C/Y2Ir2O7.nscf.in" "$T/C2/"
cp -p "$DST/PW/src/pw.x" "$T/C2/pw.x"
sed "s/#\$ -N YIOnscfC/#\$ -N YIOnscfC2/" "$T/C/send_job.sh" > "$T/C2/send_job.sh"
{ echo "NSCF symmetry test C2: as C (193145 fixed SCF save, 60 k, nbnd 180) with pw.x from $DST"
  echo "(~/q-e 934f8cd + paw_segni_lsign.patch: PAW segni_rad uses ux only if lsign)."
  echo "SCF save: $T/C/tmp.tar.gz (= test_newqe_48/193145bc889/tmp)"
  md5sum "$T/C2/pw.x"; } > "$T/C2/NOTE"
cat "$T/C2/NOTE"
diff "$T/C/Y2Ir2O7.nscf.in" "$T/C2/Y2Ir2O7.nscf.in" && echo "nscf.in identical to C"
cd "$T/C2" && qsub send_job.sh
