# Shared by the step scripts (source'd, never run alone): where the binaries, pseudopotentials and code are, and the
# helpers that write job scripts, submit them and check markers.
#
# Markers in each step folder (written by the job, or by the step script for SUBMITTED):
#   SUBMITTED  the job id, written at submission         OK    the job id, written by the job when every check passed
#   NOTE       one line per submission (date, inputs)    FAIL  the job id, written by the job on any failed check
# The cluster epilog copies the job's scratch folder back to <step folder>/<job id>bc889/; the first step that reads
# it renames it out/ (need_ok). One step folder = one calculation; a variant is another step folder. c2 tasks use the same markers with the task number appended (OK_3, FAIL_3, ...).

AX=/mnt/wkdk/bc889/axion
QE=$AX/bin/qe_db4b99b4                         # pw.x db4b99b4, pw2wannier90.x 8f5368ba, ph.x 554fabc7 (dfpt.sh)
PH_MD5=554fabc7                                # ph.x, built 2026-11-01 in the same tree (make ph)
W90=$AX/bin/wannier90_26c8b506/wannier90.x     # 26c8b506
PSEUDO=${PSEUDO:-$AX/pseudo/rel-pz}            # every production UPF (MD5SUMS there); a step's settings may name
                                               # another set of $AX/pseudo/ (SrTiO3: oncv-pbe-rel)
CODE=$AX/code/2026-10-29                       # axion_phonon/, pipeline/*.py, pythtb/ (external 5776792); MD5SUMS
PY=/mnt/wkdk/bc889/envs/axion/bin/python       # numpy 1.26, scipy 1.11; runs on the login node and the nodes
USER_TAG=bc889                                 # job folder = <step folder>/<job id>$USER_TAG
EXCL='!n138'                                   # n138 has no /mnt/wkdk
export PYTHONPATH=$CODE                        # for the Python the step scripts run on the login node

module load intel/2024 intel/ompi 2> /dev/null || true   # wannier90.x -pp runs on the login node (MKL)

die() { echo "$*" >&2; exit 1; }

check_tools() {  # binaries, pseudopotentials and code are the production ones
  [ "$(md5sum < $QE/pw.x | cut -c1-8)" = db4b99b4 ] || die "pw.x md5"
  [ "$(md5sum < $QE/pw2wannier90.x | cut -c1-8)" = 8f5368ba ] || die "pw2wannier90.x md5"
  [ "$(md5sum < $W90 | cut -c1-8)" = 26c8b506 ] || die "wannier90.x md5"
  (cd $PSEUDO && md5sum -c --quiet MD5SUMS) || die "pseudopotential md5"
  (cd $CODE && md5sum -c --quiet MD5SUMS) || die "code md5 ($CODE/MD5SUMS)"
}

ready() {  # step folder [suffix]: no OK yet and no job in flight; then clear an old FAIL
  local d=$1 m=${2:-}
  mkdir -p $d
  [ -e $d/OK$m ] && die "$d: already OK$m (job $(cat $d/OK$m)); move the folder aside to redo it"
  if [ -e $d/SUBMITTED$m ] && [ ! -e $d/FAIL$m ]; then
    die "$d: job $(cat $d/SUBMITTED$m) has no OK$m or FAIL$m yet (qstat; if it died, look at its job folder, then delete SUBMITTED$m)"
  fi
  return 0
}

need_ok() {  # step folder [suffix] -> prints its output folder out<suffix>; fails unless OK and the job left the queue.
            # The first call renames the job's folder <job id>bc889 (named by the cluster) to out<suffix>.
  local d=$1 m=${2:-} j o=$1/out${2:-}
  [ -s $d/OK$m ] || die "not OK: $d/OK$m"
  j=$(cat $d/OK$m)
  qstat -j $j > /dev/null 2>&1 && die "$d: job $j is still in the queue (copy-back not finished)"
  if [ ! -e $o ]; then
    [ -d $d/$j$USER_TAG ] || die "$d: no out$m and no job folder $j$USER_TAG"
    mv $d/$j$USER_TAG $o
  fi
  echo $o
}

header() {  # name slots h_vmem h_rt queues [suffix] -> the top of a job script
  local m=${6:-}
  printf "%s\n" "#!/bin/bash" "#\$ -N $1" "#\$ -pe orte_one $2" "#\$ -q $5" "#\$ -j y" "#\$ -l h_vmem=$3" "#\$ -l h_rt=$4" \
    "#\$ -v WIEN_DMFT_ROOT,WIENROOT,LD_LIBRARY_PATH,PATH" \
    "source \$TMPDIR/sge_init.sh" "set -u" \
    "echo \"job \$JOB_ID on \$HOSTNAME in \$PWD\"; date" \
    "export MODULEPATH=/opt/apps/modulefiles:/opt/intel/modulefiles:/opt/pgi/modulefiles:/opt/gnu/modulefiles:/opt/sw/modulefiles" \
    "module load intel/2024 intel/ompi" \
    "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$CODE" \
    "fail() { echo \"FAIL: \$*\"; echo \$JOB_ID > \$SGE_O_WORKDIR/FAIL$m; exit 1; }" \
    "ok() { echo \$JOB_ID > \$SGE_O_WORKDIR/OK$m; echo OK; date; }" \
    "get() {  # folder file: copy into the current folder, uncompressed (the epilog gzips large files)" \
    "  if [ -s \$1/\$2.gz ]; then pigz -dc \$1/\$2.gz > \$2; elif [ -s \$1/\$2 ]; then cp -p \$1/\$2 .; fi" \
    "  [ -s \$2 ] || fail \"missing \$1/\$2\"" \
    "}" \
    "get_save() {  # save folder, destination tmp folder, number of k points: copy, gunzip, check" \
    "  mkdir -p \$2 && cp -r \$1 \$2/ || fail \"copy \$1\"" \
    "  find \$2 -name '*.gz' -print0 | xargs -0 -r -P 8 -n 16 gunzip -f" \
    "  local s=\$2/\$(basename \$1)" \
    "  [ -s \$s/data-file-schema.xml ] && [ -s \$s/charge-density.dat ] || fail \"incomplete save \$1\"" \
    "  [ \$3 = 0 ] || [ \$(ls \$s | grep -c '^wfc[0-9]*\\.dat\$') = \$3 ] || fail \"\$1: not \$3 wfc files\"" \
    "}" \
    "scratch() {  # GB needed" \
    "  local a=\$(df -BG --output=avail . | tail -1 | tr -dc 0-9); echo \"scratch free \$a GB\"; [ \$a -ge \$1 ] || fail \"scratch \$a GB < \$1\"" \
    "}"
}

fill() {  # source .py, target, NAME=value ...: copy, setting each settings line 'NAME = ...' (must occur once)
  local src=$1 dst=$2 a n v
  shift 2
  cp -p $src $dst
  for a in "$@"; do
    n=${a%%=*}; v=${a#*=}
    [ "$(grep -c "^$n = " $dst)" = 1 ] || die "$src: no single '$n = ' line"
    sed -i "s|^$n = .*|$n = $v|" $dst
  done
}

submit() {  # step folder, job script, [suffix], [extra qsub arguments]; SUBMIT=false only prints the qsub line
  local d=$1 s=$2 m=${3:-} j
  if [ "$SUBMIT" != true ]; then
    echo "dry run, not submitted: (cd $d && qsub ${4:-} -l hostname='$EXCL' $s)"
    return 0
  fi
  rm -f $d/FAIL$m
  j=$(cd $d && qsub ${4:-} -l hostname="$EXCL" $s | grep -o "job [0-9]*" | grep -o "[0-9]*")
  [ -n "$j" ] || die "$d: qsub failed"
  echo $j > $d/SUBMITTED$m
  echo "$d: job $j"
}

note() {  # step folder, text: one line in its NOTE (only when submitting)
  [ "$SUBMIT" = true ] && echo "$(date +%F\ %H:%M) $2" >> $1/NOTE
  return 0
}
