#!/usr/bin/env bash
set -euo pipefail

# User defaults for cluster runs. Fill these in if you want to run this script
# without the five positional path arguments. Command-line arguments still win.
PREFIX_DEFAULT="Y2Ir2O7"
BASE="/home/bc889/mnt/w/Y2Ir2O7/phonon/base/u_3.0/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181427bc889"
M="/home/bc889/mnt/w/Y2Ir2O7/phonon/mode_1/178365bc889/trial_01_Y_d_Ir_spd_O_sp/dis_nomaxloc"
BASE_SAVE="/home/bc889/mnt/w/Y2Ir2O7/phonon/base/u_3.0/178382bc889/tmp/Y2Ir2O7.save"
M_SAVE="/home/bc889/mnt/w/Y2Ir2O7/phonon/mode_1/178365bc889/tmp/Y2Ir2O7.save"

# Optional AMN output directories. Leave blank for the default directories under
# M. OUTPUT_DIR applies to either mode; the mode-specific directories apply only
# when OUTPUT_DIR is blank. Command-line --output still wins over all of these.
OUTPUT_DIR="${OUTPUT_DIR:-}"
RIGID_SHIFT_DIR="${RIGID_SHIFT_DIR:-}"
NO_DISPLACEMENT_DIR="/home/bc889/mnt/w/Y2Ir2O7/phonon/mode_1/178365bc889/trial_02_Y_p_Ir_d_O_sp/no_displacement"

usage() {
    cat <<'EOF'
Usage:
  generate_rigid_shift_amn.sh [options] [PREFIX BASE_W90_DIR PERT_W90_DIR BASE_SAVE_DIR PERT_SAVE_DIR]

Generate a Wannier90 .amn file for a perturbed/displaced calculation using
QEWavefunctions.jl's rigid-shift AMN generator.

Positional arguments:
  PREFIX         Wannier90 seed name, e.g. Y2Ir2O7 or MnBi2Te4
  BASE_W90_DIR   Reference W90 directory containing PREFIX.chk and PREFIX.nnkp
  PERT_W90_DIR   Perturbed W90 directory containing PREFIX.nnkp
  BASE_SAVE_DIR  Reference QE save directory containing wfc*.hdf5 or wfc*.dat
  PERT_SAVE_DIR  Perturbed QE save directory containing wfc*.hdf5 or wfc*.dat

If no positional arguments are given, the script uses the editable defaults near
the top of this file:

  PREFIX_DEFAULT  Wannier90 seed name, e.g. Y2Ir2O7
  BASE            Reference W90 directory
  M               Perturbed W90 directory
  BASE_SAVE       Reference QE save directory
  M_SAVE          Perturbed QE save directory
  OUTPUT_DIR      Directory receiving PREFIX.amn for either shifted mode
  RIGID_SHIFT_DIR Directory receiving PREFIX.amn for regular rigid shift
  NO_DISPLACEMENT_DIR
                  Directory receiving PREFIX.amn for --no-displacement

If BASE_SAVE or M_SAVE is left empty, it defaults to BASE/tmp/PREFIX.save or
M/tmp/PREFIX.save.

Options:
  --output PATH          Output .amn path.
                         Default: PERT_W90_DIR/rigid_shift/PREFIX.amn
                         With --no-displacement: PERT_W90_DIR/no_displacement/PREFIX.amn
  --output-dir DIR       Directory for output PREFIX.amn in either mode.
  --rigid-shift-dir DIR  Directory for output PREFIX.amn in regular rigid-shift mode.
  --no-displacement-dir DIR
                         Directory for output PREFIX.amn with --no-displacement.
  --qewf-dir DIR         QEWavefunctions.jl checkout.
                         Default: $QEWAVEFUNCTIONS_DIR or $HOME/mnt/w/QEWavefunctions.jl
  --generator PATH       Explicit generate_amn.jl path.
                         Default: QEWAVEFUNCTIONS_DIR/scripts/generate_amn.jl
  --julia PATH           Julia executable. Default: $JULIA or julia
  --no-displacement      Pass through to generate_amn.jl; overlap unshifted
                         reference WFs with perturbed bands.
  --no-unzip             Do not auto-run gunzip -k for wfc*.dat.gz files.
  --force                Allow overwriting an existing output .amn.
  --dry-run              Print the Julia command without running it.
  -h, --help             Show this help.

Examples:
  scripts/generate_rigid_shift_amn.sh Y2Ir2O7 \
    /home/bc889/mnt/w/Y2Ir2O7/phonon/base/u_3.0/174895bc889 \
    /home/bc889/mnt/w/Y2Ir2O7/phonon/mode_1/176111bc889 \
    /home/bc889/mnt/w/Y2Ir2O7/phonon/base/u_3.0/174895bc889/tmp/Y2Ir2O7.save \
    /home/bc889/mnt/w/Y2Ir2O7/phonon/mode_1/176111bc889/tmp/Y2Ir2O7.save

  scripts/generate_rigid_shift_amn.sh --no-displacement Y2Ir2O7 BASE MODE BASE_SAVE MODE_SAVE

  # After filling the defaults at the top of this file:
  scripts/generate_rigid_shift_amn.sh
  scripts/generate_rigid_shift_amn.sh --no-displacement
EOF
}

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

info() {
    printf '# %s\n' "$*"
}

quote_cmd() {
    printf '%q ' "$@"
    printf '\n'
}

require_file() {
    local path=$1
    [[ -f "$path" ]] || die "missing file: $path"
}

require_dir() {
    local path=$1
    [[ -d "$path" ]] || die "missing directory: $path"
}

has_wfc_files() {
    local dir=$1
    compgen -G "$dir/wfc*.hdf5" >/dev/null || compgen -G "$dir/wfc*.dat" >/dev/null
}

has_gzipped_wfc_files() {
    local dir=$1
    compgen -G "$dir/wfc*.dat.gz" >/dev/null
}

prepare_wfc_dir() {
    local label=$1
    local dir=$2

    require_dir "$dir"
    if has_wfc_files "$dir"; then
        info "$label QE save directory has readable wfc files: $dir"
        return
    fi

    if has_gzipped_wfc_files "$dir"; then
        if [[ "$UNZIP_WFC" == "1" ]]; then
            info "$label QE save directory has only gzipped wfc*.dat files; running gunzip -k"
            gunzip -k "$dir"/wfc*.dat.gz
            has_wfc_files "$dir" || die "gunzip completed but no wfc*.dat files were found in $dir"
            return
        fi
        die "$label QE save directory has wfc*.dat.gz but no uncompressed wfc files; rerun without --no-unzip"
    fi

    die "$label QE save directory has no wfc*.hdf5, wfc*.dat, or wfc*.dat.gz files: $dir"
}

QEWAVEFUNCTIONS_DIR=${QEWAVEFUNCTIONS_DIR:-"$HOME/mnt/w/QEWavefunctions.jl"}
GENERATOR=${GENERATOR:-}
JULIA_BIN=${JULIA:-julia}
OUTPUT=
NO_DISPLACEMENT=0
UNZIP_WFC=1
FORCE=0
DRY_RUN=0

args=()
while (($#)); do
    case "$1" in
        --output)
            (($# >= 2)) || die "--output requires a path"
            OUTPUT=$2
            shift 2
            ;;
        --output-dir|--save-dir)
            (($# >= 2)) || die "$1 requires a directory"
            OUTPUT_DIR=$2
            shift 2
            ;;
        --rigid-shift-dir)
            (($# >= 2)) || die "--rigid-shift-dir requires a directory"
            RIGID_SHIFT_DIR=$2
            shift 2
            ;;
        --no-displacement-dir)
            (($# >= 2)) || die "--no-displacement-dir requires a directory"
            NO_DISPLACEMENT_DIR=$2
            shift 2
            ;;
        --qewf-dir)
            (($# >= 2)) || die "--qewf-dir requires a directory"
            QEWAVEFUNCTIONS_DIR=$2
            shift 2
            ;;
        --generator)
            (($# >= 2)) || die "--generator requires a path"
            GENERATOR=$2
            shift 2
            ;;
        --julia)
            (($# >= 2)) || die "--julia requires an executable path"
            JULIA_BIN=$2
            shift 2
            ;;
        --no-displacement)
            NO_DISPLACEMENT=1
            shift
            ;;
        --no-unzip)
            UNZIP_WFC=0
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            args+=("$@")
            break
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

if ((${#args[@]} == 5)); then
    PREFIX=${args[0]}
    BASE_W90_DIR=${args[1]}
    PERT_W90_DIR=${args[2]}
    BASE_SAVE_DIR=${args[3]}
    PERT_SAVE_DIR=${args[4]}
elif ((${#args[@]} == 0)); then
    PREFIX=${PREFIX_DEFAULT:-${PREFIX:-}}
    BASE_W90_DIR=$BASE
    PERT_W90_DIR=$M
    BASE_SAVE_DIR=${BASE_SAVE:-}
    PERT_SAVE_DIR=${M_SAVE:-}
    [[ -n "$PREFIX" ]] || die "set PREFIX_DEFAULT at the top of the script, or pass PREFIX as an argument"
    [[ -n "$BASE_W90_DIR" ]] || die "set BASE at the top of the script, or pass BASE_W90_DIR as an argument"
    [[ -n "$PERT_W90_DIR" ]] || die "set M at the top of the script, or pass PERT_W90_DIR as an argument"
    [[ -n "$BASE_SAVE_DIR" ]] || BASE_SAVE_DIR="$BASE_W90_DIR/tmp/$PREFIX.save"
    [[ -n "$PERT_SAVE_DIR" ]] || PERT_SAVE_DIR="$PERT_W90_DIR/tmp/$PREFIX.save"
else
    usage >&2
    exit 1
fi

[[ -n "$GENERATOR" ]] || GENERATOR="$QEWAVEFUNCTIONS_DIR/scripts/generate_amn.jl"

require_dir "$QEWAVEFUNCTIONS_DIR"
require_file "$GENERATOR"
require_dir "$BASE_W90_DIR"
require_dir "$PERT_W90_DIR"
require_file "$BASE_W90_DIR/$PREFIX.chk"
require_file "$BASE_W90_DIR/$PREFIX.nnkp"
require_file "$PERT_W90_DIR/$PREFIX.nnkp"
prepare_wfc_dir "Reference" "$BASE_SAVE_DIR"
prepare_wfc_dir "Perturbed" "$PERT_SAVE_DIR"

if [[ -z "$OUTPUT" ]]; then
    if [[ -n "$OUTPUT_DIR" ]]; then
        OUTPUT="$OUTPUT_DIR/$PREFIX.amn"
    elif [[ "$NO_DISPLACEMENT" == "1" && -n "$NO_DISPLACEMENT_DIR" ]]; then
        OUTPUT="$NO_DISPLACEMENT_DIR/$PREFIX.amn"
    elif [[ "$NO_DISPLACEMENT" != "1" && -n "$RIGID_SHIFT_DIR" ]]; then
        OUTPUT="$RIGID_SHIFT_DIR/$PREFIX.amn"
    elif [[ "$NO_DISPLACEMENT" == "1" ]]; then
        OUTPUT="$PERT_W90_DIR/no_displacement/$PREFIX.amn"
    else
        OUTPUT="$PERT_W90_DIR/rigid_shift/$PREFIX.amn"
    fi
fi

mkdir -p "$(dirname "$OUTPUT")"
if [[ -e "$OUTPUT" && "$FORCE" != "1" ]]; then
    die "output already exists: $OUTPUT (use --force to overwrite)"
fi

cmd=("$JULIA_BIN" "--project=$QEWAVEFUNCTIONS_DIR" "$GENERATOR")
if [[ "$NO_DISPLACEMENT" == "1" ]]; then
    cmd+=("--no-displacement")
fi
cmd+=("$PREFIX" "$BASE_W90_DIR" "$PERT_W90_DIR" "$BASE_SAVE_DIR" "$PERT_SAVE_DIR" "$OUTPUT")

info "QEWavefunctions.jl checkout : $QEWAVEFUNCTIONS_DIR"
info "Generator                  : $GENERATOR"
info "Output                     : $OUTPUT"
info "No-displacement mode       : $NO_DISPLACEMENT"
info "Command:"
quote_cmd "${cmd[@]}"

if [[ "$DRY_RUN" == "1" ]]; then
    exit 0
fi

"${cmd[@]}"
