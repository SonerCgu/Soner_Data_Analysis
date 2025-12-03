#!/bin/bash

###############################################################################
# BLuSH v24-fix1 — FULLY ANNOTATED VERSION (PART 1)
#
# This script has been expanded with detailed documentation intended
# to teach, clarify, and future-proof each step of the pipeline.
#
# COMMENT STRATEGY (Hybrid Mode):
#   • At the start of each major block, you will find a **conceptual
#     explanation** (why this step exists, what problem it solves).
#   • Inside the code, you will find **short inline comments** describing
#     the implementation logic.
#
# All explanations are comments only — the pipeline remains 100% runnable.
###############################################################################



###############################################################################
# SECTION 0 — SAFETY SETTINGS & GLOBAL VARIABLES
#
# WHY THIS SECTION EXISTS:
# BLuSH is a multi-stage pipeline with many calls to FSL, AFNI, Python,
# and file manipulations. Silent failures would cause corrupted outputs.
#
# Therefore we activate strict shell safeguards:
#   -e : abort if any command fails
#   -u : abort on undefined variable
#   -o pipefail : abort if ANY part of a pipeline fails
#   -E : ensures traps are inherited properly
#
# We also:
#   • Create consistent date/time stamps
#   • Set US-like floating-point behavior (important for Python/FSL)
#   • Suppress Python warnings to keep output clean
###############################################################################

set -Eeuo pipefail                   # strict mode: fail early, loudly
trap 'echo "[ERROR] Script failed at line ${LINENO}" >&2' ERR
shopt -s nullglob                   # globs that match nothing → empty, not literal
export LC_NUMERIC=C                 # ensure floats use dot notation
export PYTHONWARNINGS="ignore"      # suppress unhelpful nibabel warnings



###############################################################################
# SECTION 1 — TERMINAL COLORS & TIMESTAMP FUNCTION
#
# WHY THIS SECTION EXISTS:
# Purely cosmetic: BLuSH uses colored terminal messages to indicate stages
# and to visually separate major pipeline steps.
###############################################################################

RED="\033[38;5;160m"
PURPLE="\033[38;5;135m"
GREEN="\033[92m"
BLUE="\033[38;5;33m"
RESET="\033[0m"
BOLD="\033[1m"

timestamp(){ date +"%Y%m%d_%H%M%S"; }   # consistent timestamp helper



###############################################################################
# SECTION 2 — FOLDER STRUCTURE PREPARATION
#
# WHY THIS SECTION EXISTS:
# BLuSH writes dozens of intermediate and QC files across many stages.
# A consistent folder structure ensures:
#   • no overwriting previous runs
#   • easy navigation
#   • consistent outputs for publication & debugging
#
# DIRECTORY STRUCTURE CREATED:
#
# soner_pipeline/
#   ├── filtering/         (HP/LP/BP/Notch outputs)
#   ├── scm_outputs/       (PSC maps, baseline/signal maps, norm_func)
#   ├── detrended/         (polynomial drift-removed data)
#   ├── pca_outputs/       (PCA, ICA, PLS denoised data)
#   │      ├── raw/
#   │      ├── hp/
#   │      ├── lp/
#   │      ├── bp/
#   │      ├── notch/
#   │      ├── trim/
#   │      ├── smoothed/
#   │      └── detrended/
#   └── qc/
#         ├── global/      (global signal plots)
#         └── grid/        (component GUI grids)
###############################################################################

mkdir -p soner_pipeline/{filtering,scm_outputs,logs,detrended}
mkdir -p soner_pipeline/pca_outputs/{raw,hp,lp,bp,notch,detrended,trim,smoothed}
mkdir -p soner_pipeline/qc/{global,grid}



###############################################################################
# SECTION 3 — MASK VALIDATION
#
# WHY THIS SECTION EXISTS:
# The mask is optional, but if present:
#   • SCM operations use it for -mas
#   • Global signal measures focus on brain only
#
# If mask is missing, BLuSH automatically disables masking.
###############################################################################

MASK="mask_mean_mc_func.nii.gz"

# The cleaned input file must exist — everything depends on it
[[ ! -f cleaned_mc_func.nii.gz ]] && {
    echo "[ERR] cleaned_mc_func.nii.gz not found" >&2
    exit 1
}

# If mask missing → MASK="" so FSL math skipping -mas is safe
[[ ! -f "$MASK" ]] && MASK=""



###############################################################################
# SECTION 4 — HELPER FUNCTIONS
#
# WHY THIS SECTION EXISTS:
# BLuSH uses two central helpers:
#
# 1. center_line(text)
#    Prints a stylized ASCII header between pipeline steps.
#
# 2. get_tr(filename)
#    Extracts TR safely using nibabel.
#    Important because:
#       • 3D files do not contain TR
#       • corrupted headers are common
#       • filtering requires TR
###############################################################################

center_line(){
  printf "\n────────────────────────────────────────────────────────\n✳️  %s ✳️\n────────────────────────────────────────────────────────\n" "$1"
}

get_tr(){
python3 - "$1" <<'PY'
import nibabel as nib, sys
try:
    hdr = nib.load(sys.argv[1]).header
    zooms = list(hdr.get_zooms())
    # 3D files have no TR stored → default to 1.0
    if len(zooms) < 4:
        print("1.000000")
    else:
        print(f"{float(zooms[3]):.6f}")
except:
    # fallback for corrupted NIfTI
    print("1.000000")
PY
}



###############################################################################
# SECTION 5 — GLOBAL SIGNAL PNG GENERATOR
#
# WHY THIS SECTION EXISTS:
# Global signal QC is the backbone of BLuSH:
#   • After each preprocessing step
#   • Before & after decomposition
#   • Before SCM window selection
#
# This function handles BOTH 3D & 4D cases:
#   • 3D → produce a single global value (warn user)
#   • 4D → apply mask (if available) or flatten
#
# The output is a timestamped PNG stored in:
#       soner_pipeline/qc/global/
###############################################################################

save_global_png(){
  local nii="$1"
  local mask="$2"
  local tr="$3"
  local label="$4"

  python3 - "$nii" "$mask" "$tr" "$label" <<'PY'
import sys, os, nibabel as nib, numpy as np, matplotlib.pyplot as plt, time

f,m,tr,label=sys.argv[1:5]
tr=float(tr)

img=nib.load(f)
d=img.get_fdata()

# Handle 3D case safely
if d.ndim == 3:
    print("[WARN] 3D data → global signal is single mean value.")
    sig = np.array([d.mean()])
    x = np.array([0])

else:
    try:
        mask = nib.load(m).get_fdata()>0
        d = d[mask]               # restrict to brain voxels
    except:
        d = d.reshape(-1, d.shape[-1])

    sig = d.mean(0)
    x = np.arange(len(sig)) * tr

# Plot
plt.figure(figsize=(10,4))
plt.plot(x, sig, 'k', lw=1)
plt.title(f"Global Mean – {os.path.basename(f)}")
plt.xlabel("Time (s)")
plt.ylabel("Mean Intensity")
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save to global QC folder
out=f"soner_pipeline/qc/global/global_{label}_{time.strftime('%Y%m%d_%H%M%S')}.png"
plt.savefig(out, dpi=160)
plt.close()

print("[OK] Global PNG:", out)
PY
}
###############################################################################
# SECTION 6 — INTERACTIVE GLOBAL SIGNAL VIEW
#
# WHY THIS SECTION EXISTS:
# Before detrending, filtering, or SCM, the user must visually inspect
# the global signal to:
#   • detect scanner drifts
#   • detect motion spikes
#   • choose trimming or baseline/signal windows
#
# This viewer is interactive (matplotlib) and does NOT save anything.
###############################################################################

plot_global_signal(){
python3 - "$1" "$2" "$3" <<'PY'
import sys, nibabel as nib, numpy as np, matplotlib.pyplot as plt

f, m, tr = sys.argv[1], sys.argv[2], float(sys.argv[3])
img = nib.load(f)
d = img.get_fdata()

# 3D case
if d.ndim == 3:
    print("[WARN] 3D → only single global value exists.")
    sig = np.array([d.mean()])
    x = np.array([0])

# 4D case
else:
    try:
        mask = nib.load(m).get_fdata()>0
        d = d[mask]
    except:
        d = d.reshape(-1, d.shape[-1])

    sig = d.mean(0)
    x = np.arange(len(sig))*tr

plt.figure(figsize=(10,4))
plt.plot(x, sig, 'k')
plt.title(f"Global Mean Signal – {f}")
plt.xlabel("Time (s)")
plt.ylabel("Intensity")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
PY
}



###############################################################################
# SECTION 7 — STEP 1: TRIMMING
#
# WHY THIS STEP EXISTS:
# • The first few volumes often have unstable magnetization (T1 saturation)
# • The last few volumes may show motion due to early stop
#
# AFNI INDEXING DETAILS:
#   func[start..end]
#
# 3D CASE:
#   • 3D volumes cannot be trimmed → warn + skip
###############################################################################

TR=$(get_tr cleaned_mc_func.nii.gz)
TMP="cleaned_mc_func.nii.gz"

echo -e "${PURPLE}${BOLD}STEP 1: TRIMMING${RESET}"
center_line "Trim initial/end volumes for stable baseline"

read -p "Trim from START (vols): " TSTART
read -p "Trim from END   (vols): " TEND

# Detect 3D vs 4D
IS3D=$(python3 - <<EOF
import nibabel as nib
print(1 if nib.load("cleaned_mc_func.nii.gz").get_fdata().ndim==3 else 0)
EOF
)

if [[ "$IS3D" == "1" ]]; then
    echo "[WARN] 3D file → trimming skipped."

else
    # Only trim if values are non-zero
    if [[ "${TSTART:-0}" != "0" || "${TEND:-0}" != "0" ]]; then

        # Total volumes
        NV=$(python3 - <<EOF
import nibabel as nib
print(nib.load("cleaned_mc_func.nii.gz").shape[-1])
EOF
)

        # Compute final index for AFNI syntax
        END=$(( NV - TEND - 1 ))
        OUT="trimmed_${TSTART}-${TEND}_cleaned_mc_func_$(timestamp).nii.gz"

        # Trim with AFNI
        3dTcat -prefix "$OUT" cleaned_mc_func.nii.gz"[${TSTART}..${END}]"

        TMP="$OUT"
    fi
fi

# QC
save_global_png "$TMP" "$MASK" "$TR" "trim"



###############################################################################
# SECTION 8 — STEP 2: SMOOTHING (TEMPORAL + SPATIAL)
#
# WHY SMOOTH?
# • Temporal smoothing reduces frame-to-frame noise (NOT a replacement for
#   filtering — it's a basic denoising step).
#
# • Spatial smoothing increases SNR by reducing voxel noise.
#
# TEMPORAL SMOOTHING:
#   uniform_filter1d() does a moving-average filter of width:
#       window (in vol) = round(seconds / TR)
#
# SPATIAL SMOOTHING:
#   fslmaths -s uses σ (NOT FWHM)
###############################################################################

echo -e "${PURPLE}${BOLD}STEP 2: SMOOTHING${RESET}"
center_line "Temporal smoothing + Spatial smoothing"

AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz 2>/dev/null))

# list the available files
for i in "${!AVAILABLE[@]}"; do
    printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
done

read -p "Select file to smooth [1-${#AVAILABLE[@]}]: " idx
TMP="${AVAILABLE[$((idx-1))]}"

read -p "Temporal smoothing window (seconds, 0=skip): " TSM
read -p "Spatial σ blur (mm, 0=skip): " SSM


##############################
# TEMPORAL SMOOTHING (Python)
##############################
if [[ "${TSM:-0}" != "0" ]]; then
python3 - "$TMP" "$TSM" "$TR" <<'PY'
import sys, nibabel as nib, numpy as np
from scipy.ndimage import uniform_filter1d

f=sys.argv[1]
win_sec=float(sys.argv[2])
tr=float(sys.argv[3])

img=nib.load(f)
d=img.get_fdata()

# Skip for 3D
if d.ndim==3:
    print("[WARN] 3D file → temporal smoothing skipped.")
    sys.exit(0)

# compute window length in timepoints
win=max(1,int(round(win_sec/tr)))

# apply moving average
sm=uniform_filter1d(d,size=win,axis=-1,mode="nearest")

out=f"temporal_smoothed_{win_sec}s_{f.rsplit('/',1)[-1]}"
nib.save(nib.Nifti1Image(sm.astype(np.float32),img.affine,img.header),out)
print("[OK] Temporal smoothing →",out)
PY

    if [[ -f temporal_smoothed_${TSM}s_${TMP##*/} ]]; then
        TMP="temporal_smoothed_${TSM}s_${TMP##*/}"
        save_global_png "$TMP" "$MASK" "$TR" "tsmooth"
    fi
fi


##############################
# SPATIAL SMOOTHING (FSL)
##############################
if [[ "${SSM:-0}" != "0" ]]; then
    OUT="spatial_smoothed_${SSM}mm_${TMP%.nii.gz}_$(timestamp).nii.gz"

    # sigma-based smoothing
    fslmaths "$TMP" -s "$SSM" "$OUT"

    TMP="$OUT"
    echo "[OK] Spatial smoothing → $OUT"
    save_global_png "$TMP" "$MASK" "$TR" "ssmooth"
fi



###############################################################################
# SECTION 9 — STEP 2.5: DETRENDING
#
# WHY DETRENDING?
# Polynomial drifts appear due to:
#   • scanner heating
#   • slow head motion
#   • physiological noise
#
# AFNI 3dDetrend removes polynomial trends:
#   polort 1 = linear drift
#   polort 2 = quadratic drift
#   polort N = N-th order polynomial drift
#
# This step also includes BEFORE/AFTER global signal inspection.
###############################################################################

echo -e "${PURPLE}${BOLD}STEP 2.5: DETRENDING${RESET}"
center_line "Remove scanner drift (polynomial detrending)"

read -p "Do you want to detrend? (y/n): " DO_DET

if [[ "$DO_DET" == "y" ]]; then

    AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz temporal_smoothed_*.nii.gz spatial_smoothed_*.nii.gz 2>/dev/null))

    # list options
    for i in "${!AVAILABLE[@]}"; do
        printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
    done

    read -p "Select file for detrending [1-${#AVAILABLE[@]}]: " idx
    TMP="${AVAILABLE[$((idx-1))]}"

    # show global signal first
    plot_global_signal "$TMP" "$MASK" "$TR"

    read -p "Apply detrending? (y/n): " DETR
    if [[ "$DETR" == "y" ]]; then
        while true; do
            echo "1) Linear"
            echo "2) Quadratic"
            echo "3) Custom polynomial"
            echo "4) Cancel"
            read -p "Choice [1–4]: " CHD

            case "$CHD" in
                1) POL=1 ;;
                2) POL=2 ;;
                3) read -p "Polynomial order: " POL ;;
                4) break ;;
                *) continue ;;
            esac

            BASE=$(basename "$TMP" .nii.gz)
            OUT="soner_pipeline/detrended/detrended_${BASE}_pol${POL}_$(timestamp).nii.gz"

            # AFNI detrending
            3dDetrend -polort $POL -prefix "$OUT" "$TMP"

            # QC
            save_global_png "$OUT" "$MASK" "$TR" "detrended_pol${POL}"
            plot_global_signal "$OUT" "$MASK" "$TR"

            read -p "Use this detrended file? (y/n): " OK
            [[ "$OK" == "y" ]] && TMP="$OUT" && break
        done
    fi
fi

###############################################################################
# FUNCTION: scm_blush_fast
#
# Computes baseline, signal, PSC map, and norm_func quickly using Python.
#
# INPUTS:
#   1 = NIfTI file
#   2 = baseline_start
#   3 = baseline_end
#   4 = signal_start
#   5 = signal_end
#   6 = mask ("" allowed)
#   7 = output folder
###############################################################################
scm_blush_fast(){
  local nii="$1"
  local b1="$2"
  local b2="$3"
  local s1="$4"
  local s2="$5"
  local mask="$6"
  local outdir="$7"

python3 - "$nii" "$b1" "$b2" "$s1" "$s2" "$mask" "$outdir" <<'PY'
import sys, os, nibabel as nib, numpy as np

nii, b1, b2, s1, s2, mask, outdir = sys.argv[1:]
b1=int(b1); b2=int(b2); s1=int(s1); s2=int(s2)

img=nib.load(nii)
d=img.get_fdata()
aff=img.affine
hdr=img.header

if d.ndim!=4:
    print("[ERR] SCM requires 4D data.")
    sys.exit(1)

# Apply mask if provided
if mask and os.path.exists(mask):
    m = nib.load(mask).get_fdata()>0
    d = d * m[...,None]

# Compute baseline and signal means
baseline = d[..., b1:b2+1].mean(axis=3)
signal   = d[..., s1:s2+1].mean(axis=3)

# Percent signal change map
psc = 100 * (signal - baseline) / (baseline + 1e-8)

# norm_func
norm = (d - baseline[...,None]) / (baseline[...,None] + 1e-8)

# Save all outputs
np.set_printoptions(suppress=True)

nib.save(nib.Nifti1Image(baseline.astype(np.float32),aff,hdr),
         f"{outdir}/baseline_{b1}_{b2}.nii.gz")
nib.save(nib.Nifti1Image(signal.astype(np.float32),aff,hdr),
         f"{outdir}/signal_{b1}_{b2}.nii.gz")
nib.save(nib.Nifti1Image(psc.astype(np.float32),aff,hdr),
         f"{outdir}/signal_change_map_{b1}_{b2}_{s1}_{s2}.nii.gz")
nib.save(nib.Nifti1Image(norm.astype(np.float32),aff,hdr),
         f"{outdir}/norm_func_{b1}_{b2}.nii.gz")

print("[OK] SCM outputs written to", outdir)
PY
}


###############################################################################
# SECTION 10 — STEP 3: SCM QC (STATIC MAPS)
#
# WHY SCM?
# SCM computes:
#   baseline = mean(baseline vols)
#   signal   = mean(signal vols)
#   PSC = ((signal - baseline) / baseline) * 100
#
# SCM outputs:
#   baseline_*.nii.gz
#   signal_*.nii.gz
#   signal_change_map_*.nii.gz (PSC map)
#   norm_func_*.nii.gz (whole time-normalized series)
#
# These are viewed in FSLeyes with the correct overlay order:
#   norm_func (underlay)
#   baseline
#   signal_change_map (overlay)
###############################################################################

echo -e "${PURPLE}${BOLD}STEP 3: SCM QC${RESET}"
center_line "Baseline / Signal mean + %change maps"

read -p "Run SCM QC? (y/n): " DO_SCM

if [[ "$DO_SCM" == "y" ]]; then

    AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz spatial_smoothed_*.nii.gz temporal_smoothed_*.nii.gz 2>/dev/null))

    # list files
    for i in "${!AVAILABLE[@]}"; do
        printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
    done

    read -p "Select file [1-${#AVAILABLE[@]}]: " idx
    SCM_FILE="${AVAILABLE[$((idx-1))]}"

    # choose windows after reviewing global signal
    plot_global_signal "$SCM_FILE" "$MASK" "$TR"

    read -p "Baseline vols (start end): " B1 B2
    read -p "Signal vols (start end): "   S1 S2

    OUT="soner_pipeline/scm_outputs/$(basename "$SCM_FILE" .nii.gz)_B${B1}-${B2}_S${S1}-${S2}_$(timestamp)"
    mkdir -p "$OUT"

    # SCM computation
    scm_blush_fast "$SCM_FILE" "$B1" "$B2" "$S1" "$S2" "$MASK" "$OUT"

    # open in FSLeyes with correct layer ordering
    nohup fsleyes \
        "$OUT/norm_func_${B1}_${B2}.nii.gz" \
        "$OUT/baseline_${B1}_${B2}.nii.gz" \
        "$OUT/signal_change_map_${B1}_${B2}_${S1}_${S2}.nii.gz" >/dev/null 2>&1 &
fi



###############################################################################
# SECTION 11 — SAFE DEFAULTS FOR FILTERING
#
# WHY THIS SECTION EXISTS:
# Bash's `set -u` causes errors if variables are used before assignment.
#
# Before entering filtering logic, we define empty defaults for:
#   HPF, LPF, BPL, BPH, NC, NW
#
# This prevents "unbound variable" crashes.
###############################################################################

HPF=""; LPF=""; BPL=""; BPH=""; NC=""; NW=""



###############################################################################
# SECTION 12 — STEP 4: TEMPORAL FILTERING (HP / LP / BP / NOTCH)
#
# WHY FILTER?
# Filtering isolates specific frequency components:
#
#   High-pass: remove drift
#   Low-pass: remove high-frequency noise
#   Band-pass: keep specific band only
#   Notch: remove a narrow band (ex: respiration)
#
# INTERNALS:
#   • Butterworth filter (order=4) chosen for smooth response
#   • filtfilt() → zero-phase filtering (no lag)
#   • Only voxels with finite + non-zero variance are filtered
###############################################################################

echo -e "${PURPLE}${BOLD}STEP 4: FILTERING${RESET}"
center_line "High-pass / Low-pass / Band-pass / Notch filtering"

read -p "Do you want to apply a temporal filter? (y/n): " DO_FILT

if [[ "$DO_FILT" == "y" ]]; then

    AVAILABLE=($(ls -tr \
        cleaned_mc_func.nii.gz \
        trimmed_*.nii.gz \
        temporal_smoothed_*.nii.gz \
        spatial_smoothed_*.nii.gz \
        soner_pipeline/detrended/*.nii.gz \
        2>/dev/null))

    for i in "${!AVAILABLE[@]}"; do
        printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
    done

    read -p "Select file to filter [1-${#AVAILABLE[@]}]: " idx
    FILT_FILE="${AVAILABLE[$((idx-1))]}"

    echo -e "\nChoose filter type:"
    echo "1) High-pass"
    echo "2) Low-pass"
    echo "3) Band-pass"
    echo "4) Notch"
    echo "5) Cancel"
    read -p "Choice [1–5]: " FT

    # reset defaults safely
    HPF=""; LPF=""; BPL=""; BPH=""; NC=""; NW=""

    case "$FT" in
        1) read -p "High-pass cutoff (Hz): " HPF ;;
        2) read -p "Low-pass cutoff (Hz): " LPF ;;
        3)
            read -p "Band-pass LOW cutoff (Hz): " BPL
            read -p "Band-pass HIGH cutoff (Hz): " BPH
            ;;
        4)
            read -p "Notch center (Hz): " NC
            read -p "Notch width  (Hz): " NW
            ;;
        5)
            echo "[INFO] Filtering cancelled."
            FT=0
            ;;
        *)
            echo "[ERR] Invalid option"
            FT=0
            ;;
    esac

   #############################################  
    #
    # PYTHON FILTER ENGINE
    #
   #############################################
    if [[ "$FT" != "0" ]]; then

python3 <<PY
import numpy as np, nibabel as nib, sys
from scipy.signal import butter, filtfilt
import datetime

fname = "$FILT_FILE"
ftype = "$FT"
tr = float("$TR")

img = nib.load(fname)
d = img.get_fdata().astype(np.float32)

if d.ndim == 3:
    print("[ERR] Cannot filter 3D data.")
    sys.exit(0)

X, Y, Z, T = d.shape
flat = d.reshape(-1, T)

nyq = 0.5 / tr  # Nyquist

# Choose filter type
if ftype == "1":
    HPF = float("$HPF")
    b,a = butter(4, HPF/nyq, "highpass")

elif ftype == "2":
    LPF = float("$LPF")
    b,a = butter(4, LPF/nyq, "lowpass")

elif ftype == "3":
    BPL = float("$BPL")
    BPH = float("$BPH")
    b,a = butter(4, [BPL/nyq, BPH/nyq], "bandpass")

elif ftype == "4":
    NC = float("$NC")
    NW = float("$NW")
    low  = (NC - NW/2) / nyq
    high = (NC + NW/2) / nyq
    b,a = butter(4, [low, high], "bandstop")

else:
    print("Filter cancelled.")
    sys.exit(0)

# valid voxels
good = np.isfinite(flat).all(axis=1) & (flat.std(axis=1) > 0)

# allocate output
out = flat.copy()

# apply zero-phase filtering
if good.any():
    out[good] = filtfilt(b, a, flat[good], axis=1, method="gust")

# reshape back
out = out.reshape(X, Y, Z, T)

# keep original for NaNs
out[np.isnan(out)] = d[np.isnan(out)]

base = fname.rsplit('/',1)[-1].replace(".nii.gz","")
ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

# create output filename
if ftype == "1":
    outname = f"soner_pipeline/filtering/hp_${HPF}Hz_{base}_{ts}.nii.gz"
elif ftype == "2":
    outname = f"soner_pipeline/filtering/lp_${LPF}Hz_{base}_{ts}.nii.gz"
elif ftype == "3":
    outname = f"soner_pipeline/filtering/bp_${BPL}-${BPH}Hz_{base}_{ts}.nii.gz"
else:
    outname = f"soner_pipeline/filtering/notch_${NC}Hz_{base}_{ts}.nii.gz"

# save
nib.save(nib.Nifti1Image(out.astype(np.float32), img.affine, img.header), outname)
print("[OK] Saved:", outname)
PY

        # Update TMP
        FILT_OUT=$(ls -t soner_pipeline/filtering/*.nii.gz | head -n1)
        TMP="$FILT_OUT"

        # QC
        save_global_png "$TMP" "$MASK" "$TR" "filter"
        plot_global_signal "$TMP" "$MASK" "$TR"
    fi
fi
###############################################################################
# SECTION 13 — STEP 5: DECOMPOSITION (PCA / Temporal ICA / Spatial ICA / PLS)
#
# WHY DECOMPOSITION?
# fMRI signals contain mixtures of:
#   • neural signal
#   • physiological noise (heartbeat, respiration)
#   • motion artifacts
#   • scanner drift
#
# BLuSH provides four decomposition methods:
#
# 1) PCA  (Principal Component Analysis)
#      - finds orthogonal variance components
#      - good for catching high-variance noise
#
# 2) Temporal ICA (tICA)
#      - separates components by temporal independence
#      - ideal for removing periodic physiological noise (resp, heart)
#
# 3) Spatial ICA (sICA)
#      - separates components based on spatial independence
#      - widely used in fMRI (ICA-AROMA style logic)
#
# 4) PLS (Partial Least Squares)
#      - uses a regressor (stimulus/injection window)
#      - finds components most correlated with regressor
#
# GUI FEATURES (MATPLOTLIB):
#   • Up to 25 components shown
#   • Click a component → toggle drop
#   • CTRL-click → mark a specific timepoint
#   • ENTER → save the grid PNG
#
# OUTPUT:
#   soner_pipeline/pca_outputs/<category>/<category>_<method>_func_denoised_TIMESTAMP.nii.gz
#
###############################################################################

while true; do
  echo -e "${PURPLE}${BOLD}STEP 5: DECOMPOSITION${RESET}"
  center_line "Decompose signal via PCA/ICA/PLS and drop unwanted components"

  echo "1) PCA"
  echo "2) Temporal ICA"
  echo "3) Spatial ICA"
  echo "4) PLS"
  echo "5) Exit"
  read -p "Choice [1–5]: " CH

  [[ "$CH" == "5" ]] && break

  case "$CH" in
    1) METHOD="pca" ;;
    2) METHOD="tica" ;;
    3) METHOD="sica" ;;
    4) METHOD="pls" ;;
    *) continue ;;
  esac

  # Gather candidates for decomposition
  AVAILABLE=($(ls -tr \
    cleaned_mc_func.nii.gz \
    trimmed_*.nii.gz \
    spatial_smoothed_*.nii.gz \
    temporal_smoothed_*.nii.gz \
    soner_pipeline/detrended/*.nii.gz \
    soner_pipeline/filtering/*.nii.gz \
    2>/dev/null))

  if [[ ${#AVAILABLE[@]} -eq 0 ]]; then
    echo "[WARN] No candidate NIfTI files found for decomposition."
    break
  fi

  # Choose file
  for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
  read -p "Select file for ${METHOD^^} [1-${#AVAILABLE[@]}]: " idx
  TMP="${AVAILABLE[$((idx-1))]}"

  # Extra input for PLS
  if [[ "$METHOD" == "pls" ]]; then
    read -p "PLS regressor start/end vols (e.g. 1300 1500): " DS1 DS2
  else
    DS1=0; DS2=0
  fi

  # Infer preprocessing category
  SUB="raw"
  [[ "$TMP" == *"hp_"* ]]               && SUB="hp"
  [[ "$TMP" == *"lp_"* ]]               && SUB="lp"
  [[ "$TMP" == *"bp_"* ]]               && SUB="bp"
  [[ "$TMP" == *"notch_"* ]]            && SUB="notch"
  [[ "$TMP" == *"detrended"* ]]         && SUB="detrended"
  [[ "$TMP" == trimmed_* ]]             && SUB="trim"
  [[ "$TMP" == temporal_smoothed_* ]]   && SUB="smoothed"
  [[ "$TMP" == spatial_smoothed_* ]]    && SUB="smoothed"

  # QC before decomposition
  save_global_png "$TMP" "$MASK" "$TR" "${METHOD}_before"



###############################################################################
# SECTION 14 — PYTHON DECOMPOSITION ENGINE & INTERACTIVE GUI
#
# INTERNAL WORKFLOW:
#   1. Load 4D NIfTI
#   2. Reshape into (voxels × time)
#   3. Compute components based on selected method
#   4. Launch GUI to inspect components
#   5. Zero-out dropped components
#   6. Reconstruct denoised data
#   7. Save output in pca_outputs/<category>
#
###############################################################################

python3 <<PY
import os, sys, nibabel as nib, numpy as np, matplotlib.pyplot as plt, datetime
from sklearn.decomposition import PCA, FastICA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

fname="$TMP"
method="$METHOD"
sub="$SUB"
reg_start=int("$DS1")
reg_end=int("$DS2")
ts=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

img=nib.load(fname)
data=img.get_fdata().astype(np.float32)

if data.ndim!=4:
    print("[ERR] Decomposition requires 4D data.")
    sys.exit(0)

X,Y,Z,T=data.shape
flat=data.reshape(-1,T)

out_grid=f"soner_pipeline/qc/grid/{method}_grid_{ts}.png"
os.makedirs(f"soner_pipeline/pca_outputs/{sub}",exist_ok=True)

###############################################################################
# GUI FUNCTION — allows dropping components interactively
###############################################################################
def click_gui(comp,title,grid_path):
    n_show=min(comp.shape[1],25)
    taxis=np.arange(comp.shape[0])
    drop=set()
    fig,axes=plt.subplots(int(np.ceil(n_show/5)),5,figsize=(18,9))
    axes=axes.ravel()
    coord_text=fig.text(0.5,0.01,"",ha="center",fontsize=10)
    markers={}

    def refresh(ax,i):
        for sp in ax.spines.values():
            sp.set_edgecolor('red' if i in drop else 'black')
        ax.set_title(f"C{i+1}"+(" ✗" if i in drop else ""),fontsize=9)

    for i in range(n_show):
        ax=axes[i]
        tsig=(comp[:,i]-comp[:,i].mean())/(comp[:,i].std()+1e-6)
        ax.plot(taxis,tsig,lw=0.8)
        refresh(ax,i)

    plt.suptitle(title+" | click=drop | CTRL=marker | ENTER=save",fontsize=14)

    def on_click(event):
        if event.inaxes not in axes: return
        ax=event.inaxes
        idx=int(np.where(axes==ax)[0])

        # CTRL+click marks a time point on this component
        if event.key in ("control","ctrl"):
            if event.xdata is None: return
            t=int(round(event.xdata))
            if not(0<=t<comp.shape[0]): return
            val=comp[t,idx]
            print(f"[MARK] C{idx+1}: t={t}, val={val:.4f}")

            if idx in markers:
                for obj in markers[idx]: obj.remove()

            vline=ax.axvline(t,color='red',lw=1)
            dot=ax.plot(t,val,'ro',markersize=4)[0]
            txt=ax.text(0.02,0.95,f"t={t}, v={val:.3f}",
                        transform=ax.transAxes,color='red',va='top',fontsize=8)
            markers[idx]=(vline,dot,txt)
            fig.canvas.draw_idle()
            return

        # Simple click toggles component drop
        if idx in drop: drop.remove(idx)
        else: drop.add(idx)
        refresh(ax,idx)
        fig.canvas.draw_idle()

    def on_move(event):
        if event.inaxes not in axes:
            coord_text.set_text("")
            fig.canvas.draw_idle()
            return
        ax=event.inaxes
        idx=int(np.where(axes==ax)[0])
        if event.xdata is None:
            coord_text.set_text(f"C{idx+1}")
        else:
            t=int(round(event.xdata))
            if 0<=t<comp.shape[0]:
                v=comp[t,idx]
                coord_text.set_text(f"C{idx+1}   t={t}   v={v:.3f}")
            else:
                coord_text.set_text(f"C{idx+1}")
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key=="enter":
            plt.savefig(grid_path,dpi=160)
            print("[OK] Saved grid:",grid_path)
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event",on_click)
    fig.canvas.mpl_connect("motion_notify_event",on_move)
    fig.canvas.mpl_connect("key_press_event",on_key)

    plt.show()
    return sorted(drop)


###############################################################################
# METHOD-SPECIFIC COMPUTATION
###############################################################################

# ----- PCA -----
if method=="pca":
    ncomp=min(50,T)
    pca=PCA(n_components=ncomp)
    comps=pca.fit_transform(flat.T)
    drops=click_gui(comps,"PCA Components",out_grid)

    # zero unwanted
    for d in drops:
        if 0<=d<comps.shape[1]:
            comps[:,d]=0

    # inverse_transform → reconstruct 4D
    recon=pca.inverse_transform(comps).T.reshape(X,Y,Z,T)

# ----- TEMPORAL ICA -----
elif method=="tica":
    ncomp=min(25,T)
    ica=FastICA(n_components=ncomp,random_state=0,max_iter=500)
    S=ica.fit_transform(flat.T)
    drops=click_gui(S,"Temporal ICA Components",out_grid)

    for d in drops:
        if 0<=d<S.shape[1]:
            S[:,d]=0

    recon=(S @ ica.mixing_.T).T.reshape(X,Y,Z,T)

# ----- SPATIAL ICA -----
elif method=="sica":
    ncomp=min(25,T)
    flat2=flat-flat.mean(axis=1,keepdims=True)
    Xn=StandardScaler(with_mean=False).fit_transform(flat2.T)
    ica=FastICA(n_components=ncomp,random_state=0,max_iter=500,whiten="unit-variance")

    S=ica.fit_transform(Xn)
    drops=click_gui(S,"Spatial ICA Components",out_grid)

    for d in drops:
        if 0<=d<S.shape[1]:
            S[:,d]=0

    recon=(S @ ica.mixing_.T).T.reshape(X,Y,Z,T)

# ----- PLS -----
elif method=="pls":
    Yv=np.zeros((T,1),dtype=float)
    if 0<=reg_start<T and 0<=reg_end<T and reg_start<=reg_end:
        Yv[reg_start:reg_end+1,0]=1.0
    else:
        print(f"[WARN] Invalid PLS regressor window {reg_start}-{reg_end} for T={T}")

    ncomp=min(25,T)
    pls=PLSRegression(n_components=ncomp)
    pls.fit(flat.T,Yv)
    S=pls.x_scores_

    drops=click_gui(S,"PLS Components",out_grid)
    for d in drops:
        if 0<=d<S.shape[1]:
            S[:,d]=0

    recon=(S @ pls.x_loadings_.T).T.reshape(X,Y,Z,T)

###############################################################################
# SAVE OUTPUT
###############################################################################
outf=f"soner_pipeline/pca_outputs/{sub}/{sub}_{method}_func_denoised_{ts}.nii.gz"
nib.save(nib.Nifti1Image(recon.astype(np.float32),img.affine,img.header),outf)
print("[OK] Saved →", outf)
PY



###############################################################################
# SECTION 15 — POST-DECOMPOSITION SCM & QC
###############################################################################

DEN=$(ls -t soner_pipeline/pca_outputs/${SUB}/${SUB}_${METHOD}_func_denoised_*.nii.gz 2>/dev/null | head -n1 || true)

if [[ -n "$DEN" ]]; then

    save_global_png "$DEN" "$MASK" "$TR" "${METHOD}_after"
    plot_global_signal "$DEN" "$MASK" "$TR"

    echo -e "\n${BLUE}${BOLD}Define baseline and signal windows for denoised output:${RESET}"
    read -p "Baseline vols (start end): " XB1 XB2
    read -p "Signal vols (start end): "   XS1 XS2

    OUTD="soner_pipeline/scm_outputs/${SUB}_${METHOD}_B${XB1}-${XB2}_S${XS1}-${XS2}_$(timestamp)"
    mkdir -p "$OUTD"

    scm_blush_fast "$DEN" "$XB1" "$XB2" "$XS1" "$XS2" "$MASK" "$OUTD"

    # Open with correct overlay ordering
    nohup fsleyes "$OUTD/norm_func_${XB1}_${XB2}.nii.gz" \
                  "$OUTD/baseline_${XB1}_${XB2}.nii.gz" \
                  "$OUTD/signal_change_map_${XB1}_${XB2}_${XS1}_${XS2}.nii.gz" >/dev/null 2>&1 &

else
    echo "[WARN] No denoised file found to run SCM on."
fi

read -p "Run another decomposition? (y/n): " AGAIN
[[ "$AGAIN" != "y" ]] && break

done   # END OF DECOMPOSITION LOOP



###############################################################################
# SECTION 16 — HTML SUMMARY REPORT
#
# WHY THIS EXISTS:
# After the full BLuSH pipeline:
#   • multiple global plots are generated
#   • multiple decomposition grids exist
#
# This block produces a consolidated HTML QC page containing:
#   • All global mean plots
#   • All ICA/PCA/PLS component grids
#
# Automatically opens in browser (Linux, macOS, Windows)
###############################################################################

echo -e "${BLUE}${BOLD}Generating HTML QC summary…${RESET}"

python3 <<'PY'
import os, glob, datetime, html, platform, subprocess

ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

QC_GLOBAL=sorted(glob.glob("soner_pipeline/qc/global/*.png"))
QC_GRID=sorted(glob.glob("soner_pipeline/qc/grid/*.png"))

def make_section(title,files):
    if not files:
        return f"<h2>{title}</h2><p>No images found.</p>"
    rows="\n".join(
        f"<tr><td>{html.escape(os.path.basename(p))}</td>"
        f"<td><img src='{p}' width='360'></td></tr>"
        for p in files
    )
    return f"<h2>{title}</h2><table>{rows}</table>"

html_doc=f"""
<html>
<head>
<meta charset='utf-8'>
<title>BLuSH QC Summary</title>

<style>
body {{
    background:#0f1116;
    color:#eaeef2;
    font-family:Arial, sans-serif;
    padding:24px;
}}
h1 {{ color:#8bd3ff; }}
h2 {{ color:#62b0ff; margin-top:30px; }}
td, th {{
    border:1px solid #2a2f3a;
    padding:6px;
}}
table {{
    border-collapse:collapse;
    margin-bottom:24px;
}}
</style>

</head>
<body>

<h1>BLuSH QC Summary</h1>
<p><i>Generated: {ts}</i></p>

{make_section("Global Mean Plots", QC_GLOBAL)}
{make_section("Component Grid Plots", QC_GRID)}

</body>
</html>
"""

OUT=f"soner_pipeline/qc_summary_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.html"

with open(OUT,"w",encoding="utf-8") as f: f.write(html_doc)
print("[OK] Summary HTML →",OUT)

try:
    if platform.system()=="Linux":
        subprocess.Popen(["xdg-open",OUT])
    elif platform.system()=="Darwin":
        subprocess.Popen(["open",OUT])
    elif platform.system()=="Windows":
        subprocess.Popen(["start",OUT],shell=True)
except Exception as e:
    print("[WARN] Could not auto-open HTML:",e)

PY
