#!/bin/bash
# =====================================================================
#  BLuSH v20.1  (based on v20 + skip options)
#
#  New in v20.1:
#    ✔ Option to skip DETRENDING (Step 2.5)
#    ✔ Option to skip SCM QC (Step 3)
#    ✔ Everything else identical to v20
# =====================================================================

set -Eeuo pipefail
trap 'echo "[ERROR] line ${LINENO}" >&2' ERR
shopt -s nullglob
export LC_NUMERIC=C
export PYTHONWARNINGS="ignore"

# ---------------- Colors ----------------
RED="\033[38;5;160m"; PURPLE="\033[38;5;135m"; GREEN="\033[92m"
BLUE="\033[38;5;33m"; RESET="\033[0m"; BOLD="\033[1m"

timestamp(){ date +"%Y%m%d_%H%M%S"; }

# ---------------- Folders ----------------
mkdir -p soner_pipeline/{filtering,scm_outputs,logs,detrended,qc}
mkdir -p soner_pipeline/pca_outputs/{raw,hp,lp,bp,detrended,trim,smoothed}

MASK="mask_mean_mc_func.nii.gz"
[[ ! -f cleaned_mc_func.nii.gz ]] && { echo "[ERR] cleaned_mc_func.nii.gz missing"; exit 1; }
[[ ! -f "$MASK" ]] && MASK=""

# ---------------- Helpers ----------------
center_line(){
  printf "\n────────────────────────────────────────────────────────\n✳️  %s ✳️\n────────────────────────────────────────────────────────\n" "$1"
}

get_tr(){
python3 - "$1" <<'PY'
import sys, nibabel as nib
try:
    print(float(nib.load(sys.argv[1]).header.get_zooms()[3]))
except:
    print(1.0)
PY
}

# =====================================================================
#  GLOBAL MEAN PNG
# =====================================================================
save_global_png(){
  local nii="$1"; local mask="$2"; local tr="$3"; local label="$4"
  mkdir -p soner_pipeline/qc
  python3 - "$nii" "$mask" "$tr" "$label" <<'PY'
import sys, os, nibabel as nib, numpy as np, matplotlib.pyplot as plt, time
f,m,tr,label=sys.argv[1:5]; tr=float(tr)

img=nib.load(f); d=img.get_fdata()
try:
    mask=nib.load(m).get_fdata()>0
    d=d[mask]
except:
    d=d.reshape(-1,d.shape[-1])

sig=d.mean(0)
x=np.arange(len(sig))*tr

plt.figure(figsize=(10,4))
plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean – {os.path.basename(f)}")
plt.xlabel("Time (s)")
plt.ylabel("Mean Intensity")
plt.grid(True,alpha=0.3)
plt.tight_layout()

out=f"soner_pipeline/qc/global_{label}_{time.strftime('%Y%m%d_%H%M%S')}.png"
plt.savefig(out,dpi=160)
plt.close()
print("[OK] Global PNG:", out)
PY
}

# =====================================================================
#  SPECTRA PNG
# =====================================================================
save_spectra(){
  local nii="$1"; local mask="$2"; local tr="$3"; local label="$4"
  mkdir -p soner_pipeline/qc
  python3 - "$nii" "$mask" "$tr" "$label" <<'PY'
import sys, os, nibabel as nib, numpy as np, matplotlib.pyplot as plt, time
f,m,tr,label=sys.argv[1:5]; tr=float(tr)

img=nib.load(f); d=img.get_fdata()
try:
    mask=nib.load(m).get_fdata()>0
    d=d[mask]
except:
    d=d.reshape(-1,d.shape[-1])

sig=d.mean(0)
n=len(sig)

freq=np.fft.rfftfreq(n,d=tr)
pow=(np.abs(np.fft.rfft(sig-np.mean(sig)))**2)

for lo,hi in [(0,0.1),(0,0.5)]:
    sel=(freq>=lo)&(freq<=hi)
    plt.figure(figsize=(8,4))
    plt.plot(freq[sel],pow[sel],'k',lw=1)
    plt.title(f"Spectrum {label} ({lo}-{hi} Hz)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power")
    plt.grid(True,alpha=0.3)
    plt.tight_layout()
    out=f"soner_pipeline/qc/spectrum_{label}_{lo}-{hi}Hz_{time.strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(out,dpi=160)
    plt.close()
    print("[OK] Spectrum PNG:", out)
PY
}

# =====================================================================
#  INTERACTIVE GLOBAL SIGNAL PLOT
# =====================================================================
plot_global_signal(){
python3 - "$1" "$2" "$3" <<'PY'
import sys, nibabel as nib, numpy as np, matplotlib.pyplot as plt
f,m,tr=sys.argv[1],sys.argv[2],float(sys.argv[3])

img=nib.load(f); d=img.get_fdata()
try:
    mask=nib.load(m).get_fdata()>0
    d=d[mask]
except:
    d=d.reshape(-1,d.shape[-1])

sig=d.mean(0)
x=np.arange(len(sig))*tr

plt.figure(figsize=(10,4))
plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean Signal – {f}")
plt.xlabel("Time (s)")
plt.ylabel("Mean Intensity")
plt.grid(True,alpha=0.3)
plt.tight_layout()
plt.show()
PY
}

# =====================================================================
#  SCM CREATOR (same as v20)
# =====================================================================
scm_blush_fast(){
  local func="$1" base="$2" bend="$3" sig="$4" send="$5" mask="$6" outd="$7"
  mkdir -p "$outd"

  3dTstat -mean -prefix "$outd/baseline_${base}_${bend}.nii.gz" "$func[${base}..${bend}]"
  3dTstat -mean -prefix "$outd/signal_${sig}_${send}.nii.gz"   "$func[${sig}..${send}]"

  fslmaths "$outd/signal_${sig}_${send}" -sub "$outd/baseline_${base}_${bend}" \
           -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+ -mas "$mask"} \
           "$outd/signal_change_map_${base}_${bend}_${sig}_${send}.nii.gz"

  fslmaths "$func" -sub "$outd/baseline_${base}_${bend}" \
           -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+ -mas "$mask"} \
           "$outd/norm_func_${base}_${bend}.nii.gz"
}

# =====================================================================
#  MAIN PIPELINE — STEP 1 (TRIMMING)
# =====================================================================
TR=$(get_tr cleaned_mc_func.nii.gz)
TMP="cleaned_mc_func.nii.gz"

echo -e "${RED}${BOLD}====== BLuSH v20.1 ======${RESET}"

# --------------------
# STEP 1: TRIMMING
# --------------------
echo -e "${PURPLE}${BOLD}STEP 1: TRIMMING${RESET}"
center_line "Trim initial/end volumes for a stable baseline"

read -p "Trim from START (vols): " TSTART
read -p "Trim from END   (vols): " TEND

if [[ "${TSTART:-0}" != "0" || "${TEND:-0}" != "0" ]]; then
  END=$(python3 -c "import nibabel as nib; n=nib.load('cleaned_mc_func.nii.gz').shape[-1]; print(n-int('$TEND')-1)")
  OUT="trimmed_${TSTART}-${TEND}_cleaned_mc_func_$(timestamp).nii.gz"
  3dTcat -prefix "$OUT" cleaned_mc_func.nii.gz"[${TSTART}..${END}]"
  TMP="$OUT"
fi

save_global_png "$TMP" "$MASK" "$TR" "trim"
save_spectra     "$TMP" "$MASK" "$TR" "trim"

# --------------------
# STEP 2: SMOOTHING
# --------------------
echo -e "${PURPLE}${BOLD}STEP 2: SMOOTHING${RESET}"
center_line "Temporal smoothing reduces noise; spatial smoothing boosts SNR"

AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done

read -p "Select file to smooth [1-${#AVAILABLE[@]}]: " idx
TMP="${AVAILABLE[$((idx-1))]}"

read -p "Temporal smoothing (seconds, 0=skip): " TSM
read -p "Spatial σ blur (mm, 0=skip): " SSM

if [[ "${TSM:-0}" != "0" ]]; then
  python3 - <<PY "$TMP" "$TSM" "$TR"
import sys,nibabel as nib,numpy as np
from scipy.ndimage import uniform_filter1d

f = sys.argv[1]
sec = float(sys.argv[2])
tr = float(sys.argv[3])
win = max(1, int(round(sec/tr)))

img = nib.load(f)
d = img.get_fdata()
sm = uniform_filter1d(d, size=win, axis=-1, mode='nearest')

out = f"temporal_smoothed_{sec}s_{f.rsplit('/',1)[-1]}"
nib.save(nib.Nifti1Image(sm.astype(np.float32), img.affine, img.header), out)
print("[OK] Temporal smoothing →", out)
PY

  TMP="temporal_smoothed_${TSM}s_${TMP##*/}"
  save_global_png "$TMP" "$MASK" "$TR" "tsmooth"
  save_spectra     "$TMP" "$MASK" "$TR" "tsmooth"
fi

# --------------------
# SPATIAL SMOOTHING  (add σ naming)
# --------------------
if [[ "${SSM:-0}" != "0" ]]; then
  OUT="spatial_smoothed_${SSM}mm_${TMP%.nii.gz}_$(timestamp).nii.gz"
  fslmaths "$TMP" -s "$SSM" "$OUT"
  echo "[OK] Spatial smoothing → $OUT"
  TMP="$OUT"

  save_global_png "$TMP" "$MASK" "$TR" "ssmooth"
  save_spectra     "$TMP" "$MASK" "$TR" "ssmooth"
fi

# =====================================================================
# STEP 2.5: DETRENDING  (NEW skip option)
# =====================================================================
echo -e "${PURPLE}${BOLD}STEP 2.5: DETRENDING${RESET}"
center_line "Remove scanner drifts with polynomial fits"

read -p "Do you want to do detrending? (y/n): " DO_DET
if [[ "$DO_DET" == "y" ]]; then

  AVAILABLE=($(ls -tr \
    cleaned_mc_func.nii.gz \
    trimmed_*.nii.gz \
    temporal_smoothed_*.nii.gz \
    spatial_smoothed_*.nii.gz \
    2>/dev/null))

  for i in "${!AVAILABLE[@]}"; do
    printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
  done

  read -p "Select file for detrending [1-${#AVAILABLE[@]}]: " idx
  TMP="${AVAILABLE[$((idx-1))]}"

  plot_global_signal "$TMP" "$MASK" "$TR"

  read -p "Apply detrending? (y/n): " DETR
  if [[ "$DETR" == "y" ]]; then
    while true; do
      echo "1) Linear      – remove constant + linear drift"
      echo "2) Quadratic   – remove slow parabolic drift"
      echo "3) Polynomial  – higher-orderbaseline removal"
      echo "4) Skip"
      read -p "Select detrend type [1–4]: " DTYPE

      case "$DTYPE" in
        1) POL=1;;
        2) POL=2;;
        3) read -p "Enter polynomial order (>=3): " POL;;
        4) break;;
        *) continue;;
      esac

      BASE=$(basename "$TMP" .nii.gz)
      OUT="soner_pipeline/detrended/detrended_${BASE}_pol${POL}_$(timestamp).nii.gz"

      3dDetrend -polort $POL -prefix "$OUT" "$TMP"
      [[ -f "$OUT" ]] && echo "[OK] Detrended → $OUT"

      save_global_png "$OUT" "$MASK" "$TR" "detrended_pol${POL}"
      save_spectra     "$OUT" "$MASK" "$TR" "detrended_pol${POL}"

      plot_global_signal "$OUT" "$MASK" "$TR"

      read -p "Use this detrended file? (y/n): " OK
      [[ "$OK" == "y" ]] && TMP="$OUT" && break
    done
  fi

fi   # END skip detrending
# =====================================================================
#  STEP 3: SCM QC   (NEW skip option)
# =====================================================================
echo -e "${PURPLE}${BOLD}STEP 3: SCM QC${RESET}"
center_line "Build baseline/signal means and %change maps (QC)"

read -p "Do you want to run SCM QC? (y/n): " DO_SCM
if [[ "$DO_SCM" == "y" ]]; then

  AVAILABLE=($(ls -tr \
    cleaned_mc_func.nii.gz \
    trimmed_*.nii.gz \
    spatial_smoothed_*.nii.gz \
    temporal_smoothed_*.nii.gz \
    soner_pipeline/detrended/*.nii.gz \
    soner_pipeline/filtering/*.nii.gz \
    2>/dev/null))

  for i in "${!AVAILABLE[@]}"; do
    printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
  done

  read -p "Select file for SCM QC [1-${#AVAILABLE[@]}]: " idx
  SCM_FILE="${AVAILABLE[$((idx-1))]}"

  plot_global_signal "$SCM_FILE" "$MASK" "$TR"

  read -p "Baseline vols (start end): " B1 B2
  read -p "Signal   vols (start end): " S1 S2

  OUT="soner_pipeline/scm_outputs/$(basename "$SCM_FILE" .nii.gz)_B${B1}-${B2}_S${S1}-${S2}_$(timestamp)"
  mkdir -p "$OUT"

  scm_blush_fast "$SCM_FILE" "$B1" "$B2" "$S1" "$S2" "$MASK" "$OUT"

  # Correct FSLeyes order:
  #   bottom  → norm_func
  #   middle  → baseline
  #   top     → signal_change_map
  nohup fsleyes "$OUT/norm_func_${B1}_${B2}.nii.gz" \
                "$OUT/baseline_${B1}_${B2}.nii.gz" \
                "$OUT/signal_change_map_${B1}_${B2}_${S1}_${S2}.nii.gz" \
                >/dev/null 2>&1 &
fi   # END skip SCM QC


# =====================================================================
#  STEP 4: FILTERING
# =====================================================================
echo -e "${PURPLE}${BOLD}STEP 4: FILTERING${RESET}"
center_line "High-pass removes drifts; Low-pass suppresses fast noise"

AVAILABLE=($(ls -tr \
  cleaned_mc_func.nii.gz \
  trimmed_*.nii.gz \
  spatial_smoothed_*.nii.gz \
  temporal_smoothed_*.nii.gz \
  soner_pipeline/detrended/*.nii.gz \
  2>/dev/null))

for i in "${!AVAILABLE[@]}"; do
  printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
done

read -p "Select file to filter [1-${#AVAILABLE[@]}]: " idx
TMP="${AVAILABLE[$((idx-1))]}"

while true; do
  echo "1) High-pass   – remove slow scanner drift (< cutoff)"
  echo "2) Low-pass    – suppress fast noise (< cutoff)"
  echo "3) Band-pass   – retain band (low–high)"
  echo "4) Continue → Decomposition"
  read -p "Choose filter type [1–4]: " FTYPE

  [[ "$FTYPE" == "4" ]] && break

  stamp=$(timestamp)

  case "$FTYPE" in
    1)
      read -p "Cut-off Hz (e.g. 0.01): " FC
      OUT="soner_pipeline/filtering/hp_${FC}Hz_${stamp}.nii.gz"
      HP=$(python3 -c "print(1/(2*$TR*$FC))")
      fslmaths "$TMP" -bptf $HP -1 "$OUT"
      TMP="$OUT"
      ;;
    2)
      read -p "Cut-off Hz (e.g. 0.1): " FC
      OUT="soner_pipeline/filtering/lp_${FC}Hz_${stamp}.nii.gz"
      LP=$(python3 -c "print(1/(2*$TR*$FC))")
      fslmaths "$TMP" -bptf -1 $LP "$OUT"
      TMP="$OUT"
      ;;
    3)
      read -p "Band range Hz (e.g. 0.01 0.1): " FL FH
      OUT="soner_pipeline/filtering/bp_${FL}-${FH}Hz_${stamp}.nii.gz"
      HP=$(python3 -c "print(1/(2*$TR*$FL))")
      LP=$(python3 -c "print(1/(2*$TR*$FH))")
      fslmaths "$TMP" -bptf $HP $LP "$OUT"
      TMP="$OUT"
      ;;
    *)
      continue
      ;;
  esac

  # ----- QC for filtered output -----
  save_global_png "$TMP" "$MASK" "$TR" "filter"
  save_spectra     "$TMP" "$MASK" "$TR" "filter"

  # ----- AFNI despiking -----
  read -p "Run AFNI despiking on this output? (y/n): " DESP
  if [[ "$DESP" == "y" ]]; then
    DSP="soner_pipeline/filtering/despiked_${stamp}.nii.gz"
    3dDespike -NEW -localedit -ignore 5 -prefix "$DSP" "$TMP" || true

    if [[ -f "$DSP" ]]; then
      TMP="$DSP"
      save_global_png "$TMP" "$MASK" "$TR" "despiked"
      save_spectra     "$TMP" "$MASK" "$TR" "despiked"
    fi
  fi
done
# =====================================================================
#  STEP 5: DECOMPOSITION  (PCA / TICA / SICA / PLS)
# =====================================================================
while true; do
  echo -e "${PURPLE}${BOLD}STEP 5: DECOMPOSITION${RESET}"
  center_line "Decompose signal via PCA/ICA/PLS; drop unwanted components"

  echo "1) PCA  – variance-driven denoising"
  echo "2) Temporal ICA  – independent time-courses"
  echo "3) Spatial ICA   – independent spatial maps"
  echo "4) PLS           – regressor-correlated components"
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

  AVAILABLE=($(ls -tr \
    cleaned_mc_func.nii.gz \
    trimmed_*.nii.gz \
    spatial_smoothed_*.nii.gz \
    temporal_smoothed_*.nii.gz \
    soner_pipeline/detrended/*.nii.gz \
    soner_pipeline/filtering/*.nii.gz \
    2>/dev/null))

  for i in "${!AVAILABLE[@]}"; do
    printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"
  done

  read -p "Select file for ${METHOD^^} [1-${#AVAILABLE[@]}]: " idx
  TMP="${AVAILABLE[$((idx-1))]}"

  # PLS regressor window
  if [[ "$METHOD" == "pls" ]]; then
    read -p "PLS regressor start/end vols (e.g. 1300 1500): " DS1 DS2
  else
    DS1=0
    DS2=0
  fi

  # ----- Determine output category -----
  SUB="raw"
  [[ "$TMP" == *"hp_"* ]]          && SUB="hp"
  [[ "$TMP" == *"lp_"* ]]          && SUB="lp"
  [[ "$TMP" == *"bp_"* ]]          && SUB="bp"
  [[ "$TMP" == *"detrended"* ]]    && SUB="detrended"
  [[ "$TMP" == trimmed_* ]]         && SUB="trim"
  [[ "$TMP" == temporal_smoothed_* ]] && SUB="smoothed"
  [[ "$TMP" == spatial_smoothed_* ]]  && SUB="smoothed"

  # ----- QC before decomposition -----
  save_global_png "$TMP" "$MASK" "$TR" "${METHOD}_before"
  save_spectra     "$TMP" "$MASK" "$TR" "${METHOD}_before"


  # =====================================================================
  #  PYTHON DECOMPOSITION WITH MOUSE-TRACKING GUI
  # =====================================================================
python3 <<PY
import os, sys, nibabel as nib, numpy as np, matplotlib.pyplot as plt, datetime
from sklearn.decomposition import PCA, FastICA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

fname = "$TMP"
method = "$METHOD"
sub = "$SUB"
reg_start = int("$DS1")
reg_end   = int("$DS2")

img = nib.load(fname)
data = img.get_fdata().astype(np.float32)
X, Y, Z, T = data.shape
flat = data.reshape(-1, T)

os.makedirs(f"soner_pipeline/pca_outputs/{sub}", exist_ok=True)

# ======================================================================
#  COMPONENT VIEWER WITH MOUSE-TRACKING AND CTRL-PRECISE MARKER
# ======================================================================
def click_gui(comp, title, path):
    n_show = min(comp.shape[1], 25)
    taxis = np.arange(comp.shape[0])
    drop = set()

    fig, axes = plt.subplots(int(np.ceil(n_show/5)), 5, figsize=(18, 9))
    axes = axes.ravel()

    coord_text = fig.text(0.5, 0.01, "", ha="center", va="bottom", fontsize=10)

    markers = {}  # store (vline, dot, label) per component

    # ----- Update subplot border + title -----
    def refresh(ax, i):
        for sp in ax.spines.values():
            sp.set_edgecolor('red' if i in drop else 'black')
        ax.set_title(f"C{i+1}" + (" ✗" if i in drop else ""), fontsize=9)

    # ----- Initial plotting -----
    for i in range(n_show):
        ax = axes[i]
        ts = (comp[:,i] - np.mean(comp[:,i])) / (np.std(comp[:,i]) + 1e-8)
        ax.plot(taxis, ts, lw=0.8)
        refresh(ax, i)

    plt.suptitle(title + " | Click=drop | CTRL-click=marker | ENTER=accept", fontsize=13)

    # =============================
    # CLICK EVENTS
    # =============================
    def on_click(event):
        if event.inaxes not in axes:
            return

        ax = event.inaxes
        idx = np.where(axes == ax)[0][0]

        # ---- CTRL click: precise marker ----
        if event.key == "control" or event.key == "ctrl":
            if event.xdata is None: return
            t = int(round(event.xdata))
            if not (0 <= t < comp.shape[0]): return

            value = comp[t, idx]
            print(f"[MARKER] Component {idx+1}: t={t}, value={value:.4f}")

            # remove older markers
            if idx in markers:
                for obj in markers[idx]:
                    obj.remove()

            vline = ax.axvline(t, color='red', lw=1)
            dot   = ax.plot(t, value, 'ro', markersize=4)[0]
            label = ax.text(
                0.02, 0.95,
                f"t={t}, v={value:.3f}",
                transform=ax.transAxes, color='red',
                fontsize=8, va='top'
            )
            markers[idx] = (vline, dot, label)
            fig.canvas.draw_idle()
            return

        # ---- Normal click: toggle drop ----
        if idx in drop: drop.remove(idx)
        else: drop.add(idx)
        refresh(ax, idx)
        fig.canvas.draw_idle()

    # =============================
    # MOUSE MOVE = coordinate display
    # =============================
    def on_move(event):
        if event.inaxes not in axes:
            coord_text.set_text("")
            fig.canvas.draw_idle()
            return

        ax = event.inaxes
        idx = np.where(axes == ax)[0][0]

        t = int(round(event.xdata)) if event.xdata is not None else None
        if t is not None and 0 <= t < comp.shape[0]:
            value = comp[t, idx]
            coord_text.set_text(f"C{idx+1}   t={t}   value={value:.3f}")
        else:
            coord_text.set_text(f"C{idx+1}")

        fig.canvas.draw_idle()

    # =============================
    # ENTER = accept
    # =============================
    def on_key(event):
        if event.key == "enter":
            plt.savefig(path, dpi=180)
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("motion_notify_event", on_move)
    fig.canvas.mpl_connect("key_press_event", on_key)

    plt.show()
    return sorted(drop)


# ======================================================================
#  Perform decomposition
# ======================================================================
if method == "pca":
    pca = PCA(n_components=min(50,T))
    comps = pca.fit_transform(flat.T)
    drop = click_gui(comps, "PCA Components", "soner_pipeline/qc/pca_grid.png")
    for d in drop: comps[:,d] = 0
    recon = pca.inverse_transform(comps).T.reshape(X,Y,Z,T)

elif method == "tica":
    ica = FastICA(n_components=25, random_state=0, max_iter=500)
    S = ica.fit_transform(flat.T)
    drop = click_gui(S, "Temporal ICA Components", "soner_pipeline/qc/tica_grid.png")
    for d in drop: S[:,d] = 0
    recon = (S @ ica.mixing_.T).T.reshape(X,Y,Z,T)

elif method == "sica":
    flat2 = flat - flat.mean(axis=1, keepdims=True)
    Xn = StandardScaler(with_mean=False).fit_transform(flat2.T)
    ica = FastICA(n_components=25, random_state=0, max_iter=500, whiten="unit-variance")
    S = ica.fit_transform(Xn)
    drop = click_gui(S, "Spatial ICA Components", "soner_pipeline/qc/sica_grid.png")
    for d in drop: S[:,d] = 0
    recon = (S @ ica.mixing_.T).T.reshape(X,Y,Z,T)

elif method == "pls":
    Yv = np.zeros((T,1))
    Yv[reg_start:reg_end+1] = 1
    pls = PLSRegression(n_components=25)
    pls.fit(flat.T, Yv)
    S = pls.x_scores_
    drop = click_gui(S, "PLS Components", "soner_pipeline/qc/pls_grid.png")
    for d in drop: S[:,d] = 0
    recon = (S @ pls.x_loadings_.T).T.reshape(X,Y,Z,T)

# ======================================================================
#  Save output
# ======================================================================
outf = (
    f"soner_pipeline/pca_outputs/{sub}/"
    f"{sub}_{method}_func_denoised_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.nii.gz"
)

nib.save(nib.Nifti1Image(recon.astype(np.float32), img.affine, img.header), outf)
print("[OK] Saved →", outf)
PY

  # ----- Locate saved output -----
  DEN=$(ls -t "soner_pipeline/pca_outputs/${SUB}/${SUB}_${METHOD}_func_denoised_"*.nii.gz 2>/dev/null | head -n1)
  [[ -f "$DEN" ]] || { echo "[ERR] No denoised file found."; continue; }

  # ----- QC AFTER decomposition -----
  plot_global_signal "$DEN" "$MASK" "$TR"
  save_global_png "$DEN" "$MASK" "$TR" "${METHOD}_after"
  save_spectra     "$DEN" "$MASK" "$TR" "${METHOD}_after"

  echo -e "\n${BLUE}${BOLD}Define baseline and signal windows for denoised output:${RESET}"
  read -p "Baseline vols (start end): " XB1 XB2
  read -p "Signal vols (start end): "   XS1 XS2

  OUTD="soner_pipeline/scm_outputs/${SUB}_${METHOD}_B${XB1}-${XB2}_S${XS1}-${XS2}_$(timestamp)"
  scm_blush_fast "$DEN" "$XB1" "$XB2" "$XS1" "$XS2" "$MASK" "$OUTD"

  nohup fsleyes "$OUTD/norm_func_${XB1}_${XB2}.nii.gz" \
                "$OUTD/baseline_${XB1}_${XB2}.nii.gz" \
                "$OUTD/signal_change_map_${XB1}_${XB2}_${XS1}_${XS2}.nii.gz" \
                >/dev/null 2>&1 &

  read -p "Run another decomposition? (y/n): " AGAIN
  [[ "$AGAIN" != "y" ]] && break
done


# =====================================================================
#  HTML QC SUMMARY
# =====================================================================
echo -e "${BLUE}${BOLD}Generating HTML QC summary…${RESET}"

python3 <<PY
import os, glob, datetime, html

ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
pngs = sorted(glob.glob('soner_pipeline/qc/*.png'))

rows = "\n".join(
    f"<tr><td>{html.escape(os.path.basename(p))}</td>"
    f"<td><img src='{p}' width='360'></td></tr>"
    for p in pngs
)

doc = f"""
<html><head><meta charset='utf-8'><title>BLuSH QC Summary</title>
<style>
body{{background:#0f1116;color:#eaeef2;font-family:Arial;padding:24px}}
h1{{color:#8bd3ff}}
table{{border-collapse:collapse}}
td,th{{border:1px solid #2a2f3a;padding:6px}}
</style>
</head><body>
<h1>BLuSH QC Summary</h1>
<p><i>{ts}</i></p>
<table>{rows}</table>
</body></html>
"""

out = f"soner_pipeline/qc_summary_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
open(out, "w", encoding="utf-8").write(doc)
print("[OK] Summary HTML →", out)

try:
    import platform, subprocess
    if platform.system()=="Darwin":
        subprocess.Popen(["open", out])
    elif platform.system()=="Linux":
        subprocess.Popen(["xdg-open", out])
    elif platform.system()=="Windows":
        s



