#!/bin/bash
# =====================================================================
#  filtering_decomposition_v18_corrected.sh
#
#  v18 oriented to v17 (same control flow):
#    1) TRIMMING → 2) SMOOTHING → 2.5) DETRENDING → 3) SCM QC
#    → 4) FILTERING → 5) DECOMPOSITION (with loop / returns)
#
#  Fixes:
#   - Correct FSLeyes order: norm_func (bottom) → baseline → signal_change_map (top)
#   - Timestamps added to outputs to prevent overwrites
#   - Save global-mean PNG (‘blob up’) + spectra PNGs (0–0.1 Hz, 0–0.5 Hz)
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

# -------------- Helpers -----------------
center_line(){ printf "\n────────────────────────────────────────────────────────\n✳️  %s ✳️\n────────────────────────────────────────────────────────\n" "$1"; }

get_tr(){ python3 - "$1" <<'PY'
import sys,nibabel as nib
try: print(float(nib.load(sys.argv[1]).header.get_zooms()[3]))
except: print(1.0)
PY
}

# --- Save global mean as PNG (non-interactive) ---
save_global_png(){
  local nii="$1"; local mask="$2"; local tr="$3"; local label="$4"
  mkdir -p soner_pipeline/qc
  python3 - "$nii" "$mask" "$tr" "$label" <<'PY'
import sys,os,nibabel as nib,numpy as np,matplotlib.pyplot as plt, time
f,m,tr,label=sys.argv[1:5]; tr=float(tr)
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); x=np.arange(len(sig))*tr
plt.figure(figsize=(10,4)); plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean – {os.path.basename(f)}"); plt.xlabel("Time (s)"); plt.ylabel("Mean Intensity")
plt.grid(True,alpha=0.3); plt.tight_layout()
out=f"soner_pipeline/qc/global_{label}_{time.strftime('%Y%m%d_%H%M%S')}.png"
plt.savefig(out,dpi=160); plt.close(); print("[OK] Global PNG:",out)
PY
}
# --- Save spectra in two bands (0–0.1 Hz, 0–0.5 Hz) ---
save_spectra(){
  local nii="$1"; local mask="$2"; local tr="$3"; local label="$4"
  mkdir -p soner_pipeline/qc
  python3 - "$nii" "$mask" "$tr" "$label" <<'PY'
import sys,os,nibabel as nib,numpy as np,matplotlib.pyplot as plt, time
f,m,tr,label=sys.argv[1:5]; tr=float(tr)
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); n=len(sig)
freq=np.fft.rfftfreq(n,d=tr); pow=(np.abs(np.fft.rfft(sig-np.mean(sig)))**2)
for lo,hi in [(0,0.1),(0,0.5)]:
    sel=(freq>=lo)&(freq<=hi)
    plt.figure(figsize=(8,4)); plt.plot(freq[sel],pow[sel],'k',lw=1)
    plt.title(f"Spectrum {label} ({lo}-{hi} Hz)"); plt.xlabel("Frequency (Hz)"); plt.ylabel("Power")
    plt.grid(True,alpha=0.3); plt.tight_layout()
    out=f"soner_pipeline/qc/spectrum_{label}_{lo}-{hi}Hz_{time.strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(out,dpi=160); plt.close(); print("[OK] Spectrum PNG:",out)
PY
}

# --- Interactive global plot (as in v17) ---
plot_global_signal(){
python3 - "$1" "$2" "$3" <<'PY'
import sys,nibabel as nib,numpy as np,matplotlib.pyplot as plt,os
f,m,tr=sys.argv[1],sys.argv[2],float(sys.argv[3])
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); x=np.arange(len(sig))*tr
plt.figure(figsize=(10,4)); plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean Signal – {os.path.basename(f)}")
plt.xlabel("Time (s)"); plt.ylabel("Mean Intensity")
plt.grid(True,alpha=0.3); plt.tight_layout(); plt.show()
PY
}
# --- SCM creator (unchanged math) ---
scm_blush_fast(){
  local func="$1" base="$2" bend="$3" sig="$4" send="$5" mask="$6" outd="$7"
  mkdir -p "$outd"
  3dTstat -mean -prefix "$outd/baseline_${base}_${bend}.nii.gz" "$func[${base}..${bend}]"
  3dTstat -mean -prefix "$outd/signal_${sig}_${send}.nii.gz"   "$func[${sig}..${send}]"
  fslmaths "$outd/signal_${sig}_${send}" -sub "$outd/baseline_${base}_${bend}" \
    -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
    "$outd/signal_change_map_${base}_${bend}_${sig}_${send}.nii.gz"
  fslmaths "$func" -sub "$outd/baseline_${base}_${bend}" \
    -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
    "$outd/norm_func_${base}_${bend}.nii.gz"
}

# --- Filtering wrapper (hp/lp/bp with FSL -bptf) ---
apply_freq_filter(){
  local inp="$1"; local outp="$2"; local fl="$3"; local fh="$4"; local tr="$5"; local ftype="$6"
  case "$ftype" in
    1) HP=$(python3 -c "print(1/(2*$tr*$fl))"); fslmaths "$inp" -bptf $HP -1 "$outp";;
    2) LP=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf -1 $LP "$outp";;
    3) HP=$(python3 -c "print(1/(2*$tr*$fl))"); LP=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf $HP $LP "$outp";;
  esac
}

# --- DECOMPOSITION (interactive component selection + live global plot) ---
run_decomposition(){
  local infile="$1" method="$2" subdir="$3" reg_start="$4" reg_end="$5"
  python3 - "$infile" "$method" "$subdir" "$reg_start" "$reg_end" <<'PY'
import sys,os,nibabel as nib,numpy as np,matplotlib.pyplot as plt,datetime
from sklearn.decomposition import PCA,FastICA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

fname,method,sub,reg_start,reg_end=sys.argv[1:6]
img=nib.load(fname); data=img.get_fdata().astype(np.float32)
X,Y,Z,T=data.shape; flat=data.reshape(-1,T)
os.makedirs(f"soner_pipeline/pca_outputs/{sub}",exist_ok=True)
reg_start,reg_end=int(reg_start),int(reg_end)

def click_gui(comp,title,path):
    n_show=min(comp.shape[1],25)
    taxis=np.arange(comp.shape[0])
    drop=set()
    fig,axes=plt.subplots(int(np.ceil(n_show/5)),5,figsize=(18,9))
    axes=axes.ravel()
    def refresh(ax,i):
        for sp in ax.spines.values():
            sp.set_edgecolor('red' if i in drop else 'black')
        ax.set_title(f"C{i+1}" + (" ✗" if i in drop else ""),fontsize=9)
    for i in range(n_show):
        ax=axes[i]
        ts=(comp[:,i]-np.mean(comp[:,i]))/(np.std(comp[:,i])+1e-8)
        ax.plot(taxis,ts,lw=0.8)
        refresh(ax,i)
    plt.suptitle(title+" | Click=drop, Enter=confirm",fontsize=13)
    def on_click(e):
        if e.inaxes in axes:
            i=np.where(axes==e.inaxes)[0][0]
            if i in drop: drop.remove(i)
            else: drop.add(i)
            refresh(e.inaxes,i); fig.canvas.draw_idle()
    def on_key(e):
        if e.key=="enter":
            plt.savefig(path,dpi=180)
            plt.close(fig)
    fig.canvas.mpl_connect('button_press_event',on_click)
    fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
    return sorted(drop)

# --- Decomposition logic ---
if method=="pca":
    pca=PCA(n_components=min(50,T))
    comps=pca.fit_transform(flat.T)
    drop=click_gui(comps,"PCA Components","soner_pipeline/qc/pca_grid.png")
    for d in drop: comps[:,d]=0
    recon=pca.inverse_transform(comps).T.reshape(X,Y,Z,T)

elif method=="tica":
    ica=FastICA(n_components=25,random_state=0,max_iter=500)
    S=ica.fit_transform(flat.T)
    drop=click_gui(S,"Temporal ICA Components","soner_pipeline/qc/tica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@ica.mixing_.T).T.reshape(X,Y,Z,T)

elif method=="sica":
    flat-=flat.mean(axis=1,keepdims=True)
    Xn=StandardScaler(with_mean=False).fit_transform(flat.T)
    ica=FastICA(n_components=25,random_state=0,max_iter=500,whiten='unit-variance')
    S=ica.fit_transform(Xn)
    drop=click_gui(S,"Spatial ICA Components","soner_pipeline/qc/sica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@ica.mixing_.T).T.reshape(X,Y,Z,T)

elif method=="pls":
    Yv=np.zeros((T,1)); Yv[reg_start:reg_end+1]=1
    pls=PLSRegression(n_components=25)
    pls.fit(flat.T,Yv)
    S=pls.x_scores_
    drop=click_gui(S,"PLS Components","soner_pipeline/qc/pls_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@pls.x_loadings_.T).T.reshape(X,Y,Z,T)

outf=f"soner_pipeline/pca_outputs/{sub}/{sub}_{method}_func_denoised_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.nii.gz"
nib.save(nib.Nifti1Image(recon.astype(np.float32),img.affine,img.header),outf)
print("[OK] Saved →",outf)
PY
}

# ====================== MAIN (v17 order preserved) ======================
TR=$(get_tr cleaned_mc_func.nii.gz)
TMP="cleaned_mc_func.nii.gz"

echo -e "${RED}${BOLD}====== BLuSH v18 (corrected, v17 flow) ======${RESET}"

# --- STEP 1: TRIMMING ---------------------------------------------------
echo -e "${PURPLE}${BOLD}STEP 1: TRIMMING${RESET}"
center_line "Trim initial/end volumes for a stable baseline"
read -p "Trim from START (vols, e.g. 0 or 10): " TSTART
read -p "Trim from END   (vols, e.g. 0 or 5):  " TEND
if [[ "${TSTART:-0}" != "0" || "${TEND:-0}" != "0" ]]; then
  END=$(python3 -c "import nibabel as nib; n=nib.load('cleaned_mc_func.nii.gz').shape[-1]; print(n-int('$TEND')-1)")
  OUT="trimmed_${TSTART}-${TEND}_cleaned_mc_func_$(timestamp).nii.gz"
  3dTcat -prefix "$OUT" cleaned_mc_func.nii.gz"[${TSTART}..${END}]"
  TMP="$OUT"
fi
# Save QC on current TMP
save_global_png "$TMP" "$MASK" "$TR" "trim"
save_spectra     "$TMP" "$MASK" "$TR" "trim"

# --- STEP 2: SMOOTHING --------------------------------------------------
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
f,tr=sys.argv[1],float(sys.argv[3]); win=max(1,int(round(float(sys.argv[2])/tr)))
img=nib.load(f); d=img.get_fdata(); sm=uniform_filter1d(d,size=win,axis=-1,mode='nearest')
out="temporal_smoothed_"+f.rsplit('/',1)[-1]
nib.save(nib.Nifti1Image(sm.astype(np.float32),img.affine,img.header),out)
print("[OK] Temporal smoothing →",out)
PY
  TMP="temporal_smoothed_${TMP##*/}"
  save_global_png "$TMP" "$MASK" "$TR" "tsmooth"
  save_spectra     "$TMP" "$MASK" "$TR" "tsmooth"
fi
if [[ "${SSM:-0}" != "0" ]]; then
  OUT="spatial_smoothed_${TMP%.nii.gz}_$(timestamp).nii.gz"
  fslmaths "$TMP" -s "$SSM" "$OUT"
  TMP="$OUT"
  save_global_png "$TMP" "$MASK" "$TR" "ssmooth"
  save_spectra     "$TMP" "$MASK" "$TR" "ssmooth"
fi

# --- STEP 2.5: DETRENDING ----------------------------------------------
echo -e "${PURPLE}${BOLD}STEP 2.5: DETRENDING${RESET}"
center_line "Remove scanner drifts with polynomial fits"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz temporal_smoothed_*.nii.gz spatial_smoothed_*.nii.gz spatial_smoothed_*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file for detrending [1-${#AVAILABLE[@]}]: " idx
TMP="${AVAILABLE[$((idx-1))]}"
plot_global_signal "$TMP" "$MASK" "$TR"
read -p "Apply detrending? (y/n): " DETR
if [[ "$DETR" == "y" ]]; then
  while true; do
    echo "1) Linear      – remove constant + linear drift"
    echo "2) Quadratic   – remove slow parabolic drift"
    echo "3) Polynomial  – higher-order baseline removal"
    echo "4) Skip        – keep as-is"
    read -p "Select detrend type [1–4]: " DTYPE
    case "$DTYPE" in
      1) POL=1;;
      2) POL=2;;
      3) read -p "Enter polynomial order (≥3): " POL;;
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
    if [[ "$OK" == "y" ]]; then TMP="$OUT"; break; fi
  done
fi

# --- STEP 3: SCM QC -----------------------------------------------------
echo -e "${PURPLE}${BOLD}STEP 3: SCM QC${RESET}"
center_line "Build baseline/signal means and %change maps (QC)"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz spatial_smoothed_*.nii.gz temporal_smoothed_*.nii.gz soner_pipeline/detrended/*.nii.gz soner_pipeline/filtering/*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file for SCM QC [1-${#AVAILABLE[@]}]: " idx
SCM_FILE="${AVAILABLE[$((idx-1))]}"
plot_global_signal "$SCM_FILE" "$MASK" "$TR"
read -p "Baseline vols (start end): " B1 B2
read -p "Signal   vols (start end): " S1 S2
OUT="soner_pipeline/scm_outputs/$(basename "$SCM_FILE" .nii.gz)_B${B1}-${B2}_S${S1}-${S2}_$(timestamp)"
scm_blush_fast "$SCM_FILE" "$B1" "$B2" "$S1" "$S2" "$MASK" "$OUT"
# Correct FSLeyes order: norm_func (bottom) → baseline → signal_change_map (top)
nohup fsleyes "$OUT/norm_func_${B1}_${B2}.nii.gz" \
              "$OUT/baseline_${B1}_${B2}.nii.gz" \
              "$OUT/signal_change_map_${B1}_${B2}_${S1}_${S2}.nii.gz" >/dev/null 2>&1 &

# --- STEP 4: FILTERING --------------------------------------------------
echo -e "${PURPLE}${BOLD}STEP 4: FILTERING${RESET}"
center_line "High-pass removes drifts; Low-pass keeps slow oscillations"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz spatial_smoothed_*.nii.gz temporal_smoothed_*.nii.gz soner_pipeline/detrended/*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
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
    1) read -p "Cut-off Hz (e.g. 0.01): " FC; OUT="soner_pipeline/filtering/hp_${FC}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FC" "0" "$TR" "1"; TMP="$OUT";;
    2) read -p "Cut-off Hz (e.g. 0.1): " FC; OUT="soner_pipeline/filtering/lp_${FC}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "0" "$FC" "$TR" "2"; TMP="$OUT";;
    3) read -p "Band range Hz (e.g. 0.01 0.1): " FL FH; OUT="soner_pipeline/filtering/bp_${FL}-${FH}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FL" "$FH" "$TR" "3"; TMP="$OUT";;
    *) continue;;
  esac
  # Save QC for the filtered output
  save_global_png "$TMP" "$MASK" "$TR" "filter"
  save_spectra     "$TMP" "$MASK" "$TR" "filter"
  # Optional despike (as in v17 enhancement)
  read -p "Run AFNI despiking on this output? (y/n): " DESP
  if [[ "$DESP" == "y" ]]; then
    DSP="soner_pipeline/filtering/despiked_${stamp}.nii.gz"
    3dDespike -NEW -localedit -ignore 5 -prefix "$DSP" "$TMP" || true
    [[ -f "$DSP" ]] && TMP="$DSP" && save_global_png "$TMP" "$MASK" "$TR" "despiked" && save_spectra "$TMP" "$MASK" "$TR" "despiked"
  fi
done

# --- STEP 5: DECOMPOSITION (interactive + live global plot) ---
while true; do
  echo -e "${PURPLE}${BOLD}STEP 5: DECOMPOSITION${RESET}"
  center_line "Decompose signal via PCA/ICA/PLS; drop unwanted components"
  echo "1) PCA (variance-based denoising)"
  echo "2) Temporal ICA (independent time courses)"
  echo "3) Spatial ICA (independent spatial maps)"
  echo "4) PLS (stimulus-correlated components)"
  echo "5) Exit"
  read -p "Choice [1–5]: " CH
  [[ "$CH" == "5" ]] && break
  case "$CH" in
    1) METHOD="pca";;
    2) METHOD="tica";;
    3) METHOD="sica";;
    4) METHOD="pls";;
    *) continue;;
  esac

  AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz spatial_smoothed_*.nii.gz temporal_smoothed_*.nii.gz soner_pipeline/detrended/*.nii.gz soner_pipeline/filtering/*.nii.gz 2>/dev/null))
  for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
  read -p "Select file for ${METHOD^^} [1-${#AVAILABLE[@]}]: " idx
  TMP="${AVAILABLE[$((idx-1))]}"

  if [[ "$METHOD" == "pls" ]]; then
    read -p "PLS regressor start/end vols (e.g. 1300 1500): " DS1 DS2
  else
    DS1=0; DS2=0
  fi

  SUB="raw"
  [[ "$TMP" == *"hp_"* ]] && SUB="hp"
  [[ "$TMP" == *"lp_"* ]] && SUB="lp"
  [[ "$TMP" == *"bp_"* ]] && SUB="bp"
  [[ "$TMP" == *"detrended"* ]] && SUB="detrended"
  [[ "$TMP" == trimmed_* ]] && SUB="trim"
  [[ "$TMP" == spatial_smoothed_* || "$TMP" == temporal_smoothed_* ]] && SUB="smoothed"

  # QC before decomp
  save_global_png "$TMP" "$MASK" "$TR" "${METHOD}_before"
  save_spectra     "$TMP" "$MASK" "$TR" "${METHOD}_before"

  run_decomposition "$TMP" "$METHOD" "$SUB" "$DS1" "$DS2"
  DEN=$(ls -t "soner_pipeline/pca_outputs/${SUB}/${SUB}_${METHOD}_func_denoised_"*.nii.gz 2>/dev/null | head -n1)
  [[ -f "$DEN" ]] || { echo "[ERR] No denoised file found."; continue; }

  # --- NEW: show global mean interactively after decomposition ---
  plot_global_signal "$DEN" "$MASK" "$TR"
  save_global_png "$DEN" "$MASK" "$TR" "${METHOD}_after"
  save_spectra     "$DEN" "$MASK" "$TR" "${METHOD}_after"

  echo -e "\n${BLUE}${BOLD}Define baseline and signal windows for denoised global mean:${RESET}"
  read -p "Baseline vols (start end): " XB1 XB2
  read -p "Signal vols (start end): " XS1 XS2

  OUTD="soner_pipeline/scm_outputs/${SUB}_${METHOD}_B${XB1}-${XB2}_S${XS1}-${XS2}_$(timestamp)"
  scm_blush_fast "$DEN" "$XB1" "$XB2" "$XS1" "$XS2" "$MASK" "$OUTD"

  # Correct FSLeyes order: norm_func (bottom) → baseline → signal_change_map (top)
  nohup fsleyes "$OUTD/norm_func_${XB1}_${XB2}.nii.gz" \
                "$OUTD/baseline_${XB1}_${XB2}.nii.gz" \
                "$OUTD/signal_change_map_${XB1}_${XB2}_${XS1}_${XS2}.nii.gz" >/dev/null 2>&1 &

  read -p "Run another decomposition? (y/n): " AGAIN
  [[ "$AGAIN" != "y" ]] && break
done
# --- HTML QC Summary (lists all QC PNGs) -------------------------------
echo -e "${BLUE}${BOLD}Generating HTML QC summary…${RESET}"
python3 <<'PY'
import os,glob,datetime,html
ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
pngs=sorted(glob.glob('soner_pipeline/qc/*.png'))
rows='\n'.join(f"<tr><td>{html.escape(os.path.basename(p))}</td><td><img src='{p}' width='360'></td></tr>" for p in pngs)
doc=f"""<html><head><meta charset='utf-8'><title>BLuSH QC Summary</title>
<style>body{{background:#0f1116;color:#eaeef2;font-family:Arial;padding:24px}}
h1{{color:#8bd3ff}} table{{border-collapse:collapse}} td,th{{border:1px solid #2a2f3a;padding:6px}}</style>
</head><body><h1>BLuSH QC Summary</h1><p><i>{ts}</i></p>
<table>{rows}</table></body></html>"""
out=f"soner_pipeline/qc_summary_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
open(out,'w',encoding='utf-8').write(doc)
print("[OK] Summary HTML →",out)
try:
    import platform,subprocess
    if platform.system()=="Darwin": subprocess.Popen(["open",out])
    elif platform.system()=="Linux": subprocess.Popen(["xdg-open",out])
    elif platform.system()=="Windows": subprocess.Popen(["start",out],shell=True)
except Exception as e:
    print("[WARN] Could not auto-open:",e)
PY

echo -e "\n${GREEN}${BOLD}All steps complete — BLuSH v18 (corrected) finished.${RESET}"
