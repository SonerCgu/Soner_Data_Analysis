#!/bin/bash
# =====================================================================
#  filtering_decomposition_v12_final.sh
#
#  BLuSH-style preclinical fMRI preprocessing + decomposition pipeline
#  Based on v12, adds bold, centered warm-orange explanations and
#  replaces "Cubic" detrending with generic "Polynomial (user-defined)".
#
#  Author: Naman Jain (MPI for Biological Cybernetics)
#  Date:   November 2025
# =====================================================================

set -Eeuo pipefail
trap 'echo "[ERROR] line ${LINENO}" >&2' ERR
shopt -s nullglob
export LC_NUMERIC=C
export PYTHONWARNINGS="ignore"

# ---------------------------------------------------------------
#  COLORS
# ---------------------------------------------------------------
RED="\033[38;5;160m"
BLUE="\033[38;5;33m"
GREEN="\033[92m"
YELLOW="\033[93m"
ORANGE_BOLD="\033[1;38;5;208m"
RESET="\033[0m"

# ---------------------------------------------------------------
#  FOLDER INITIALIZATION
# ---------------------------------------------------------------
mkdir -p soner_pipeline/{filtering,scm_outputs,logs,detrended,qc}
mkdir -p soner_pipeline/pca_outputs/{raw,hp,lp,bp,despiked}

MASK="mask_mean_mc_func.nii.gz"
[[ ! -f cleaned_mc_func.nii.gz ]] && { echo "[ERR] cleaned_mc_func.nii.gz missing"; exit 1; }
[[ ! -f "$MASK" ]] && MASK=""

timestamp(){ date +"%Y%m%d_%H%M%S"; }
LOG="soner_pipeline/logs/blush_run_$(timestamp).csv"
echo "Step,Parameter,Value" > "$LOG"

# ---------------------------------------------------------------
#  UTILITIES
# ---------------------------------------------------------------
get_tr(){ python3 - "$1" <<'PY'
import sys,nibabel as nib
try: print(float(nib.load(sys.argv[1]).header.get_zooms()[3]))
except: print(1.0)
PY
}

plot_global_signal(){ python3 - "$1" "$2" "$3" <<'PY'
import sys,nibabel as nib,numpy as np,matplotlib.pyplot as plt,os
f,m,tr=sys.argv[1],sys.argv[2],float(sys.argv[3])
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); x=np.arange(len(sig))*tr
plt.figure(figsize=(10,4)); plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean Signal – {os.path.basename(f)}")
plt.xlabel("Time (s)"); plt.ylabel("Mean Intensity"); plt.grid(True,alpha=0.3)
plt.tight_layout(); plt.show()
PY
}

# --- Safe SCM ---------------------------------------------------
scm_blush_fast(){
local func="$1" base="$2" bend="$3" sig="$4" send="$5" mask="$6" outd="$7"
mkdir -p "$outd"
nvol=$(3dinfo -nv "$func" 2>/dev/null || echo 0)
[[ "$nvol" == "0" ]] && { echo "[ERR] Cannot read $func"; return; }

# Clamp indices
for v in base bend sig send; do
  eval val=\$$v
  (( val < 0 )) && eval $v=0
  (( val >= nvol )) && eval $v=$((nvol-1))
done
(( bend <= base )) && bend=$((base+1))
(( send <= sig )) && send=$((sig+1))
echo "[INFO] Using baseline ${base}-${bend}, signal ${sig}-${send} (max=$((nvol-1)))"

3dTstat -mean -prefix "$outd/baseline_${base}_${bend}.nii.gz" "$func[${base}..${bend}]"
3dTstat -mean -prefix "$outd/signal_${sig}_${send}.nii.gz" "$func[${sig}..${send}]"
fslmaths "$outd/signal_${sig}_${send}" -sub "$outd/baseline_${base}_${bend}" \
 -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
 "$outd/signal_change_map_${base}_${bend}_${sig}_${send}.nii.gz"
fslmaths "$func" -sub "$outd/baseline_${base}_${bend}" \
 -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
 "$outd/norm_func_${base}_${bend}.nii.gz"
}

# --- Frequency filter helper -----------------------------------
apply_freq_filter(){
inp="$1"; outp="$2"; fl="$3"; fh="$4"; tr="$5"; ftype="$6"
mkdir -p soner_pipeline/filtering
stamp=$(timestamp)
case "$ftype" in
  1) HP=$(python3 -c "print(1/(2*$tr*$fl))"); fslmaths "$inp" -bptf $HP -1 "$outp"; TYPE="High-pass (${fl} Hz removes slow drifts)";;
  2) LP=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf -1 $LP "$outp"; TYPE="Low-pass (${fh} Hz removes fast noise)";;
  3) HP=$(python3 -c "print(1/(2*$tr*$fl))"); LP=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf $HP $LP "$outp"; TYPE="Band-pass (${fl}-${fh} Hz isolates target band)";;
esac
echo "Filtering,$TYPE,$outp" >> "$LOG"
plot_global_signal "$outp" "$MASK" "$tr"
read -p "Run AFNI despiking? (y/n): " DESP
if [[ "$DESP" == "y" ]]; then
  DSP="soner_pipeline/filtering/despiked_${stamp}.nii.gz"
  3dDespike -NEW -localedit -ignore 5 -prefix "$DSP" "$outp" || true
  [[ -f "$DSP" ]] && outp="$DSP"
  echo "Despiking,Enabled,$DSP" >> "$LOG"
  plot_global_signal "$outp" "$MASK" "$tr"
fi
}

# --- Decomposition helper --------------------------------------
run_decomposition(){
local infile="$1"; local method="$2"; local subdir="$3"; local reg_start="$4"; local reg_end="$5"
python3 - "$infile" "$method" "$subdir" "$reg_start" "$reg_end" <<'PY'
import sys,os,nibabel as nib,numpy as np,matplotlib.pyplot as plt
from sklearn.decomposition import PCA,FastICA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
fname,method,sub,reg_start,reg_end=sys.argv[1:6]
img=nib.load(fname); data=img.get_fdata().astype(np.float32)
X,Y,Z,T=data.shape; flat=data.reshape(-1,T)
os.makedirs(f"soner_pipeline/pca_outputs/{sub}",exist_ok=True)
reg_start,reg_end=int(reg_start),int(reg_end)

def click_gui(comp,title,path):
    n_show=min(25,comp.shape[1]); taxis=np.arange(comp.shape[0]); drop=set()
    fig,axes=plt.subplots(int(np.ceil(n_show/5)),5,figsize=(18,9)); axes=axes.ravel()
    def refresh(ax,i):
        ax.set_title(f"C{i+1}"+(" (DROP)" if i in drop else ""),fontsize=9)
        for sp in ax.spines.values(): sp.set_edgecolor('red' if i in drop else 'black')
    for i in range(n_show):
        ax=axes[i]; ts=(comp[:,i]-np.mean(comp[:,i]))/(np.std(comp[:,i])+1e-8)
        ax.plot(taxis,ts,lw=0.8); refresh(ax,i)
    plt.suptitle(title+" | Click=drop | ENTER=confirm",fontsize=13)
    def on_click(e):
        if e.inaxes in axes:
            i=np.where(axes==e.inaxes)[0][0]
            if i in drop: drop.remove(i)
            else: drop.add(i)
            refresh(e.inaxes,i); fig.canvas.draw_idle()
    def on_key(e):
        if e.key=="enter": plt.savefig(path,dpi=200); plt.close(fig)
    fig.canvas.mpl_connect('button_press_event',on_click)
    fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show(); return sorted(drop)

print(f"[INFO] Running {method.upper()} on {data.shape}")

if method=="pca":
    pca=PCA(n_components=min(50,T)); pcs=pca.fit_transform(flat.T)
    drop=click_gui(pcs,"PCA","soner_pipeline/qc/pca_grid.png")
    for d in drop: pcs[:,d]=0
    recon=pca.inverse_transform(pcs).T.reshape(X,Y,Z,T)
elif method=="tica":
    ica=FastICA(n_components=25,random_state=0,max_iter=500)
    S=ica.fit_transform(flat.T); A=ica.mixing_
    drop=click_gui(S,"tICA","soner_pipeline/qc/tica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@A.T).T.reshape(X,Y,Z,T)
elif method=="sica":
    flat-=flat.mean(axis=1,keepdims=True)
    Xn=StandardScaler(with_mean=False).fit_transform(flat.T)
    ica=FastICA(n_components=25,random_state=0,max_iter=500,whiten='unit-variance')
    S=ica.fit_transform(Xn); A=ica.mixing_.T
    drop=click_gui(S,"sICA","soner_pipeline/qc/sica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@A).T.reshape(X,Y,Z,T)
elif method=="pls":
    Yv=np.zeros((T,1)); Yv[reg_start:reg_end+1]=1
    pls=PLSRegression(n_components=25)
    pls.fit(flat.T,Yv); S=pls.x_scores_
    drop=click_gui(S,"PLS","soner_pipeline/qc/pls_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@pls.x_loadings_.T).T.reshape(X,Y,Z,T)

outf=f"soner_pipeline/pca_outputs/{sub}/{sub}_{method}_func_denoised.nii.gz"
nib.save(nib.Nifti1Image(recon.astype(np.float32),img.affine,img.header),outf)
print(f"[OK] Saved → {outf}")
PY
}
# ===============================================================
#  MAIN PIPELINE
# ===============================================================
TR=$(get_tr cleaned_mc_func.nii.gz)
TMP="cleaned_mc_func.nii.gz"
echo -e "${RED}========== filtering_decomposition_v12_final ==========${RESET}"

# --- STEP 1: TRIMMING ------------------------------------------
echo -e "${RED}STEP 1: TRIMMING${RESET}"
echo -e "${ORANGE_BOLD}      Remove unstable initial/final volumes. Improves baseline stability.      ${RESET}"
read -p "Trim from START (vols, e.g. 10): " TSTART
read -p "Trim from END (vols, e.g. 5): " TEND
if [[ "$TSTART" != "0" ]] || [[ "$TEND" != "0" ]]; then
  END=$(python3 -c "import nibabel as nib; n=nib.load('cleaned_mc_func.nii.gz').shape[-1]; print(n-int('$TEND')-1)")
  3dTcat -prefix trimmed_cleaned_mc_func.nii.gz cleaned_mc_func.nii.gz"[${TSTART}..${END}]"
  TMP="trimmed_cleaned_mc_func.nii.gz"
fi

# --- STEP 2: SMOOTHING -----------------------------------------
echo -e "${RED}STEP 2: SMOOTHING${RESET}"
echo -e "${ORANGE_BOLD}      Reduce noise. Temporal = time smooth; spatial = voxel blur.      ${RESET}"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file to smooth [1-${#AVAILABLE[@]}]: " idx; TMP="${AVAILABLE[$((idx-1))]}"
read -p "Temporal smoothing (sec, 0 = skip): " TSM
read -p "Spatial σ blur (e.g. 0.2, 0 = skip): " SSM
if [[ "$TSM" != "0" ]]; then
  python3 - <<PY "$TMP" "$TSM" "$TR"
import sys,nibabel as nib,numpy as np; from scipy.ndimage import uniform_filter1d
f,tr=sys.argv[1],float(sys.argv[3]); win=max(1,int(round(float(sys.argv[2])/tr)))
img=nib.load(f); d=img.get_fdata(); sm=uniform_filter1d(d,size=win,axis=-1,mode='nearest')
out="smoothed_"+f; nib.save(nib.Nifti1Image(sm.astype(np.float32),img.affine,img.header),out)
print(f"[OK] Smoothed → {out}")
PY
  TMP="smoothed_${TMP}"
fi
if [[ "$SSM" != "0" ]]; then fslmaths "$TMP" -s "$SSM" "spatial_${TMP}"; TMP="spatial_${TMP}"; fi

# --- STEP 2.5: DETRENDING --------------------------------------
echo -e "${RED}STEP 2.5: DETRENDING${RESET}"
echo -e "${ORANGE_BOLD}      Remove scanner drift using AFNI polynomial detrending.      ${RESET}"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz smoothed_*.nii.gz spatial_*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file for detrending [1-${#AVAILABLE[@]}]: " idx; TMP="${AVAILABLE[$((idx-1))]}"
plot_global_signal "$TMP" "$MASK" "$TR"
read -p "Apply detrending? (y/n): " DETR
if [[ "$DETR" == "y" ]]; then
  while true; do
    echo "1) Linear (remove drift)"
    echo "2) Quadratic (correct curvature)"
    echo "3) Polynomial (custom order)"
    echo "4) Skip"
    read -p "Select detrend type [1–4]: " DTYPE
    case "$DTYPE" in
      1) POL=1 ;;
      2) POL=2 ;;
      3) read -p "Enter polynomial order (e.g. 4, 6, 10): " POL ;;
      4) break ;;
    esac
    BASE=$(basename "$TMP" .nii.gz)
    OUT="soner_pipeline/detrended/detrended_${BASE}_pol${POL}_$(timestamp).nii.gz"
    3dDetrend -polort $POL -prefix "$OUT" "$TMP"
    [[ -f "$OUT" ]] && echo "[OK] Saved → $OUT"
    plot_global_signal "$OUT" "$MASK" "$TR"
    read -p "Satisfied? (y/n): " OK; [[ "$OK" == "y" ]] && TMP="$OUT" && break
  done
fi

# --- STEP 3: SCM QC --------------------------------------------
echo -e "${RED}STEP 3: SCM QC${RESET}"
echo -e "${ORANGE_BOLD}      Compute percent-signal-change maps for visual QC.      ${RESET}"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz smoothed_*.nii.gz spatial_*.nii.gz soner_pipeline/detrended/*.nii.gz soner_pipeline/filtering/*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file for SCM QC [1-${#AVAILABLE[@]}]: " idx; SCM_FILE="${AVAILABLE[$((idx-1))]}"
plot_global_signal "$SCM_FILE" "$MASK" "$TR"
read -p "Baseline (start end): " B1 B2; read -p "Signal (start end): " S1 S2
(( S2 <= S1 )) && { echo "[WARN] Swapping signal indices"; tmp=$S1; S1=$S2; S2=$tmp; }
OUT="soner_pipeline/scm_outputs/$(basename "$SCM_FILE" .nii.gz)_B${B1}-${B2}_S${S1}-${S2}"
scm_blush_fast "$SCM_FILE" "$B1" "$B2" "$S1" "$S2" "$MASK" "$OUT"
nohup fsleyes "$OUT/norm_func_${B1}_${B2}.nii.gz" "$OUT/baseline_${B1}_${B2}.nii.gz" "$OUT/signal_change_map_${B1}_${B2}_${S1}_${S2}.nii.gz" >/dev/null 2>&1 &

# --- STEP 4: FILTERING ------------------------------------------
echo -e "${RED}STEP 4: FILTERING${RESET}"
echo -e "${ORANGE_BOLD}      Isolate desired frequency band (HP, LP, BP).      ${RESET}"
AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz smoothed_*.nii.gz spatial_*.nii.gz soner_pipeline/detrended/*.nii.gz soner_pipeline/filtering/*.nii.gz 2>/dev/null))
for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
read -p "Select file to filter [1-${#AVAILABLE[@]}]: " idx; TMP="${AVAILABLE[$((idx-1))]}"
while true; do
  echo "1) High-pass"
  echo "2) Low-pass"
  echo "3) Band-pass"
  echo "4) Continue → Decomposition"
  read -p "Choose filter type [1–4]: " FTYPE; [[ "$FTYPE" == "4" ]] && break
  stamp=$(timestamp)
  case "$FTYPE" in
    1)read -p "Cut-off Hz (e.g. 0.01): " FC; OUT="soner_pipeline/filtering/hp_${FC}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FC" "0" "$TR" "1"; TMP="$OUT";;
    2)read -p "Cut-off Hz (e.g. 0.1): " FC; OUT="soner_pipeline/filtering/lp_${FC}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "0" "$FC" "$TR" "2"; TMP="$OUT";;
    3)read -p "Band range Hz (e.g. 0.01 0.1): " FL FH; OUT="soner_pipeline/filtering/bp_${FL}-${FH}Hz_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FL" "$FH" "$TR" "3"; TMP="$OUT";;
  esac
done

# --- STEP 5: DECOMPOSITION --------------------------------------
echo -e "${RED}STEP 5: DECOMPOSITION${RESET}"
echo -e "${ORANGE_BOLD}      Extract independent / orthogonal components for denoising.      ${RESET}"
while true; do
  echo "--------------------------------------------------------------"
  echo "1) PCA  –  Remove correlated noise via orthogonal axes"
  echo "2) Temporal ICA  –  Separate independent time courses"
  echo "3) Spatial ICA   –  Separate independent spatial maps"
  echo "4) PLS  –  Regressor-based components linked to stimulus"
  echo "5) End  –  Finish and generate summary"
  echo "--------------------------------------------------------------"
  read -p "Choice [1–5]: " CH
  [[ "$CH" == "5" ]] && break

  case "$CH" in
    1) METHOD="pca";;
    2) METHOD="tica";;
    3) METHOD="sica";;
    4) METHOD="pls";;
    *) echo "[WARN] Invalid choice."; continue;;
  esac

  echo -e "${YELLOW}Tip:${RESET} PLS needs injection window; others auto-learn."
  AVAILABLE=($(ls -tr cleaned_mc_func.nii.gz trimmed_*.nii.gz smoothed_*.nii.gz spatial_*.nii.gz soner_pipeline/detrended/*.nii.gz soner_pipeline/filtering/*.nii.gz 2>/dev/null))
  for i in "${!AVAILABLE[@]}"; do printf "%3d) %s\n" $((i+1)) "${AVAILABLE[$i]}"; done
  read -p "Select input file for ${METHOD^^}: " idx
  TMP="${AVAILABLE[$((idx-1))]}"

  SRC=$(basename "$TMP")
  if [[ "$SRC" == *"hp_"* ]]; then SUB="hp"
  elif [[ "$SRC" == *"lp_"* ]]; then SUB="lp"
  elif [[ "$SRC" == *"bp_"* ]]; then SUB="bp"
  elif [[ "$SRC" == *"detrended"* ]]; then SUB="detrended"
  elif [[ "$SRC" == *"trimmed"* ]]; then SUB="trimmed"
  else SUB="raw"; fi

  if [[ "$METHOD" == "pls" ]]; then
    read -p "Signal (injection) start/end (e.g. 1500 2000): " DS1 DS2
  else
    DS1=0; DS2=0
  fi

  echo -e "${BLUE}[INFO] Running ${METHOD^^} on $TMP → subdir=${SUB}${RESET}"
  stamp=$(timestamp)
  save_png="soner_pipeline/qc/global_before_${METHOD}_${stamp}.png"
  python3 - "$TMP" "$MASK" "$TR" "$save_png" <<'PY'
import sys,nibabel as nib,numpy as np,matplotlib.pyplot as plt,os
f,m,tr,out=sys.argv[1:5]; tr=float(tr)
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); x=np.arange(len(sig))*tr
plt.figure(figsize=(8,4)); plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean Before {os.path.basename(f)}"); plt.xlabel("Time (s)")
plt.ylabel("Mean Intensity"); plt.grid(True,alpha=0.3)
plt.tight_layout(); plt.savefig(out,dpi=150); plt.close()
PY

  run_decomposition "$TMP" "$METHOD" "$SUB" "$DS1" "$DS2"

  DEN="soner_pipeline/pca_outputs/${SUB}/${SUB}_${METHOD}_func_denoised.nii.gz"
  if [[ -f "$DEN" ]]; then
    echo -e "${GREEN}[OK] Denoised output saved → $DEN${RESET}"
    save_png2="soner_pipeline/qc/global_after_${METHOD}_${stamp}.png"
    python3 - "$DEN" "$MASK" "$TR" "$save_png2" <<'PY'
import sys,nibabel as nib,numpy as np,matplotlib.pyplot as plt,os
f,m,tr,out=sys.argv[1:5]; tr=float(tr)
img=nib.load(f); d=img.get_fdata()
try: mask=nib.load(m).get_fdata()>0; d=d[mask]
except: d=d.reshape(-1,d.shape[-1])
sig=d.mean(0); x=np.arange(len(sig))*tr
plt.figure(figsize=(8,4)); plt.plot(x,sig,'k',lw=1)
plt.title(f"Global Mean After {os.path.basename(f)}"); plt.xlabel("Time (s)")
plt.ylabel("Mean Intensity"); plt.grid(True,alpha=0.3)
plt.tight_layout(); plt.savefig(out,dpi=150); plt.close()
PY

    echo -e "${YELLOW}Verify baseline/signal of denoised data...${RESET}"
    plot_global_signal "$DEN" "$MASK" "$TR"
    read -p "Baseline start/end: " XB1 XB2
    read -p "Signal start/end: " XS1 XS2
    (( XS2 <= XS1 )) && { tmp=$XS1; XS1=$XS2; XS2=$tmp; }
    OUTD="soner_pipeline/scm_outputs/${SUB}_${METHOD}_B${XB1}-${XB2}_S${XS1}-${XS2}"
    scm_blush_fast "$DEN" "$XB1" "$XB2" "$XS1" "$XS2" "$MASK" "$OUTD"
    nohup fsleyes "$OUTD/norm_func_${XB1}_${XB2}.nii.gz" \
                  "$OUTD/baseline_${XB1}_${XB2}.nii.gz" \
                  "$OUTD/signal_change_map_${XB1}_${XB2}_${XS1}_${XS2}.nii.gz" >/dev/null 2>&1 &
  fi
done

# --- HTML SUMMARY -----------------------------------------------
python3 <<'PY' "$LOG"
import subprocess, sys, os, datetime
try: import pandas as pd
except ImportError: subprocess.run([sys.executable,"-m","pip","install","pandas","--quiet"])
import pandas as pd
log=sys.argv[1]; df=pd.read_csv(log)
ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
html=f"<html><head><title>BLuSH Summary</title><style>body{{font-family:Arial;background:#111;color:#eee;padding:2em}}h2{{color:#8ef}}td,th{{padding:4px;border:1px solid #333}}</style></head><body><h1>BLuSH Summary</h1><p><i>{ts}</i></p>{df.to_html(index=False)}</body></html>"
out=log.replace('blush_run_','blush_summary_').replace('.csv','.html')
open(out,'w').write(html)
print(f"[OK] Summary HTML → {out}")
os.system(f"open '{out}' || xdg-open '{out}'" if os.name=='posix' else f"start {out}")
PY

echo -e "${GREEN}All steps complete — BLuSH v12_final finished successfully.${RESET}"
