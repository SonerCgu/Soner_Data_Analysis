#!/bin/bash
# =====================================================================
#  filtering_decomposition_v4.sh
#  BLuSH-style fMRI preprocessing + decomposition + HTML summary
#
#  Author: Soner Caner Cagun + ChatGPT co-dev
#  Version: v4 (2025-11-07)
# =====================================================================

set -Eeuo pipefail
trap 'echo "[ERROR] line ${LINENO}" >&2' ERR
export LC_NUMERIC=C
export PYTHONWARNINGS="ignore"

# ---------------------------------------------------------------
# --- Terminal Color Palette
# ---------------------------------------------------------------
DARKRED="\033[38;5;160m"
DARKBLUE="\033[38;5;33m"
GREEN="\033[92m"
YELLOW="\033[93m"
RESET="\033[0m"

# ---------------------------------------------------------------
# --- Folder Initialization
# ---------------------------------------------------------------
mkdir -p soner_pipeline/{filtering,scm_outputs,logs}
mkdir -p soner_pipeline/pca_outputs/{raw,hp,lp,bp,despiked}

MASK="mask_mean_mc_func.nii.gz"
[[ ! -f cleaned_mc_func.nii.gz ]] && { echo "[ERR] cleaned_mc_func.nii.gz missing"; exit 1; }
[[ ! -f "$MASK" ]] && MASK=""

timestamp() { date +"%Y%m%d_%H%M%S"; }
LOG="soner_pipeline/logs/blush_run_$(timestamp).csv"
echo "Step,Parameter,Value" > "$LOG"

# ===============================================================
#  UTILITIES
# ===============================================================
get_tr() { python3 - "$1" <<'PY'
import sys,nibabel as nib
try: print(float(nib.load(sys.argv[1]).header.get_zooms()[3]))
except: print(1.0)
PY
}

plot_global_signal() { python3 - "$1" "$2" "$3" <<'PY'
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

scm_blush_fast() {
# Static Contrast Map: compute baseline vs. signal mean maps and %
# change to visualize response amplitude.
local func="$1" base="$2" bend="$3" sig="$4" send="$5" mask="$6" outd="$7"
mkdir -p "$outd"
3dTstat -mean -prefix "$outd/baseline_${base}_${bend}.nii.gz" "$func[${base}..${bend}]"
3dTstat -mean -prefix "$outd/signal_${sig}_${send}.nii.gz" "$func[${sig}..${send}]"
fslmaths "$outd/signal_${sig}_${send}" -sub "$outd/baseline_${base}_${bend}" \
 -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
 "$outd/signal_change_map_${base}_${bend}_${sig}_${send}.nii.gz"
fslmaths "$func" -sub "$outd/baseline_${base}_${bend}" \
 -div "$outd/baseline_${base}_${bend}" -mul 100 ${mask:+-mas "$mask"} \
 "$outd/norm_func_${base}_${bend}.nii.gz"
}

# ===============================================================
#  FILTERING + DESPIKING
# ===============================================================
apply_freq_filter() {
inp="$1"; outp="$2"; fl="$3"; fh="$4"; tr="$5"; ftype="$6"
mkdir -p soner_pipeline/filtering; stamp=$(timestamp)
case "$ftype" in
  1) HP_SIG=$(python3 -c "print(1/(2*$tr*$fl))"); fslmaths "$inp" -bptf $HP_SIG -1 "$outp"; TYPE="High-pass (${fl} Hz)";;
  2) LP_SIG=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf -1 $LP_SIG "$outp"; TYPE="Low-pass (${fh} Hz)";;
  3) HP_SIG=$(python3 -c "print(1/(2*$tr*$fl))"); LP_SIG=$(python3 -c "print(1/(2*$tr*$fh))"); fslmaths "$inp" -bptf $HP_SIG $LP_SIG "$outp"; TYPE="Band-pass (${fl}–${fh} Hz)";;
esac
echo "Filtering,$TYPE,$outp" >> "$LOG"
plot_global_signal "$outp" "$MASK" "$tr"

read -p "Run AFNI despiking? (y/n): " DESP
if [[ "$DESP" == "y" ]]; then
  DSP_OUT="soner_pipeline/filtering/despiked_${stamp}.nii.gz"
  3dDespike -NEW -localedit -ignore 5 -prefix "$DSP_OUT" "$outp" || true
  [[ -f "$DSP_OUT" ]] && outp="$DSP_OUT"
  echo "Despiking,Enabled,$DSP_OUT" >> "$LOG"
else
  echo "Despiking,Skipped," >> "$LOG"
fi
plot_global_signal "$outp" "$MASK" "$tr"
}

# ===============================================================
#  DECOMPOSITION (PCA / tICA / sICA / PLS)
# ===============================================================
run_decomposition() {
local infile="$1"; local method="$2"; local subdir="$3"
python3 - "$infile" "$method" "$subdir" <<'PY'
import sys,os,nibabel as nib,numpy as np,matplotlib.pyplot as plt
from sklearn.decomposition import PCA,FastICA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
fname,method,sub=sys.argv[1:4]
img=nib.load(fname); data=img.get_fdata().astype(np.float32)
X,Y,Z,T=data.shape; flat=data.reshape(-1,T)
qc_dir=f"soner_pipeline/pca_outputs/{sub}/qc"; os.makedirs(qc_dir,exist_ok=True)

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
    drop=click_gui(pcs,"PCA",f"{qc_dir}/pca_grid.png")
    for d in drop: pcs[:,d]=0
    recon=pca.inverse_transform(pcs).T.reshape(X,Y,Z,T)
elif method=="tica":
    ica=FastICA(n_components=25,random_state=0,max_iter=500)
    S=ica.fit_transform(flat.T); A=ica.mixing_
    drop=click_gui(S,"tICA",f"{qc_dir}/tica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@A.T).T.reshape(X,Y,Z,T)
elif method=="sica":
    flat-=flat.mean(axis=1,keepdims=True)
    Xn=StandardScaler(with_mean=False).fit_transform(flat.T)
    ica=FastICA(n_components=25,random_state=0,max_iter=500,whiten='unit-variance')
    S=ica.fit_transform(Xn); A=ica.mixing_.T
    drop=click_gui(S,"sICA",f"{qc_dir}/sica_grid.png")
    for d in drop: S[:,d]=0
    recon=(S@A).T.reshape(X,Y,Z,T)
elif method=="pls":
    Yv=np.linspace(0,1,T).reshape(-1,1)
    pls=PLSRegression(n_components=25)
    pls.fit(flat.T,Yv); S=pls.x_scores_
    drop=click_gui(S,"PLS",f"{qc_dir}/pls_grid.png")
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
echo -e "${DARKRED}========== filtering_decomposition_v4 ==========${RESET}"
echo "TR,$TR," >> "$LOG"

# --- STEP 1: Trimming ---
echo -e "${DARKRED}STEP 1: TRIMMING${RESET}"
read -p "Trim from START (vols, e.g. 10): " TSTART
read -p "Trim from END (vols, e.g. 5): " TEND
if [[ "$TSTART" != "0" ]] || [[ "$TEND" != "0" ]]; then
  END=$(python3 -c "import nibabel as nib; n=nib.load('cleaned_mc_func.nii.gz').shape[-1]; print(n-int('$TEND')-1)")
  3dTcat -prefix trimmed_cleaned_mc_func.nii.gz cleaned_mc_func.nii.gz"[${TSTART}..${END}]"
  TMP="trimmed_cleaned_mc_func.nii.gz"
fi
echo "Trim,Start,$TSTART" >> "$LOG"; echo "Trim,End,$TEND" >> "$LOG"

# --- STEP 2: Smoothing ---
echo -e "${DARKRED}STEP 2: SMOOTHING${RESET}"
read -p "Temporal smoothing (sec, e.g. 2): " TSM
read -p "Spatial σ blur (e.g. 0.2): " SSM
[[ "$TSM" != "0" ]] && python3 - <<PY "$TMP" "$TSM" "$TR"
import sys,nibabel as nib,numpy as np; from scipy.ndimage import uniform_filter1d
f,tr=sys.argv[1],float(sys.argv[3]); win=max(1,int(round(float(sys.argv[2])/tr)))
img=nib.load(f); d=img.get_fdata(); sm=uniform_filter1d(d,size=win,axis=-1,mode='nearest')
out="smoothed_"+f; nib.save(nib.Nifti1Image(sm.astype(np.float32),img.affine,img.header),out)
print("[INFO] Smoothed →",out)
PY
[[ "$TSM" != "0" ]] && TMP="smoothed_${TMP}"
[[ "$SSM" != "0" ]] && fslmaths "$TMP" -s "$SSM" "spatial_${TMP}" && TMP="spatial_${TMP}"
echo "Temporal smoothing,$TSM," >> "$LOG"; echo "Spatial blur,$SSM," >> "$LOG"

# --- STEP 3: SCM QC ---
echo -e "${DARKRED}STEP 3: SCM QC${RESET}"
plot_global_signal "$TMP" "$MASK" "$TR"
read -p "Baseline (start end e.g. 1000 1050): " B1 B2
read -p "Signal (start end e.g. 1500 2000): " S1 S2
OUT="soner_pipeline/scm_outputs/raw_B${B1}-${B2}_S${S1}-${S2}"
scm_blush_fast "$TMP" "$B1" "$B2" "$S1" "$S2" "$MASK" "$OUT"
fsleyes "$OUT/norm_func_${B1}_${B2}.nii.gz" "$OUT/baseline_${B1}_${B2}.nii.gz" "$OUT/signal_change_map_${B1}_${B2}_${S1}_${S2}.nii.gz" &
echo "Baseline,$B1-$B2," >> "$LOG"; echo "Signal,$S1-$S2," >> "$LOG"

# --- STEP 4: FILTERING ---
declare -a FILTERED_FILES=()
echo -e "${DARKRED}STEP 4: FILTERING${RESET}"
echo -e "${DARKBLUE}HP removes drifts; LP removes fast noise. Ultra-slow (~0.1 Hz) BOLD oscillations reflect resting-state vasoneuronal coupling.${RESET}"
while true; do
  echo "1) High-pass"; echo "2) Low-pass"; echo "3) Band-pass"; echo "4) Continue → Decomposition"
  read -p "Choose [1–4]: " FTYPE; [[ "$FTYPE" == "4" ]] && break
  stamp=$(timestamp)
  case "$FTYPE" in
    1) read -p "Cut-off Hz (e.g. 0.01): " FC; OUT="soner_pipeline/filtering/hp_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FC" "0" "$TR" "1"; TMP="$OUT";;
    2) read -p "Cut-off Hz (e.g. 0.1): " FC; OUT="soner_pipeline/filtering/lp_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "0" "$FC" "$TR" "2"; TMP="$OUT";;
    3) read -p "Low–High Hz (e.g. 0.01 0.1): " FL FH; OUT="soner_pipeline/filtering/bp_${stamp}.nii.gz"; apply_freq_filter "$TMP" "$OUT" "$FL" "$FH" "$TR" "3"; TMP="$OUT";;
  esac
  FILTERED_FILES+=("$TMP")
done

# --- STEP 5: DECOMPOSITION ---
while true; do
  echo -e "${DARKRED}STEP 5: DECOMPOSITION${RESET}"
  echo "1) PCA"; echo "2) Temporal ICA"; echo "3) Spatial ICA"; echo "4) PLS"; echo "5) Re-run SCM"; echo "6) End"
  read -p "Choice [1–6]: " CH
  case "$CH" in
    6) break;;
    5)
      plot_global_signal "cleaned_mc_func.nii.gz" "$MASK" "$TR"
      read -p "Baseline (e.g. 1000 1050): " RB1 RB2; read -p "Signal (e.g. 1500 2000): " RS1 RS2
      OUTR="soner_pipeline/scm_outputs/ReRun_RAW_B${RB1}-${RB2}_S${RS1}-${RS2}"
      scm_blush_fast "cleaned_mc_func.nii.gz" "$RB1" "$RB2" "$RS1" "$RS2" "$MASK" "$OUTR"
      fsleyes "$OUTR/norm_func_${RB1}_${RB2}.nii.gz" "$OUTR/baseline_${RB1}_${RB2}.nii.gz" "$OUTR/signal_change_map_${RB1}_${RB2}_${RS1}_${RS2}.nii.gz" &
      continue;;
  esac
  case "$CH" in 1)METHOD="pca";; 2)METHOD="tica";; 3)METHOD="sica";; 4)METHOD="pls";; esac
  echo "Decomposition,$METHOD," >> "$LOG"

  ALL=("cleaned_mc_func.nii.gz" "${FILTERED_FILES[@]}")
  select FILE_CHOICE in "${ALL[@]}"; do TMP="$FILE_CHOICE"; break; done
  SRC=$(basename "$TMP")
  if [[ "$SRC" == *"hp_"* ]]; then SUB="hp"
  elif [[ "$SRC" == *"lp_"* ]]; then SUB="lp"
  elif [[ "$SRC" == *"bp_"* ]]; then SUB="bp"
  elif [[ "$SRC" == *"despiked"* ]]; then SUB="despiked"
  else SUB="raw"; fi

  run_decomposition "$TMP" "$METHOD" "$SUB"
  DEN="soner_pipeline/pca_outputs/${SUB}/${SUB}_${METHOD}_func_denoised.nii.gz"
  if [[ -f "$DEN" ]]; then
    plot_global_signal "$DEN" "$MASK" "$TR"
    read -p "Baseline (e.g. 1000 1050): " DB1 DB2; read -p "Signal (e.g. 1500 2000): " DS1 DS2
    OUTD="soner_pipeline/scm_outputs/${SUB}_${METHOD}_B${DB1}-${DB2}_S${DS1}-${DS2}"
    scm_blush_fast "$DEN" "$DB1" "$DB2" "$DS1" "$DS2" "$MASK" "$OUTD"
    fsleyes "$OUTD/norm_func_${DB1}_${DB2}.nii.gz" "$OUTD/baseline_${DB1}_${DB2}.nii.gz" "$OUTD/signal_change_map_${DB1}_${DB2}_${DS1}_${DS2}.nii.gz" &
  fi
done

# ===============================================================
#  HTML SUMMARY GENERATOR
# ===============================================================
python3 - <<PY "$LOG"
import sys,os,pandas as pd,datetime,glob
log=sys.argv[1]
df=pd.read_csv(log)
ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
html=f"<html><head><title>BLuSH Summary</title><style>body{{font-family:Arial;background:#111;color:#eee;padding:2em}}h2{{color:#8ef}}td,th{{padding:4px;border:1px solid #333}}</style></head><body><h1>BLuSH Summary Report</h1><p><i>Generated {ts}</i></p><h2>Parameter Log</h2>{df.to_html(index=False)}</body></html>"
out=log.replace("blush_run_","blush_summary_").replace(".csv",".html")
open(out,"w").write(html)
print(f"[OK] Summary HTML → {out}")
os.system(f"open {out}" if os.name=="posix" else f"start {out}")
PY

echo -e "${GREEN}All steps complete. HTML summary created.${RESET}"
