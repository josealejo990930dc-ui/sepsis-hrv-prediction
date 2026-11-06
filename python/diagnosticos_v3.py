# -*- coding: utf-8 -*-
"""
diagnosticos_v3.py  (read-only diagnostic cells; paste into Colab)
==================================================================
D1. Quality control of HRV extraction v3 (size, quality, deciles, flow).
D2. Check of R-peak placement (Pan-Tompkins vs NeuroKit vs NeuroKit + parabolic
    refinement) on ~20 real 5-min segments. This is the test that exposed the
    detector error of the previous extraction. To reproduce it on the OLD
    extraction, `hrv` must be the old CSV (hrv_completo_reauditado.csv).
D3. Comparison of SOFA windows with and without a new laboratory measurement in
    the horizon (basis of "Decision A": the new-lab criterion is a sensitivity
    analysis only).
No cell writes files or modifies tables.
"""
# %% ===== D1. QUALITY CONTROL OF HRV EXTRACTION v3 =====
from google.colab import drive
drive.mount('/content/drive')
import pandas as pd
import numpy as np

RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
v3  = pd.read_csv(f'{RUTA_BASE}/hrv_completo_v3.csv')
err = pd.read_csv(f'{RUTA_BASE}/errores_v3.csv')
rep = pd.read_csv(f'{RUTA_BASE}/df_reporte_final_unificado.csv')

# ---------- 1. Size ----------
print('=== 1. Tamaño ===')
print('Pacientes con HRV :', v3['icustay_id'].nunique())
print('Horas válidas     :', len(v3))
print(v3.groupby('icustay_id').size().describe().round(1))

# ---------- 2. Signal quality ----------
print('\n=== 2. Calidad de la señal ===')
rho1 = 1 - (v3['RMSSD_T0'] / v3['SDNN_T0']) ** 2 / 2
samp = v3['SampEn_T0']
print(f"rho1 mediana        : {rho1.median():.2f}   (antes -0.28 | esperado > 0)")
print(f"SD1 > SD2           : {(v3['SD1_T0'] > v3['SD2_T0']).mean() * 100:.1f} %   (antes 68.9 % | esperado bajo)")
print(f"DFA a1 mediana      : {v3['DFA_a1_T0'].median():.2f}   (antes ~0.30 | esperado > 0.6)")
print(f"SDNN mediana (ms)   : {v3['SDNN_T0'].median():.1f}")
print(f"RMSSD mediana (ms)  : {v3['RMSSD_T0'].median():.1f}")
print(f"SampEn inf/negativa : {(~np.isfinite(samp) | (samp < 0)).sum()}   (esperado 0)")

# ---------- 3. By SDNN decile ----------
print('\n=== 3. Por decil de SDNN ===')
t = v3[['SDNN_T0', 'RMSSD_T0', 'SD1_T0', 'SD2_T0', 'DFA_a1_T0']].dropna().copy()
t['decil'] = pd.qcut(t['SDNN_T0'], 10, labels=False) + 1
t['sd1_mayor'] = t['SD1_T0'] > t['SD2_T0']
t['rho1'] = 1 - (t['RMSSD_T0'] / t['SDNN_T0']) ** 2 / 2
print(t.groupby('decil').agg(
    SDNN=('SDNN_T0', 'median'),
    pct_SD1_mayor_SD2=('sd1_mayor', lambda x: round(100 * x.mean(), 1)),
    rho1=('rho1', 'median'),
    DFA_a1=('DFA_a1_T0', 'median')).round(2))

# ---------- 4. Patient flow ----------
print('\n=== 4. Flujo de pacientes ===')
con_hrv = set(v3['icustay_id'])
con_error = set(err['icustay_id'])
print('Cohorte sincronizada            :', rep['icustay_id'].nunique())
print('Aptos para extracción           :', int(rep['apto_extraccion'].sum()))
print('Con HRV extraída                :', len(con_hrv))
print('Con error registrado            :', len(con_error))
print('En ambos (reintentados con éxito):', len(con_hrv & con_error))
print('Total procesados                :', len(con_hrv | con_error))
print('\nPacientes únicos por etapa de error:')
print(err.groupby('etapa')['icustay_id'].nunique().sort_values(ascending=False))
print('Detalle de los fallos técnicos:')
print(err[err['etapa'].isin(['runtime_error', 'extraccion_ecg'])][['icustay_id', 'etapa', 'detalle']])


# %% ===== D2. R-PEAK DETECTOR CHECK =====
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA_REPORTE = f'{RUTA_BASE}/df_reporte_final_unificado.csv'
hrv = pd.read_csv(f'{RUTA_BASE}/hrv_completo_reauditado.csv')  # OLD extraction (Pan-Tompkins)
# Read-only. Downloads ~20 five-minute segments from PhysioNet (the same ones behind each
# *_T0 value: epoch 12 of the hour) and compares three R-peak detectors:
#   1. Pan-Tompkins, exactly as in the old extraction
#   2. NeuroKit default detector
#   3. NeuroKit + parabolic peak refinement (sub-sample precision)
import subprocess, re, warnings
subprocess.run(['pip', 'install', 'wfdb', 'neurokit2', '-q'], check=True)
import wfdb, neurokit2 as nk
warnings.filterwarnings('ignore')

N_POR_DECIL = 2
SEMILLA = 42

# --- Old extraction functions, copied verbatim to reproduce the stored values ---
def _mascara_ectopicos(rr_ms):
    n = len(rr_ms); fuera = (rr_ms < 300) | (rr_ms > 2000); dev = np.zeros(n, bool)
    for i in range(n):
        v = np.concatenate([rr_ms[max(0, i - 2):i], rr_ms[i + 1:min(n, i + 3)]])
        if len(v) >= 2:
            med = np.median(v)
            if med > 0: dev[i] = abs(rr_ms[i] - med) / med > 0.20
    return fuera | dev

def _corregir_pipeline(rpeaks_idx, fs):
    rr = np.diff(rpeaks_idx).astype(float) / fs * 1000; n = len(rr); rc = rr.copy()
    ect = _mascara_ectopicos(rr)
    for i in range(n):
        if ect[i]:
            izq = i - 1
            while izq >= 0 and ect[izq]: izq -= 1
            der = i + 1
            while der < n and ect[der]: der -= 1  # as in the old extraction (bug fixed in v3: der += 1)
            a = rr[izq] if izq >= 0 else np.nan; b = rr[der] if der < n else np.nan
            if not np.isnan(a) and not np.isnan(b): rc[i] = (a + b) / 2
            elif not np.isnan(a): rc[i] = a
            elif not np.isnan(b): rc[i] = b
    out = np.empty(len(rpeaks_idx), dtype=int); out[0] = rpeaks_idx[0]
    for j in range(n): out[j + 1] = out[j] + int(round(rc[j] / 1000 * fs))
    return out

# --- RR statistics (excluding intervals flagged as ectopic) ---
def stats_rr(t_s):
    rr = np.diff(t_s) * 1000
    ok = ~_mascara_ectopicos(rr)
    par = ok[:-1] & ok[1:]
    d = np.diff(rr)[par]
    return {'n_picos': len(t_s), 'carga_ect_%': round(100 * (~ok).mean(), 1),
            'SDNN': rr[ok].std(ddof=1), 'RMSSD': np.sqrt(np.mean(d ** 2)),
            'rho1': np.corrcoef(rr[:-1][par], rr[1:][par])[0, 1], 'FC': 60000 / rr[ok].mean()}

def refinar_parabolico(sig, peaks):
    signo = np.sign(np.median(sig[peaks])) or 1.0; s = sig * signo; out = []
    for p in peaks:
        k = int(np.clip(p, 1, len(s) - 2))
        k = k - 1 + int(np.argmax(s[k - 1:k + 2]))  # ensure local maximum
        k = int(np.clip(k, 1, len(s) - 2))
        y0, y1, y2 = s[k - 1], s[k], s[k + 1]; den = y0 - 2 * y1 + y2
        out.append(k + (np.clip(0.5 * (y0 - y2) / den, -1, 1) if den != 0 else 0.0))
    return np.array(out)

def analizar_segmento(ecg, fs):
    res = {}
    # 1. Pan-Tompkins (identical to the old extraction)
    cl_pt = nk.ecg_clean(ecg, sampling_rate=fs, method='pantompkins1985')
    _, i_pt = nk.ecg_peaks(cl_pt, sampling_rate=fs, method='pantompkins1985')
    pk_pt = np.asarray(i_pt['ECG_R_Peaks'])
    res['PanTompkins'] = stats_rr(pk_pt / fs)
    # RMSSD reproduced with the full old pipeline (alignment check)
    hrv_t = nk.hrv_time({'ECG_R_Peaks': _corregir_pipeline(pk_pt, fs)}, sampling_rate=fs, show=False)
    res['RMSSD_reproducido'] = float(hrv_t['HRV_RMSSD'].iloc[0])
    # 2. NeuroKit
    cl_nk = nk.ecg_clean(ecg, sampling_rate=fs, method='neurokit')
    _, i_nk = nk.ecg_peaks(cl_nk, sampling_rate=fs, method='neurokit')
    pk_nk = np.asarray(i_nk['ECG_R_Peaks'])
    res['NeuroKit'] = stats_rr(pk_nk / fs)
    # 3. NeuroKit + parabolic refinement
    res['NeuroKit_refinado'] = stats_rr(refinar_parabolico(cl_nk, pk_nk) / fs)
    return res, cl_nk, pk_pt, pk_nk

# --- Sample stratified by SDNN decile ---
rep = pd.read_csv(RUTA_REPORTE)
rep['icustay_id'] = rep['icustay_id'].astype('int64')
for c in ['inicio_analisis_real', 'waveform_start']:
    rep[c] = pd.to_datetime(rep[c]).dt.tz_localize(None)

m = hrv[['subject_id', 'icustay_id', 'hour_index', 'SDNN_T0', 'RMSSD_T0', 'DFA_a1_T0']].copy()
m['decil'] = pd.qcut(m['SDNN_T0'], 10, labels=False) + 1
muestra = (m.groupby('decil', group_keys=False).apply(lambda g: g.sample(min(N_POR_DECIL, len(g)), random_state=SEMILLA))
             .merge(rep[['icustay_id', 'ruta_registro', 'inicio_analisis_real', 'waveform_start']], on='icustay_id', how='left'))

filas, ejemplos = [], {}
for _, r in muestra.iterrows():
    try:
        partes = r['ruta_registro'].split('/'); rec = partes[-1]
        pn = f"mimic3wdb-matched/1.0/{'/'.join(partes[:-1])}"
        fs = int(wfdb.rdheader(rec, pn_dir=pn).fs)
        base = max(0, int((r['inicio_analisis_real'] - r['waveform_start']).total_seconds() * fs))
        ini = base + int(r['hour_index']) * 3600 * fs + 11 * 300 * fs  # epoch 12 = the one behind *_T0
        rc = wfdb.rdrecord(rec, pn_dir=pn, sampfrom=ini, sampto=ini + 300 * fs)
        idx = next(i for i, nm in enumerate(rc.sig_name) if re.search(r'\bII\b', nm.upper()))
        ecg = np.nan_to_num(rc.p_signal[:, idx])
        res, cl, pk_pt, pk_nk = analizar_segmento(ecg, fs)
        fila = {'icustay_id': r['icustay_id'], 'hora': int(r['hour_index']), 'decil_SDNN': int(r['decil']), 'fs': fs,
                'RMSSD_guardado': r['RMSSD_T0'], 'RMSSD_reproducido': res['RMSSD_reproducido']}
        for met in ['PanTompkins', 'NeuroKit', 'NeuroKit_refinado']:
            for k, v in res[met].items():
                fila[f'{met}_{k}'] = v
        filas.append(fila)
        ejemplos[int(r['decil'])] = (cl, pk_pt, pk_nk, fs)
        print(f"  ✓ decil {int(r['decil'])} | icustay {r['icustay_id']} hora {int(r['hour_index'])} | fs={fs}")
    except Exception as e:
        print(f"  ✗ icustay {r['icustay_id']} hora {int(r['hour_index'])}: {e}")

diag = pd.DataFrame(filas)

# --- 1. Are we reading the same segment as the extraction? ---
print("\n1. ALINEACIÓN: RMSSD guardado vs reproducido (deberían coincidir)")
print(diag[['decil_SDNN', 'RMSSD_guardado', 'RMSSD_reproducido']].round(1))
print(f"   Diferencia absoluta mediana: {(diag.RMSSD_guardado - diag.RMSSD_reproducido).abs().median():.2f} ms")

# --- 2. Detector comparison ---
print("\n2. COMPARACIÓN DE DETECTORES (medianas sobre la muestra)")
comp = pd.DataFrame({met: {k: diag[f'{met}_{k}'].median() for k in ['n_picos', 'carga_ect_%', 'FC', 'SDNN', 'RMSSD', 'rho1']}
                     for met in ['PanTompkins', 'NeuroKit', 'NeuroKit_refinado']}).T
comp['RMSSD/SDNN'] = comp['RMSSD'] / comp['SDNN']
print(comp.round(2))

print("\n   rho1 por decil de SDNN:")
print(diag.groupby('decil_SDNN')[['PanTompkins_rho1', 'NeuroKit_rho1', 'NeuroKit_refinado_rho1']].median().round(2))

# --- 3. Plot: 10 s of ECG with the peaks of each detector ---
for dec in [d for d in [1, 5, 10] if d in ejemplos]:
    cl, pk_pt, pk_nk, fs = ejemplos[dec]
    a, b = 60 * fs, 70 * fs
    t = np.arange(a, b) / fs
    plt.figure(figsize=(15, 3))
    plt.plot(t, cl[a:b], lw=.8, c='gray')
    p1 = pk_pt[(pk_pt >= a) & (pk_pt < b)]; p2 = pk_nk[(pk_nk >= a) & (pk_nk < b)]
    plt.scatter(p1 / fs, cl[p1], c='red', marker='x', s=60, label='Pan-Tompkins (extracción)', zorder=3)
    plt.scatter(p2 / fs, cl[p2], facecolors='none', edgecolors='green', s=80, label='NeuroKit', zorder=3)
    plt.title(f'Decil de SDNN {dec}: ¿dónde cae cada pico?'); plt.xlabel('s'); plt.legend(loc='upper right')
    plt.tight_layout(); plt.show()

print('''
>> CÓMO LEERLO
   1. Si RMSSD guardado ≈ reproducido, estamos analizando exactamente el segmento que generó los datos.
   2. Si Pan-Tompkins da rho1 negativa y RMSSD/SDNN ~1,5-1,7, mientras NeuroKit da rho1 positiva
      y RMSSD/SDNN < 1,2 en el mismo segmento, el problema está en el detector.
   3. En la gráfica, las X rojas deberían caer sobre la punta de la onda R. Si caen desplazadas
      o saltan entre la R y otro punto del complejo, eso es el error de ubicación.
''')


# %% ===== D3. SOFA WINDOWS WITH AND WITHOUT A NEW LAB VALUE IN THE HORIZON =====
import pandas as pd
df = pd.read_csv('/content/drive/MyDrive/tesis_sepsis/ultimate_dataset_v3.csv')

s = df[(df['incl_hrv_6h'] == 1) & df['label_sofa_responder'].notna() & (df['muere_en_horizonte'] == 0)].copy()
s['grupo'] = s['new_lab_6h'].fillna(0).map({1: 'Con laboratorio nuevo', 0: 'Sin laboratorio nuevo'})

print(s.groupby('grupo').agg(
    ventanas=('icustay_id', 'size'),
    pacientes=('icustay_id', 'nunique'),
    pct_delta_cero=('delta_sofa_pareado', lambda x: round(100 * (x == 0).mean(), 1)),
    pct_responder=('label_sofa_responder', lambda x: round(100 * x.mean(), 1)),
    pct_deterioro=('label_sofa_deterioro_2pts', lambda x: round(100 * x.mean(), 1)),
    sofa_t0_mediana=('sofa_t0_parcial', 'median'),
    pct_con_vasopresor=('ned', lambda x: round(100 * (x > 0).mean(), 1)),
    hora_t0_mediana=('hour_index', 'median')))
