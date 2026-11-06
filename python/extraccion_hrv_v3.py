"""
extraccion_hrv_v3.py  (HRV extraction v3, all available signal, capped at 96 h)
================================================================================
One extraction per patient, up to 96 h of ECG (or whatever is available).
Written as a Google Colab script: it runs top to bottom when executed (no main()).

Input : Sepsis-3 cohort and HRV anchor time from BigQuery
        (sepsis3_final_cohort, hrv_analysis_mapping) and the ECG waveforms of the
        MIMIC-III Waveform Database Matched Subset (mimic3wdb-matched/1.0, via wfdb).
Output: hrv_completo_v3.csv (one row per patient-hour with a valid HRV hour:
        10 metrics at T0 plus 110 deltas), df_reporte_final_unificado.csv (synchronised
        cohort and time anchor), errores_v3.csv, checkpoint_v3.csv, reporte_pacientes_v3.txt.

Method:
  - R peaks: NeuroKit (method='neurokit') + parabolic refinement (sub-sample
    precision; at 125 Hz one sample = 8 ms). Pan-Tompkins was rejected: at 125 Hz
    it misplaced the R peak, inflating SDNN ~2x and RMSSD ~5x (checked on 20 real
    segments; see diagnosticos_v3.py, D2).
  - Ectopic beats flagged and corrected on RR intervals in ms; epochs with
    > 5 % flagged beats are invalid; infinite or negative SampEn -> invalid epoch.
  - Valid hour: >= 8/11 valid previous epochs + valid epoch 12 (T0 of the hour);
    output = 10 metrics at T0 + 110 deltas; >= 3 valid hours per patient.
  - Patients with < 12 h of signal from the anchor are skipped (no 12-h window possible).
  - The >= 6 valid hours per window filter is applied later, on the window grid.
In the code a 5-min segment is called "ventana"; the README calls it an "epoch"
to avoid confusion with the 12-h prediction windows.
"""

# Installation and record index (equivalent to !pip and !wget in Colab)
import subprocess
subprocess.run(['pip', 'install', '-q', 'wfdb', 'neurokit2', 'pandas', 'numpy', 'google-cloud-bigquery'], check=True)
subprocess.run(['wget', '-q', '-N', '-c', 'https://physionet.org/files/mimic3wdb-matched/1.0/RECORDS-waveforms'], check=True)

# ─── 0. SETUP ───────────────────────────────────────────────────────────────
import os
from google.colab import drive

try:
    drive.flush_and_unmount()
    print("Drive desmontado correctamente.")
except Exception as e:
    print(f"Aviso al desmontar: {e}")

mountpoint_path = '/content/drive'
if os.path.exists(mountpoint_path) and os.listdir(mountpoint_path):
    print(f"Mountpoint '{mountpoint_path}' contiene residuos — limpiando.")
    os.system(f"rm -rf {mountpoint_path}/*")
drive.mount(mountpoint_path, force_remount=True)

RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
os.makedirs(RUTA_BASE, exist_ok=True)

import pandas as pd
import numpy as np
import re
import time
import warnings
from datetime import datetime
from tqdm.auto import tqdm
import wfdb
import neurokit2 as nk
warnings.filterwarnings('ignore')

from google.colab import auth
auth.authenticate_user()
from google.cloud import bigquery
client = bigquery.Client(project='sepsis-mimic-def')

HORIZONTE_UNIFICADO_HORAS = 12.0  # buscar_inicio_laxo: single horizon
HORAS_OBJETIVO_TOPE = 96          # maximum hours extracted per patient
HORAS_MINIMAS = 3                 # minimum valid hours per patient
HORAS_MINIMAS_SENAL = 12          # v3: minimum signal from the anchor to form >= 1 window of 12 h

# v3 outputs (new files: do not reuse those of the previous extraction)
RUTA_HRV_V3        = f'{RUTA_BASE}/hrv_completo_v3.csv'
RUTA_ERRORES_V3    = f'{RUTA_BASE}/errores_v3.csv'
RUTA_CHECKPOINT_V3 = f'{RUTA_BASE}/checkpoint_v3.csv'
RUTA_REPORTE_V3    = f'{RUTA_BASE}/reporte_pacientes_v3.txt'

# ─── 1. PHYSIONET RECORD INDEX ──────────────────────────────────────────────
df_wfdb = pd.read_csv('RECORDS-waveforms', header=None, names=['ruta_registro'])
df_wfdb['subject_id'] = df_wfdb['ruta_registro'].str.extract(r'p(\d{6})').astype(float)
df_wfdb = df_wfdb.dropna(subset=['subject_id'])
print(f"Estructura indexada. Se encontraron {len(df_wfdb)} registros vinculables.")

# ─── 2. CLINICAL COHORT FROM BIGQUERY ───────────────────────────────────────
# tipo_t0 classifies the sepsis onset relative to ICU admission: before admission,
# within the first hour, or later during the stay.
query_cohort = '''
SELECT
  cohorte.subject_id,
  cohorte.icustay_id,
  ie.intime,
  cohorte.t0_sepsis,
  map.hrv_analysis_start_time,
  CASE
    WHEN DATETIME_DIFF(cohorte.t0_sepsis, ie.intime, MINUTE) < 0 THEN 'sepsis_previa_al_ingreso'
    WHEN DATETIME_DIFF(cohorte.t0_sepsis, ie.intime, MINUTE) BETWEEN 0 AND 60 THEN 'sepsis_primera_hora'
    ELSE 'sepsis_durante_estancia'
  END AS tipo_t0
FROM `sepsis-mimic-def.sepsis3_workspace_def.sepsis3_final_cohort` cohorte
INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON cohorte.icustay_id = ie.icustay_id
INNER JOIN `sepsis-mimic-def.sepsis3_workspace_def.hrv_analysis_mapping` map ON cohorte.icustay_id = map.icustay_id
'''
df_mapping = client.query(query_cohort).to_dataframe()
df_mapping['intime'] = pd.to_datetime(df_mapping['intime'])
df_mapping['t0_sepsis'] = pd.to_datetime(df_mapping['t0_sepsis'])
df_mapping['hrv_analysis_start_time'] = pd.to_datetime(df_mapping['hrv_analysis_start_time'])
print(f"Cohorte cargada: {len(df_mapping)} entradas")
display(df_mapping['tipo_t0'].value_counts())


# ─── 3. buscar_inicio_laxo (unchanged, single 12-h horizon) ─────────────────
def buscar_inicio_laxo(row, horizonte_prediccion_horas=HORIZONTE_UNIFICADO_HORAS):
    """Find the actual start of the HRV analysis inside the waveform record.

    - If the theoretical anchor (hrv_analysis_start_time) lies inside the record,
      use it.
    - If sepsis started before ICU admission and the record starts later, use the
      record start, provided the shift is < horizon (12 h).
    - Otherwise try the anchor + 1, 2, ... (horizon - 1) hours.
    Returns NaT if no valid start is found (patient discarded)."""
    t_teorico = row['hrv_analysis_start_time']
    w_start = row['waveform_start']
    w_end = row['waveform_end']

    if pd.isna(t_teorico):
        return pd.NaT

    limite = pd.to_timedelta(horizonte_prediccion_horas, unit='h')

    if w_start <= t_teorico <= w_end:
        return t_teorico

    if row['tipo_t0'] == 'sepsis_previa_al_ingreso' and w_start > t_teorico:
        corrimiento = w_start - t_teorico
        if corrimiento < limite:
            return w_start
        return pd.NaT

    max_saltos = int(horizonte_prediccion_horas)
    for salto_hora in range(1, max(max_saltos, 1)):
        t_alternativo = t_teorico + pd.to_timedelta(salto_hora, unit='h')
        if w_start <= t_alternativo <= w_end:
            return t_alternativo

    return pd.NaT

print("Función 'buscar_inicio_laxo' definida (horizonte único = 12h).")


# ─── 4. WAVEFORM SYNCHRONISATION (unchanged) ────────────────────────────────
# Link each ICU stay to the waveform record whose start is 0-48 h before the
# HRV anchor, keeping the closest one per icustay_id.
df_cohorte_final = df_mapping.copy()
df_cohorte_final['hrv_analysis_start_time_dt'] = df_cohorte_final['hrv_analysis_start_time'].dt.tz_localize(None)

patron_wfdb = r'p(\d{6})/p\d{6}-(\d{4}-\d{2}-\d{2}-\d{2}-\d{2})'
extraccion = df_wfdb['ruta_registro'].str.extract(patron_wfdb)
df_records_parsed = df_wfdb.copy()
df_records_parsed['subject_id'] = extraccion[0].astype(int)
df_records_parsed['waveform_start_orig'] = pd.to_datetime(
    extraccion[1], format='%Y-%m-%d-%H-%M', errors='coerce'
).dt.tz_localize(None)
df_records_parsed = df_records_parsed.dropna(subset=['subject_id', 'waveform_start_orig'])

df_merged = pd.merge(df_cohorte_final, df_records_parsed, on='subject_id', how='inner')
df_merged['desfase_horas'] = (
    df_merged['hrv_analysis_start_time_dt'] - df_merged['waveform_start_orig']
).dt.total_seconds() / 3600

df_procesamiento = (
    df_merged[(df_merged['desfase_horas'] >= 0) & (df_merged['desfase_horas'] <= 48)]
    .sort_values('desfase_horas')
    .drop_duplicates(subset=['icustay_id'], keep='first')
    .reset_index(drop=True)
)
print(f"Pacientes tras el filtro de 48h y selección del waveform más cercano: {len(df_procesamiento)}")

print("Leyendo encabezados de waveform (metadatos reales)...")
waveform_details = []
for _, row in tqdm(df_procesamiento.iterrows(), total=len(df_procesamiento), desc="Leyendo encabezados"):
    try:
        partes = row['ruta_registro'].split('/')
        header = wfdb.rdheader(partes[-1], pn_dir=f"mimic3wdb-matched/1.0/{'/'.join(partes[:-1])}")
        duration_s = header.sig_len / header.fs
        waveform_details.append({
            'icustay_id': row['icustay_id'],
            'waveform_start': row['waveform_start_orig'],
            'waveform_end': row['waveform_start_orig'] + pd.to_timedelta(duration_s, unit='s'),
            'duration_hours_real': duration_s / 3600
        })
    except Exception as e:
        print(f"  [WARN] icustay_id={row.get('icustay_id', 'NA')}: {e}")

df_reporte_base = pd.merge(df_procesamiento, pd.DataFrame(waveform_details), on='icustay_id', how='inner')
print(f"Pacientes con metadatos de waveform leídos: {len(df_reporte_base)}")

df_reporte_final = df_reporte_base.copy()
df_reporte_final['inicio_analisis_real'] = df_reporte_final.apply(
    lambda row: buscar_inicio_laxo(row, horizonte_prediccion_horas=HORIZONTE_UNIFICADO_HORAS), axis=1)

n_total = len(df_reporte_final)
n_descartados = df_reporte_final['inicio_analisis_real'].isna().sum()
df_reporte_final = df_reporte_final.dropna(subset=['inicio_analisis_real']).reset_index(drop=True)
print(f"Total antes del filtro de corrimiento: {n_total}")
print(f"Descartados por corrimiento >= {HORIZONTE_UNIFICADO_HORAS}h: {n_descartados} ({n_descartados/n_total*100:.1f}%)")
print(f"Cohorte final sincronizada (única): {len(df_reporte_final)}")

df_reporte_final['horas_disponibles'] = (
    df_reporte_final['waveform_end'] - df_reporte_final['inicio_analisis_real']
).dt.total_seconds() / 3600.0
df_reporte_final['apto_extraccion'] = df_reporte_final['horas_disponibles'] >= HORAS_MINIMAS

df_reporte_final['corrimiento_horas'] = (
    pd.to_datetime(df_reporte_final['inicio_analisis_real']).dt.tz_localize(None)
    - pd.to_datetime(df_reporte_final['hrv_analysis_start_time']).dt.tz_localize(None)
).dt.total_seconds() / 3600.0
n_invade = (df_reporte_final['corrimiento_horas'] >= HORIZONTE_UNIFICADO_HORAS).sum()
print(f"Pacientes con corrimiento >= horizonte: {n_invade} (debería ser 0)")
assert n_invade == 0, "Aún hay corrimientos que invaden la ventana de predicción."
print("OK — ningún corrimiento invade la ventana de predicción.\n")

RUTA_COHORTE = f'{RUTA_BASE}/df_reporte_final_unificado.csv'
df_reporte_final.to_csv(RUTA_COHORTE, index=False)
print(f"✅ Cohorte unificada guardada: {RUTA_COHORTE} ({len(df_reporte_final)} filas)")
print("\nDistribución de horas_disponibles (antes del tope de 96h):")
print(df_reporte_final['horas_disponibles'].describe())
print(f"Pacientes con ≥{HORAS_MINIMAS_SENAL}h de señal (pueden formar ventanas): "
      f"{(df_reporte_final['horas_disponibles'] >= HORAS_MINIMAS_SENAL).sum()}")


# ═══════════════════════════════════════════════════════════════════════════
# 5. HRV FUNCTIONS — v3
# ═══════════════════════════════════════════════════════════════════════════
# NOTE: these two lists are only used here to enumerate the 10 metrics
# (TODAS_LAS_METRICAS); they are an older grouping. The domain grouping used in
# the dataset and the models (SD1 and SD2 in the linear/time domain, Brennan et
# al. 2001) is defined in construir_ultimate_dataset_v3.py.
METRICAS_LINEALES    = ['SDNN', 'RMSSD', 'LF', 'HF', 'LFHF']
METRICAS_NO_LINEALES = ['SD1', 'SD2', 'SampEn', 'ApEn', 'DFA_a1']
TODAS_LAS_METRICAS   = METRICAS_LINEALES + METRICAS_NO_LINEALES  # 10 in total

# One representative metric per domain to decide whether an epoch is valid
METRICAS_VALIDACION = ['SDNN', 'RMSSD', 'SampEn']


def _refinar_parabolico(sig: np.ndarray, peaks: np.ndarray) -> np.ndarray:
    """Place each peak at the local maximum of the R wave with sub-sample
    precision (parabolic interpolation over the 3 samples around the maximum)."""
    signo = np.sign(np.median(sig[peaks])) or 1.0      # supports an inverted QRS
    s = sig * signo
    out = np.empty(len(peaks))
    for j, p in enumerate(peaks):
        k = int(np.clip(p, 1, len(s) - 2))
        k = int(np.clip(k - 1 + int(np.argmax(s[k - 1:k + 2])), 1, len(s) - 2))
        y0, y1, y2 = s[k - 1], s[k], s[k + 1]
        den = y0 - 2 * y1 + y2
        out[j] = k + (np.clip(0.5 * (y0 - y2) / den, -1, 1) if den != 0 else 0.0)
    return out


def _detectar_rr_ms(fragmento: np.ndarray, fs: int) -> np.ndarray:
    """ECG (5 min) -> RR intervals in ms (float), from refined peaks."""
    ecg_limpio = nk.ecg_clean(fragmento, sampling_rate=fs, method='neurokit')
    _, info = nk.ecg_peaks(ecg_limpio, sampling_rate=fs, method='neurokit')
    picos = np.asarray(info['ECG_R_Peaks'])
    if len(picos) < 3:
        return np.array([])
    return np.diff(_refinar_parabolico(ecg_limpio, picos)) / fs * 1000


def _mascara_ectopicos_rr(rr_ms: np.ndarray) -> np.ndarray:
    """Flag ectopic/artefactual beats using two physiological criteria:
      1. Out of range: RR < 300 ms or > 2000 ms.
      2. Deviation > 20 % from the local median (up to 2 neighbours on each side)."""
    n = len(rr_ms)
    fuera = (rr_ms < 300) | (rr_ms > 2000)
    dev = np.zeros(n, dtype=bool)
    for i in range(n):
        v = np.concatenate([rr_ms[max(0, i - 2):i], rr_ms[i + 1:min(n, i + 3)]])
        if len(v) >= 2:
            med = np.median(v)
            if med > 0:
                dev[i] = abs(rr_ms[i] - med) / med > 0.20
    return fuera | dev


def _corregir_rr(rr_ms: np.ndarray, ect: np.ndarray) -> np.ndarray:
    """Replace each flagged RR by the mean of its nearest valid neighbours on the
    left and on the right (or by the only one available)."""
    rc = rr_ms.copy()
    n = len(rr_ms)
    for i in np.where(ect)[0]:
        izq = i - 1
        while izq >= 0 and ect[izq]:
            izq -= 1
        der = i + 1
        while der < n and ect[der]:
            der += 1                                      # v3: fixed (previously der -= 1)
        vecinos = [rr_ms[k] for k in (izq, der) if 0 <= k < n]
        if vecinos:
            rc[i] = np.mean(vecinos)
    return rc


def _metricas_desde_rr(rr_ms: np.ndarray) -> dict:
    """10 HRV metrics in three domains, computed from the corrected RR intervals
    with NeuroKit2 defaults (hrv_time, hrv_frequency, hrv_nonlinear)."""
    rri = {'RRI': rr_ms, 'RRI_Time': np.cumsum(rr_ms) / 1000}
    hrv_t  = nk.hrv_time(rri, show=False)
    hrv_f  = nk.hrv_frequency(rri, show=False)
    hrv_nl = nk.hrv_nonlinear(rri, show=False)

    def _val(df, col):
        return float(df[col].iloc[0]) if col in df.columns else np.nan

    return {
        'SDNN'  : _val(hrv_t,  'HRV_SDNN'),
        'RMSSD' : _val(hrv_t,  'HRV_RMSSD'),
        'LF'    : _val(hrv_f,  'HRV_LF'),
        'HF'    : _val(hrv_f,  'HRV_HF'),
        'LFHF'  : _val(hrv_f,  'HRV_LFHF'),
        'SD1'   : _val(hrv_nl, 'HRV_SD1'),
        'SD2'   : _val(hrv_nl, 'HRV_SD2'),
        'SampEn': _val(hrv_nl, 'HRV_SampEn'),
        'ApEn'  : _val(hrv_nl, 'HRV_ApEn'),
        'DFA_a1': _val(hrv_nl, 'HRV_DFA_alpha1'),
    }


def _metricas_ventana(fragmento: np.ndarray, fs: int) -> dict | None:
    """One 5-min epoch -> metrics, or None if it fails quality control."""
    rr = _detectar_rr_ms(fragmento, fs)
    if len(rr) < 2:
        return None
    ect = _mascara_ectopicos_rr(rr)
    if ect.mean() > 0.05:                                 # ectopic load > 5 % -> epoch excluded
        return None
    m = _metricas_desde_rr(_corregir_rr(rr, ect))
    if not np.isfinite(m['SampEn']) or m['SampEn'] < 0:   # v3: impossible value -> invalid epoch
        m['SampEn'] = np.nan
    return m


def procesar_hora_hrv(senal_ecg: np.ndarray, fs: int = 125, verbose: bool = False) -> dict | None:
    """
    Process 1 hour of continuous ECG (60 x 60 x fs samples).
      1. Split the signal into 12 non-overlapping 5-min epochs.
      2. Detect R peaks (NeuroKit) and refine them to sub-sample precision.
      3. Exclude the epoch if the ectopic load exceeds 5 %.
      4. Correct the remaining flagged RR intervals (interpolation with valid neighbours).
      5. Extract the 10 HRV metrics.
      6. Eligibility rule: >= 8 of the 11 previous epochs valid and the anchor
         epoch (epoch 12) valid, evaluated on METRICAS_VALIDACION (three domains).
      7. Return the 10 metrics at T0 (epoch 12) + 110 deltas, or None.
    Delta_{m}_minus_{N}min = value at epoch 12 - value of the epoch N minutes earlier.
    """
    if senal_ecg is None or len(senal_ecg) < (60 * 60 * fs):
        return None

    muestras_5min = 5 * 60 * fs
    nan_row = {m: np.nan for m in TODAS_LAS_METRICAS}
    ventanas = []
    for i in range(12):
        frag = senal_ecg[i * muestras_5min:(i + 1) * muestras_5min]
        try:
            m = _metricas_ventana(frag, fs)
        except Exception as e:
            if verbose:
                print(f"    ventana {i + 1}: error {e}")
            m = None
        ventanas.append(m if m is not None else nan_row.copy())
    df_ventanas = pd.DataFrame(ventanas)[TODAS_LAS_METRICAS]

    previas_validas = df_ventanas.iloc[:11][METRICAS_VALIDACION].notna().all(axis=1).sum()
    anclaje_valido = df_ventanas.iloc[11][METRICAS_VALIDACION].notna().all()
    if verbose:
        print(f"    previas válidas: {previas_validas}/11 | ancla válida: {anclaje_valido}")
    if previas_validas < 8 or not anclaje_valido:
        return None

    features = {}
    ancla = df_ventanas.iloc[11]
    for met in TODAS_LAS_METRICAS:
        features[f'{met}_T0'] = ancla[met]
    for iv in range(11):
        minutos = (11 - iv) * 5                           # 55, 50, ..., 5 minutes
        previa = df_ventanas.iloc[iv]
        for met in TODAS_LAS_METRICAS:
            features[f'Delta_{met}_minus_{minutos}min'] = ancla[met] - previa[met]
    return features                                       # 10 + 110 = 120 features


def procesar_bloque_longitudinal(senal_completa, fs, horas_objetivo=6, horas_minimas=3, verbose=False):
    """Process up to `horas_objetivo` hours of signal, one row per valid hour.
    Discards the patient if there are fewer than `horas_minimas` valid hours.
    Returns (DataFrame | None, number of valid hours). Unchanged in v3."""
    muestras_por_hora = 3600 * fs
    total_horas_disponibles = len(senal_completa) // muestras_por_hora
    horas_a_procesar = min(horas_objetivo, total_horas_disponibles)

    all_hourly_features = []
    iterador = tqdm(range(horas_a_procesar), desc="Horas", leave=False) if not verbose else range(horas_a_procesar)
    for h in iterador:
        inicio = h * muestras_por_hora
        fragmento_1h = senal_completa[inicio: inicio + muestras_por_hora]
        features_1h = procesar_hora_hrv(fragmento_1h, fs=fs, verbose=verbose)
        if features_1h:
            features_1h['hour_index'] = h
            all_hourly_features.append(features_1h)

    horas_validas = len(all_hourly_features)
    if horas_validas < horas_minimas:
        return None, horas_validas

    df_resultado = pd.DataFrame(all_hourly_features)
    df_resultado['horas_validas_total'] = horas_validas
    df_resultado['horas_objetivo'] = horas_objetivo
    return df_resultado, horas_validas

print("✅ Funciones HRV v3 cargadas (NeuroKit + refinamiento parabólico).")


# ═══════════════════════════════════════════════════════════════════════════
# 6. UNIFIED EXTRACTION LOOP
# ═══════════════════════════════════════════════════════════════════════════
def _wfdb_con_reintentos(func, *args, max_intentos=4, espera_base=5, **kwargs):
    """Call a wfdb download function with exponential back-off retries
    (no retry on 404 / Not Found)."""
    ultimo_error = None
    for intento in range(1, max_intentos + 1):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            ultimo_error = e
            mensaje = str(e)
            if '404' in mensaje or 'Not Found' in mensaje:
                raise
            if intento < max_intentos:
                time.sleep(espera_base * (2 ** (intento - 1)))
    raise ultimo_error


def correr_extraccion_hrv_unificada(
    df_reporte_final, horas_objetivo_tope=HORAS_OBJETIVO_TOPE, horas_minimas=HORAS_MINIMAS,
    horas_minimas_senal=HORAS_MINIMAS_SENAL,
    ruta_salida_csv=RUTA_HRV_V3,
    ruta_log_errores=RUTA_ERRORES_V3,
    ruta_checkpoint=RUTA_CHECKPOINT_V3,
    ruta_reporte_txt=RUTA_REPORTE_V3,
    flush_cada_n_pacientes=5, flush_cada_n_segundos=180,
):
    """One run per patient, up to `horas_objetivo_tope` hours of signal (lead II).
    Checkpoint per patient: if Colab disconnects, the run resumes from the next
    pending patient without reprocessing completed ones. Results are buffered and
    flushed to disk every few patients or seconds."""
    RUTA_GUARDADO = ruta_salida_csv
    RUTA_LOG_ERRORES = ruta_log_errores
    RUTA_CHECKPOINT = ruta_checkpoint
    RUTA_REPORTE_TXT = ruta_reporte_txt

    def registrar_error(sub_id, icu_id, etapa, detalle):
        hay_header = not os.path.exists(RUTA_LOG_ERRORES)
        fila = pd.DataFrame([{'fecha': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                              'subject_id': sub_id, 'icustay_id': icu_id, 'etapa': etapa, 'detalle': detalle}])
        with open(RUTA_LOG_ERRORES, 'a', newline='', encoding='utf-8') as f:
            fila.to_csv(f, header=hay_header, index=False)
            f.flush(); os.fsync(f.fileno())

    def registrar_checkpoint(sub_id, estado):
        hay_header = not os.path.exists(RUTA_CHECKPOINT)
        fila = pd.DataFrame([{'fecha': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                              'subject_id': sub_id, 'estado': estado}])
        with open(RUTA_CHECKPOINT, 'a', newline='', encoding='utf-8') as f:
            fila.to_csv(f, header=hay_header, index=False)
            f.flush(); os.fsync(f.fileno())

    def registrar_reporte_paciente(row, horas_reales, horas_procesadas):
        t0_teorico = row.get('hrv_analysis_start_time', None)
        inicio_real = row.get('inicio_analisis_real', None)
        desfase_str = 'N/A'
        if t0_teorico is not None and inicio_real is not None:
            try:
                desfase_h = (pd.to_datetime(inicio_real) - pd.to_datetime(t0_teorico)).total_seconds() / 3600
                desfase_str = f"{desfase_h:+.2f}h"
            except Exception:
                pass
        lineas = [
            "=======================================================================",
            f"Fecha de Procesamiento : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"  subject_id        : {row['subject_id']}",
            f"  icustay_id        : {row['icustay_id']}",
            f"  tipo_t0           : {row.get('tipo_t0', 'N/A')}",
            f"  t0_sepsis_teorico : {t0_teorico if t0_teorico is not None else 'N/A'}",
            f"  inicio_analisis   : {inicio_real if inicio_real is not None else 'N/A'}",
            f"  desfase_vs_t0     : {desfase_str}",
            f"  horas_disponibles : {horas_reales:.2f}h",
            f"  horas_procesadas  : {horas_procesadas}h (tope={horas_objetivo_tope}h)",
            "========================================================================\n"
        ]
        with open(RUTA_REPORTE_TXT, 'a', encoding='utf-8') as f:
            f.write('\n'.join(lineas))
            f.flush(); os.fsync(f.fileno())

    def vaciar_buffer(buffer, ruta, es_primer_lote):
        if not buffer:
            return es_primer_lote
        modo = 'w' if es_primer_lote else 'a'
        with open(ruta, modo, newline='', encoding='utf-8') as f:
            pd.DataFrame(buffer).to_csv(f, header=es_primer_lote, index=False)
            f.flush(); os.fsync(f.fileno())
        return False

    def flush_todo(buffer, primer_lote, motivo=''):
        if buffer:
            tqdm.write(f"  [FLUSH - {motivo}] {len(buffer)} filas")
        primer_lote = vaciar_buffer(buffer, RUTA_GUARDADO, primer_lote)
        buffer.clear()
        return primer_lote

    if os.path.exists(RUTA_CHECKPOINT):
        df_ckpt = pd.read_csv(RUTA_CHECKPOINT)
        sujetos_procesados = set(df_ckpt[df_ckpt['estado'].isin({'exitoso', 'fallido'})]['subject_id'].tolist())
        tqdm.write(f"Retomando: {len(sujetos_procesados)} pacientes ya procesados.")
        tqdm.write(f"  Exitosos: {df_ckpt['estado'].value_counts().get('exitoso', 0)} | "
                   f"Fallidos: {df_ckpt['estado'].value_counts().get('fallido', 0)}")
    else:
        sujetos_procesados = set()
        tqdm.write("Sin checkpoint previo. Iniciando desde cero.")

    df_pendiente = df_reporte_final[
        (df_reporte_final['apto_extraccion']) &
        (~df_reporte_final['subject_id'].isin(sujetos_procesados))
    ].copy()
    tqdm.write(f"Pacientes pendientes: {len(df_pendiente)}")
    for tipo, n in df_pendiente['tipo_t0'].value_counts().items():
        tqdm.write(f"  {tipo}: {n}")

    buffer, n_exitosos, n_fallidos = [], 0, 0
    primer_lote = not os.path.exists(RUTA_GUARDADO)
    pacientes_en_buffer = 0
    ultimo_flush_time = time.time()

    barra = tqdm(df_pendiente.iterrows(), total=len(df_pendiente), desc="HRV v3 (hasta 96h)")

    for _, row in barra:
        sub_id = int(row['subject_id'])
        icu_id = int(row['icustay_id'])
        tipo_t0 = row.get('tipo_t0', 'N/A')
        estado_final = 'fallido'

        try:
            partes = row['ruta_registro'].split('/')
            record_name = partes[-1]
            pn_dir = f"mimic3wdb-matched/1.0/{'/'.join(partes[:-1])}"
            header = _wfdb_con_reintentos(wfdb.rdheader, record_name, pn_dir=pn_dir)
            fs = int(header.fs)

            offset_segundos = (row['inicio_analisis_real'] - row['waveform_start']).total_seconds()
            offset = int(offset_segundos * fs)

            if offset >= header.sig_len:
                registrar_error(sub_id, icu_id, 'offset_invalido', f'offset={offset} >= sig_len={header.sig_len}')
                registrar_checkpoint(sub_id, 'fallido')
                tqdm.write(f"  [{sub_id}] ({tipo_t0}) Offset supera señal disponible → SKIP")
                n_fallidos += 1
                barra.set_postfix(ok=n_exitosos, err=n_fallidos)
                continue

            # v3: actual signal from the anchor (if offset < 0 the download starts at 0)
            horas_reales = (header.sig_len - max(0, offset)) / (fs * 3600)
            if horas_reales < horas_minimas_senal:
                registrar_error(sub_id, icu_id, 'menos_12h_senal',
                                f'horas_reales={horas_reales:.1f} < {horas_minimas_senal} '
                                f'(no puede formar ninguna ventana de 12h)')
                registrar_checkpoint(sub_id, 'fallido')
                tqdm.write(f"  [{sub_id}] ({tipo_t0}) Solo {horas_reales:.1f}h de señal → SKIP")
                n_fallidos += 1
                barra.set_postfix(ok=n_exitosos, err=n_fallidos)
                continue

            muestras_objetivo = int(horas_objetivo_tope * 3600 + 900)      # cap + 15 min margin (in seconds; x fs below)
            sampto_real = min(max(0, offset) + muestras_objetivo * fs, header.sig_len)

            record = _wfdb_con_reintentos(wfdb.rdrecord, record_name, pn_dir=pn_dir,
                                          sampfrom=max(0, offset), sampto=sampto_real)
            horas_descargadas = (sampto_real - max(0, offset)) / (fs * 3600)

            # ECG lead II only
            ecg = None
            for i, name in enumerate(record.sig_name):
                if re.search(r'\bII\b', name.upper()):
                    ecg = record.p_signal[:, i]
                    break

            if ecg is None:
                registrar_error(sub_id, icu_id, 'extraccion_ecg', 'Canal II no encontrado')
                registrar_checkpoint(sub_id, 'fallido')
                tqdm.write(f"  [{sub_id}] ({tipo_t0}) Canal II no encontrado")
                n_fallidos += 1
                barra.set_postfix(ok=n_exitosos, err=n_fallidos)
                continue

            resultado, total_valid_hours = procesar_bloque_longitudinal(
                ecg, fs=fs, horas_objetivo=horas_objetivo_tope, horas_minimas=horas_minimas, verbose=False)

            if resultado is None or resultado.empty:
                registrar_error(sub_id, icu_id, 'hrv_total',
                                f'Sin HRV válido (descargadas={horas_descargadas:.1f}h). Valid hours: {total_valid_hours}')
                registrar_checkpoint(sub_id, 'fallido')
                tqdm.write(f"  [{sub_id}] ({tipo_t0}) Sin HRV válido | descargadas={horas_descargadas:.1f}h | "
                           f"valid={total_valid_hours}")
                n_fallidos += 1
                barra.set_postfix(ok=n_exitosos, err=n_fallidos)
                continue

            resultado['subject_id'] = sub_id
            resultado['icustay_id'] = icu_id
            n_filas = len(resultado)

            if n_filas >= horas_minimas:
                buffer.extend(resultado.to_dict('records'))
                n_exitosos += 1
                estado_final = 'exitoso'
                registrar_reporte_paciente(row, horas_reales, n_filas)
                pacientes_en_buffer += 1
            else:
                registrar_error(sub_id, icu_id, 'hrv_criterio_horas',
                                f'Solo {n_filas}h válidas, mínimo requerido {horas_minimas}h.')
                n_fallidos += 1

            registrar_checkpoint(sub_id, estado_final)
            tqdm.write(f"  [{sub_id}] ({tipo_t0}): "
                       f"{'OK(' + str(n_filas) + 'h)' if estado_final == 'exitoso' else 'FAIL(' + str(n_filas) + 'h)'}")

        except Exception as e:
            registrar_error(sub_id, icu_id, 'runtime_error', str(e))
            n_fallidos += 1
            registrar_checkpoint(sub_id, 'fallido')
            tqdm.write(f"  [{sub_id}] ({tipo_t0}) EXCEPCIÓN: {str(e)[:120]}")

        tiempo_desde_flush = time.time() - ultimo_flush_time
        if pacientes_en_buffer >= flush_cada_n_pacientes or tiempo_desde_flush >= flush_cada_n_segundos:
            motivo = (f"lote={pacientes_en_buffer}" if pacientes_en_buffer >= flush_cada_n_pacientes
                      else f"tiempo={int(tiempo_desde_flush)}s")
            primer_lote = flush_todo(buffer, primer_lote, motivo)
            pacientes_en_buffer = 0
            ultimo_flush_time = time.time()

        barra.set_postfix(ok=n_exitosos, err=n_fallidos)

    tqdm.write("Vaciando buffer final...")
    flush_todo(buffer, primer_lote, 'final')

    tqdm.write(f"\n{'=' * 60}")
    tqdm.write("PROCESO FINALIZADO (extracción v3 hasta 96h)")
    tqdm.write(f"  Exitosos: {n_exitosos} | Fallidos: {n_fallidos}")

    if os.path.exists(RUTA_GUARDADO):
        df_res = pd.read_csv(RUTA_GUARDADO)
        tqdm.write("\nHoras por paciente (CSV final):")
        tqdm.write(df_res.groupby('subject_id')['hour_index'].count().describe().to_string())
    if os.path.exists(RUTA_LOG_ERRORES):
        df_err = pd.read_csv(RUTA_LOG_ERRORES)
        tqdm.write("\nErrores por etapa:")
        tqdm.write(df_err['etapa'].value_counts().to_string())


# ═══════════════════════════════════════════════════════════════════════════
# 7. QUALITY CONTROL v3
#    Run once 10-20 patients have been processed (stop the extraction, run this
#    function and resume: the checkpoint continues where it stopped).
# ═══════════════════════════════════════════════════════════════════════════
def control_calidad_v3(ruta=RUTA_HRV_V3):
    """Sanity checks against the known failure of the old detector:
    rho1 = 1 - (RMSSD/SDNN)^2 / 2 should be positive, SD1 > SD2 should be rare,
    DFA alpha1 should be > 0.6 and no SampEn should be infinite or negative."""
    v3 = pd.read_csv(ruta)
    rho1 = 1 - (v3.RMSSD_T0 / v3.SDNN_T0) ** 2 / 2
    print(f"pacientes: {v3.icustay_id.nunique()} | horas: {len(v3)}")
    print(f"SD1 > SD2      : {(v3.SD1_T0 > v3.SD2_T0).mean() * 100:.1f} %   (antes 68.9 %; esperado: bajo)")
    print(f"rho1 mediana   : {rho1.median():.2f}      (antes ≈ -0.28; esperado: positiva)")
    print(f"DFA a1 mediana : {v3.DFA_a1_T0.median():.2f}      (antes ≈ 0.30; esperado: > 0.6)")
    print(f"RMSSD mediana  : {v3.RMSSD_T0.median():.1f} ms  | SDNN mediana: {v3.SDNN_T0.median():.1f} ms")
    print(f"SampEn inf/neg : {(~np.isfinite(v3.SampEn_T0) | (v3.SampEn_T0 < 0)).sum()} (esperado: 0)")


# ─── EXECUTION ──────────────────────────────────────────────────────────────
correr_extraccion_hrv_unificada(df_reporte_final)

# control_calidad_v3()      # <- uncomment to run the quality control
