# -*- coding: utf-8 -*-
"""
construir_ultimate_dataset_v3.py  (dataset builder)
===================================================
Assembles ultimate_dataset_v3_1.csv (v3.1: LF/HF/LFHF on the natural-log scale)
following the Mentasti et al. sliding-window design and HRV extraction v3.

Identifiers are in Spanish; English terms are in docs/data_dictionary.md
(label_sofa_responder = SOFA improvement, label_cardio = haemodynamic stability
or recovery, label_sofa_deterioro_2pts = SOFA deterioration).

Design:
  - Unit of analysis: window (icustay_id, T0), T0 every 6 h from hour 12;
    12-h observation window (features) and 6-h horizon (labels).
  - Hourly HRV value = median of the valid 5-min epochs of the hour, rebuilt
    from the epoch-12 value and its 11 deltas. LF, HF and LFHF are
    log-transformed per epoch, before aggregation.
  - Linear HRV = SDNN, RMSSD, SD1, SD2, LF, HF, LFHF (SD1/SD2 in the time domain,
    Brennan et al. 2001); non-linear = SampEn, ApEn, DFA_a1.
  - Delta SOFA is paired per component; complete case is a sensitivity analysis.
    *_nl labels keep only windows with a new creatinine, bilirubin or platelet
    value (or death) in the horizon: sensitivity analysis.
  - Death in the horizon = no improvement and deterioration; discharge alive
    inside the horizon = excluded.
  - MAP < 20 or > 200 mmHg -> NaN; lactato_medido_12h flags a lactate in the
    12 h before T0.
  - sofa_sepsis_onset is not used (it includes data up to +24 h after the
    suspicion of infection: leakage).
Requires the tables of 03_queries_vigentes_v4.sql.

Outputs (earlier versions are not overwritten):
  - ultimate_dataset_v3_1.csv: ALL candidate windows with flags.
  - columnas_v3_1.json      : column groups and columns forbidden as features.
  - flujo_cohorte_v3_1.csv  : counts for the cohort flow diagram.

Project rule: this is the ONLY modelling script that queries BigQuery.
"""

import json
import numpy as np
import pandas as pd

# ─── Configuration ──────────────────────────────────────────────────────────
PROYECTO = 'sepsis-mimic-def'
DS = 'sepsis-mimic-def.sepsis3_workspace_def'
RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA_HRV = f'{RUTA_BASE}/hrv_completo_v3.csv'
RUTA_SALIDA = f'{RUTA_BASE}/ultimate_dataset_v3_1.csv'
RUTA_COLUMNAS = f'{RUTA_BASE}/columnas_v3_1.json'
RUTA_FLUJO = f'{RUTA_BASE}/flujo_cohorte_v3_1.csv'

METRICAS = ['SDNN', 'RMSSD', 'LF', 'HF', 'LFHF', 'SD1', 'SD2', 'SampEn', 'ApEn', 'DFA_a1']
METRICAS_LINEALES = ['SDNN', 'RMSSD', 'SD1', 'SD2', 'LF', 'HF', 'LFHF']   # SD1 and SD2: linear (Brennan et al., 2001)
METRICAS_NOLINEALES = ['SampEn', 'ApEn', 'DFA_a1']
VENTANA_OBS_H = 12
HORIZONTE_H = 6
UMBRALES_HRV = [4, 6, 8]          # 6 = primary (fixed a priori); 4 and 8 = sensitivity
MIN_PUNTOS_SLOPE = 3              # always met with the >= 6 h filter; applies to the >= 4 h sensitivity analysis

# ─── 1. Master query (one row per candidate window) ─────────────────────────
SQL_MASTER = f"""
WITH p AS (
  SELECT subject_id, icustay_id, hour_index, t0_prediccion, n_horas_validas
  FROM `{DS}.hrv_prediction_points`
),
dem AS (
  SELECT icustay_id, age, gender          -- sofa_sepsis_onset removed (leakage)
  FROM `{DS}.sepsis3_final_cohort`
),
fin AS (
  SELECT ie.icustay_id,
         CAST(ie.outtime AS DATETIME)    AS outtime,
         CAST(adm.deathtime AS DATETIME) AS deathtime
  FROM `physionet-data.mimiciii_clinical.icustays` ie
  JOIN `physionet-data.mimiciii_clinical.admissions` adm USING (hadm_id)
),
newlab AS (
  -- >= 1 new creatinine, bilirubin or platelet measurement in the horizon
  SELECT p.icustay_id, p.hour_index,
         MAX(IF(l.charttime IS NOT NULL, 1, 0)) AS new_lab_6h
  FROM p
  LEFT JOIN `{DS}._raw_labs` l
    ON p.icustay_id = l.icustay_id
    AND l.charttime > p.t0_prediccion
    AND l.charttime <= DATETIME_ADD(p.t0_prediccion, INTERVAL 6 HOUR)
    AND (l.creatinine IS NOT NULL OR l.bilirubin IS NOT NULL OR l.platelet IS NOT NULL)
  GROUP BY p.icustay_id, p.hour_index
)
SELECT
  p.subject_id, p.icustay_id, p.hour_index, p.t0_prediccion, p.n_horas_validas,
  dem.age, dem.gender,
  fin.outtime, fin.deathtime, nl.new_lab_6h,

  -- SOFA (v3.3)
  s.sofa_t0, s.sofa_t0_valid, s.sofa_t0_6h, s.sofa_t0_6h_valid,
  s.sofa_resp_t0, s.sofa_coag_t0, s.sofa_liver_t0,
  s.sofa_cardio_t0, s.sofa_neuro_t0, s.sofa_renal_t0,
  s.sofa_t0_parcial, s.n_comp_validos_t0,
  s.delta_pareado_valid, s.delta_sofa_pareado,

  -- Vital signs and drugs at T0
  v.map_t0, v.rr_t0, v.spo2_t0, v.temp_t0,
  d.ned, d.vis, d.cumvp, d.sedante_bin, d.opioide_bin,

  -- MAP and lactate (current 12 h window = feature; next 6 h window = outcome)
  c.mbp_mean_actual12h, c.n_mbp_actual12h, c.lactate_last_actual12h,
  c.mbp_mean_next6h, c.n_mbp_next6h, c.mbp_pct_next6h_gt65,
  c.mbp_incremento_pct, c.lactate_last_next6h

FROM p
LEFT JOIN dem USING (icustay_id)
LEFT JOIN fin USING (icustay_id)
LEFT JOIN `{DS}.sofa_longitudinal`            s USING (icustay_id, hour_index)
LEFT JOIN `{DS}.vitals_longitudinal`          v USING (icustay_id, hour_index)
LEFT JOIN `{DS}.drugs_longitudinal`           d USING (icustay_id, hour_index)
LEFT JOIN `{DS}.cardio_lactato_longitudinal`  c USING (icustay_id, hour_index)
LEFT JOIN newlab                               nl USING (icustay_id, hour_index)
"""


# ─── 2a. Hourly value = median of the valid epochs ──────────────────────────
def valor_horario_mediana(hrv: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the 12 five-minute epochs of each hour and replace {m}_T0 by the
    median of the valid epochs.

    Last epoch = {m}_T0; epoch N minutes earlier = {m}_T0 - Delta_{m}_minus_{N}min.
    Invalid epochs have a NaN delta and therefore drop out of the median."""
    out = hrv.copy()
    for m in METRICAS:
        ancla = hrv[f'{m}_T0'].replace([np.inf, -np.inf], np.nan)
        epocas = [ancla - hrv[f'Delta_{m}_minus_{k}min'].replace([np.inf, -np.inf], np.nan)
                  for k in range(5, 60, 5)]
        epocas.append(ancla)
        out[f'{m}_T0'] = pd.concat(epocas, axis=1).median(axis=1, skipna=True)
    return out


# ─── 2. HRV aggregation per window ──────────────────────────────────────────
def agregar_hrv(hrv: pd.DataFrame, puntos: pd.DataFrame) -> pd.DataFrame:
    """Summarise the hourly HRV rows of each 12-h observation window.

    hrv    : one row per (icustay_id, hour_index_hrv) with the {m}_T0 columns.
    puntos : one row per window (icustay_id, hour_index = hour of T0).
    The window of a T0 at hour t covers HRV hours t-12 .. t-1 (each HRV row
    represents the END of its hour; t0_prediccion = anchor + t hours).

    Per metric: mean, SD, last value and slope (OLS against the actual hour,
    in units per hour, so gaps do not distort it). Also returns
    n_horas_hrv_ventana (valid HRV hours in the window) and
    horas_desde_ultima_hrv (hours between the last valid HRV hour and T0).
    """
    cols_t0 = [f'{m}_T0' for m in METRICAS]
    h = hrv[['icustay_id', 'hour_index'] + cols_t0].rename(columns={'hour_index': 'h'})
    w = puntos[['icustay_id', 'hour_index']].drop_duplicates()

    x = w.merge(h, on='icustay_id', how='inner')
    x = x[(x['h'] >= x['hour_index'] - VENTANA_OBS_H) & (x['h'] < x['hour_index'])].copy()
    x['t'] = x['h'] - (x['hour_index'] - VENTANA_OBS_H)          # 0..11 within the window
    x = x.sort_values(['icustay_id', 'hour_index', 't'])
    g = x.groupby(['icustay_id', 'hour_index'], sort=False)

    out = pd.DataFrame(index=g.size().index)
    out['n_horas_hrv_ventana'] = g.size()
    out['horas_desde_ultima_hrv'] = (VENTANA_OBS_H - 1) - g['t'].max()

    for m, col in zip(METRICAS, cols_t0):
        y = x[col]
        ok = y.notna()
        xt = x['t'].where(ok)
        tmp = pd.DataFrame({
            'icustay_id': x['icustay_id'], 'hour_index': x['hour_index'],
            'n': ok.astype(float), 'sx': xt, 'sy': y.where(ok),
            'sxy': xt * y.where(ok), 'sxx': xt ** 2,
        }).groupby(['icustay_id', 'hour_index'], sort=False).sum(min_count=1)
        den = tmp['n'] * tmp['sxx'] - tmp['sx'] ** 2
        slope = (tmp['n'] * tmp['sxy'] - tmp['sx'] * tmp['sy']) / den.replace(0, np.nan)
        slope[tmp['n'] < MIN_PUNTOS_SLOPE] = np.nan

        out[f'{m}_mean'] = g[col].mean()
        out[f'{m}_sd'] = g[col].std(ddof=1)
        out[f'{m}_last'] = g[col].last()        # last non-null value in the window
        out[f'{m}_slope'] = slope

    return out.reset_index()


# ─── 3. Labels, death and discharge in the horizon ──────────────────────────
def construir_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Create the outcome labels, the horizon flags and the inclusion filters.

    Label values are 1.0 / 0.0, or NaN when the window is not evaluable.
      label_sofa_responder      : SOFA improvement (paired delta SOFA < 0)      [main]
      label_cardio              : haemodynamic stability or recovery            [main]
      label_sofa_deterioro_2pts : SOFA deterioration (paired delta SOFA >= 2)     [secondary]
      *_cc                      : complete-case SOFA versions                   [sensitivity]
      *_nl                      : SOFA labels restricted to windows with a new lab [sensitivity]
    """
    df = df.copy()
    t0 = pd.to_datetime(df['t0_prediccion'])
    fin_h = t0 + pd.Timedelta(hours=HORIZONTE_H)
    deathtime = pd.to_datetime(df['deathtime'])
    outtime = pd.to_datetime(df['outtime'])

    df['muere_en_horizonte'] = (deathtime.notna() & (deathtime > t0) & (deathtime <= fin_h)).astype(int)
    df['alta_viva_en_horizonte'] = ((outtime < fin_h) & (df['muere_en_horizonte'] == 0)).astype(int)
    muere = df['muere_en_horizonte'] == 1
    alta = df['alta_viva_en_horizonte'] == 1

    # SOFA - paired delta (main)
    d = df['delta_sofa_pareado'].where(df['delta_pareado_valid'] == 1)
    df['label_sofa_responder'] = np.where(d.notna(), (d < 0).astype(float), np.nan)
    df['label_sofa_deterioro_2pts'] = np.where(d.notna(), (d >= 2).astype(float), np.nan)

    # SOFA - complete case (sensitivity)
    cc_ok = (df['sofa_t0_valid'] == 1) & (df['sofa_t0_6h_valid'] == 1)
    dcc = (df['sofa_t0_6h'] - df['sofa_t0']).where(cc_ok)
    df['label_sofa_responder_cc'] = np.where(dcc.notna(), (dcc < 0).astype(float), np.nan)
    df['label_sofa_deterioro_2pts_cc'] = np.where(dcc.notna(), (dcc >= 2).astype(float), np.nan)

    # Cardiovascular / haemodynamic (Mentasti). Evaluable only if MAP and lactate
    # exist in both the current (12 h) and the next (6 h) window.
    # MAP criterion  : > 50 % of next-window MAP values above 65 mmHg, or MAP up >= 10 %.
    # Lactate criterion: < 2 in both windows; 2-4 -> < 2; > 4 -> <= 4 (mmol/L).
    la, ln = df['lactate_last_actual12h'], df['lactate_last_next6h']
    cardio_eval = (df['n_mbp_actual12h'].fillna(0) > 0) & (df['n_mbp_next6h'].fillna(0) > 0) \
                  & la.notna() & ln.notna()
    mbp_ok = (df['mbp_pct_next6h_gt65'] > 0.5) | (df['mbp_incremento_pct'] >= 0.10)
    lac_ok = ((la < 2) & (ln < 2)) | (la.between(2, 4) & (ln < 2)) | ((la > 4) & (ln <= 4))
    df['label_cardio'] = np.where(cardio_eval, (mbp_ok & lac_ok).astype(float), np.nan)

    # Death in the horizon = worst outcome: no improvement and deterioration (all versions)
    for c in ['label_sofa_responder', 'label_sofa_responder_cc', 'label_cardio']:
        df.loc[muere, c] = 0.0
    for c in ['label_sofa_deterioro_2pts', 'label_sofa_deterioro_2pts_cc']:
        df.loc[muere, c] = 1.0

    # Discharge alive in the horizon: the label does not reflect those 6 h -> exclude
    for c in ['label_sofa_responder', 'label_sofa_deterioro_2pts', 'label_sofa_responder_cc',
              'label_sofa_deterioro_2pts_cc', 'label_cardio']:
        df.loc[alta, c] = np.nan

    # HRV inclusion filters (primary and sensitivity)
    for u in UMBRALES_HRV:
        df[f'incl_hrv_{u}h'] = (df['n_horas_validas'] >= u).astype(int)

    # Sensitivity: only windows with a new laboratory measurement in the horizon
    # (or death, whose outcome is observed)
    nl_ok = (df['new_lab_6h'].fillna(0) == 1) | muere
    for c in ['label_sofa_responder', 'label_sofa_deterioro_2pts']:
        df[f'{c}_nl'] = df[c].where(nl_ok)

    df['eval_comun'] = (df['label_sofa_responder'].notna() & df['label_cardio'].notna()).astype(int)
    return df


# ─── 4. Column groups for modelling ─────────────────────────────────────────
def grupos_columnas() -> dict:
    """Column groups written to columnas_v3_1.json and read by the modelling
    pipeline. 'desenlace_no_usar_como_feature' lists columns that are measured
    after T0 or that define a label; they must never be predictors."""
    hrv = [f'{m}_{s}' for m in METRICAS for s in ['mean', 'sd', 'last', 'slope']] + ['horas_desde_ultima_hrv']
    return {
        'id': ['subject_id', 'icustay_id', 'hour_index', 't0_prediccion'],
        'features_hrv': hrv,
        'features_hrv_lineal': [f'{m}_{s}' for m in METRICAS_LINEALES for s in ['mean', 'sd', 'last', 'slope']],
        'features_hrv_nolineal': [f'{m}_{s}' for m in METRICAS_NOLINEALES for s in ['mean', 'sd', 'last', 'slope']],
        'features_clinicas': ['age', 'gender', 'map_t0', 'rr_t0', 'spo2_t0',
                              'ned', 'vis', 'cumvp', 'sedante_bin', 'opioide_bin'],
        'features_sofa': ['sofa_t0_parcial', 'n_comp_validos_t0'],
        'features_pam_lactato_actual': ['mbp_mean_actual12h', 'lactate_last_actual12h', 'lactato_medido_12h'],
        'labels_principales': ['label_sofa_responder', 'label_cardio'],
        'labels_secundarios': ['label_sofa_deterioro_2pts'],
        'labels_sensibilidad': ['label_sofa_responder_cc', 'label_sofa_deterioro_2pts_cc',
                                'label_sofa_responder_nl', 'label_sofa_deterioro_2pts_nl'],
        'filtros': [f'incl_hrv_{u}h' for u in UMBRALES_HRV] + ['eval_comun'],
        # NEVER as a feature: measured after T0 or defining the label
        'desenlace_no_usar_como_feature': [
            'sofa_t0_6h', 'sofa_t0_6h_valid', 'delta_sofa_pareado', 'delta_pareado_valid',
            'mbp_mean_next6h', 'n_mbp_next6h', 'mbp_pct_next6h_gt65', 'mbp_incremento_pct',
            'lactate_last_next6h', 'muere_en_horizonte', 'alta_viva_en_horizonte',
            'outtime', 'deathtime', 'new_lab_6h'],
    }


# ─── 5. Cohort summary ──────────────────────────────────────────────────────
def resumen(df: pd.DataFrame) -> pd.DataFrame:
    """Window and patient counts at each step, for the cohort flow diagram."""
    filas = []
    def add(paso, sub):
        filas.append({'paso': paso, 'ventanas': len(sub), 'pacientes': sub['icustay_id'].nunique()})
    add('Ventanas candidatas (≥1 ventana de 12h)', df)
    p = df[df['incl_hrv_6h'] == 1]
    add('Con ≥6h de HRV válida', p)
    add('  Muerte en horizonte (forzadas)', p[p['muere_en_horizonte'] == 1])
    add('  Alta viva en horizonte (excluidas)', p[p['alta_viva_en_horizonte'] == 1])
    for lab in ['label_sofa_responder', 'label_cardio', 'label_sofa_deterioro_2pts',
                'label_sofa_responder_cc', 'label_sofa_responder_nl']:
        s = p[p[lab].notna()]
        add(f'Evaluable {lab}', s)
        filas[-1]['prevalencia_%'] = round(s[lab].mean() * 100, 1) if len(s) else np.nan
    add('Cohorte común (SOFA + cardio)', p[p['eval_comun'] == 1])
    s = p[p['label_sofa_responder'].notna()]
    add('  SOFA evaluable sin medición nueva de laboratorio', s[s['new_lab_6h'].fillna(0) == 0])
    return pd.DataFrame(filas)


# ─── 6. Execution in Colab ──────────────────────────────────────────────────
def main():
    """Query BigQuery, aggregate HRV, build labels and write the three outputs."""
    from google.colab import auth, drive
    drive.mount('/content/drive')
    auth.authenticate_user()
    from google.cloud import bigquery
    client = bigquery.Client(project=PROYECTO)

    print('Consultando BigQuery...')
    master = client.query(SQL_MASTER).to_dataframe()
    for c in ['subject_id', 'icustay_id', 'hour_index', 'n_horas_validas']:
        master[c] = master[c].astype('int64')
    assert not master.duplicated(['icustay_id', 'hour_index']).any(), 'Ventanas duplicadas en master'
    print(f'  master: {len(master)} ventanas, {master.icustay_id.nunique()} pacientes')

    print('Agregando HRV por ventana...')
    hrv = pd.read_csv(RUTA_HRV)
    hrv['icustay_id'] = hrv['icustay_id'].astype('int64')
    hrv['hour_index'] = hrv['hour_index'].astype('int64')
    hrv = valor_horario_mediana(hrv)
    # Change 7: natural-log scale for the spectral metrics (strongly skewed distribution)
    for m in ['LF', 'HF', 'LFHF']:
        hrv[f'{m}_T0'] = np.log(hrv[f'{m}_T0'].clip(lower=1e-6))
    agg = agregar_hrv(hrv, master)

    df = master.merge(agg, on=['icustay_id', 'hour_index'], how='left', validate='one_to_one')
    # Consistency: the hours counted here must match the hours in the uploaded grid
    dif = (df['n_horas_hrv_ventana'].fillna(0).astype(int) != df['n_horas_validas']).sum()
    assert dif == 0, f'{dif} ventanas con n_horas distinto a la grilla (revisar convención de horas)'

    df = construir_labels(df)
    # Change 8: physiologically implausible MAP values -> missing
    for c in ['map_t0', 'mbp_mean_actual12h']:
        df.loc[(df[c] < 20) | (df[c] > 200), c] = np.nan
    # Change 9: lactate-measured indicator (absence of a measurement is informative)
    df['lactato_medido_12h'] = df['lactate_last_actual12h'].notna().astype(int)
    # v3.1 control: LF/HF/LFHF must be on the ln scale (negative values exist)
    for _m in ['LF', 'HF', 'LFHF']:
        assert (df[f'{_m}_mean'].dropna() < 0).any(), f'{_m}_mean no está en ln'
    print('✅ LF, HF y LFHF en ln:', {m: round(df[f"{m}_mean"].median(), 3) for m in ['LF', 'HF', 'LFHF']})
    df.to_csv(RUTA_SALIDA, index=False)
    with open(RUTA_COLUMNAS, 'w') as f:
        json.dump(grupos_columnas(), f, indent=2, ensure_ascii=False)
    flujo = resumen(df)
    flujo.to_csv(RUTA_FLUJO, index=False)

    print(f'\n✅ {RUTA_SALIDA}: {df.shape}')
    print(flujo.to_string(index=False))


if __name__ == '__main__':
    main()
