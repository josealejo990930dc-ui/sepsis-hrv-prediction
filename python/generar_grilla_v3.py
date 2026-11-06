# -*- coding: utf-8 -*-
"""
generar_grilla_v3.py  (prediction-window grid)
==============================================
Builds the grid of prediction points from hrv_completo_v3.csv and uploads it to
BigQuery as hrv_prediction_points (replacing the previous grid).

Grid (Mentasti design):
  - T0 at hours 12, 18, 24... (<= 96 and <= last valid HRV hour + 1).
  - Observation window of a T0 at hour t = HRV rows with hour_index in [t-12, t-1]
    (each HRV row represents the END of its hour).
  - t0_prediccion = inicio_analisis_real + t hours.
  - ALL candidate windows are uploaded with n_horas_validas; the >= 6 h filter
    (and the >= 4 h / >= 8 h sensitivity filters) is applied in the dataset builder.

Before uploading, the previous grid is backed up as hrv_prediction_points_v2_backup
(only the first time; an existing backup is never overwritten).

After running this script: queries_vigentes_v4.sql, part B.
"""

import pandas as pd

PROYECTO = 'sepsis-mimic-def'
DS = 'sepsis-mimic-def.sepsis3_workspace_def'
RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA_HRV = f'{RUTA_BASE}/hrv_completo_v3.csv'
RUTA_REPORTE = f'{RUTA_BASE}/df_reporte_final_unificado.csv'
RUTA_VENTANAS = f'{RUTA_BASE}/ventanas_candidatas_v3.csv'
TABLA = f'{DS}.hrv_prediction_points'
RESPALDO = f'{DS}.hrv_prediction_points_v2_backup'

PRIMER_T0, PASO_H, VENTANA_OBS_H, TOPE_H = 12, 6, 12, 96   # first T0 (h), step (h), observation window (h), cap (h)


def construir_grilla(hrv: pd.DataFrame) -> pd.DataFrame:
    """One row per candidate window: icustay_id, t0_h, n_horas_validas
    (number of valid HRV hours inside the 12-h observation window)."""
    filas = []
    for icu, h in hrv.groupby('icustay_id')['hour_index']:
        horas = set(int(x) for x in h)
        ultima = max(horas)
        for t0 in range(PRIMER_T0, min(ultima + 2, TOPE_H + 1), PASO_H):
            n = sum(t0 - VENTANA_OBS_H <= x < t0 for x in horas)
            filas.append((int(icu), t0, n))
    return pd.DataFrame(filas, columns=['icustay_id', 't0_h', 'n_horas_validas'])


def resumen_grilla(hrv: pd.DataFrame, win: pd.DataFrame) -> None:
    """Print patient and window counts for each HRV-hours threshold."""
    n_pac = hrv['icustay_id'].nunique()
    print(f"Pacientes en HRV v3            : {n_pac}")
    print(f"Pacientes sin ninguna ventana  : {n_pac - win['icustay_id'].nunique()}")
    print(f"Ventanas candidatas            : {len(win)}")
    for u in [4, 6, 8]:
        w = win[win['n_horas_validas'] >= u]
        print(f"  ≥{u} h válidas: {len(w)} ventanas, {w['icustay_id'].nunique()} pacientes")


def main():
    """Build the grid, attach the time anchor and upload it to BigQuery."""
    from google.colab import auth, drive
    drive.mount('/content/drive')
    auth.authenticate_user()
    from google.cloud import bigquery
    client = bigquery.Client(project=PROYECTO)

    # ── 1. HRV v3 and structure checks ──────────────────────────────────────
    hrv = pd.read_csv(RUTA_HRV)
    hrv['icustay_id'] = hrv['icustay_id'].astype('int64')
    hrv['hour_index'] = hrv['hour_index'].astype('int64')
    assert not hrv.duplicated(['icustay_id', 'hour_index']).any(), 'Horas duplicadas en el CSV v3'
    g = hrv.groupby('icustay_id')
    assert (g.size() == g['horas_validas_total'].first()).all(), 'Filas ≠ horas_validas_total'

    # ── 2. Grid ─────────────────────────────────────────────────────────────
    win = construir_grilla(hrv)
    resumen_grilla(hrv, win)
    win.to_csv(RUTA_VENTANAS, index=False)

    # ── 3. Time anchor ──────────────────────────────────────────────────────
    meta = pd.read_csv(RUTA_REPORTE, usecols=['subject_id', 'icustay_id', 'inicio_analisis_real'])
    meta = meta.dropna(subset=['icustay_id', 'inicio_analisis_real'])
    meta['icustay_id'] = meta['icustay_id'].astype('int64')
    meta['subject_id'] = meta['subject_id'].astype('int64')
    meta['inicio_analisis_real'] = pd.to_datetime(meta['inicio_analisis_real']).dt.tz_localize(None)
    assert meta['icustay_id'].is_unique, 'Ancla duplicada por icustay_id'

    puntos = win.merge(meta, on='icustay_id', how='left', validate='many_to_one')
    sin_ancla = puntos['inicio_analisis_real'].isna().sum()
    assert sin_ancla == 0, f'{sin_ancla} ventanas sin ancla: revisar df_reporte_final_unificado.csv'
    puntos['t0_prediccion'] = puntos['inicio_analisis_real'] + pd.to_timedelta(puntos['t0_h'], unit='h')
    puntos = puntos.rename(columns={'t0_h': 'hour_index'})
    puntos = puntos[['subject_id', 'icustay_id', 'hour_index', 't0_prediccion', 'n_horas_validas']]

    # ── 4. Back up the previous grid (only if no backup exists) ─────────────
    client.query(f"CREATE TABLE IF NOT EXISTS `{RESPALDO}` AS SELECT * FROM `{TABLA}`").result()
    print(f"Respaldo: {RESPALDO} ({client.get_table(RESPALDO).num_rows} filas)")

    # ── 5. Upload ───────────────────────────────────────────────────────────
    schema = [
        bigquery.SchemaField('subject_id', 'INTEGER'),
        bigquery.SchemaField('icustay_id', 'INTEGER'),
        bigquery.SchemaField('hour_index', 'INTEGER'),
        bigquery.SchemaField('t0_prediccion', 'DATETIME'),
        bigquery.SchemaField('n_horas_validas', 'INTEGER'),
    ]
    client.load_table_from_dataframe(
        puntos, TABLA,
        job_config=bigquery.LoadJobConfig(schema=schema, write_disposition='WRITE_TRUNCATE'),
    ).result()
    subidas = client.get_table(TABLA).num_rows
    assert subidas == len(puntos), f'Se subieron {subidas} filas de {len(puntos)}'
    print(f"✅ {TABLA}: {subidas} ventanas, {puntos['icustay_id'].nunique()} pacientes")
    print("Siguiente paso: queries_vigentes_v4.sql, parte B (empezando por B0).")


if __name__ == '__main__':
    main()
