# -*- coding: utf-8 -*-
"""
bootstrap_v5.py  (patient-level bootstrap comparisons)
======================================================
Compares models from the out-of-fold predictions of pipeline_modelado_v5.py,
using a patient-level bootstrap.

Identifiers are in Spanish; English terms are in docs/data_dictionary.md.

Contrasts (new model - reference model, same model family)
----------------------------------------------------------
- H1 (main outcomes, AUROC):  integrada_sofa - clinico_sofa
        (integrated HRV + SOFA minus clinical + SOFA).
- H2 (main outcomes, AUROC):  hrv_nolineal - hrv_lineal
        (non-linear HRV minus linear HRV).
- H3 (progression, AUPRC):    integrada_sofa - clinico_sofa.
- Complementary (the outcome's own metric): hrv_integrada - clinico_basal
        (integrated HRV + clinical minus clinical, both without SOFA).
Sensitivity analyses repeat H1 and H3 with their own predictions.
The one-sided p-values written here are re-expressed as two-sided p-values for
H2 by h2_bilateral_v5_1.py (H2 is bilateral in the research proposal).

Method
------
In every replicate PATIENTS are resampled with replacement, keeping all their
windows, and the metric of both models is recomputed on the same rows.
- Point difference: on the original data.
- 95 % CI: 2.5th and 97.5th percentiles of the bootstrap differences.
- One-sided p: share of replicates WITHOUT improvement (delta <= 0 for AUROC and
  AUPRC; delta >= 0 for Brier, where lower is better).
- Benjamini-Hochberg correction across families, within each
  (analysis, hypothesis, outcome), on the main metric.
- Significant improvement: 95 % CI excluding 0 in the direction of improvement
  and corrected p < 0.05.

Output: modelado_v5_1/comparaciones_bootstrap.csv
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA_SALIDA = f'{RUTA_BASE}/modelado_v5_1'
N_REPLICAS = 2000
SEMILLA = 2026

CONTRASTES = [
    # hypothesis, outcomes, reference, new model, main metric (None = the outcome's own)
    ('H1', ['sofa_respuesta', 'hemodinamico'], 'clinico_sofa', 'integrada_sofa', 'auroc'),
    ('H2', ['sofa_respuesta', 'hemodinamico'], 'hrv_lineal', 'hrv_nolineal', 'auroc'),
    ('H3', ['progresion'], 'clinico_sofa', 'integrada_sofa', 'auprc'),
    ('Complementario', ['sofa_respuesta', 'hemodinamico', 'progresion'], 'clinico_basal', 'hrv_integrada', None),
]
METRICA_DESENLACE = {'sofa_respuesta': 'auroc', 'hemodinamico': 'auroc', 'progresion': 'auprc'}


def metrica(nombre, y, p_raw, p_cal):
    """AUROC / AUPRC on the uncalibrated probability, Brier on the calibrated one.
    NaN if the sample contains a single class."""
    if len(np.unique(y)) < 2:
        return np.nan
    if nombre == 'auroc':
        return roc_auc_score(y, p_raw)
    if nombre == 'auprc':
        return average_precision_score(y, p_raw)
    return brier_score_loss(y, p_cal)


def benjamini_hochberg(p):
    """Benjamini-Hochberg adjusted p-values (q-values); NaNs are left as NaN."""
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    if ok.sum() == 0:
        return q
    pv = p[ok]; n = len(pv)
    orden = np.argsort(pv)
    ajust = pv[orden] * n / np.arange(1, n + 1)
    ajust = np.minimum.accumulate(ajust[::-1])[::-1]
    res = np.empty(n); res[orden] = np.minimum(ajust, 1.0)
    q[ok] = res
    return q


def bootstrap_par(a, b, rng):
    """Patient-level paired bootstrap of the metric difference (a - b).

    a, b: predictions of the two models, aligned by window (a = new model,
    b = reference). Returns (point values {metric: (a, b)}, bootstrap deltas {metric: array})."""
    pacientes = a['icustay_id'].to_numpy()
    unicos = np.unique(pacientes)
    filas_de = {pid: np.flatnonzero(pacientes == pid) for pid in unicos}
    y = a['y'].to_numpy()
    pa_r, pa_c = a['p_sin_calibrar'].to_numpy(), a['p_calibrada'].to_numpy()
    pb_r, pb_c = b['p_sin_calibrar'].to_numpy(), b['p_calibrada'].to_numpy()

    puntual = {m: (metrica(m, y, pa_r, pa_c), metrica(m, y, pb_r, pb_c)) for m in ['auroc', 'auprc', 'brier']}
    deltas = {m: [] for m in puntual}
    for _ in range(N_REPLICAS):
        idx = np.concatenate([filas_de[pid] for pid in rng.choice(unicos, size=len(unicos), replace=True)])
        for m in deltas:
            deltas[m].append(metrica(m, y[idx], pa_r[idx], pa_c[idx]) - metrica(m, y[idx], pb_r[idx], pb_c[idx]))
    return puntual, {m: np.asarray(v, dtype=float) for m, v in deltas.items()}


def main():
    """Run every contrast for every analysis present in predicciones_oof.csv and
    write comparaciones_bootstrap.csv (one row per contrast x family x metric)."""
    try:
        from google.colab import drive
        drive.mount('/content/drive')
    except ImportError:
        pass
    pred = pd.read_csv(f'{RUTA_SALIDA}/predicciones_oof.csv')
    rng = np.random.default_rng(SEMILLA)
    clave = ['icustay_id', 'hour_index']
    filas = []

    for an in pred['analisis'].unique():
        pa = pred[pred['analisis'] == an]
        for hip, desenlaces, ref, nuevo, met_principal in CONTRASTES:
            if an != 'principal' and hip not in ('H1', 'H3'):
                continue
            for des in desenlaces:
                met_p = met_principal or METRICA_DESENLACE[des]
                pd_des = pa[pa['desenlace'] == des]
                familias = sorted(set(pd_des.loc[pd_des.conjunto == ref, 'familia'])
                                  & set(pd_des.loc[pd_des.conjunto == nuevo, 'familia']))
                for fam in familias:
                    b = pd_des[(pd_des.conjunto == ref) & (pd_des.familia == fam)].sort_values(clave)
                    a = pd_des[(pd_des.conjunto == nuevo) & (pd_des.familia == fam)].sort_values(clave)
                    a = a.merge(b[clave], on=clave)            # same windows in both models
                    b = b.merge(a[clave], on=clave)
                    assert (a['y'].to_numpy() == b['y'].to_numpy()).all(), 'Etiquetas no coinciden'
                    puntual, deltas = bootstrap_par(a.reset_index(drop=True), b.reset_index(drop=True), rng)
                    for m in ['auroc', 'auprc', 'brier']:
                        d = deltas[m][~np.isnan(deltas[m])]
                        mejora = d >= 0 if m == 'brier' else d <= 0   # replicates WITHOUT improvement
                        filas.append({
                            'analisis': an, 'hipotesis': hip, 'desenlace': des, 'familia': fam,
                            'referencia': ref, 'nuevo': nuevo, 'metrica': m, 'principal': m == met_p,
                            'valor_referencia': puntual[m][1], 'valor_nuevo': puntual[m][0],
                            'delta': puntual[m][0] - puntual[m][1],
                            'ic95_inf': np.percentile(d, 2.5), 'ic95_sup': np.percentile(d, 97.5),
                            'p_unilateral': mejora.mean() if len(d) else np.nan,
                            'replicas_validas': len(d), 'ventanas': len(a), 'pacientes': a.icustay_id.nunique()})
                print(f'{an} · {hip} · {des}: {len(familias)} familias')

    res = pd.DataFrame(filas)
    res['p_bh'] = np.nan
    prin = res['principal']
    for _, idx in res[prin].groupby(['analisis', 'hipotesis', 'desenlace']).groups.items():
        res.loc[idx, 'p_bh'] = benjamini_hochberg(res.loc[idx, 'p_unilateral'])
    excluye_cero = np.where(res['metrica'] == 'brier', res['ic95_sup'] < 0, res['ic95_inf'] > 0)
    res['mejora_significativa'] = prin & excluye_cero & (res['p_bh'] < 0.05)
    res = res.round(4)
    res.to_csv(f'{RUTA_SALIDA}/comparaciones_bootstrap.csv', index=False)

    vista = res[(res.analisis == 'principal') & res.principal][
        ['hipotesis', 'desenlace', 'familia', 'metrica', 'valor_referencia', 'valor_nuevo',
         'delta', 'ic95_inf', 'ic95_sup', 'p_bh', 'mejora_significativa']]
    print('\nAnálisis principal — métrica de cada hipótesis:')
    print(vista.to_string(index=False))
    return res


if __name__ == '__main__':
    main()
