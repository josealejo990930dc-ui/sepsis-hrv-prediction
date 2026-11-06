# -*- coding: utf-8 -*-
"""
exploratorio_proxy_fragmentacion.py  (exploratory, decision D38; not part of the main analysis)
==============================================================================================
Question: does adding the RMSSD/SDNN ratio (a proxy of heart-rate fragmentation) to the
clinical + SOFA model improve (a) in-hospital mortality (first window per patient, AUROC)
or (b) SOFA deterioration (all evaluable windows, AUPRC)?

Decision rule written before running: re-extract fragmentation indices only if the
improvement exceeds 0.01. Note: the code flags `supera_regla` as delta > 0.01; the
confidence interval is reported but not used in that flag.

Fixed hyperparameters (no search), 5-fold CV repeated 5 times (grouped by patient for
deterioration), patient-level bootstrap (2,000 replicates, seed 2026).
Input : ultimate_dataset_v3_1.csv, columnas_v3_1.json. Does not touch modelado_v5_1.
Output: paso0_proxy_fragmentacion.csv
Originally run as a Colab cell after pipeline_modelado_v5.py.
"""
import json, warnings, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
warnings.filterwarnings('ignore')

RUTA = '/content/drive/MyDrive/tesis_sepsis'
d = pd.read_csv(f'{RUTA}/ultimate_dataset_v3_1.csv'); J = json.load(open(f'{RUTA}/columnas_v3_1.json'))
d = d[d['incl_hrv_6h'] == 1].copy()
d['sexo_masculino'] = (d['gender'].astype(str).str.upper().str[0] == 'M').astype(int) \
    if not pd.api.types.is_numeric_dtype(d['gender']) else d['gender']
if 'lactato_medido_12h' not in d: d['lactato_medido_12h'] = d['lactate_last_actual12h'].notna().astype(int)
for c in ['map_t0', 'mbp_mean_actual12h']:
    d.loc[(d[c] < 20) | (d[c] > 200), c] = np.nan

# The dataset has window aggregates only, so the ratio is computed on the _mean and _last values
for s in ['mean', 'last']:
    d[f'frag_ratio_{s}'] = d[f'RMSSD_{s}'] / d[f'SDNN_{s}'].replace(0, np.nan)
PROXY = ['frag_ratio_mean', 'frag_ratio_last']
BASE = [c if c != 'gender' else 'sexo_masculino' for c in J['features_clinicas']] \
       + J['features_pam_lactato_actual'] + J['features_sofa']
print('Proxy RMSSD/SDNN (mediana [p25–p75]):',
      round(d['frag_ratio_mean'].median(), 2), d['frag_ratio_mean'].quantile([.25, .75]).round(2).tolist())

def rf(): return make_pipeline(SimpleImputer(strategy='median'),
                               RandomForestClassifier(500, min_samples_leaf=5, n_jobs=-1, random_state=0))
def en(): return make_pipeline(SimpleImputer(strategy='median'), StandardScaler(),
                               LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.1, max_iter=5000))

def oof(x, y, cols, modelo, grupos=None, reps=5):
    pr = np.zeros(len(x))
    for r in range(reps):
        sp = (StratifiedGroupKFold(5, shuffle=True, random_state=r).split(x, y, grupos) if grupos is not None
              else StratifiedKFold(5, shuffle=True, random_state=r).split(x, y))
        for tr, te in sp:
            w = None if grupos is None else 1 / pd.Series(grupos[tr]).map(pd.Series(grupos[tr]).value_counts()).values
            m = modelo(); m.fit(x.iloc[tr][cols], y[tr], **{m.steps[-1][0] + '__sample_weight': w} if w is not None else {})
            pr[te] += m.predict_proba(x.iloc[te][cols])[:, 1] / reps
    return pr

def boot_delta(y, p_ref, p_new, met, grupos, n=2000, seed=2026):
    rng = np.random.default_rng(seed); u = np.unique(grupos)
    filas = {g: np.flatnonzero(grupos == g) for g in u}; dl = []
    for _ in range(n):
        i = np.concatenate([filas[g] for g in rng.choice(u, len(u))])
        if 0 < y[i].sum() < len(i): dl.append(met(y[i], p_new[i]) - met(y[i], p_ref[i]))
    dl = np.array(dl)
    return met(y, p_new) - met(y, p_ref), np.percentile(dl, 2.5), np.percentile(dl, 97.5), (dl <= 0).mean()

res = []
# (a) In-hospital mortality, first window per patient
p = d.sort_values(['icustay_id', 'hour_index']).groupby('icustay_id').head(1).reset_index(drop=True)
yp = p['deathtime'].notna().astype(int).values
assert 0.05 < yp.mean() < 0.35, 'Prevalencia de mortalidad implausible'
# (b) SOFA deterioration, evaluable windows
v = d[d['label_sofa_deterioro_2pts'].notna()].reset_index(drop=True)
yv = v['label_sofa_deterioro_2pts'].astype(int).values; gv = v['icustay_id'].values

for nombre, x, y, g, met, mname in [('mortalidad', p, yp, None, roc_auc_score, 'AUROC'),
                                     ('progresion', v, yv, gv, average_precision_score, 'AUPRC')]:
    for fam, mod in [('ElasticNet', en), ('RandomForest', rf)]:
        pr_ref = oof(x, y, BASE, mod, g); pr_new = oof(x, y, BASE + PROXY, mod, g)
        grupos = g if g is not None else np.arange(len(x))
        dlt, lo, hi, pu = boot_delta(y, pr_ref, pr_new, met, grupos)
        res.append({'desenlace': nombre, 'familia': fam, 'metrica': mname, 'eventos': int(y.sum()), 'n': len(y),
                    'clinico_sofa': round(met(y, pr_ref), 3), '+proxy': round(met(y, pr_new), 3),
                    'delta': round(dlt, 3), 'IC95': f'[{lo:+.3f}, {hi:+.3f}]', 'p_unilateral': round(pu, 3),
                    'supera_regla': dlt > 0.01})
R = pd.DataFrame(res); print(R.to_string(index=False))
R.to_csv(f'{RUTA}/paso0_proxy_fragmentacion.csv', index=False)
print('\nDECISIÓN (regla preespecificada):',
      'proceder a re-extracción' if R['supera_regla'].any() else 'NO re-extraer; fragmentación como línea futura')
