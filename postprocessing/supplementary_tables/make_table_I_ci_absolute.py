"""Patient-level bootstrap 95% CIs for absolute AUROC and AUPRC (primary analysis, modelado_v5_1).

Post-processing of the saved out-of-fold predictions only: no model is refitted.
Same design as bootstrap_v5.py: resample patients (icustay_id) with replacement,
keep all their windows, 2000 replicates, seed 2026, percentile CI.
AUROC/AUPRC on the uncalibrated probability (as in resumen_metricas.csv).

Usage: python ci_absolute.py <tesis_sepsis dir> <output csv>
"""
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

N_REP, SEED = 2000, 2026
base, out = sys.argv[1], sys.argv[2]
pred = pd.read_csv(f'{base}/modelado_v5_1/predicciones_oof.csv')
pred = pred[pred.analisis == 'principal']
ref = pd.read_csv(f'{base}/modelado_v5_1/resumen_metricas.csv')
ref = ref[ref.analisis == 'principal'].set_index(['desenlace', 'conjunto', 'familia'])

rows = []
for (des, cj, fam), g in pred.groupby(['desenlace', 'conjunto', 'familia']):
    y, p = g.y.to_numpy(), g.p_sin_calibrar.to_numpy()
    ids = g.icustay_id.to_numpy()
    uniq, inv = np.unique(ids, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(len(uniq))]
    rng = np.random.default_rng(SEED)
    auc, ap = [], []
    for _ in range(N_REP):
        idx = np.concatenate([members[k] for k in rng.integers(0, len(uniq), len(uniq))])
        yy = y[idx]
        if yy.min() == yy.max():
            continue
        auc.append(roc_auc_score(yy, p[idx]))
        ap.append(average_precision_score(yy, p[idx]))
    a0, ap0 = roc_auc_score(y, p), average_precision_score(y, p)
    r = ref.loc[(des, cj, fam)]
    rows.append(dict(desenlace=des, conjunto=cj, familia=fam, n=len(y), pacientes=len(uniq), eventos=int(y.sum()),
                     auroc=a0, auroc_lo=np.percentile(auc, 2.5), auroc_hi=np.percentile(auc, 97.5),
                     auprc=ap0, auprc_lo=np.percentile(ap, 2.5), auprc_hi=np.percentile(ap, 97.5),
                     replicas=len(auc), check_auroc=abs(a0 - r.auroc) < 5e-4, check_auprc=abs(ap0 - r.auprc) < 5e-4))
    print(des, cj, fam, f'{a0:.3f} ({rows[-1]["auroc_lo"]:.3f}-{rows[-1]["auroc_hi"]:.3f})', flush=True)
res = pd.DataFrame(rows)
res.to_csv(out, index=False)
print('all point estimates match resumen_metricas:', bool(res.check_auroc.all() and res.check_auprc.all()))
