"""Post hoc alignment check: univariable AUROC of each of the 41 HRV predictors for each outcome.

Manuscript: Results > Post hoc analyses of the HRV signal > Alignment check.
Purpose: confirm that HRV predictors and outcome labels are correctly aligned (the expected
univariable signal is present), after the null primary result was known. Descriptive only.

Windows: primary analytic windows (incl_hrv_6h == 1) evaluable for each outcome (label not missing).
AUROC is reported both raw (higher value -> outcome) and oriented, max(AUROC, 1 - AUROC), so that
values above 0.5 indicate separation in either direction; 'direction' says which.

Usage:
  python 01_alignment_check.py --data-dir /path/to/tesis_sepsis --out-dir results/post_hoc
Inputs : ultimate_dataset_v3_1.csv, columnas_v3_1.json
Outputs: alignment_check_auroc.csv (all 41 x 3), alignment_check_top.csv (top 5 per outcome)
Expected (manuscript): haemodynamic SD2_last 0.613, SDNN_last 0.606; deterioration SD2_mean 0.639,
DFA_a1_mean 0.607 (lower values -> deterioration); SOFA improvement maximum 0.544 (HF_sd).
"""
import argparse
import json
import os

import pandas as pd
from sklearn.metrics import roc_auc_score

OUTCOMES = {'label_sofa_responder': 'SOFA improvement',
            'label_cardio': 'Haemodynamic stability or recovery',
            'label_sofa_deterioro_2pts': 'SOFA deterioration'}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data-dir', required=True)
    ap.add_argument('--out-dir', default='results/post_hoc')
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    d = pd.read_csv(os.path.join(a.data_dir, 'ultimate_dataset_v3_1.csv'))
    hrv = json.load(open(os.path.join(a.data_dir, 'columnas_v3_1.json'), encoding='utf8'))['features_hrv']
    assert len(hrv) == 41
    d = d[d.incl_hrv_6h == 1]

    rows = []
    for lab, name in OUTCOMES.items():
        e = d[d[lab].notna()]
        for c in hrv:
            ok = e[c].notna()
            auc = roc_auc_score(e.loc[ok, lab].astype(int), e.loc[ok, c])
            rows.append({'outcome': name, 'predictor': c, 'windows': int(ok.sum()),
                         'patients': e.loc[ok, 'icustay_id'].nunique(), 'events': int(e.loc[ok, lab].sum()),
                         'missing_fraction': round(1 - ok.mean(), 4), 'auroc_raw': round(auc, 4),
                         'auroc_oriented': round(max(auc, 1 - auc), 4),
                         'direction': 'higher value -> outcome' if auc >= 0.5 else 'lower value -> outcome'})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(a.out_dir, 'alignment_check_auroc.csv'), index=False)
    top = (res.sort_values('auroc_oriented', ascending=False).groupby('outcome', sort=False).head(5)
           .sort_values(['outcome', 'auroc_oriented'], ascending=[True, False]))
    top.to_csv(os.path.join(a.out_dir, 'alignment_check_top.csv'), index=False)
    print(top[['outcome', 'predictor', 'auroc_oriented', 'direction', 'windows', 'events']].to_string(index=False))


if __name__ == '__main__':
    main()
