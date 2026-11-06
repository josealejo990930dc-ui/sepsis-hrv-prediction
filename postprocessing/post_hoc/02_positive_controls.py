"""Post hoc positive controls: HRV and in-hospital mortality, and HRV and age.

Manuscript: Results > Post hoc analyses of the HRV signal > Positive controls; S1 Appendix Table L.
Unit: the first eligible window of each patient (incl_hrv_6h == 1, lowest hour_index), n = 487.
In-hospital death: deathtime not missing (MIMIC-III admissions.deathtime, in-hospital deaths only).

1. Mortality: univariable AUROC of window means, oriented as stated in the manuscript
   (low HRV value -> death; high SOFA -> death), with 95% CI from 2000 bootstrap resamples of
   patients (one window per patient, so a simple patient bootstrap), seed 2026.
2. Age: Spearman correlation of age with mean RMSSD, SDNN and SD2.

Usage:
  python 02_positive_controls.py --data-dir /path/to/tesis_sepsis --out-dir results/post_hoc
Inputs : ultimate_dataset_v3_1.csv
Outputs: positive_control_mortality.csv, positive_control_age.csv
Expected (manuscript): DFA_a1_mean 0.693 (0.628-0.754); partial SOFA 0.679 (0.602-0.750);
LFHF_mean 0.636; SD2_mean 0.616; SDNN_mean 0.587; RMSSD_mean 0.463; SampEn_mean 0.415;
age vs RMSSD rho 0.29 (p 4.9e-11), SDNN 0.14, SD2 0.09.
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

N_REP, SEED = 2000, 2026
# variable -> sign (-1: low value predicts death; +1: high value predicts death)
MORTALITY_VARS = {'DFA_a1_mean': -1, 'LFHF_mean': -1, 'SD2_mean': -1, 'SDNN_mean': -1, 'LF_mean': -1,
                  'RMSSD_mean': -1, 'SampEn_mean': -1, 'sofa_t0_parcial': 1, 'sofa_t0': 1}
AGE_VARS = ['RMSSD_mean', 'SDNN_mean', 'SD2_mean']


def first_windows(path):
    d = pd.read_csv(path)
    f = d[d.incl_hrv_6h == 1].sort_values(['icustay_id', 'hour_index']).groupby('icustay_id').head(1)
    assert len(f) == 487, len(f)
    return f.assign(death=f.deathtime.notna().astype(int))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data-dir', required=True)
    ap.add_argument('--out-dir', default='results/post_hoc')
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    f = first_windows(os.path.join(a.data_dir, 'ultimate_dataset_v3_1.csv'))
    death = f.death.to_numpy()
    print(f'patients {len(f)}, in-hospital deaths {death.sum()}')

    rng = np.random.default_rng(SEED)
    boot = [rng.integers(0, len(f), len(f)) for _ in range(N_REP)]
    rows = []
    for v, sign in MORTALITY_VARS.items():
        x = sign * f[v].to_numpy(dtype=float)
        ok = ~np.isnan(x)
        auc = roc_auc_score(death[ok], x[ok])
        reps = [roc_auc_score(death[i], x[i]) for i in (b[ok[b]] for b in boot) if death[i].min() != death[i].max()]
        rows.append({'variable': v, 'orientation': 'low value -> death' if sign < 0 else 'high value -> death',
                     'n': int(ok.sum()), 'deaths': int(death[ok].sum()), 'auroc': round(auc, 4),
                     'ci95_low': round(np.percentile(reps, 2.5), 4), 'ci95_high': round(np.percentile(reps, 97.5), 4),
                     'replicates': len(reps)})
    mort = pd.DataFrame(rows)
    mort.to_csv(os.path.join(a.out_dir, 'positive_control_mortality.csv'), index=False)
    print(mort.to_string(index=False))

    age = []
    for v in AGE_VARS:
        r = spearmanr(f.age, f[v], nan_policy='omit')
        age.append({'variable': v, 'n': int(f[v].notna().sum()), 'spearman_rho': round(r.statistic, 4),
                    'p_value': float(f'{r.pvalue:.3g}')})
    age = pd.DataFrame(age)
    age.to_csv(os.path.join(a.out_dir, 'positive_control_age.csv'), index=False)
    print(age.to_string(index=False))


if __name__ == '__main__':
    main()
