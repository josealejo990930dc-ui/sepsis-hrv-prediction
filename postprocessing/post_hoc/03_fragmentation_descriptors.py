"""Post hoc fragmentation descriptors by age tertile and by in-hospital death.

Manuscript: Results > Post hoc analyses of the HRV signal > Fragmentation descriptors;
Discussion > Beat-to-beat HRV and mortality: a hypothesis.
Unit: first eligible window of each patient (incl_hrv_6h == 1, lowest hour_index), n = 487.
Descriptive only; formal fragmentation indices (PIP, IALS, PSS, PAS) were not computed.

Outputs
  fragmentation_by_age_tertile.csv : median age, RMSSD, SDNN, RMSSD/SDNN ratio (median of the
                                     per-patient ratio) and DFA alpha1 by age tertile (pd.qcut, 3 groups)
  fragmentation_by_death.csv       : median RMSSD/SDNN ratio in patients who died vs survivors
  fragmentation_correlation.csv    : Spearman correlation between mean RMSSD and mean DFA alpha1
Usage:
  python 03_fragmentation_descriptors.py --data-dir /path/to/tesis_sepsis --out-dir results/post_hoc
Expected (manuscript): tertile medians of age 46 / 63 / 79 years; RMSSD/SDNN 0.55 / 0.70 / 0.82;
DFA alpha1 1.13 / 1.00 / 0.84; ratio 0.95 in deaths vs 0.64 in survivors; RMSSD vs DFA rho -0.28.
"""
import argparse
import os

import pandas as pd
from scipy.stats import spearmanr


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data-dir', required=True)
    ap.add_argument('--out-dir', default='results/post_hoc')
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    d = pd.read_csv(os.path.join(a.data_dir, 'ultimate_dataset_v3_1.csv'))
    f = d[d.incl_hrv_6h == 1].sort_values(['icustay_id', 'hour_index']).groupby('icustay_id').head(1)
    assert len(f) == 487, len(f)
    f = f.assign(ratio=f.RMSSD_mean / f.SDNN_mean, death=f.deathtime.notna(),
                 age_tertile=pd.qcut(f.age, 3, labels=['lower', 'middle', 'upper']))

    terc = f.groupby('age_tertile', observed=True).agg(
        patients=('icustay_id', 'size'), age_median=('age', 'median'), rmssd_median_ms=('RMSSD_mean', 'median'),
        sdnn_median_ms=('SDNN_mean', 'median'), rmssd_sdnn_ratio_median=('ratio', 'median'),
        dfa_a1_median=('DFA_a1_mean', 'median')).round(2).reset_index()
    terc.to_csv(os.path.join(a.out_dir, 'fragmentation_by_age_tertile.csv'), index=False)
    print(terc.to_string(index=False))

    death = (f.groupby('death').agg(patients=('icustay_id', 'size'), rmssd_sdnn_ratio_median=('ratio', 'median'))
             .round(2).rename(index={True: 'died', False: 'survived'}).reset_index())
    death.to_csv(os.path.join(a.out_dir, 'fragmentation_by_death.csv'), index=False)
    print(death.to_string(index=False))

    r = spearmanr(f.RMSSD_mean, f.DFA_a1_mean, nan_policy='omit')
    cor = pd.DataFrame([{'x': 'RMSSD_mean', 'y': 'DFA_a1_mean', 'n': len(f),
                         'spearman_rho': round(r.statistic, 4), 'p_value': float(f'{r.pvalue:.3g}')}])
    cor.to_csv(os.path.join(a.out_dir, 'fragmentation_correlation.csv'), index=False)
    print(cor.to_string(index=False))


if __name__ == '__main__':
    main()
