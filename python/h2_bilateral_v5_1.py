# -*- coding: utf-8 -*-
"""
h2_bilateral_v5_1.py  (two-sided p-values for H2)
=================================================
Recomputes H2 with a TWO-SIDED test, as pre-specified in the research proposal
(delta AUROC != 0 between non-linear HRV and linear HRV; numbered H1 in the
proposal and H2 in the code and the results files).

bootstrap_v5.py (already run) computed a one-sided p for every hypothesis. H2
needs two tails. The bootstrap does not have to be repeated:
  p_bilateral = min(1, 2 * min(p_unilateral, 1 - p_unilateral))
where p_unilateral is the share of replicates with delta <= 0 (2,000 replicates,
continuous AUROC, so exact ties at 0 are negligible). The 95 % percentile CI was
already two-sided.

Benjamini-Hochberg correction across the 4 families, within each outcome (as in
bootstrap_v5.py). Significant difference (in either direction): 95 % CI excluding
0 and BH-corrected two-sided p < 0.05.

Input : modelado_v5_1/comparaciones_bootstrap.csv  (output of bootstrap_v5.py; not modified)
Output: modelado_v5_1/comparaciones_h2_bilateral.csv
"""

import numpy as np
import pandas as pd

RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA = f'{RUTA_BASE}/modelado_v5_1'


def benjamini_hochberg(p):
    """Benjamini-Hochberg adjusted p-values (no NaNs expected here)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    orden = np.argsort(p)
    ajust = p[orden] * n / np.arange(1, n + 1)
    ajust = np.minimum.accumulate(ajust[::-1])[::-1]
    res = np.empty(n)
    res[orden] = np.minimum(ajust, 1.0)
    return res


def main(ruta=RUTA):
    """Read the H2 rows (AUROC, main analysis) and add two-sided p-values."""
    c = pd.read_csv(f'{ruta}/comparaciones_bootstrap.csv')
    h2 = c[(c.analisis == 'principal') & (c.hipotesis == 'H2') & (c.metrica == 'auroc')].copy()
    assert len(h2) == 8, f'Se esperaban 8 filas de H2 (2 desenlaces × 4 familias); hay {len(h2)}'
    p = h2['p_unilateral'].to_numpy()
    h2['p_bilateral'] = np.minimum(1.0, 2 * np.minimum(p, 1 - p))
    h2['p_bh_bilateral'] = np.nan
    for _, idx in h2.groupby('desenlace').groups.items():
        h2.loc[idx, 'p_bh_bilateral'] = benjamini_hochberg(h2.loc[idx, 'p_bilateral'])
    excluye_cero = (h2['ic95_inf'] > 0) | (h2['ic95_sup'] < 0)
    h2['diferencia_significativa'] = excluye_cero & (h2['p_bh_bilateral'] < 0.05)
    cols = ['desenlace', 'familia', 'valor_referencia', 'valor_nuevo', 'delta', 'ic95_inf', 'ic95_sup',
            'p_unilateral', 'p_bh', 'p_bilateral', 'p_bh_bilateral', 'diferencia_significativa',
            'ventanas', 'pacientes']
    out = h2[cols].round(4)
    out.to_csv(f'{ruta}/comparaciones_h2_bilateral.csv', index=False)
    print(out.to_string(index=False))
    print(f"\nDiferencias significativas (bilateral, BH): {int(out['diferencia_significativa'].sum())} de {len(out)}")
    return out


if __name__ == '__main__':
    try:
        from google.colab import drive
        drive.mount('/content/drive')
    except ImportError:
        pass
    main()
