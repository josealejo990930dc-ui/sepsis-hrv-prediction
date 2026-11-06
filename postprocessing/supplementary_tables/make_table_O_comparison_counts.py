"""S1 Appendix Table O: number of prespecified bootstrap comparisons by analysis block.

Counts the rows of comparaciones_bootstrap.csv that use the primary metric of each hypothesis
(column 'principal' == True), by analysis, hypothesis, outcome and model family, and checks the
totals (35 primary + 48 sensitivity = 83; x 3 metrics = 249 rows).

Usage:
  python make_table_O_comparison_counts.py --results-dir /path/to/tesis_sepsis/modelado_v5_1 --out results/tables
Outputs: table_O_comparison_counts.csv and .tex
"""
import argparse
import os

import pandas as pd

BS = chr(92)
NL = ' ' + BS + BS
FAM = [('logistica', 'LR'), ('elastic_net', 'EN'), ('random_forest', 'RF'), ('xgboost', 'XGB'), ('pca_logistica', 'PCA-LR')]
OUTC = {'sofa_respuesta': 'SOFA improvement', 'hemodinamico': 'Haemodynamic stability or recovery',
        'progresion': 'SOFA deterioration'}
ANAL = [('hrv_4h', 'HRV $' + BS + 'geq$4~h'), ('hrv_8h', 'HRV $' + BS + 'geq$8~h'),
        ('sofa_caso_completo', 'SOFA complete case'), ('laboratorio_nuevo', 'New laboratory'),
        ('cohorte_comun', 'Common cohort')]
SENS_LABEL = {('sofa_respuesta', 'hemodinamico', 'progresion'): 'SOFA, haemodynamic, deterioration',
              ('sofa_respuesta', 'progresion'): 'SOFA improvement, deterioration',
              ('sofa_respuesta', 'hemodinamico'): 'SOFA, haemodynamic'}
HYP = [('H1', 'H1'), ('H2', 'H2 (two-sided)'), ('H3', 'H3'), ('Complementario', 'Complementary')]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--results-dir', required=True)
    ap.add_argument('--out', default='results/tables')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    b = pd.read_csv(os.path.join(a.results_dir, 'comparaciones_bootstrap.csv'))
    n_all = len(b)
    b = b[b.principal.astype(str) == 'True']

    rows = []
    for an in ['principal'] + [k for k, _ in ANAL]:
        for hk, _ in HYP:
            for dk in OUTC:
                g = b[(b.analisis == an) & (b.hipotesis == hk) & (b.desenlace == dk)]
                if len(g):
                    rows.append({'analysis': an, 'hypothesis': hk, 'outcome': dk,
                                 **{fl: int((g.familia == fk).sum()) for fk, fl in FAM}, 'comparisons': len(g)})
    counts = pd.DataFrame(rows)
    counts.to_csv(os.path.join(a.out, 'table_O_comparison_counts.csv'), index=False)
    prim = counts[counts.analysis == 'principal'].comparisons.sum()
    sens = counts[counts.analysis != 'principal'].comparisons.sum()
    assert prim + sens == len(b) and len(b) * 3 == n_all, (prim, sens, len(b), n_all)

    L = [BS + 'begin{table}[htbp]', BS + 'centering',
         BS + 'caption{Number of patient-level bootstrap comparisons by analysis block.}',
         BS + 'label{tab:comparison_counts}', BS + 'small', BS + 'begin{tabular}{@{}llcccccc@{}}', BS + 'toprule',
         'Block & Outcome(s) & ' + ' & '.join(fl for _, fl in FAM) + ' & Comparisons' + NL, BS + 'midrule',
         BS + 'multicolumn{8}{@{}l}{' + BS + 'textit{Primary analysis}}' + NL]
    p = counts[counts.analysis == 'principal']
    for hk, hl in HYP:
        g = p[p.hypothesis == hk]
        if hk == 'Complementario':
            cells = [str(g[fl].sum()) if g[fl].sum() else '---' for _, fl in FAM]
            L.append(f'{BS}quad {hl} & {len(g)} outcomes & ' + ' & '.join(cells) + f' & {g.comparisons.sum()}' + NL)
            continue
        for _, r in g.iterrows():
            cells = [str(r[fl]) if r[fl] else '---' for _, fl in FAM]
            L.append(f'{BS}quad {hl} & {OUTC[r.outcome]} & ' + ' & '.join(cells) + f' & {r.comparisons}' + NL)
    L += [f'{BS}quad {BS}textit{{Subtotal}} & & & & & & & {BS}textit{{{prim}}}' + NL, BS + 'addlinespace',
          BS + 'multicolumn{8}{@{}l}{' + BS + 'textit{Sensitivity analyses (H1 and H3 only)}}' + NL]
    for ak, al in ANAL:
        g = counts[counts.analysis == ak]
        outs = SENS_LABEL[tuple(dict.fromkeys(g.outcome))]
        cells = [str(g[fl].sum()) if g[fl].sum() else '---' for _, fl in FAM]
        L.append(f'{BS}quad {al} & {outs} & ' + ' & '.join(cells) + f' & {g.comparisons.sum()}' + NL)
    L += [f'{BS}quad {BS}textit{{Subtotal}} & & & & & & & {BS}textit{{{sens}}}' + NL, BS + 'midrule',
          BS + 'multicolumn{7}{@{}l}{Total on the primary metric of each hypothesis} & ' + str(prim + sens) + NL,
          BS + 'multicolumn{7}{@{}l}{Total including the two secondary metrics (' + str(prim + sens) + ' $' + BS +
          'times$ 3 metrics)} & ' + str(n_all) + NL,
          BS + 'bottomrule', BS + 'end{tabular}',
          BS + 'par' + BS + 'smallskip' + BS + 'footnotesize The prespecified criterion (95' + BS + '% CI excluding '
          'zero and Benjamini--Hochberg-corrected ' + BS + 'textit{p} $<$0.05) applies to the ' + str(prim + sens) +
          ' primary-metric comparisons; the ' + str(n_all - prim - sens) + ' secondary-metric contrasts were not '
          'corrected for multiplicity. LR, logistic regression; EN, elastic net; RF, random forest; XGB, XGBoost; '
          'PCA-LR, principal components + logistic regression.', BS + 'end{table}']
    open(os.path.join(a.out, 'table_O_comparison_counts.tex'), 'w', encoding='utf8').write('\n'.join(L) + '\n')
    print(counts.to_string(index=False))
    print(f'primary {prim} + sensitivity {sens} = {prim + sens}; all rows {n_all}')


if __name__ == '__main__':
    main()
