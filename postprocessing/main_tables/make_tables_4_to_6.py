"""Generate LaTeX results tables directly from the package CSVs (no hand transcription).

Usage (from the manuscript folder):
    python make_tables.py ../../fuentes/paquete_v1_5/paquete_tesis_hrv_sepsis_v1_5/resultados
Writes tables/tab_performance.tex, tables/tab_hypotheses.tex, tables/tab_calibration.tex
"""
import csv
import os
import sys

BS = chr(92)
RES = sys.argv[1]
OUT = 'tables'
os.makedirs(OUT, exist_ok=True)

def read(name):
    path = os.path.join(RES, name)
    if not os.path.exists(path):  # calibration results are written to the calibracion_shap subfolder
        path = os.path.join(RES, 'calibracion_shap', name)
    return list(csv.DictReader(open(path, encoding='utf-8')))

SETS = [('hrv_completa', 'HRV, all'), ('hrv_lineal', 'HRV, linear'), ('hrv_nolineal', 'HRV, nonlinear'),
        ('clinico_basal', 'Clinical baseline'), ('clinico_sofa', 'Clinical + SOFA'),
        ('hrv_integrada', 'HRV + clinical'), ('integrada_sofa', 'Integrated + SOFA')]
FAMS = [('logistica', 'LR'), ('elastic_net', 'EN'), ('random_forest', 'RF'), ('xgboost', 'XGB'), ('pca_logistica', 'PCA-LR')]
OUTCOMES = [('sofa_respuesta', 'SOFA improvement', 'auroc'), ('hemodinamico', 'Haemodynamic stability or recovery', 'auroc'),
            ('progresion', 'SOFA deterioration', 'auprc')]
MINUS = '$-$'

def f3(x):
    # round half away from zero from the exact decimal string in the CSV (float formatting would
    # turn -0.0225 into -0.022 while the text, rounding half up, says -0.023)
    from decimal import Decimal, ROUND_HALF_UP
    d = Decimal(str(x).strip()).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)
    s = f'{d:.3f}'
    return s.replace('-', MINUS) if d < 0 else s

def fp(x):
    v = float(x)
    return '$<$0.001' if v < 0.001 else f'{v:.3f}'

# ---------- performance ----------
R = [r for r in read('resumen_metricas.csv') if r['analisis'] == 'principal']
L = [BS + 'begin{table}[htbp]', BS + 'centering',
     BS + 'caption{Out-of-fold discrimination in the primary analysis. AUROC for SOFA improvement and haemodynamic stability or recovery; '
     'AUPRC for SOFA deterioration (chance level = prevalence, 0.056).}',
     BS + 'label{tab:performance}', BS + 'small',
     BS + 'begin{tabular}{@{}l' + 'c' * len(FAMS) + '@{}}', BS + 'toprule',
     'Predictor set & ' + ' & '.join(f for _, f in FAMS) + ' ' + BS + BS, BS + 'midrule']
for key, label, met in OUTCOMES:
    rows = [r for r in R if r['desenlace'] == key]
    n, pts, prev = rows[0]['n'], rows[0]['pacientes'], float(rows[0]['prevalencia'])
    L.append(BS + 'multicolumn{%d}{@{}l}{%s{%s} (%s windows, %s patients; prevalence %.1f%s; %s)} %s%s'
             % (len(FAMS) + 1, BS + 'textit', label, n, pts, prev * 100, BS + '%', met.upper(), BS, BS))
    for sk, sl in SETS:
        cells = []
        for fk, _ in FAMS:
            m = [r for r in rows if r['conjunto'] == sk and r['familia'] == fk]
            cells.append(f3(m[0][met]) if m else '---')
        L.append(f'{BS}quad {sl} & ' + ' & '.join(cells) + f' {BS}{BS}')
    L.append(BS + 'addlinespace')
L += [BS + 'bottomrule', BS + 'end{tabular}',
      BS + 'par' + BS + 'smallskip' + BS + 'footnotesize LR, logistic regression with Welch filter; EN, elastic net; '
      'RF, random forest; XGB, XGBoost; PCA-LR, principal components + logistic regression (fitted on the two '
      'central sets only). 95' + BS + '% CIs for every configuration are given in eTable~2.', BS + 'end{table}']
open(os.path.join(OUT, 'tab_performance.tex'), 'w', encoding='utf8').write('\n'.join(L) + '\n')

# ---------- hypotheses ----------
B = [b for b in read('comparaciones_bootstrap.csv') if b['analisis'] == 'principal' and b['principal'] == 'True']
H2 = {(h['desenlace'], h['familia']): h for h in read('comparaciones_h2_bilateral.csv')}
HYP = [('H1', 'sofa_respuesta', 'H1, SOFA improvement (AUROC): integrated + SOFA vs clinical + SOFA'),
       ('H1', 'hemodinamico', 'H1, haemodynamic stability or recovery (AUROC): integrated + SOFA vs clinical + SOFA'),
       ('H2', 'sofa_respuesta', 'H2, SOFA improvement (AUROC, two-sided): nonlinear vs linear HRV'),
       ('H2', 'hemodinamico', 'H2, haemodynamic stability or recovery (AUROC, two-sided): nonlinear vs linear HRV'),
       ('H3', 'progresion', 'H3, SOFA deterioration (AUPRC): integrated + SOFA vs clinical + SOFA'),
       ('Complementario', 'sofa_respuesta', 'Complementary, SOFA improvement (AUROC): HRV + clinical vs clinical baseline'),
       ('Complementario', 'hemodinamico', 'Complementary, haemodynamic stability or recovery (AUROC): HRV + clinical vs clinical baseline'),
       ('Complementario', 'progresion', 'Complementary, SOFA deterioration (AUPRC): HRV + clinical vs clinical baseline')]
L = [BS + 'begin{table}[htbp]', BS + 'centering',
     BS + 'caption{Prespecified comparisons (patient-level bootstrap, 2000 replicates). $' + BS +
     'Delta$ = new $-$ reference model, same family and windows. ' + BS + 'textit{p} values are one-sided '
     '(improvement) except for H2 (two-sided), Benjamini--Hochberg-corrected across families.}',
     BS + 'label{tab:hypotheses}', BS + 'small',
     BS + 'begin{tabular}{@{}lccccc@{}}', BS + 'toprule',
     'Family & Reference & New & $' + BS + 'Delta$ (95' + BS + '% CI) & Corrected ' + BS + 'textit{p} & Criterion met ' + BS + BS,
     BS + 'midrule']
for h, d, title in HYP:
    L.append(BS + 'multicolumn{6}{@{}l}{' + BS + 'textit{' + title + '}} ' + BS + BS)
    for fk, fl in FAMS:
        m = [b for b in B if b['hipotesis'] == h and b['desenlace'] == d and b['familia'] == fk]
        if not m:
            continue
        b = m[0]
        p = H2[(d, fk)]['p_bh_bilateral'] if h == 'H2' else b['p_bh']
        sig = H2[(d, fk)]['diferencia_significativa'] if h == 'H2' else b['mejora_significativa']
        L.append(f"{BS}quad {fl} & {f3(b['valor_referencia'])} & {f3(b['valor_nuevo'])} & "
                 f"{f3(b['delta'])} ({f3(b['ic95_inf'])} to {f3(b['ic95_sup'])}) & {fp(p)} & "
                 f"{'Yes' if sig == 'True' else 'No'} {BS}{BS}")
    L.append(BS + 'addlinespace')
L += [BS + 'bottomrule', BS + 'end{tabular}',
      BS + 'par' + BS + 'smallskip' + BS + 'footnotesize Criterion met: for H1, H3 and the complementary comparison, improvement (95' + BS + '% CI excluding zero '
      'in the direction of improvement and corrected ' + BS + 'textit{p} $<$ 0.05); for H2 (two-sided), a difference in either '
      'direction (95' + BS + '% CI excluding zero and corrected ' + BS + 'textit{p} $<$ 0.05). Abbreviations as in Table~' + BS + 'ref{tab:performance}.', BS + 'end{table}']
open(os.path.join(OUT, 'tab_hypotheses.tex'), 'w', encoding='utf8').write('\n'.join(L) + '\n')

# ---------- calibration ----------
C = read('calibracion_metricas.csv')
def ci(s):
    lo, hi = s.replace('–', '|').replace('–', '|').split('|')
    return f'{f3(lo)} to {f3(hi)}'
L = [BS + 'begin{table}[htbp]', BS + 'centering',
     BS + 'caption{Calibration of Platt-calibrated out-of-fold predictions (95' + BS + '% CI from 1000 '
     'patient-level bootstrap replicates).}', BS + 'label{tab:calibration}', BS + 'scriptsize',
     BS + 'begin{tabular}{@{}llcccc@{}}', BS + 'toprule',
     'Set & Family & Intercept & Slope & Scaled Brier & ICI ' + BS + BS, BS + 'midrule']
for key, label, _ in OUTCOMES:
    L.append(BS + 'multicolumn{6}{@{}l}{' + BS + 'textit{' + label + '}} ' + BS + BS)
    for sk, sl in [('clinico_sofa', 'Clinical + SOFA'), ('integrada_sofa', 'Integrated + SOFA')]:
        for fk, fl in FAMS:
            m = [r for r in C if r['desenlace'] == key and r['conjunto'] == sk and r['familia'] == fk]
            if not m:
                continue
            r = m[0]
            L.append(f"{BS}quad {sl} & {fl} & {f3(r['intercepto'])} ({ci(r['intercepto_ic'])}) & "
                     f"{f3(r['pendiente'])} ({ci(r['pendiente_ic'])}) & {f3(r['brier_escalado'])} "
                     f"({ci(r['brier_escalado_ic'])}) & {f3(r['ICI'])} {BS}{BS}")
    L.append(BS + 'addlinespace')
L += [BS + 'bottomrule', BS + 'end{tabular}',
      BS + 'par' + BS + 'smallskip' + BS + 'footnotesize Ideal values: intercept 0, slope 1; scaled Brier $>$0 '
      'indicates better than predicting the prevalence. ICI, integrated calibration index. Expected-to-observed '
      'ratios and ICI CIs are in the shared results files.', BS + 'end{table}']
open(os.path.join(OUT, 'tab_calibration.tex'), 'w', encoding='utf8').write('\n'.join(L) + '\n')
print('written:', os.listdir(OUT))
