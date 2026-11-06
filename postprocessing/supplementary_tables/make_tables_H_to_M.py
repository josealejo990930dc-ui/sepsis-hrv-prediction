"""Supplementary eTables generated from the results files (no hand transcription).

Usage: python make_tables_H_to_M.py <tesis_sepsis dir> <results dir> <tables out dir>
  (set H2_CSV to the path of comparaciones_h2_bilateral.csv if it is not in <tesis_sepsis>/modelado_v5_1)
Writes etables.tex = S1 Appendix Tables H-M:
  H missing data; I performance with 95% CIs (needs <results dir>/ci_absolute_v5_1.csv from
  make_table_I_ci_absolute.py); J all prespecified comparisons; K SHAP group shares;
  L positive control (needs <results dir>/positive_control_mortality.csv from post_hoc/02_positive_controls.py);
  M clinical SMD by outcome.
"""
import os
import sys

import pandas as pd

BASE, AN, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
BS = chr(92)
M = f'{BASE}/modelado_v5_1'
NL = ' ' + BS + BS

FAM = {'logistica': 'LR', 'elastic_net': 'EN', 'random_forest': 'RF', 'xgboost': 'XGB', 'pca_logistica': 'PCA-LR'}
FAM_ORDER = list(FAM)
SET = {'hrv_completa': 'HRV, all', 'hrv_lineal': 'HRV, linear', 'hrv_nolineal': 'HRV, nonlinear',
       'clinico_basal': 'Clinical baseline', 'clinico_sofa': 'Clinical + SOFA', 'hrv_integrada': 'HRV + clinical',
       'integrada_sofa': 'Integrated + SOFA'}
OUTC = {'sofa_respuesta': 'SOFA improvement', 'hemodinamico': 'Haemodynamic stability or recovery', 'progresion': 'SOFA deterioration'}
ANAL = {'principal': 'Primary', 'hrv_4h': 'HRV $' + BS + 'geq$4~h', 'hrv_8h': 'HRV $' + BS + 'geq$8~h',
        'sofa_caso_completo': 'SOFA complete case', 'laboratorio_nuevo': 'New laboratory', 'cohorte_comun': 'Common cohort'}
HYP = {'H1': 'H1', 'H2': 'H2', 'H3': 'H3', 'Complementario': 'Compl.'}
VAR = {'lactate_last_actual12h': 'Last lactate, 12~h', 'sofa_resp_t0': 'SOFA respiratory component',
       'sofa_liver_t0': 'SOFA hepatic component', 'temp_t0': 'Temperature (not a predictor)',
       'sofa_neuro_t0': 'SOFA neurological component', 'map_t0': 'MAP at T0', 'sofa_cardio_t0': 'SOFA cardiovascular component',
       'sofa_coag_t0': 'SOFA coagulation component', 'sofa_renal_t0': 'SOFA renal component', 'rr_t0': 'Respiratory rate',
       'spo2_t0': 'SpO' + BS + 'textsubscript{2}', 'mbp_mean_actual12h': 'Mean MAP, 12~h', 'age': 'Age',
       'ned': 'Norepinephrine equivalent dose', 'sexo_masculino': 'Sex', 'sofa_t0_parcial': 'Partial SOFA',
       'opioide_bin': 'Opioids', 'cumvp': 'Cumulative vasopressor load', 'vis': 'Vasoactive-inotropic score',
       'sedante_bin': 'Sedation', 'n_comp_validos_t0': 'Number of SOFA components'}

def num(x, d=3):
    # round half away from zero from the decimal string (consistent with make_tables.py and the text)
    from decimal import Decimal, ROUND_HALF_UP
    v = Decimal(repr(float(x)) if not isinstance(x, str) else x.strip())
    v = v.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    s = f'{v:.{d}f}'
    return s.replace('-', '$-$') if v < 0 else s

def table(caption, label, spec, header, body, note=None, size='scriptsize', long=False):
    if long:
        # longtable cannot sit inside a table float: caption, header and note go inside it
        L = ['{' + BS + size, BS + 'begin{longtable}{' + spec + '}',
             BS + 'caption{' + caption + '}' + BS + 'label{' + label + '}' + NL,
             BS + 'toprule', header + NL, BS + 'midrule', BS + 'endfirsthead',
             BS + 'toprule', header + NL, BS + 'midrule', BS + 'endhead',
             BS + 'bottomrule', BS + 'endfoot']
        L += body
        L += [BS + 'end{longtable}']
        if note:
            L.append(BS + 'footnotesize ' + note + BS + 'par')
        L.append('}')
        return '\n'.join(L)
    L = [BS + 'begin{table}[htbp]', BS + 'centering', BS + 'caption{' + caption + '}', BS + 'label{' + label + '}',
         BS + size, BS + 'begin{tabular}{' + spec + '}', BS + 'toprule', header + NL, BS + 'midrule']
    L += body
    L += [BS + 'bottomrule', BS + 'end{tabular}']
    if note:
        L.append(BS + 'par' + BS + 'smallskip' + BS + 'footnotesize ' + note)
    L.append(BS + 'end{table}')
    return '\n'.join(L)

parts = []  # numbering is set by the including document (PLOS S1 Appendix: Tables A, B, ...)

# eTable 1 - missing data
c2 = pd.read_csv(f'{BASE}/eda_clinico_v3/C2_faltantes.csv')
body = [f"{VAR.get(r.columna, r.columna)} & {r.pct_faltante:.1f}" + NL for r in c2.itertuples()]
parts.append(table('Missing data by variable across the 2555 analytic windows.', 'etab:missing', '@{}lr@{}',
                    'Variable & Missing (' + BS + '%)', body,
                    'Missing clinical predictors were imputed with the training-fold median. HRV predictors had no '
                    'missing values in eligible windows. Source: clinical EDA, C2.', size='small'))

# eTable 2 - full performance with CIs
ci_path = f'{AN}/ci_absolute_v5_1.csv'
if os.path.exists(ci_path):
    ci = pd.read_csv(ci_path)
    rs = pd.read_csv(f'{M}/resumen_metricas.csv')
    rs = rs[rs.analisis == 'principal'][['desenlace', 'conjunto', 'familia', 'brier']]
    ci = ci.merge(rs, on=['desenlace', 'conjunto', 'familia'])
    body = []
    for d in OUTC:
        body.append(BS + 'multicolumn{5}{@{}l}{' + BS + 'textit{' + OUTC[d] + '}}' + NL)
        for s in SET:
            for f in FAM_ORDER:
                r = ci[(ci.desenlace == d) & (ci.conjunto == s) & (ci.familia == f)]
                if r.empty:
                    continue
                r = r.iloc[0]
                body.append(f"{BS}quad {SET[s]} & {FAM[f]} & {num(r.auroc)} ({num(r.auroc_lo)}--{num(r.auroc_hi)}) & "
                            f"{num(r.auprc)} ({num(r.auprc_lo)}--{num(r.auprc_hi)}) & {num(r.brier)}" + NL)
    parts.append(table('Out-of-fold performance of all primary configurations with patient-level bootstrap 95' + BS +
                       '% CIs (2000 replicates).', 'etab:performance_full', '@{}llccc@{}',
                       'Set & Family & AUROC (95' + BS + '% CI) & AUPRC (95' + BS + '% CI) & Brier', body,
                       'AUROC and AUPRC from uncalibrated predictions; Brier score from calibrated predictions. CIs '
                       'computed post hoc from the saved out-of-fold predictions (no refitting).', long=True))

# eTable 3 - all bootstrap comparisons (primary metric of each hypothesis, all analyses)
b = pd.read_csv(f'{M}/comparaciones_bootstrap.csv')
b = b[b.principal]
h2 = pd.read_csv(os.environ.get('H2_CSV', f'{M}/comparaciones_h2_bilateral.csv')).set_index(['desenlace', 'familia'])
body = []
for a in ANAL:
    for h in HYP:
        for d in OUTC:
            g = b[(b.analisis == a) & (b.hipotesis == h) & (b.desenlace == d)]
            if g.empty:
                continue
            body.append(BS + 'multicolumn{6}{@{}l}{' + BS + 'textit{' + f'{ANAL[a]} · {HYP[h]} · {OUTC[d]} '
                        f'({g.ventanas.iloc[0]} windows, {g.pacientes.iloc[0]} patients; {g.metrica.iloc[0].upper()})' + '}}' + NL)
            for f in FAM_ORDER:
                r = g[g.familia == f]
                if r.empty:
                    continue
                r = r.iloc[0]
                p = h2.loc[(d, f), 'p_bh_bilateral'] if h == 'H2' else r.p_bh
                body.append(f"{BS}quad {FAM[f]} & {num(r.valor_referencia)} & {num(r.valor_nuevo)} & "
                            f"{num(r.delta)} ({num(r.ic95_inf)} to {num(r.ic95_sup)}) & {num(p)} & No" + NL)
parts.append(table('All prespecified comparisons on the primary metric of each hypothesis, primary and sensitivity '
                   'analyses.', 'etab:comparisons', '@{}lccccc@{}',
                   'Family & Reference & New & $' + BS + 'Delta$ (95' + BS + '% CI) & Corrected ' + BS +
                   'textit{p} & Improvement', body,
                   'One-sided ' + BS + 'textit{p} values for improvement (two-sided for H2), Benjamini--Hochberg-corrected '
                   'across families. No comparison met the criterion for improvement.', long=True))

# eTable 4 - SHAP group shares
sg = pd.read_csv(f'{M}/calibracion_shap/shap_pct_por_grupo.csv')
G = [('Clínicas en T0', 'Clinical'), ('PAM y lactato (12 h)', 'MAP/lactate'), ('SOFA', 'SOFA'),
     ('HRV temporal', 'HRV time'), ('HRV frecuencial', 'HRV freq.'), ('HRV no lineal', 'HRV nonlin.'), ('HRV total', 'HRV total')]
body = [f"{OUTC[r.desenlace]} & {FAM[r.familia]} & " + ' & '.join(f'{float(getattr(r, "_" + str(i + 3))):.1f}'
        for i in range(len(G))) + NL for r in sg.itertuples()] if False else []
for _, r in sg.iterrows():
    body.append(f"{OUTC[r.desenlace]} & {FAM[r.familia]} & " + ' & '.join(f'{float(r[c]):.1f}' for c, _ in G) + NL)
parts.append(table('Share of mean absolute SHAP value (' + BS + '%) by predictor group, integrated + SOFA models.',
                   'etab:shap', '@{}ll' + 'r' * len(G) + '@{}', 'Outcome & Family & ' + ' & '.join(l for _, l in G), body,
                   'Out-of-fold exact SHAP values on the log-odds scale. A large group can accumulate a large share '
                   'even when it does not improve performance.'))

# eTable 5 - positive control
pc = pd.read_csv(f'{AN}/positive_control_mortality.csv')  # output of post_hoc/02_positive_controls.py
pc = pc.rename(columns={'ci95_low': 'lo', 'ci95_high': 'hi'})
NAMES = {'DFA_a1_mean': 'Mean DFA $' + BS + 'alpha_1$', 'LFHF_mean': 'Mean ln LF/HF', 'SD2_mean': 'Mean SD2',
         'SDNN_mean': 'Mean SDNN', 'LF_mean': 'Mean ln LF', 'RMSSD_mean': 'Mean RMSSD', 'SampEn_mean': 'Mean SampEn',
         'sofa_t0_parcial': 'Partial SOFA at T0', 'sofa_t0': 'Complete SOFA at T0'}
body = [f"{NAMES[r.variable]} & {'Low' if r.orientation.startswith('low') else 'High'} & {r.n} & "
        f"{num(r.auroc)} ({num(r.lo)}--{num(r.hi)})" + NL for r in pc.itertuples()]
parts.append(table('Post hoc positive control: discrimination of in-hospital mortality in the first eligible window of '
                   'each patient (69 deaths among 487 patients).', 'etab:positive_control', '@{}lccc@{}',
                   'Variable & Value predicting death & $n$ & AUROC (95' + BS + '% CI)', body,
                   'Bootstrap over patients, 2000 replicates. Complete SOFA exists only in windows with all six components.',
                   size='small'))

# eTable 6 - clinical SMD by outcome
c6 = pd.read_csv(f'{BASE}/eda_clinico_v3/C6_smd_clinica_por_desenlace.csv')
c6.columns = ['variable', 'sofa', 'hemo', 'prog']
ES = {'Edad (años)': 'Age', 'PAM en T0 (mmHg)': 'MAP at T0', 'Frecuencia respiratoria (rpm)': 'Respiratory rate',
      'NED en T0 (µg/kg/min)': 'Norepinephrine equivalent dose', 'VIS en T0': 'Vasoactive-inotropic score',
      'cumVP (12 h)': 'Cumulative vasopressor load', 'Sedación (12 h)': 'Sedation', 'Opioides (12 h)': 'Opioids',
      'Sexo masculino': 'Male sex', 'PAM media 12 h (mmHg)': 'Mean MAP, 12~h', 'Último lactato 12 h (mmol/L)': 'Last lactate, 12~h',
      'SOFA en T0 (componentes disponibles)': 'Partial SOFA', 'Componentes SOFA medidos': 'Number of SOFA components',
      'Vasopresor activo en T0': 'Vasopressor at T0', 'Lactato medido en 12 h': 'Lactate measured'}
body = []
for r in c6.itertuples():
    name = ES.get(r.variable, 'SpO' + BS + 'textsubscript{2}' if r.variable.startswith('SpO') else r.variable)
    body.append(f"{name} & {num(r.sofa, 2)} & {num(r.hemo, 2)} & {num(r.prog, 2)}" + NL)
parts.append(table('Standardised mean differences of clinical predictors between outcome classes (class 1 $-$ class 0).',
                   'etab:smd_clinical', '@{}lrrr@{}', 'Variable & SOFA improvement & Haemodynamic stability or recovery & SOFA deterioration',
                   body, 'Descriptive (windows are not independent). Class 1 = SOFA improvement, haemodynamic stability or recovery and '
                   'deterioration for SOFA deterioration.', size='small'))

open(f'{OUT}/etables.tex', 'w', encoding='utf8').write('\n\n'.join(parts) + '\n')
print('etables written; eTable 2 included:', os.path.exists(ci_path))
