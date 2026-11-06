"""English-language manuscript figures generated from the saved results (no model refitting).

Usage: python make_figures.py <tesis_sepsis dir> <figures dir>
Outputs (PNG, 300 dpi):
  fig2_forest.png             prespecified comparisons, delta with 95% CI (main text)
  efig_decision_curves.png    decision curves, 5 families x 3 outcomes
  efig_calibration.png        calibration curves from calibrated OOF predictions
  efig_shap_groups.png        share of mean |SHAP| by predictor group (integrated + SOFA)
Palette: validated categorical slots (dataviz reference palette), checked with validate_palette.js.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

BASE, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
M = f'{BASE}/modelado_v5_1'
BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, GRID = '#1f1f1e', '#5c5b55', '#d9d8d2'
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.edgecolor': INK2,
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 100})

FAM = {'logistica': 'Logistic', 'elastic_net': 'Elastic net', 'random_forest': 'Random forest',
       'xgboost': 'XGBoost', 'pca_logistica': 'PCA + logistic'}
FAM_ORDER = ['logistica', 'elastic_net', 'random_forest', 'xgboost', 'pca_logistica']
OUTC = {'sofa_respuesta': 'SOFA improvement', 'hemodinamico': 'Haemodynamic stability or recovery',
        'progresion': 'SOFA deterioration'}
OUTC_ORDER = ['sofa_respuesta', 'hemodinamico', 'progresion']
SETS = {'clinico_sofa': 'Clinical + SOFA', 'integrada_sofa': 'Integrated + SOFA'}

# ---------------- Figure 2: forest plot ----------------
b = pd.read_csv(f'{M}/comparaciones_bootstrap.csv')
b = b[(b.analisis == 'principal') & (b.principal)]
panels = [('H1', 'sofa_respuesta', 'H1 · SOFA improvement (ΔAUROC)'),
          ('H1', 'hemodinamico', 'H1 · Haemodynamic stability/recovery (ΔAUROC)'),
          ('H2', 'sofa_respuesta', 'H2 · SOFA improvement (ΔAUROC, NL − L)'),
          ('H2', 'hemodinamico', 'H2 · Haemodynamic stability/recovery (ΔAUROC, NL − L)'),
          ('H3', 'progresion', 'H3 · SOFA deterioration (ΔAUPRC)'),
          ('Complementario', 'sofa_respuesta', 'Compl. · SOFA improvement (ΔAUROC)'),
          ('Complementario', 'hemodinamico', 'Compl. · Haemodynamic stability/recovery (ΔAUROC)'),
          ('Complementario', 'progresion', 'Compl. · SOFA deterioration (ΔAUPRC)')]
fig, axes = plt.subplots(4, 2, figsize=(7.2, 8.2), sharex=False)
for ax, (h, d, title) in zip(axes.T.flatten(), panels):
    g = b[(b.hipotesis == h) & (b.desenlace == d)].set_index('familia')
    fams = [f for f in FAM_ORDER if f in g.index][::-1]
    for i, f in enumerate(fams):
        r = g.loc[f]
        ax.plot([r.ic95_inf, r.ic95_sup], [i, i], color=INK, lw=1.5, solid_capstyle='round')
        ax.plot(r.delta, i, 'o', color=INK, ms=5, mec='white', mew=1)
    ax.axvline(0, color=INK2, lw=1, ls=(0, (3, 2)))
    ax.set_yticks(range(len(fams)), [FAM[f] for f in fams])
    lim = max(0.06, np.ceil(max(abs(g.ic95_inf).max(), abs(g.ic95_sup).max()) * 50) / 50)
    ax.set_xlim(-lim, lim)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, symmetric=True))
    ax.set_title(title, fontsize=8, loc='left', color=INK)
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
for ax in axes[-1]:
    ax.set_xlabel('Difference (new − reference), 95% CI')
fig.text(0.01, 0.005, 'Values left of the dashed line favour the reference model. H2: NL − L, nonlinear minus linear HRV. Compl., complementary comparison. '
         'Patient-level bootstrap, 2000 replicates.', fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.015, 1, 1))
fig.savefig(f'{OUT}/fig2_forest.png', dpi=300)
plt.close(fig)

# ---------------- eFigure: decision curves ----------------
dc = pd.read_csv(f'{M}/calibracion_shap/curvas_decision.csv')
fig, axes = plt.subplots(5, 3, figsize=(7.2, 9.6))
for j, d in enumerate(OUTC_ORDER):
    for i, f in enumerate(FAM_ORDER):
        ax = axes[i, j]
        g = dc[(dc.desenlace == d) & (dc.familia == f)]
        for cj, col in [('clinico_sofa', BLUE), ('integrada_sofa', ORANGE)]:
            s = g[g.conjunto == cj].sort_values('umbral')
            ax.plot(s.umbral, s.beneficio_neto, color=col, lw=1.5, label=SETS[cj])
        s = g[g.conjunto == 'clinico_sofa'].sort_values('umbral')
        ax.plot(s.umbral, s.intervenir_a_todos, color=INK2, lw=1, ls=(0, (3, 2)), label='Treat all')
        ax.axhline(0, color=INK, lw=0.8, label='Treat none')
        top = g.beneficio_neto.max()
        ax.set_ylim(-0.1 * top, top * 1.08)
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        if i == 0:
            ev = {'sofa_respuesta': 'no improvement', 'hemodinamico': 'neither stable nor recovered'}.get(d, 'deterioration')
            ax.set_title(f'{OUTC[d]}\n(event = {ev})', fontsize=8)
        if j == 0:
            ax.set_ylabel(f'{FAM[f]}\nNet benefit')
        if i == 4:
            ax.set_xlabel('Threshold probability')
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(f'{OUT}/efig_decision_curves.png', dpi=300)
plt.close(fig)

# ---------------- eFigure: calibration curves ----------------
pred = pd.read_csv(f'{M}/predicciones_oof.csv')
pred = pred[(pred.analisis == 'principal') & pred.conjunto.isin(SETS)]

def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    hw = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return c - hw, c + hw

fig, axes = plt.subplots(5, 3, figsize=(7.2, 10.2))
for j, d in enumerate(OUTC_ORDER):
    q = 5 if d == 'progresion' else 10
    for i, f in enumerate(FAM_ORDER):
        ax = axes[i, j]
        g = pred[(pred.desenlace == d) & (pred.familia == f)]
        mx = 0
        for k, (cj, col) in enumerate([('clinico_sofa', BLUE), ('integrada_sofa', ORANGE)]):
            s = g[g.conjunto == cj]
            if s.empty:
                continue
            p, y = s.p_calibrada.to_numpy(), s.y.to_numpy()
            bins = pd.qcut(p, q, duplicates='drop')
            t = pd.DataFrame({'p': p, 'y': y, 'b': bins}).groupby('b', observed=True).agg(
                p=('p', 'mean'), k=('y', 'sum'), n=('y', 'size'))
            lo, hi = wilson(t.k.to_numpy(), t.n.to_numpy())
            off = (k - 0.5) * 0.004
            ax.errorbar(t.p + off, t.k / t.n, yerr=[t.k / t.n - lo, hi - t.k / t.n], fmt='o', ms=3.5,
                        color=col, mec='white', mew=0.6, elinewidth=1, capsize=0, label=SETS[cj])
            sm = lowess(y, p, frac=0.6, it=0, return_sorted=True)  # it=0: robustness iterations break on binary y
            mx = max(mx, t.p.max(), hi.max())
            sm = sm[sm[:, 0] <= t.p.max()]
            ax.plot(sm[:, 0], sm[:, 1], color=col, lw=1.3)
        mx = min(1, mx * 1.1)
        ax.plot([0, mx], [0, mx], color=INK2, lw=0.8, ls=(0, (3, 2)))
        ax.set_xlim(0, mx)
        ax.set_ylim(0, mx)
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        if i == 0:
            ax.set_title(OUTC[d], fontsize=8)
        if j == 0:
            ax.set_ylabel(f'{FAM[f]}\nObserved proportion')
        if i == 4:
            ax.set_xlabel('Predicted probability (calibrated)')
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=2, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.025, 1, 1))
fig.savefig(f'{OUT}/efig_calibration.png', dpi=300)
plt.close(fig)

# ---------------- eFigure: SHAP group shares ----------------
sg = pd.read_csv(f'{M}/calibracion_shap/shap_pct_por_grupo.csv')
GROUPS = [('Clínicas en T0', 'Clinical at T0'), ('PAM y lactato (12 h)', 'MAP and lactate (12 h)'),
          ('SOFA', 'SOFA'), ('HRV temporal', 'HRV time domain'), ('HRV frecuencial', 'HRV frequency domain'),
          ('HRV no lineal', 'HRV nonlinear')]
rows = [(d, f) for d in OUTC_ORDER for f in ['elastic_net', 'xgboost']]
fig, ax = plt.subplots(figsize=(7.2, 3.2))
for yi, (d, f) in enumerate(rows[::-1]):
    r = sg[(sg.desenlace == d) & (sg.familia == f)].iloc[0]
    left = 0
    for (col_es, lab), c in zip(GROUPS, SLOTS):
        w = float(r[col_es])
        ax.barh(yi, w, left=left, color=c, edgecolor='white', linewidth=1.2, height=0.62,
                label=lab if yi == 0 else None)
        if w >= 6:
            ax.text(left + w / 2, yi, f'{w:.0f}', ha='center', va='center', fontsize=8, color=INK)
        left += w
ax.set_yticks(range(len(rows)), [f'{OUTC[d]} · {FAM[f]}' for d, f in rows[::-1]])
ax.set_xlim(0, 100)
ax.set_xlabel('Share of mean |SHAP| (%), integrated + SOFA models')
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.22), ncol=3, frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f'{OUT}/efig_shap_groups.png', dpi=300)
plt.close(fig)
print('figures written to', OUT)
