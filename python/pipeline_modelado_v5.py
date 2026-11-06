# -*- coding: utf-8 -*-
"""
pipeline_modelado_v5.py  (modelling pipeline)
=============================================
Trains and evaluates the thesis models (Mentasti design, dataset v3.1) with
nested cross-validation grouped by patient, and saves the OUT-OF-FOLD
PREDICTIONS of every window. Model comparisons (patient-level bootstrap) are
done afterwards by bootstrap_v5.py from those predictions.

Identifiers are in Spanish; English terms are in docs/data_dictionary.md.

Design
------
- Outcomes: SOFA improvement and haemodynamic stability or recovery (main);
  SOFA progression, delta SOFA >= 2 (secondary).
- Feature sets (7): full HRV, linear HRV, non-linear HRV, clinical, clinical +
  SOFA, integrated HRV, integrated HRV + SOFA. PCA + logistic regression uses
  only its two specific sets (clinical + SOFA and integrated + SOFA).
- Model families (5): logistic regression (with Welch t-test filter, p < 0.20),
  Elastic Net, Random Forest, XGBoost, PCA + logistic regression.
- Validation: StratifiedGroupKFold by patient, 5 outer and 3 inner folds. The
  outer folds are THE SAME for every model of one outcome, so comparisons are paired.
- Hyperparameters: grid search in the inner folds, optimising the metric of the
  outcome's hypothesis (AUROC for the main outcomes, AUPRC for progression).
- Calibration: Platt scaling fitted on inner-fold predictions. It does not change
  AUROC or AUPRC (both are computed on the uncalibrated probabilities); it does
  change the Brier score (computed on the calibrated probabilities).
- Weighting: each window weighs 1 / (number of windows of its patient),
  computed inside the training set.
- No imbalance correction (neither SMOTE nor class weights): van den
  Goorbergh et al., JAMIA 2022.
- Imputation (median), scaling, univariate filter and calibration are fitted
  only on the training data of each fold (no leakage).

Analyses
--------
- principal: windows with >= 6 h of HRV; all sets and families.
- Sensitivity (core comparison only, clinical + SOFA versus integrated HRV + SOFA,
  with the 4 non-PCA families): hrv_4h, hrv_8h, sofa_caso_completo
  (complete-case SOFA), laboratorio_nuevo (new lab in horizon), cohorte_comun
  (common cohort evaluable for both main outcomes).

Outputs (folder modelado_v5_1 on Drive)
---------------------------------------
- predicciones_oof.csv : one row per window x configuration, with the
  uncalibrated and calibrated probability and the outer fold.
- metricas_por_particion.csv : AUROC, AUPRC and Brier per outer fold.
- registro_configuraciones.csv : finished configurations, best hyperparameters
  per fold and run time. Doubles as a checkpoint: if Colab disconnects,
  finished configurations are skipped when the script is run again.
- resumen_metricas.csv : metrics on the pooled out-of-fold predictions
  (the ones that are reported).

Environment: CPU is enough (about 2,500 windows x <= 60 variables). A GPU does
not speed anything up at this size.
"""

import os, json, time, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from xgboost import XGBClassifier

warnings.filterwarnings('ignore')

# ─── Configuration ──────────────────────────────────────────────────────────
RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'
RUTA_DATASET = f'{RUTA_BASE}/ultimate_dataset_v3_1.csv'
RUTA_COLUMNAS = f'{RUTA_BASE}/columnas_v3_1.json'
RUTA_SALIDA = f'{RUTA_BASE}/modelado_v5_1'

SEMILLA = 42
N_EXTERNAS = 5
N_INTERNAS = 3
N_JOBS = -1                    # parallelism of the hyperparameter search
RAPIDO = False                 # True = minimal grids, only to check that everything runs

DESENLACES = {
    'sofa_respuesta': {'label': 'label_sofa_responder',      'metrica': 'roc_auc'},
    'hemodinamico':   {'label': 'label_cardio',              'metrica': 'roc_auc'},
    'progresion':     {'label': 'label_sofa_deterioro_2pts', 'metrica': 'average_precision'},
}

CENTRALES = ['clinico_sofa', 'integrada_sofa']
FAMILIAS_SENSIBILIDAD = ['logistica', 'elastic_net', 'random_forest', 'xgboost']

# analysis: window filter, label per outcome (None = does not apply), feature sets and families
ANALISIS = {
    'principal':          {'filtro': 'incl_hrv_6h', 'labels': {}, 'conjuntos': None, 'familias': None},
    'hrv_4h':             {'filtro': 'incl_hrv_4h', 'labels': {}, 'conjuntos': CENTRALES, 'familias': FAMILIAS_SENSIBILIDAD},
    'hrv_8h':             {'filtro': 'incl_hrv_8h', 'labels': {}, 'conjuntos': CENTRALES, 'familias': FAMILIAS_SENSIBILIDAD},
    'sofa_caso_completo': {'filtro': 'incl_hrv_6h',
                           'labels': {'sofa_respuesta': 'label_sofa_responder_cc',
                                      'progresion': 'label_sofa_deterioro_2pts_cc', 'hemodinamico': None},
                           'conjuntos': CENTRALES, 'familias': FAMILIAS_SENSIBILIDAD},
    'laboratorio_nuevo':  {'filtro': 'incl_hrv_6h',
                           'labels': {'sofa_respuesta': 'label_sofa_responder_nl',
                                      'progresion': 'label_sofa_deterioro_2pts_nl', 'hemodinamico': None},
                           'conjuntos': CENTRALES, 'familias': FAMILIAS_SENSIBILIDAD},
    'cohorte_comun':      {'filtro': 'incl_hrv_6h', 'extra': 'eval_comun',
                           'labels': {'progresion': None},
                           'conjuntos': CENTRALES, 'familias': FAMILIAS_SENSIBILIDAD},
}


# ─── Univariate filter (logistic regression only) ───────────────────────────
class FiltroWelch(BaseEstimator, TransformerMixin):
    """Keep the variables with p < threshold in the Welch t-test between classes,
    computed on the training data only. If none passes, keep the one with the
    smallest p."""
    def __init__(self, umbral=0.20):
        self.umbral = umbral

    def fit(self, X, y):
        X = np.asarray(X, dtype=float); y = np.asarray(y)
        _, p = stats.ttest_ind(X[y == 1], X[y == 0], equal_var=False, nan_policy='omit')
        p = np.where(np.isfinite(p), p, 1.0)
        self.mascara_ = p < self.umbral
        if not self.mascara_.any():
            self.mascara_[np.argmin(p)] = True
        return self

    def transform(self, X):
        return np.asarray(X, dtype=float)[:, self.mascara_]


# ─── Model families and grids ───────────────────────────────────────────────
def construir_familias(rapido=False):
    """Return {family name: (sklearn Pipeline, hyperparameter grid)}.

    Every pipeline starts with median imputation; scaling is applied to the
    linear families only (trees are scale-invariant). The step named 'model'
    receives the per-patient sample weights."""
    lr_sin_pen = LogisticRegression(penalty=None, solver='lbfgs', max_iter=5000)
    en = LogisticRegression(penalty='elasticnet', solver='saga', max_iter=5000, random_state=SEMILLA)
    rf = RandomForestClassifier(n_estimators=100 if rapido else 300, n_jobs=1, random_state=SEMILLA)
    xgb = XGBClassifier(n_estimators=100 if rapido else 300, subsample=0.8, colsample_bytree=0.8,
                        tree_method='hist', eval_metric='logloss', n_jobs=1, random_state=SEMILLA,
                        verbosity=0)
    pca_lr = LogisticRegression(max_iter=5000)

    fam = {
        'logistica': (Pipeline([('imp', SimpleImputer(strategy='median')), ('esc', StandardScaler()),
                                ('filtro', FiltroWelch(0.20)), ('model', lr_sin_pen)]),
                      {}),
        'elastic_net': (Pipeline([('imp', SimpleImputer(strategy='median')), ('esc', StandardScaler()),
                                  ('model', en)]),
                        {'model__C': [0.01, 0.1, 1.0, 10.0], 'model__l1_ratio': [0.1, 0.5, 0.9]}),
        'random_forest': (Pipeline([('imp', SimpleImputer(strategy='median')), ('model', rf)]),
                          {'model__max_depth': [3, 6, None], 'model__min_samples_leaf': [5, 20],
                           'model__max_features': ['sqrt', 0.3]}),
        'xgboost': (Pipeline([('imp', SimpleImputer(strategy='median')), ('model', xgb)]),
                    {'model__max_depth': [2, 3, 4], 'model__learning_rate': [0.03, 0.1],
                     'model__min_child_weight': [1, 5]}),
        'pca_logistica': (Pipeline([('imp', SimpleImputer(strategy='median')), ('esc', StandardScaler()),
                                    ('pca', PCA(svd_solver='full', random_state=SEMILLA)), ('model', pca_lr)]),
                          {'pca__n_components': [0.80, 0.90, 0.95]}),
    }
    if rapido:
        fam['elastic_net'] = (fam['elastic_net'][0], {'model__C': [0.1], 'model__l1_ratio': [0.5]})
        fam['random_forest'] = (fam['random_forest'][0], {'model__max_depth': [3]})
        fam['xgboost'] = (fam['xgboost'][0], {'model__max_depth': [2]})
        fam['pca_logistica'] = (fam['pca_logistica'][0], {'pca__n_components': [0.9]})
    return fam


# ─── Data and feature sets ──────────────────────────────────────────────────
def cargar_datos():
    """Load the dataset and build the 7 feature sets from columnas_v3_1.json.

    Includes the leakage barrier: no outcome column, label, inclusion flag or
    identifier may appear in any feature set (checked with assert)."""
    df = pd.read_csv(RUTA_DATASET)
    with open(RUTA_COLUMNAS) as f:
        cols = json.load(f)

    # Sex as a numeric variable (1 = male)
    df['sexo_masculino'] = df['gender'].astype(str).str.upper().str.startswith('M').astype(int)
    # Idempotent: if the dataset already carries them, nothing changes
    for c in ['map_t0', 'mbp_mean_actual12h']:
        df.loc[(df[c] < 20) | (df[c] > 200), c] = np.nan
    if 'lactato_medido_12h' not in df.columns:
        df['lactato_medido_12h'] = df['lactate_last_actual12h'].notna().astype(int)

    clin = [c if c != 'gender' else 'sexo_masculino' for c in cols['features_clinicas']]
    pamlac = list(dict.fromkeys(cols['features_pam_lactato_actual'] + ['lactato_medido_12h']))
    sofa = cols['features_sofa']
    hrv = cols['features_hrv']
    basal = clin + pamlac

    conjuntos = {
        'hrv_completa': hrv,
        'hrv_lineal': cols['features_hrv_lineal'],
        'hrv_nolineal': cols['features_hrv_nolineal'],
        'clinico_basal': basal,
        'clinico_sofa': basal + sofa,
        'hrv_integrada': hrv + basal,
        'integrada_sofa': hrv + basal + sofa,
    }

    # Leakage barrier: no outcome column or label as a predictor
    prohibidas = set(cols['desenlace_no_usar_como_feature'])
    prohibidas |= {c for c in df.columns if c.startswith('label_') or c.startswith('incl_')}
    prohibidas |= {'eval_comun', 'icustay_id', 'subject_id', 'hour_index', 't0_prediccion'}
    for nombre, vs in conjuntos.items():
        fuga = prohibidas.intersection(vs)
        assert not fuga, f'Fuga de información en {nombre}: {fuga}'
        faltan = [v for v in vs if v not in df.columns]
        assert not faltan, f'Columnas inexistentes en {nombre}: {faltan}'
    return df, conjuntos


def familias_por_conjunto(familias_permitidas=None):
    """(feature set, family) combinations. PCA only on its two sets."""
    combos = []
    for cj in ['hrv_completa', 'hrv_lineal', 'hrv_nolineal', 'clinico_basal',
               'clinico_sofa', 'hrv_integrada', 'integrada_sofa']:
        for fam in ['logistica', 'elastic_net', 'random_forest', 'xgboost']:
            combos.append((cj, fam))
    combos += [('clinico_sofa', 'pca_logistica'), ('integrada_sofa', 'pca_logistica')]
    if familias_permitidas is not None:
        combos = [c for c in combos if c[1] in familias_permitidas]
    return combos


def pesos_por_paciente(grupos):
    """Sample weight of each window = 1 / (number of windows of its patient)."""
    g = pd.Series(grupos)
    return (1.0 / g.map(g.value_counts())).to_numpy()


# ─── Platt calibration with inner-fold predictions ──────────────────────────
def _logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def ajustar_platt(estimador, X, y, grupos, pesos, cv_interna):
    """Fit a Platt scaler (logistic regression on the logit of the raw
    probability) using out-of-fold predictions from the INNER folds of the
    training data, so the outer test fold is never used for calibration."""
    oof = np.full(len(y), np.nan)
    for tr, te in cv_interna.split(X, y, grupos):
        m = clone(estimador)
        m.fit(X[tr], y[tr], model__sample_weight=pesos[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    platt = LogisticRegression(C=1e6, max_iter=1000)
    platt.fit(_logit(oof).reshape(-1, 1), y)
    return platt


# ─── One configuration (analysis x outcome x feature set x family) ──────────
def correr_configuracion(datos, label, columnas, familia, metrica, particiones, familias):
    """Nested cross-validation of one model configuration.

    For each outer fold: grid search on the inner folds -> refit the best model
    on the whole training fold -> Platt calibration from inner-fold predictions
    -> predict the outer test fold.
    Returns (out-of-fold predictions, metrics per fold, best hyperparameters per fold)."""
    X = datos[columnas].to_numpy(dtype=float)
    y = datos[label].to_numpy(dtype=int)
    grupos = datos['icustay_id'].to_numpy()
    base, grilla = familias[familia]

    filas_pred, filas_met, mejores = [], [], []
    for k, (tr, te) in enumerate(particiones):
        pesos = pesos_por_paciente(grupos[tr])
        cv_int = StratifiedGroupKFold(n_splits=N_INTERNAS, shuffle=True, random_state=SEMILLA)
        busqueda = GridSearchCV(clone(base), grilla, scoring=metrica, cv=cv_int, n_jobs=N_JOBS, refit=True)
        busqueda.fit(X[tr], y[tr], groups=grupos[tr], model__sample_weight=pesos)
        mejor = busqueda.best_estimator_
        platt = ajustar_platt(busqueda.best_estimator_, X[tr], y[tr], grupos[tr], pesos, cv_int)

        p_raw = mejor.predict_proba(X[te])[:, 1]
        p_cal = platt.predict_proba(_logit(p_raw).reshape(-1, 1))[:, 1]
        y_te = y[te]
        filas_pred.append(pd.DataFrame({
            'icustay_id': grupos[te], 'hour_index': datos['hour_index'].to_numpy()[te],
            'particion': k, 'y': y_te, 'p_sin_calibrar': p_raw, 'p_calibrada': p_cal}))
        filas_met.append({'particion': k, 'n': len(te), 'eventos': int(y_te.sum()),
                          'auroc': roc_auc_score(y_te, p_raw) if len(np.unique(y_te)) > 1 else np.nan,
                          'auprc': average_precision_score(y_te, p_raw) if y_te.sum() > 0 else np.nan,
                          'brier': brier_score_loss(y_te, p_cal)})
        mejores.append(busqueda.best_params_)
    return pd.concat(filas_pred, ignore_index=True), pd.DataFrame(filas_met), mejores


def _anexar(df, ruta):
    """Append a DataFrame to a CSV (writes the header only if the file is new)."""
    df.to_csv(ruta, mode='a', header=not os.path.exists(ruta), index=False)


# ─── Execution ──────────────────────────────────────────────────────────────
def main(analisis_a_correr=None):
    """Run every analysis x outcome x (feature set, family) configuration.

    analisis_a_correr: optional list of analysis names (keys of ANALISIS) to run;
    None runs all of them. Finished configurations are skipped (checkpoint)."""
    try:
        from google.colab import drive
        drive.mount('/content/drive')
    except ImportError:
        pass
    os.makedirs(RUTA_SALIDA, exist_ok=True)
    r_pred = f'{RUTA_SALIDA}/predicciones_oof.csv'
    r_met = f'{RUTA_SALIDA}/metricas_por_particion.csv'
    r_reg = f'{RUTA_SALIDA}/registro_configuraciones.csv'

    df, conjuntos = cargar_datos()
    familias = construir_familias(RAPIDO)
    hechos = set()
    if os.path.exists(r_reg):
        reg = pd.read_csv(r_reg)
        hechos = set(zip(reg['analisis'], reg['desenlace'], reg['conjunto'], reg['familia']))
        print(f'Retomando: {len(hechos)} configuraciones ya terminadas.')

    for an, cfg in ANALISIS.items():
        if analisis_a_correr is not None and an not in analisis_a_correr:
            continue
        for des, dcfg in DESENLACES.items():
            label = cfg['labels'].get(des, dcfg['label'])
            if label is None:
                continue
            datos = df[(df[cfg['filtro']] == 1) & df[label].notna()]
            if cfg.get('extra'):
                datos = datos[datos[cfg['extra']] == 1]
            datos = datos.reset_index(drop=True)
            y = datos[label].astype(int)
            cv_ext = StratifiedGroupKFold(n_splits=N_EXTERNAS, shuffle=True, random_state=SEMILLA)
            particiones = list(cv_ext.split(datos, y, datos['icustay_id']))
            print(f'\n=== {an} · {des} ({label}): {len(datos)} ventanas, '
                  f'{datos.icustay_id.nunique()} pacientes, prevalencia {100 * y.mean():.1f} % ===')

            combos = familias_por_conjunto(cfg['familias'])
            if cfg['conjuntos'] is not None:
                combos = [c for c in combos if c[0] in cfg['conjuntos']]
            for cj, fam in combos:
                if (an, des, cj, fam) in hechos:
                    continue
                t = time.time()
                pred, met, mejores = correr_configuracion(datos, label, conjuntos[cj], fam,
                                                          dcfg['metrica'], particiones, familias)
                ident = {'analisis': an, 'desenlace': des, 'conjunto': cj, 'familia': fam}
                _anexar(pred.assign(**ident), r_pred)
                _anexar(met.assign(**ident), r_met)
                _anexar(pd.DataFrame([{**ident, 'label': label, 'n_variables': len(conjuntos[cj]),
                                       'mejores_hiperparametros': json.dumps(mejores, default=str),
                                       'segundos': round(time.time() - t, 1)}]), r_reg)
                print(f'  {cj:15s} {fam:14s} AUROC {met.auroc.mean():.3f} ± {met.auroc.std():.3f} | '
                      f'AUPRC {met.auprc.mean():.3f} | Brier {met.brier.mean():.3f} | {time.time() - t:.0f}s')

    resumir()


def resumir():
    """Metrics on the pooled out-of-fold predictions (the ones that are reported).
    Writes resumen_metricas.csv and prints the AUROC matrix of the main analysis."""
    r_pred = f'{RUTA_SALIDA}/predicciones_oof.csv'
    pred = pd.read_csv(r_pred)
    filas = []
    for clave, g in pred.groupby(['analisis', 'desenlace', 'conjunto', 'familia']):
        filas.append(dict(zip(['analisis', 'desenlace', 'conjunto', 'familia'], clave),
                          n=len(g), pacientes=g.icustay_id.nunique(), prevalencia=round(g.y.mean(), 3),
                          auroc=roc_auc_score(g.y, g.p_sin_calibrar),
                          auprc=average_precision_score(g.y, g.p_sin_calibrar),
                          brier=brier_score_loss(g.y, g.p_calibrada)))
    res = pd.DataFrame(filas).round(4)
    res.to_csv(f'{RUTA_SALIDA}/resumen_metricas.csv', index=False)
    prin = res[res.analisis == 'principal']
    for des in prin.desenlace.unique():
        t = prin[prin.desenlace == des].pivot(index='conjunto', columns='familia', values='auroc')
        print(f'\nAUROC fuera de muestra — {des}:')
        print(t.round(3).to_string())
    return res


if __name__ == '__main__':
    main()
