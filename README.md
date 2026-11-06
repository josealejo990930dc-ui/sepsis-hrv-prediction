# sepsis-hrv-prediction

Code and analysis materials for *"Incremental value of heart rate variability for predicting 6-hour clinical improvement in sepsis: development and internal validation of prediction models"* (manuscript submitted to PLOS ONE; master's thesis in Biomedical Informatics, Universidad El Bosque). No patient-level data or credentials are included.

**Question.** In septic ICU patients, does heart rate variability (HRV) from the raw ECG improve prediction of the clinical course over the next 6 hours beyond routine clinical variables and the SOFA score?

**Design.** Retrospective cohort from MIMIC-III (Sepsis-3), with sliding prediction windows adapted from Mentasti et al. (2025): a 12-h observation window, a 6-h prediction horizon and a 6-h step. Each window is identified by an ICU stay and a prediction time T0. There are 3,100 candidate windows from 521 patients; the main analysis uses 2,555 windows from 487 patients with at least 6 valid HRV hours.

**Main result.** No hypothesis was supported. HRV did not add predictive value to the clinical + SOFA model for any outcome (all 249 bootstrap comparisons non-significant).

> **Terminology.** Code, column names and result files keep their original Spanish identifiers so they match the archived results. In particular, the files say `respuesta` / `responder` (response) where the manuscript says **improvement** or **stability or recovery**. See [`docs/data_dictionary.md`](docs/data_dictionary.md) for every mapping.

| Identifier in code and results | Manuscript term |
|---|---|
| `sofa_respuesta` / `label_sofa_responder` | SOFA improvement |
| `hemodinamico` / `label_cardio` | Haemodynamic stability or recovery |
| `progresion` / `label_sofa_deterioro_2pts` | SOFA deterioration (delta SOFA ≥ 2) |

---

## Repository structure

```
sql/                 BigQuery SQL: Sepsis-3 cohort, raw tables, longitudinal tables, diagnostics
python/              Pipeline scripts (HRV extraction, window grid, dataset builder, modelling, bootstrap)
notebooks/           Exploratory analysis and calibration / decision curves / SHAP (outputs cleared)
resultados/          Aggregated result tables of the main analysis ("resultados" = results)
  eda_hrv_v3_1/      HRV exploratory analysis tables
  eda_clinico_v3/    Clinical exploratory analysis tables (incl. Table 1, EPV)
  archivo_v5_sin_ln/ Superseded run without log-transformed LF/HF (reference only)
columnas_v3_1.json   Column groups used by the modelling pipeline
docs/                Data dictionary (Markdown and CSV)
requirements.txt     Python dependencies
```

## Data access

This repository contains **no patient-level data**. Under the MIMIC-III data use agreement, neither the derived dataset (`ultimate_dataset_v3_1.csv`) nor the out-of-fold predictions can be shared. To reproduce the analysis you need:

- **MIMIC-III Clinical Database v1.4** (credentialed access through PhysioNet). The queries use the BigQuery copy `physionet-data.mimiciii_clinical`.
- **MIMIC-III Waveform Database Matched Subset v1.0** (`mimic3wdb-matched/1.0`), read directly from PhysioNet with `wfdb`.

## Environment and dependencies

| Component | Version |
|---|---|
| Python | 3.12.13 (Google Colab, CPU only) |
| NeuroKit2 | 0.2.13 (R-peak detection and HRV metrics) |
| NumPy | 2.0.2 |
| Other Python packages | `requirements.txt` |
| SQL | Google BigQuery (Standard SQL) |

The scripts were written for Google Colab with Google Drive as storage. Paths and the BigQuery project are constants at the top of each script; edit them to run elsewhere:

- `RUTA_BASE = '/content/drive/MyDrive/tesis_sepsis'` (working folder)
- `PROYECTO = 'sepsis-mimic-def'`, `DS = 'sepsis-mimic-def.sepsis3_workspace_def'` (BigQuery project and dataset where the intermediate tables are written)

Scripts that read from BigQuery authenticate through `google.colab.auth`. Outside Colab, replace that call with your own Google Cloud credentials.

## Execution order

Run the steps in this order. Each step writes new files or tables and never overwrites earlier versions.

| # | File | Runs in | Input | Output |
|---|---|---|---|---|
| 1 | `sql/01_cohorte_sepsis3.sql` (blocks 1.1 → 1.8) | BigQuery | MIMIC-III clinical | Sepsis-3 cohort (`sepsis3_final_cohort`) and HRV analysis anchor (`hrv_analysis_mapping`) |
| 2 | `sql/02_tablas_crudas.sql` | BigQuery | MIMIC-III clinical | Raw tables `_raw_labs`, `_raw_gcs`, `_raw_vitals`, `_raw_fio2`, `_raw_lactate` |
| 3 | `python/extraccion_hrv_v3.py` | Colab (runs for days; resumable through a checkpoint) | Step 1 tables + ECG waveforms | `hrv_completo_v3.csv` (one row per valid patient-hour) and `df_reporte_final_unificado.csv` (synchronised cohort, time anchor), plus logs |
| 4 | `sql/03_queries_vigentes_v4.sql`, **part A** | BigQuery | MIMIC-III clinical | Drug infusions as active intervals: `_raw_vaso`, `_raw_sedantes`, `_raw_opioides` |
| 5 | `python/generar_grilla_v3.py` | Colab | Step 3 outputs | Prediction-window grid `hrv_prediction_points` (BigQuery) |
| 6 | `sql/03_queries_vigentes_v4.sql`, **part B** (B0 → B6) | BigQuery | Steps 2, 4, 5 | `vitals_longitudinal`, `drugs_longitudinal`, `sofa_longitudinal`, `cardio_lactato_longitudinal` (one row per window); B0, B5 and B6 are checks |
| 7 | `python/construir_ultimate_dataset_v3.py` | Colab | Step 6 tables + `hrv_completo_v3.csv` | `ultimate_dataset_v3_1.csv` (3,100 × 97), `columnas_v3_1.json`, `flujo_cohorte_v3_1.csv` |
| 8 | `notebooks/EDA_HRV_v3_ventanas.ipynb`, then `notebooks/EDA_clinico_v3.ipynb` | Colab | Step 7 | Descriptive analysis (standardised mean differences, no hypothesis tests) |
| 9 | `python/pipeline_modelado_v5.py` (`RAPIDO = False`) | Colab (CPU, about 1–2 h) | Step 7 | `modelado_v5_1/`: `predicciones_oof.csv`, `metricas_por_particion.csv`, `registro_configuraciones.csv`, `resumen_metricas.csv` |
| 10 | `python/bootstrap_v5.py` | Colab | `predicciones_oof.csv` | `comparaciones_bootstrap.csv` (H1, H2, H3, complementary) |
| 11 | `python/h2_bilateral_v5_1.py` | Colab | `comparaciones_bootstrap.csv` | `comparaciones_h2_bilateral.csv` (two-sided H2; does not rerun the bootstrap) |
| 12 | `notebooks/Calibracion_SHAP_v5.ipynb` | Colab | `predicciones_oof.csv`, dataset | `calibracion_metricas.csv`, `curvas_decision.csv`, `shap_*.csv` |

Read-only checks, which never write tables or files: `sql/04_diagnosticos.sql` and `python/diagnosticos_v3.py` (detector validation, extraction QC).

Note on step 8: `EDA_clinico_v3.ipynb` was run on `ultimate_dataset_v3.csv`, the superseded dataset without log-transformed LF/HF. All clinical columns are identical in v3 and v3.1, so the clinical tables are unaffected. Only the HRV comparison in section C6 reads HRV results from that earlier run.

Rules the pipeline relies on:
- If the grid `hrv_prediction_points` changes, rerun all of part B and the dataset builder.
- Steps 9 onward read only `ultimate_dataset_v3_1.csv` and `columnas_v3_1.json`, never BigQuery.
- Set `RAPIDO = True` in step 9 only to test that the pipeline runs (minimal grids).

## What each key step does

### HRV extraction (`extraccion_hrv_v3.py`)
- ECG lead II at 125 Hz. Each patient is linked to the closest waveform record starting 0–48 h before the HRV anchor. Signal is processed up to 96 h.
- Each hour is split into 12 non-overlapping 5-min epochs. R peaks are detected with NeuroKit2 (`method='neurokit'`) and refined to sub-sample precision by parabolic interpolation.
- An RR interval is flagged as ectopic if it is < 300 ms or > 2000 ms, or deviates > 20 % from the median of up to 2 neighbours on each side. Epochs with > 5 % flagged beats are discarded. Flagged beats in the remaining epochs are replaced by the mean of the nearest valid neighbours.
- 10 metrics per epoch: SDNN, RMSSD, LF, HF, LF/HF, SD1, SD2, SampEn, ApEn and DFA α1 (NeuroKit2 defaults; infinite or negative SampEn makes the epoch invalid).
- An hour is valid if ≥ 8 of its first 11 epochs and the 12th epoch are valid (checked on SDNN, RMSSD and SampEn). Patients need ≥ 3 valid hours.
- Detector validation: Pan-Tompkins at 125 Hz misplaced R peaks and inflated short-term HRV, which is why it was replaced. The check is reproduced by `diagnosticos_v3.py` (D2) and `control_calidad_v3()`.

### Dataset builder (`construir_ultimate_dataset_v3.py`)
- **Hourly HRV value**: median of the valid 5-min epochs. LF, HF and LF/HF are natural-log transformed.
- **Window features**: for each metric, the mean, SD, last value and OLS slope over the 12-h observation window, plus hours since the last valid HRV hour (41 HRV features).
- **Clinical features at T0**: age, sex, MAP, respiratory rate, SpO₂, norepinephrine-equivalent dose, vasoactive-inotropic score, cumulative vasopressor exposure, sedative and opioid infusion flags, mean MAP and last lactate in the previous 12 h, a lactate-measured flag, partial SOFA and the number of SOFA components available.
- **Outcomes** (horizon (T0, T0 + 6 h]):
  - *SOFA improvement*: paired delta SOFA < 0. The delta is computed component by component; coagulation, cardiovascular, CNS and renal are required at both times.
  - *Haemodynamic stability or recovery*: requires MAP and lactate in both windows, and both of:
    - MAP criterion: > 50 % of MAP values > 65 mmHg, or mean MAP up ≥ 10 %;
    - lactate criterion: < 2 mmol/L in both windows, or 2–4 → < 2, or > 4 → ≤ 4.
  - *SOFA deterioration*: paired delta SOFA ≥ 2.
  - Death in the horizon = no improvement / no stability and deterioration. Discharge alive before T0 + 6 h = excluded.
- Leakage control: variables measured after T0 are listed in `desenlace_no_usar_como_feature` and blocked by an `assert` in the modelling pipeline.

### Modelling pipeline (`pipeline_modelado_v5.py`)
- 7 feature sets × 5 model families (logistic regression with a Welch filter, Elastic Net, random forest, XGBoost, PCA + logistic regression). That gives 30 configurations per outcome: 90 in the main analysis and 96 in the sensitivity analyses.
- Nested cross-validation grouped by patient: `StratifiedGroupKFold`, 5 outer × 3 inner folds, seed 42. Every model of an outcome uses the same outer folds, so comparisons are paired.
- Grid search on AUROC (main outcomes) or AUPRC (deterioration). Sample weight = 1 / windows per patient. No resampling or class weights.
- Platt calibration fitted on inner-fold predictions.
- Imputation, scaling, filtering and PCA are fitted inside each training fold only.

### Bootstrap comparisons (`bootstrap_v5.py`, `h2_bilateral_v5_1.py`)
- Patient-level bootstrap: 2,000 replicates resampling patients with all their windows, seed 2026.
- 95 % percentile CIs. One-sided p for H1, H3 and the complementary contrast; two-sided p for H2, as in the research proposal.
- Benjamini–Hochberg correction across model families within each hypothesis and outcome.
- A comparison counts as significant when the CI excludes 0 and the adjusted p < 0.05.

## Results

`resultados/` contains the aggregated tables of the main analysis (`modelado_v5_1`, dataset v3.1). Column definitions are in [`docs/data_dictionary.md`](docs/data_dictionary.md), section 4. The out-of-fold prediction files (`predicciones_oof.csv`, `shap_oof.csv`) are not included because they are window-level data.

`resultados/eda_hrv_v3_1/` and `resultados/eda_clinico_v3/` contain the descriptive tables of step 8 (all aggregated; the figures are not included). `resultados/flujo_cohorte_v3_1.csv` gives the counts for the cohort flow diagram.

`resultados/archivo_v5_sin_ln/` is an earlier run in which LF, HF and LF/HF were not log-transformed, contrary to the protocol. It was superseded by the current run and is kept only for transparency.

## Reporting

The study is reported following TRIPOD+AI. Risk of bias was self-assessed with PROBAST+AI.

## Citation

If you use this code, please cite the article: *Incremental value of heart rate variability for predicting 6-hour clinical improvement in sepsis: development and internal validation of prediction models* (PLOS ONE, under review; full reference to be added on publication).

## License

Code released under the MIT License (see `LICENSE`). MIMIC-III data are not covered by this license and remain subject to the PhysioNet credentialed data use agreement.
