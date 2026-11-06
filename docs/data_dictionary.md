# Data dictionary

This dictionary maps the Spanish identifiers used in the code, in `columnas_v3_1.json` and in the results files (`resultados/`) to the English terms used in the manuscript. The identifiers were left in Spanish on purpose: renaming them would break the link between the code and the archived results.

A machine-readable version of section 3 is in [`data_dictionary_columns.csv`](data_dictionary_columns.csv).

English terms follow the manuscript (outcomes; Table 2 for predictor sets; Tables 4–6).

## 1. Outcomes (most important for reading Tables 4–6)

Wherever the results files say *respuesta* / *responder* (response), the manuscript says **improvement** or **stability or recovery**.

| Value in `desenlace` column | Label column in the dataset | Manuscript term | Definition (6-h horizon after T0) | Primary metric |
|---|---|---|---|---|
| `sofa_respuesta` | `label_sofa_responder` | **SOFA improvement** | Paired delta SOFA < 0 | AUROC |
| `hemodinamico` | `label_cardio` | **Haemodynamic stability or recovery** | MAP criterion **and** lactate criterion (README, *Outcomes*) | AUROC |
| `progresion` | `label_sofa_deterioro_2pts` | **SOFA deterioration** (delta SOFA ≥ 2) | Paired delta SOFA ≥ 2 | AUPRC |

Rules shared by all three: death within the horizon counts as no improvement / no stability (0) and as deterioration (1). A window is excluded (label missing) if the patient was discharged alive from the ICU before T0 + 6 h.

**Event direction in the decision curves.** In `curvas_decision.csv`, for `sofa_respuesta` and `hemodinamico` the event is **no** improvement / **no** stability (the notebook uses `1 - y` and `1 - p`). For `progresion` the event is deterioration. Everywhere else, class 1 is the favourable outcome for the two main outcomes.

## 2. Values used in the results files

### 2.1 Feature sets (`conjunto`, `referencia`, `nuevo`)

| Value | English name | No. of variables | Content |
|---|---|---|---|
| `hrv_completa` | HRV, all | 41 | 10 HRV metrics × 4 summaries + `horas_desde_ultima_hrv` |
| `hrv_lineal` | HRV, linear | 28 | SDNN, RMSSD, SD1, SD2, LF, HF, LF/HF × 4 summaries |
| `hrv_nolineal` | HRV, nonlinear | 12 | SampEn, ApEn, DFA α1 × 4 summaries |
| `clinico_basal` | Clinical baseline | 13 | Clinical variables at T0 + MAP/lactate in the previous 12 h (no SOFA) |
| `clinico_sofa` | Clinical + SOFA | 15 | `clinico_basal` + partial SOFA + number of SOFA components available |
| `hrv_integrada` | HRV + clinical | 54 | `hrv_completa` + `clinico_basal` |
| `integrada_sofa` | Integrated + SOFA | 56 | `hrv_completa` + `clinico_sofa` |

### 2.2 Model families (`familia`)

| Value | English name |
|---|---|
| `logistica` | LR: logistic regression (univariate Welch t-test filter, p < 0.20, inside each training fold) |
| `elastic_net` | EN: elastic net logistic regression |
| `random_forest` | RF: random forest |
| `xgboost` | XGB: XGBoost |
| `pca_logistica` | PCA-LR: principal components + logistic regression (only with `clinico_sofa` and `integrada_sofa`) |

### 2.3 Analyses (`analisis`)

| Value | English name | Change with respect to the main analysis |
|---|---|---|
| `principal` | Main analysis | Windows with ≥ 6 valid HRV hours (threshold fixed a priori) |
| `hrv_4h` | Sensitivity: ≥ 4 h HRV | HRV threshold ≥ 4 h |
| `hrv_8h` | Sensitivity: ≥ 8 h HRV | HRV threshold ≥ 8 h |
| `sofa_caso_completo` | Sensitivity: complete-case SOFA | SOFA labels require all six components at both times (`*_cc` labels) |
| `laboratorio_nuevo` | Sensitivity: new laboratory value | SOFA labels only in windows with a new creatinine, bilirubin or platelet value, or death, in the horizon (`*_nl` labels) |
| `cohorte_comun` | Sensitivity: common cohort | Windows evaluable for both main outcomes (`eval_comun = 1`) |

### 2.4 Hypotheses (`hipotesis`)

| Value | Contrast (new − reference, same family) | Outcomes | Metric | Test |
|---|---|---|---|---|
| `H1` | `integrada_sofa` − `clinico_sofa` | SOFA improvement, haemodynamic stability or recovery | AUROC | one-sided |
| `H2` | `hrv_nolineal` − `hrv_lineal` | same two | AUROC | two-sided (`comparaciones_h2_bilateral.csv`) |
| `H3` | `integrada_sofa` − `clinico_sofa` | SOFA deterioration | AUPRC | one-sided |
| `Complementario` (complementary) | `hrv_integrada` − `clinico_basal` | all three | the outcome's metric | one-sided |

### 2.5 SHAP groups and HRV summary types

| Value (`grupo` column or column header) | English |
|---|---|
| `Clínicas en T0` | Clinical variables at T0 |
| `PAM y lactato (12 h)` | MAP and lactate (previous 12 h) |
| `SOFA` | SOFA |
| `HRV temporal` | Time-domain HRV |
| `HRV frecuencial` | Frequency-domain HRV |
| `HRV no lineal` | Non-linear HRV |
| `HRV total` | All HRV groups combined |
| `Nivel (media)` | Level (mean) |
| `Nivel (último)` | Level (last value) |
| `Tendencia (pendiente)` | Trend (slope) |
| `Variabilidad (DE)` | Variability (SD) |

## 3. Columns of `ultimate_dataset_v3_1.csv` (3,100 windows × 97 columns)

One row per prediction window (`icustay_id`, T0). The dataset is **not** in this repository (MIMIC-III data use agreement); this table documents it so the code can be audited. Roles: *predictor* = used in at least one feature set; *outcome (never a predictor)* = listed in `desenlace_no_usar_como_feature` and blocked by an `assert` in the pipeline.

| Column | Group in `columnas_v3_1.json` | Role | Unit | Description |
|---|---|---|---|---|
| `subject_id` | `id` | identifier |  | MIMIC-III patient identifier |
| `icustay_id` | `id` | identifier / CV group |  | MIMIC-III ICU stay identifier; grouping unit for cross-validation and bootstrap ("patient") |
| `hour_index` | `id` | identifier | h | Hour of the prediction time T0, counted from the start of HRV analysis (12, 18, 24 ... <= 96) |
| `t0_prediccion` | `id` | identifier | datetime | Prediction time T0 = inicio_analisis_real + hour_index hours |
| `n_horas_validas` | — | inclusion | h | Valid HRV hours in the observation window [T0-12 h, T0-1 h] (0-12), from the window grid |
| `n_horas_hrv_ventana` | — | QC | h | Same count recomputed by the dataset builder; asserted equal to n_horas_validas |
| `incl_hrv_4h` | `filtros` | inclusion | 0/1 | >= 4 valid HRV hours in the window (sensitivity) |
| `incl_hrv_6h` | `filtros` | inclusion | 0/1 | >= 6 valid HRV hours in the window (main analysis, fixed a priori) |
| `incl_hrv_8h` | `filtros` | inclusion | 0/1 | >= 8 valid HRV hours in the window (sensitivity) |
| `eval_comun` | `filtros` | inclusion | 0/1 | 1 if both main labels are evaluable (common-cohort sensitivity analysis) |
| `age` | `features_clinicas` | predictor | years | Age at ICU admission (>89 years recoded to 91.4 as in MIMIC-III) |
| `gender` | `features_clinicas` | predictor | M/F | Sex; recoded to sexo_masculino (1 = male) in the modelling pipeline |
| `map_t0` | `features_clinicas` | predictor | mmHg | Mean arterial pressure at T0 (last value carried forward up to 1 h); < 20 or > 200 set to missing |
| `rr_t0` | `features_clinicas` | predictor | breaths/min | Respiratory rate at T0 (LOCF up to 6 h) |
| `spo2_t0` | `features_clinicas` | predictor | % | Peripheral oxygen saturation at T0 (LOCF up to 6 h) |
| `ned` | `features_clinicas` | predictor | µg/kg/min | Norepinephrine-equivalent dose at T0 = norepinephrine + epinephrine + dopamine/100 (maximum rate of active infusions) |
| `vis` | `features_clinicas` | predictor | score | Vasoactive-inotropic score at T0 = dopamine + dobutamine + 100 x epinephrine + 100 x norepinephrine |
| `cumvp` | `features_clinicas` | predictor | µg/kg/min (summed) | Cumulative vasopressor exposure = sum of NED over the 12 hourly points T0 ... T0-11 h |
| `sedante_bin` | `features_clinicas` | predictor | 0/1 | Any sedative infusion (propofol, midazolam, lorazepam, dexmedetomidine) overlapping (T0-12 h, T0] |
| `opioide_bin` | `features_clinicas` | predictor | 0/1 | Any opioid infusion (fentanyl, morphine, hydromorphone) overlapping (T0-12 h, T0] |
| `mbp_mean_actual12h` | `features_pam_lactato_actual` | predictor | mmHg | Mean of MAP measurements in (T0-12 h, T0]; < 20 or > 200 set to missing |
| `lactate_last_actual12h` | `features_pam_lactato_actual` | predictor | mmol/L | Last lactate in (T0-12 h, T0] |
| `lactato_medido_12h` | `features_pam_lactato_actual` | predictor | 0/1 | 1 if at least one lactate was measured in (T0-12 h, T0] (missingness indicator) |
| `sofa_t0_parcial` | `features_sofa` | predictor | points | Partial SOFA at T0 = sum of the components available at T0 |
| `n_comp_validos_t0` | `features_sofa` | predictor | count (0-6) | Number of SOFA components available at T0 |
| `horas_desde_ultima_hrv` | `features_hrv` | predictor | h | Hours between the last valid HRV hour of the window and T0 (0 = the hour just before T0) |
| `SDNN_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ms | SDNN, time domain (linear set). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `SDNN_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ms | SDNN, time domain (linear set). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SDNN_last` | `features_hrv`, `features_hrv_lineal` | predictor | ms | SDNN, time domain (linear set). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SDNN_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ms per hour | SDNN, time domain (linear set). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `RMSSD_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ms | RMSSD, time domain (linear set). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `RMSSD_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ms | RMSSD, time domain (linear set). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `RMSSD_last` | `features_hrv`, `features_hrv_lineal` | predictor | ms | RMSSD, time domain (linear set). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `RMSSD_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ms per hour | RMSSD, time domain (linear set). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `LF_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln LF power, frequency domain (linear set). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `LF_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln LF power, frequency domain (linear set). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `LF_last` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln LF power, frequency domain (linear set). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `LF_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) per hour | ln LF power, frequency domain (linear set). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `HF_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln HF power, frequency domain (linear set). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `HF_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln HF power, frequency domain (linear set). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `HF_last` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) | ln HF power, frequency domain (linear set). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `HF_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ln(NeuroKit2 units) per hour | ln HF power, frequency domain (linear set). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `LFHF_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ln(ratio) | ln LF/HF ratio, frequency domain (linear set). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `LFHF_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ln(ratio) | ln LF/HF ratio, frequency domain (linear set). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `LFHF_last` | `features_hrv`, `features_hrv_lineal` | predictor | ln(ratio) | ln LF/HF ratio, frequency domain (linear set). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `LFHF_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ln(ratio) per hour | ln LF/HF ratio, frequency domain (linear set). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `SD1_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD1, time domain (linear set; Brennan 2001). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD1_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD1, time domain (linear set; Brennan 2001). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD1_last` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD1, time domain (linear set; Brennan 2001). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD1_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ms per hour | Poincaré SD1, time domain (linear set; Brennan 2001). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `SD2_mean` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD2, time domain (linear set; Brennan 2001). Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD2_sd` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD2, time domain (linear set; Brennan 2001). Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD2_last` | `features_hrv`, `features_hrv_lineal` | predictor | ms | Poincaré SD2, time domain (linear set; Brennan 2001). Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SD2_slope` | `features_hrv`, `features_hrv_lineal` | predictor | ms per hour | Poincaré SD2, time domain (linear set; Brennan 2001). OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `SampEn_mean` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Sample entropy, non-linear set. Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `SampEn_sd` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Sample entropy, non-linear set. Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SampEn_last` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Sample entropy, non-linear set. Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `SampEn_slope` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless per hour | Sample entropy, non-linear set. OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `ApEn_mean` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Approximate entropy, non-linear set. Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `ApEn_sd` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Approximate entropy, non-linear set. Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `ApEn_last` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | Approximate entropy, non-linear set. Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `ApEn_slope` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless per hour | Approximate entropy, non-linear set. OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `DFA_a1_mean` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | DFA alpha1 (short-term scaling exponent), non-linear set. Mean of the hourly values in the 12-h observation window. Hourly value = median of the valid 5-min epochs of the hour |
| `DFA_a1_sd` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | DFA alpha1 (short-term scaling exponent), non-linear set. Standard deviation of the hourly values in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `DFA_a1_last` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless | DFA alpha1 (short-term scaling exponent), non-linear set. Last available hourly value in the window. Hourly value = median of the valid 5-min epochs of the hour |
| `DFA_a1_slope` | `features_hrv`, `features_hrv_nolineal` | predictor | dimensionless per hour | DFA alpha1 (short-term scaling exponent), non-linear set. OLS slope of the hourly values against hour (>= 3 points). Hourly value = median of the valid 5-min epochs of the hour |
| `label_sofa_responder` | `labels_principales` | label (main) | 1/0/missing | SOFA improvement: 1 if paired delta SOFA < 0; death in horizon = 0 |
| `label_cardio` | `labels_principales` | label (main) | 1/0/missing | Haemodynamic stability or recovery: 1 if the MAP criterion AND the lactate criterion are met (see README); death = 0 |
| `label_sofa_deterioro_2pts` | `labels_secundarios` | label (secondary) | 1/0/missing | SOFA deterioration: 1 if paired delta SOFA >= 2; death in horizon = 1 |
| `label_sofa_responder_cc` | `labels_sensibilidad` | label (sensitivity) | 1/0/missing | SOFA improvement, complete-case SOFA at both times |
| `label_sofa_deterioro_2pts_cc` | `labels_sensibilidad` | label (sensitivity) | 1/0/missing | SOFA deterioration, complete-case SOFA at both times |
| `label_sofa_responder_nl` | `labels_sensibilidad` | label (sensitivity) | 1/0/missing | SOFA improvement restricted to windows with a new laboratory value (or death) in the horizon |
| `label_sofa_deterioro_2pts_nl` | `labels_sensibilidad` | label (sensitivity) | 1/0/missing | SOFA deterioration restricted to windows with a new laboratory value (or death) in the horizon |
| `n_mbp_actual12h` | — | outcome definition | count | Number of MAP measurements in (T0-12 h, T0]; used only to decide whether label_cardio is evaluable |
| `sofa_t0_valid` | — | outcome definition | 0/1 | 1 if all six SOFA components are available at T0 |
| `sofa_t0_6h` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | points | Total SOFA at T0 + 6 h |
| `sofa_t0_6h_valid` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | 0/1 | 1 if all six SOFA components are available at T0 + 6 h |
| `delta_sofa_pareado` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | points | Paired delta SOFA (T0 + 6 h minus T0), component by component (pairing option C) |
| `delta_pareado_valid` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | 0/1 | 1 if the paired delta is evaluable (coagulation, cardiovascular, CNS and renal present at both times) |
| `mbp_mean_next6h` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | mmHg | Mean MAP in (T0, T0 + 6 h] |
| `n_mbp_next6h` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | count | Number of MAP measurements in (T0, T0 + 6 h] |
| `mbp_pct_next6h_gt65` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | fraction | Share of MAP measurements in (T0, T0 + 6 h] above 65 mmHg |
| `mbp_incremento_pct` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | fraction | Relative change of mean MAP: (mbp_mean_next6h - mbp_mean_actual12h) / mbp_mean_actual12h |
| `lactate_last_next6h` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | mmol/L | Last lactate in (T0, T0 + 6 h] |
| `muere_en_horizonte` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | 0/1 | Death in (T0, T0 + 6 h]; forces no improvement and deterioration |
| `alta_viva_en_horizonte` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | 0/1 | Discharged alive from the ICU before T0 + 6 h; window excluded (labels set to missing) |
| `outtime` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | datetime | ICU discharge time |
| `deathtime` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | datetime | In-hospital death time |
| `new_lab_6h` | `desenlace_no_usar_como_feature` | outcome (never a predictor) | 0/1 | At least one creatinine, bilirubin or platelet measurement in (T0, T0 + 6 h] |
| `temp_t0` | — | not used | °C | Temperature at T0 (LOCF up to 6 h); excluded from the models |
| `sofa_t0` | — | not used | points | Total SOFA at T0 (v3.3; component-specific LOCF windows); excluded as predictor, used for complete-case labels |
| `sofa_resp_t0` | — | not used | points | SOFA respiratory component at T0 |
| `sofa_coag_t0` | — | not used | points | SOFA coagulation component at T0 |
| `sofa_liver_t0` | — | not used | points | SOFA liver component at T0 |
| `sofa_cardio_t0` | — | not used | points | SOFA cardiovascular component at T0 |
| `sofa_neuro_t0` | — | not used | points | SOFA central nervous system component at T0 |
| `sofa_renal_t0` | — | not used | points | SOFA renal component at T0 (creatinine only) |

Notes:
- **HRV hourly value.** Each HRV hour has 12 five-minute epochs. The hourly value of a metric is the median of its valid epochs. LF, HF and LF/HF are log-transformed per epoch before aggregation.
- **LF and HF units.** These are NeuroKit2 0.2.13 `hrv_frequency` defaults. The project documentation describes them as normalised to the PSD maximum, so absolute LF and HF are not comparable with ms²; only the LF/HF ratio is interpretable in physiological terms.
- **SD1/SD2 grouping.** SD1 and SD2 belong to the linear (time-domain) set (Brennan et al., 2001). A legacy list in `extraccion_hrv_v3.py` puts them in the non-linear group. That list is only used to enumerate the metrics and does not affect the models.

## 4. Columns of the results files (`resultados/`)

All files describe the main analysis run (`modelado_v5_1`, dataset v3.1). `resultados/archivo_v5_sin_ln/` holds the same files from the superseded run without the log transform of LF, HF and LF/HF, kept for reference only.

Common columns: `analisis` (2.3), `hipotesis` (2.4), `desenlace` (1), `conjunto` (2.1), `familia` (2.2). AUROC and AUPRC are always computed on uncalibrated probabilities. Brier and all calibration measures use Platt-calibrated probabilities.

### `resumen_metricas.csv` — pooled out-of-fold performance (one row per configuration)
| Column | Meaning |
|---|---|
| `n` | Windows evaluated |
| `pacientes` | Patients (ICU stays) |
| `prevalencia` | Proportion of class 1 (favourable outcome for improvement/stability; deterioration for `progresion`) |
| `auroc`, `auprc` | Area under the ROC / precision-recall curve |
| `brier` | Brier score (calibrated probabilities) |

### `metricas_por_particion.csv` — performance per outer fold
| Column | Meaning |
|---|---|
| `particion` | Outer fold (0–4) |
| `n`, `eventos` | Windows in the test fold; windows with label = 1 |
| `auroc`, `auprc`, `brier` | As above, in that fold |

### `registro_configuraciones.csv` — configuration log
| Column | Meaning |
|---|---|
| `label` | Label column used |
| `n_variables` | Number of input variables |
| `mejores_hiperparametros` | Best hyperparameters, one dict per outer fold (JSON) |
| `segundos` | Run time in seconds |

### `comparaciones_bootstrap.csv` — patient-level bootstrap comparisons (2,000 replicates)
| Column | Meaning |
|---|---|
| `referencia`, `nuevo` | Reference and new feature set |
| `metrica` | `auroc`, `auprc` or `brier` |
| `principal` | True if this is the hypothesis's primary metric |
| `valor_referencia`, `valor_nuevo` | Metric of each model on the original data |
| `delta` | `valor_nuevo` − `valor_referencia` |
| `ic95_inf`, `ic95_sup` | 95 % percentile bootstrap CI of the delta |
| `p_unilateral` | One-sided p: share of replicates without improvement |
| `replicas_validas` | Replicates with both classes present |
| `ventanas`, `pacientes` | Windows and patients compared |
| `p_bh` | Benjamini–Hochberg-adjusted p across families (primary metric only) |
| `mejora_significativa` | **Significant improvement**: CI excludes 0 in the direction of improvement and `p_bh` < 0.05 |

### `comparaciones_h2_bilateral.csv` — H2 with two-sided p
Same columns as above, plus `p_bilateral` (= 2 · min(p, 1 − p)), `p_bh_bilateral` (BH within outcome) and `diferencia_significativa` (**significant difference** in either direction: CI excludes 0 and `p_bh_bilateral` < 0.05).

### `calibracion_metricas.csv` — calibration (`clinico_sofa` and `integrada_sofa`)
Columns ending in `_ic` hold the patient-level bootstrap 95 % CI as `low–high` text.
| Column | Meaning |
|---|---|
| `n`, `eventos` | Windows; windows with label = 1 |
| `auroc`, `auprc` | Discrimination |
| `intercepto` | Calibration intercept (calibration-in-the-large; logistic model with the logit as offset) |
| `pendiente` | Calibration slope |
| `E/O` | Expected / observed events |
| `brier`, `brier_escalado` | Brier score; scaled Brier = 1 − Brier / [prevalence × (1 − prevalence)] |
| `ICI` | Integrated calibration index (LOWESS, frac = 0.75) |

### `curvas_decision.csv` — decision curve analysis
| Column | Meaning |
|---|---|
| `umbral` | Threshold probability |
| `beneficio_neto` | Net benefit of the model |
| `intervenir_a_todos` | Net benefit of treating all |

Event direction: see section 1.

### `shap_importancia_variables.csv` — SHAP importance per variable (Elastic Net and XGBoost, log-odds)
| Column | Meaning |
|---|---|
| `variable` | Input variable (section 3) |
| `grupo` | Variable group (2.5) |
| `resumen_hrv` | HRV summary type (2.5); empty for non-HRV variables |
| `media_abs_shap` | Mean absolute out-of-fold SHAP value |
| `pct` | Share of the model's total mean absolute SHAP (%) |

### `shap_pct_por_grupo.csv`
Share (%) of total SHAP importance per group in the `integrada_sofa` model, plus `HRV total`.

### `shap_media_por_variable_grupo.csv`
Mean importance per variable within each group, rescaled so that the group with the highest mean in each model = 100.

### `shap_hrv_nivel_vs_tendencia.csv`
Within the HRV variables of `integrada_sofa`: share (%) of HRV importance by summary type (level, variability, trend).

### Descriptive tables
`resultados/eda_hrv_v3_1/` and `resultados/eda_clinico_v3/` keep the Spanish column headers produced by the notebooks. Their content is described in the markdown cells of the corresponding notebook (section codes A1–B8 and C1–C10 match the file names).
