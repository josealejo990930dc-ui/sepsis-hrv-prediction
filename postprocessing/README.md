# Post-processing code — HRV and 6-hour clinical improvement in sepsis

Code for the analyses reported in *"Incremental value of heart rate variability for predicting
6-hour clinical improvement in sepsis: development and internal validation of prediction models"*
(Domínguez Camacho JA, Cely A), beyond the main modelling pipeline: the post hoc analyses, the
supplementary tables (S1 Appendix), the main results tables, the figures.

The cohort definition (SQL), HRV extraction, dataset builder, modelling pipeline
(`pipeline_modelado_v5.py`), patient-level bootstrap (`bootstrap_v5.py`), two-sided H2 correction
(`h2_bilateral_v5_1.py`) and the calibration/decision-curve/SHAP notebook are in the companion
folders of this repository: [`../sql`](../sql), [`../python`](../python) and [`../notebooks`](../notebooks).

## Data access

No patient-level data are included. MIMIC-III (v1.4) and the MIMIC-III Waveform Database Matched
Subset (v1.0) are available to credentialed users through PhysioNet (https://physionet.org).
Every script reads files produced by the main pipeline in its output folder (called
`tesis_sepsis/` below):

| File | Produced by |
|---|---|
| `ultimate_dataset_v3_1.csv`, `columnas_v3_1.json` | `construir_ultimate_dataset_v3.py` |
| `modelado_v5_1/predicciones_oof.csv`, `resumen_metricas.csv`, `comparaciones_bootstrap.csv` | `pipeline_modelado_v5.py`, `bootstrap_v5.py` |
| `modelado_v5_1/comparaciones_h2_bilateral.csv` | `h2_bilateral_v5_1.py` |
| `modelado_v5_1/calibracion_shap/*.csv` | `Calibracion_SHAP_v5.ipynb` |
| `eda_clinico_v3/C2_faltantes.csv`, `C6_smd_clinica_por_desenlace.csv` | `EDA_clinico_v3.ipynb` |

## Installation

```bash
pip install -r requirements.txt
```

## What each script reproduces

| Script | Manuscript location | Main outputs |
|---|---|---|
| `post_hoc/01_alignment_check.py` | Results > Post hoc analyses > Alignment check | `alignment_check_auroc.csv`, `alignment_check_top.csv` |
| `post_hoc/02_positive_controls.py` | Results > Positive controls; S1 Table L | `positive_control_mortality.csv`, `positive_control_age.csv` |
| `post_hoc/03_fragmentation_descriptors.py` | Results > Fragmentation descriptors; Discussion | `fragmentation_by_age_tertile.csv`, `fragmentation_by_death.csv`, `fragmentation_correlation.csv` |
| `supplementary_tables/make_table_I_ci_absolute.py` | S1 Table I (95% CIs of all 90 configurations) | `ci_absolute_v5_1.csv` |
| `supplementary_tables/make_tables_H_to_M.py` | S1 Tables H–M (missing data, performance, all comparisons, SHAP, positive control, SMD) | `etables.tex` |
| `supplementary_tables/make_table_O_comparison_counts.py` | S1 Table O (83 comparisons; 249 with secondary metrics) | `table_O_comparison_counts.csv/.tex` |
| `main_tables/make_tables_4_to_6.py` | Tables 4–6 (discrimination, prespecified comparisons, calibration) | `tables/tab_*.tex` |
| `figures/make_fig1_flow.py`, `figures/make_figures.py`, `figures/to_tiff.py` | Fig 1, Fig 2, S1 Figs A–C; TIFF export for PLOS ONE | `*.png`, `Fig1.tif`, `Fig2.tif` |

## How to run (from this folder)

```bash
D=/path/to/tesis_sepsis

# Post hoc analyses (seconds)
python post_hoc/01_alignment_check.py --data-dir $D --out-dir results/post_hoc
python post_hoc/02_positive_controls.py --data-dir $D --out-dir results/post_hoc
python post_hoc/03_fragmentation_descriptors.py --data-dir $D --out-dir results/post_hoc

# S1 Appendix tables
python supplementary_tables/make_table_I_ci_absolute.py $D results/post_hoc/ci_absolute_v5_1.csv   # ~15 min
H2_CSV=$D/modelado_v5_1/comparaciones_h2_bilateral.csv \
  python supplementary_tables/make_tables_H_to_M.py $D results/post_hoc results/tables
python supplementary_tables/make_table_O_comparison_counts.py --results-dir $D/modelado_v5_1 --out results/tables

# Main tables 4-6 (written to ./tables)
python main_tables/make_tables_4_to_6.py $D/modelado_v5_1

# Figures
python figures/make_fig1_flow.py results/figures
python figures/make_figures.py $D results/figures
cd figures && python to_tiff.py ../results/figures ../results/figures && cd ..
```

## Expected values (as reported in the manuscript)

- **Alignment check** (oriented AUROC): haemodynamic stability or recovery, SD2_last 0.613 and
  SDNN_last 0.606; SOFA deterioration, SD2_mean 0.639 and DFA_a1_mean 0.607 (lower values predict
  deterioration); SOFA improvement, maximum 0.544 (HF_sd).
- **Positive controls** (first window per patient, n = 487, 69 deaths): DFA α1 0.693 (0.628–0.754),
  partial SOFA 0.679 (0.602–0.750), LF/HF 0.636, SD2 0.616, SDNN 0.587; inverted for RMSSD 0.463 and
  SampEn 0.415. Age: Spearman ρ RMSSD 0.29 (p = 4.9 × 10⁻¹¹), SDNN 0.14, SD2 0.09.
- **Fragmentation descriptors**: age tertile medians 46 / 63 / 79 years; RMSSD/SDNN 0.55 / 0.70 / 0.82;
  DFA α1 1.13 / 1.00 / 0.84; RMSSD/SDNN 0.95 in deaths vs 0.64 in survivors; RMSSD vs DFA α1 ρ −0.28.
- **Table O**: 35 primary + 48 sensitivity = 83 comparisons on the primary metric; 249 rows in total.

## Reproducibility note

The post-processing scripts are deterministic (fixed seed 2026 for all bootstraps) and reproduce
the manuscript values exactly from the pipeline outputs. Refitting the models is not bit-for-bit
reproducible across library versions: refitting the central configurations with scikit-learn 1.9.1
and XGBoost 3.4.1 reproduced the conclusions (no comparison met its criterion) but changed
individual AUROCs by up to 0.05 relative to the original Colab run (Python 3.12.13, NumPy 2.0.2;
other packages as listed in [`../requirements.txt`](../requirements.txt)). 
