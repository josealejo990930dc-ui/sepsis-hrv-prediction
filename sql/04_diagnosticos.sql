-- ============================================================================
-- 04_diagnosticos.sql: read-only checks (no table is created or modified)
-- ============================================================================

-- 4.1 Drug itemids in d_items (led to the fix 222315 vasopressin -> 225150 dexmedetomidine)
SELECT itemid, label, dbsource
FROM `physionet-data.mimiciii_clinical.d_items`
WHERE itemid IN (222168, 221668, 221385, 222315, 225150,   -- sedatives (+ vasopressin)
                 221744, 225154, 221833,   -- opioids
                 221906, 221289, 221662, 221653, 221749);   -- vasopressors/inotropes

-- 4.2 Time range of the grid and of each longitudinal table
--   (detects tables built on an old grid)
SELECT 'grilla' AS tabla, COUNT(*) AS filas, COUNT(DISTINCT icustay_id) AS pacientes,
       MIN(hour_index) AS h_min, MAX(hour_index) AS h_max
FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
UNION ALL
SELECT 'sofa', COUNT(*), COUNT(DISTINCT icustay_id), MIN(hour_index), MAX(hour_index)
FROM `sepsis-mimic-def.sepsis3_workspace_def.sofa_longitudinal`
UNION ALL
SELECT 'drugs', COUNT(*), COUNT(DISTINCT icustay_id), MIN(hour_index), MAX(hour_index)
FROM `sepsis-mimic-def.sepsis3_workspace_def.drugs_longitudinal`
UNION ALL
SELECT 'vitals', COUNT(*), COUNT(DISTINCT icustay_id), MIN(hour_index), MAX(hour_index)
FROM `sepsis-mimic-def.sepsis3_workspace_def.vitals_longitudinal`
UNION ALL
SELECT 'cardio_lactato', COUNT(*), COUNT(DISTINCT icustay_id), MIN(hour_index), MAX(hour_index)
FROM `sepsis-mimic-def.sepsis3_workspace_def.cardio_lactato_longitudinal`;

-- 4.3 Columns of a table (detects old versions with missing columns)
SELECT column_name, data_type
FROM `sepsis-mimic-def.sepsis3_workspace_def.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name = 'sofa_longitudinal'
ORDER BY ordinal_position;

-- 4.4 Duration of the infusion segments (check of part A)
SELECT 'vaso' AS tabla, COUNT(*) AS tramos,
       COUNTIF(endtime <= starttime) AS sin_duracion,
       APPROX_QUANTILES(DATETIME_DIFF(endtime, starttime, MINUTE), 4) AS cuartiles_min
FROM `sepsis-mimic-def.sepsis3_workspace_def._raw_vaso`
UNION ALL
SELECT 'sedantes', COUNT(*), COUNTIF(endtime <= starttime),
       APPROX_QUANTILES(DATETIME_DIFF(endtime, starttime, MINUTE), 4)
FROM `sepsis-mimic-def.sepsis3_workspace_def._raw_sedantes`
UNION ALL
SELECT 'opioides', COUNT(*), COUNTIF(endtime <= starttime),
       APPROX_QUANTILES(DATETIME_DIFF(endtime, starttime, MINUTE), 4)
FROM `sepsis-mimic-def.sepsis3_workspace_def._raw_opioides`;
