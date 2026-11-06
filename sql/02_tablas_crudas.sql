-- ============================================================================
-- 02_tablas_crudas.sql: cleaned raw event tables (all ICU stays)
-- The longitudinal tables (03_queries_vigentes_v4.sql) apply LOCF to these events.
--   2.1 _raw_labs     creatinine, bilirubin, platelets
--   2.2 _raw_gcs      GCS per event
--   2.3 _raw_vitals   MAP, respiratory rate, SpO2, temperature
--   2.4 _raw_fio2     FiO2 normalised to %
--   2.5 _raw_lactate  lactate (itemid 50813)
-- Drug infusion tables are in 03_queries_vigentes_v4.sql, part A.
-- ============================================================================


-- 2.1 Laboratory values (same itemids and cleaning limits as labs_si)
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_labs` AS
SELECT
  ie.icustay_id,
  CAST(le.charttime AS DATETIME) AS charttime,
  CASE WHEN le.itemid = 50912 AND le.valuenum <= 150   THEN le.valuenum END AS creatinine,
  CASE WHEN le.itemid = 50885 AND le.valuenum <= 150   THEN le.valuenum END AS bilirubin,
  CASE WHEN le.itemid = 51265 AND le.valuenum <= 10000 THEN le.valuenum END AS platelet
FROM `physionet-data.mimiciii_clinical.icustays` ie
INNER JOIN `physionet-data.mimiciii_clinical.labevents` le ON ie.hadm_id = le.hadm_id
WHERE le.itemid IN (50912, 50885, 51265)
  AND le.valuenum IS NOT NULL AND le.valuenum > 0;


-- 2.2 GCS per event
--   The verbal score of an intubated patient ('1.0 ET/Trach' / 'No Response-ETT')
--   is coded as 15, so the sum is >= 15 and the SOFA CNS component is 0 (Johnson et al.
--   convention). The gcs column can exceed 15: do not use it as a continuous variable.
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_gcs` AS
WITH base AS (
  SELECT
    l.icustay_id,
    CAST(l.charttime AS DATETIME) AS charttime,
    MAX(CASE WHEN l.itemid IN (454, 223901) THEN l.valuenum END) AS gcs_motor,
    MAX(CASE
          WHEN l.itemid = 723    AND l.value = '1.0 ET/Trach'    THEN 15
          WHEN l.itemid = 223900 AND l.value = 'No Response-ETT' THEN 15
          WHEN l.itemid IN (723, 223900) THEN l.valuenum
        END) AS gcs_verbal,
    MAX(CASE WHEN l.itemid IN (184, 220739) THEN l.valuenum END) AS gcs_eyes
  FROM `physionet-data.mimiciii_clinical.chartevents` l
  WHERE l.itemid IN (184, 454, 723, 223900, 223901, 220739)
    AND COALESCE(l.error, 0) = 0
  GROUP BY l.icustay_id, CAST(l.charttime AS DATETIME)
)
SELECT icustay_id, charttime,
  COALESCE(gcs_motor, 6) + COALESCE(gcs_verbal, 5) + COALESCE(gcs_eyes, 4) AS gcs
FROM base;


-- 2.3 Vital signs (temperature °F -> °C)
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` AS
SELECT
  ce.icustay_id,
  CAST(ce.charttime AS DATETIME) AS charttime,
  CASE WHEN ce.itemid IN (456, 52, 6702, 443, 220052, 220181, 225312)
       AND ce.valuenum > 0 AND ce.valuenum < 300 THEN ce.valuenum END AS map_val,
  CASE WHEN ce.itemid IN (615, 618, 220210, 224690)
       AND ce.valuenum > 0 AND ce.valuenum < 70 THEN ce.valuenum END AS rr_val,
  CASE WHEN ce.itemid IN (646, 220277)
       AND ce.valuenum > 0 AND ce.valuenum <= 100 THEN ce.valuenum END AS spo2_val,
  CASE
    WHEN ce.itemid IN (676, 223762) AND ce.valuenum > 25 AND ce.valuenum < 45 THEN ce.valuenum
    WHEN ce.itemid IN (678, 223761) AND ce.valuenum > 70 AND ce.valuenum < 113 THEN (ce.valuenum - 32) / 1.8
  END AS temp_val
FROM `physionet-data.mimiciii_clinical.chartevents` ce
WHERE ce.itemid IN (456, 52, 6702, 443, 220052, 220181, 225312,   -- MAP
                    615, 618, 220210, 224690,   -- respiratory rate
                    646, 220277,   -- SpO2
                    676, 223762, 678, 223761)   -- temperature
  AND COALESCE(ce.error, 0) = 0;


-- 2.4 FiO2 (same normalisation as stg_fio2 in bloodgasarterial_si)
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_fio2` AS
SELECT
  ce.icustay_id,
  CAST(ce.charttime AS DATETIME) AS charttime,
  CASE
    WHEN ce.itemid = 223835 THEN CASE
      WHEN ce.valuenum > 0 AND ce.valuenum <= 1 THEN ce.valuenum * 100
      WHEN ce.valuenum > 1 AND ce.valuenum < 21 THEN NULL
      WHEN ce.valuenum >= 21 AND ce.valuenum <= 100 THEN ce.valuenum ELSE NULL END
    WHEN ce.itemid IN (3420, 3422) THEN ce.valuenum
    WHEN ce.itemid = 190 AND ce.valuenum > 0.20 AND ce.valuenum < 1 THEN ce.valuenum * 100
    ELSE NULL
  END AS fio2_val
FROM `physionet-data.mimiciii_clinical.chartevents` ce
WHERE ce.itemid IN (3420, 190, 223835, 3422)
  AND COALESCE(ce.error, 0) = 0;


-- 2.5 Lactate (haemodynamic outcome and lactate predictors)
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_lactate` AS
SELECT
  ie.icustay_id,
  CAST(le.charttime AS DATETIME) AS charttime,
  le.valuenum AS lactate
FROM `physionet-data.mimiciii_clinical.icustays` ie
INNER JOIN `physionet-data.mimiciii_clinical.labevents` le ON ie.hadm_id = le.hadm_id
WHERE le.itemid = 50813
  AND le.valuenum IS NOT NULL AND le.valuenum > 0 AND le.valuenum <= 50;   -- same cap as labs_si
