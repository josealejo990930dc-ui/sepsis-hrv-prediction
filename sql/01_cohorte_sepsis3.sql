-- ============================================================================
-- 01_cohorte_sepsis3.sql: Sepsis-3 cohort and HRV time anchor
-- BigQuery reimplementation of Johnson et al. (2018), "A Comparative Analysis of
-- Sepsis Identification Methods in an Electronic Database".
-- Source: physionet-data.mimiciii_clinical (MIMIC-III v1.4).
--
-- Run in order (each table depends on the previous ones):
--   1.1 suspinfect_poe        suspected infection (antibiotic + culture)
--   1.2 labs_si               laboratory values in the window [-48 h, +24 h]
--   1.3 gcs_si                worst GCS in the window
--   1.4 vitals_si             minimum MAP in the window (reconstructed, see 1.4)
--   1.5 bloodgasarterial_si   arterial blood gases and PaO2/FiO2
--   1.6 sofa_si               SOFA in the suspicion window
--   1.7 sepsis3_final_cohort  final cohort (inclusion/exclusion criteria)
--   1.8 hrv_analysis_mapping  time anchor for the HRV extraction
-- ============================================================================


-- ============================================================================
-- 1.1 suspinfect_poe: first time of suspected infection
--   Culture -> antibiotic within 72 h, or antibiotic -> culture within 24 h;
--   suspicion time = the earlier of the two events.
--   Clinical window: [suspicion - 48 h, suspicion + 24 h].
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` AS
WITH
abx_poe_list AS (
  SELECT DISTINCT drug
  FROM `physionet-data.mimiciii_clinical.prescriptions`
  WHERE drug_type IN ('MAIN','ADDITIVE')
    AND route NOT IN ('OU','OS','OD','AU','AS','AD','TP')
    AND NOT REGEXP_CONTAINS(LOWER(route), r'ear|eye')
    AND NOT REGEXP_CONTAINS(LOWER(drug), r'cream|desensitization|ophth oint|gel')
    AND REGEXP_CONTAINS(LOWER(drug), r'amikacin|amoxicillin|ampicillin|augmentin|azithromycin|aztreonam|cefazolin|cefoxitin|ceftazidime|cefepime|cefotaxime|cefuroxime|cephalexin|ciprofloxacin|clindamycin|doxycycline|erythromycin|gentamicin|levofloxacin|metronidazole|nafcillin|nitrofurantoin|oxacillin|penicillin|piperacillin|tazobactam|rifampin|rocephin|smz-tmp|sulfamethoxazole|tobramycin|vancomycin|zosyn|zyvox')
),
abx AS (
  SELECT pr.hadm_id, pr.drug AS antibiotic_name,
         CAST(pr.startdate AS DATETIME) AS antibiotic_time,
         CAST(pr.enddate AS DATETIME)   AS antibiotic_endtime
  FROM `physionet-data.mimiciii_clinical.prescriptions` pr
  INNER JOIN abx_poe_list ab ON pr.drug = ab.drug
),
ab_tbl AS (
  SELECT ie.subject_id, ie.hadm_id, ie.icustay_id, ie.intime, ie.outtime,
         abx.antibiotic_name, abx.antibiotic_time, abx.antibiotic_endtime
  FROM `physionet-data.mimiciii_clinical.icustays` ie
  LEFT JOIN abx ON ie.hadm_id = abx.hadm_id
),
me AS (
  SELECT hadm_id,
         CAST(chartdate AS DATETIME) AS chartdate,
         CAST(charttime AS DATETIME) AS charttime,
         spec_type_desc,
         MAX(CASE WHEN org_name IS NOT NULL AND org_name != '' THEN 1 ELSE 0 END) AS PositiveCulture
  FROM `physionet-data.mimiciii_clinical.microbiologyevents`
  GROUP BY hadm_id, chartdate, charttime, spec_type_desc
),
ab_fnl AS (
  SELECT ab_tbl.icustay_id, ab_tbl.intime, ab_tbl.outtime, ab_tbl.antibiotic_name, ab_tbl.antibiotic_time,
         COALESCE(me72.charttime, me72.chartdate) AS last72_charttime,
         COALESCE(me24.charttime, me24.chartdate) AS next24_charttime,
         me72.PositiveCulture AS last72_positiveculture, me72.spec_type_desc AS last72_specimen,
         me24.PositiveCulture AS next24_positiveculture, me24.spec_type_desc AS next24_specimen
  FROM ab_tbl
  -- Scenario A: culture up to 72 h BEFORE the antibiotic (96 h if only the date is known)
  LEFT JOIN me me72
    ON ab_tbl.hadm_id = me72.hadm_id
    AND ab_tbl.antibiotic_time IS NOT NULL
    AND (
      (ab_tbl.antibiotic_time > me72.charttime AND ab_tbl.antibiotic_time <= DATETIME_ADD(me72.charttime, INTERVAL 72 HOUR))
      OR
      (me72.charttime IS NULL AND ab_tbl.antibiotic_time > me72.chartdate AND ab_tbl.antibiotic_time < DATETIME_ADD(me72.chartdate, INTERVAL 96 HOUR))
    )
  -- Scenario B: culture up to 24 h AFTER the antibiotic
  LEFT JOIN me me24
    ON ab_tbl.hadm_id = me24.hadm_id
    AND ab_tbl.antibiotic_time IS NOT NULL
    AND me24.charttime IS NOT NULL
    AND (
      (ab_tbl.antibiotic_time > DATETIME_SUB(me24.charttime, INTERVAL 24 HOUR) AND ab_tbl.antibiotic_time <= me24.charttime)
      OR
      (me24.charttime IS NULL AND ab_tbl.antibiotic_time > me24.chartdate AND ab_tbl.antibiotic_time <= DATETIME_ADD(me24.chartdate, INTERVAL 24 HOUR))
    )
),
ab_laststg AS (
  SELECT icustay_id, antibiotic_name, antibiotic_time, last72_charttime, next24_charttime,
         COALESCE(last72_charttime, antibiotic_time) AS suspected_infection_time,
         CASE WHEN last72_charttime IS NOT NULL THEN last72_specimen
              WHEN next24_charttime IS NOT NULL THEN next24_specimen END AS specimen,
         CASE WHEN last72_charttime IS NOT NULL THEN last72_positiveculture
              WHEN next24_charttime IS NOT NULL THEN next24_positiveculture END AS positiveculture
  FROM ab_fnl
),
abx_ranked AS (
  SELECT icustay_id, suspected_infection_time, specimen, positiveculture, antibiotic_name, antibiotic_time,
         ROW_NUMBER() OVER (PARTITION BY icustay_id ORDER BY suspected_infection_time) AS rn
  FROM ab_laststg
  WHERE suspected_infection_time IS NOT NULL
)
SELECT ie.icustay_id, abx.antibiotic_name, abx.antibiotic_time, abx.suspected_infection_time,
       DATETIME_SUB(abx.suspected_infection_time, INTERVAL 48 HOUR) AS si_starttime,
       DATETIME_ADD(abx.suspected_infection_time, INTERVAL 24 HOUR) AS si_endtime,
       abx.specimen, abx.positiveculture
FROM `physionet-data.mimiciii_clinical.icustays` ie
LEFT JOIN abx_ranked abx ON ie.icustay_id = abx.icustay_id AND abx.rn = 1
WHERE abx.suspected_infection_time IS NOT NULL;


-- ============================================================================
-- 1.2 labs_si: minimum and maximum of each laboratory value in the suspicion window
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.labs_si` AS
WITH cohort AS (
  SELECT s.icustay_id, s.si_starttime, s.si_endtime, ie.hadm_id
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
  INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON s.icustay_id = ie.icustay_id
),
pvt AS (
  SELECT ie.icustay_id,
    CASE
      WHEN itemid = 50868 THEN 'ANION GAP'   WHEN itemid = 50862 THEN 'ALBUMIN'
      WHEN itemid = 50882 THEN 'BICARBONATE' WHEN itemid = 50885 THEN 'BILIRUBIN'
      WHEN itemid = 50912 THEN 'CREATININE'  WHEN itemid IN (50806, 50902) THEN 'CHLORIDE'
      WHEN itemid IN (50809, 50931) THEN 'GLUCOSE'    WHEN itemid IN (50810, 51221) THEN 'HEMATOCRIT'
      WHEN itemid IN (50811, 51222) THEN 'HEMOGLOBIN' WHEN itemid = 50813 THEN 'LACTATE'
      WHEN itemid = 51265 THEN 'PLATELET'    WHEN itemid IN (50822, 50971) THEN 'POTASSIUM'
      WHEN itemid = 51275 THEN 'PTT'         WHEN itemid = 51237 THEN 'INR'
      WHEN itemid = 51274 THEN 'PT'          WHEN itemid IN (50824, 50983) THEN 'SODIUM'
      WHEN itemid = 51006 THEN 'BUN'         WHEN itemid IN (51300, 51301) THEN 'WBC'
      WHEN itemid = 51144 THEN 'BANDS'
      ELSE NULL
    END AS label,
    CASE
      WHEN itemid = 50862 AND valuenum > 10 THEN NULL                 -- ALBUMIN
      WHEN itemid = 50868 AND valuenum > 10000 THEN NULL              -- ANION GAP
      WHEN itemid = 50882 AND valuenum > 10000 THEN NULL              -- BICARBONATE
      WHEN itemid = 50885 AND valuenum > 150 THEN NULL                -- BILIRUBIN
      WHEN itemid IN (50806, 50902) AND valuenum > 10000 THEN NULL    -- CHLORIDE
      WHEN itemid = 50912 AND valuenum > 150 THEN NULL                -- CREATININE
      WHEN itemid IN (50809, 50931) AND valuenum > 10000 THEN NULL    -- GLUCOSE
      WHEN itemid IN (50810, 51221) AND valuenum > 100 THEN NULL      -- HEMATOCRIT
      WHEN itemid IN (50811, 51222) AND valuenum > 50 THEN NULL       -- HEMOGLOBIN
      WHEN itemid = 50813 AND valuenum > 50 THEN NULL                 -- LACTATE
      WHEN itemid = 51265 AND valuenum > 10000 THEN NULL              -- PLATELET
      WHEN itemid IN (50822, 50971) AND valuenum > 30 THEN NULL       -- POTASSIUM
      WHEN itemid = 51275 AND valuenum > 150 THEN NULL                -- PTT
      WHEN itemid = 51237 AND valuenum > 50 THEN NULL                 -- INR
      WHEN itemid = 51274 AND valuenum > 150 THEN NULL                -- PT
      WHEN itemid IN (50824, 50983) AND valuenum > 200 THEN NULL      -- SODIUM
      WHEN itemid = 51006 AND valuenum > 300 THEN NULL                -- BUN
      WHEN itemid IN (51300, 51301) AND valuenum > 1000 THEN NULL     -- WBC
      WHEN itemid = 51144 AND (valuenum < 0 OR valuenum > 100) THEN NULL -- BANDS
      ELSE le.valuenum
    END AS valuenum
  FROM cohort ie
  LEFT JOIN `physionet-data.mimiciii_clinical.labevents` le
    ON le.hadm_id = ie.hadm_id
    AND CAST(le.charttime AS DATETIME) BETWEEN ie.si_starttime AND ie.si_endtime
    AND le.itemid IN (50868, 50862, 51144, 50882, 50885, 50912, 50902, 50806, 50931, 50809, 51221, 50810,
                      51222, 50811, 50813, 51265, 50971, 50822, 51275, 51237, 51274, 50983, 50824, 51006, 51301, 51300)
    AND le.valuenum IS NOT NULL AND le.valuenum > 0
)
SELECT icustay_id,
  MIN(CASE WHEN label = 'ANION GAP'   THEN valuenum END) AS ANIONGAP_min,
  MAX(CASE WHEN label = 'ANION GAP'   THEN valuenum END) AS ANIONGAP_max,
  MIN(CASE WHEN label = 'ALBUMIN'     THEN valuenum END) AS ALBUMIN_min,
  MAX(CASE WHEN label = 'ALBUMIN'     THEN valuenum END) AS ALBUMIN_max,
  MIN(CASE WHEN label = 'BANDS'       THEN valuenum END) AS BANDS_min,
  MAX(CASE WHEN label = 'BANDS'       THEN valuenum END) AS BANDS_max,
  MIN(CASE WHEN label = 'BICARBONATE' THEN valuenum END) AS BICARBONATE_min,
  MAX(CASE WHEN label = 'BICARBONATE' THEN valuenum END) AS BICARBONATE_max,
  MIN(CASE WHEN label = 'BILIRUBIN'   THEN valuenum END) AS BILIRUBIN_min,
  MAX(CASE WHEN label = 'BILIRUBIN'   THEN valuenum END) AS BILIRUBIN_max,
  MIN(CASE WHEN label = 'CREATININE'  THEN valuenum END) AS CREATININE_min,
  MAX(CASE WHEN label = 'CREATININE'  THEN valuenum END) AS CREATININE_max,
  MIN(CASE WHEN label = 'CHLORIDE'    THEN valuenum END) AS CHLORIDE_min,
  MAX(CASE WHEN label = 'CHLORIDE'    THEN valuenum END) AS CHLORIDE_max,
  MIN(CASE WHEN label = 'GLUCOSE'     THEN valuenum END) AS GLUCOSE_min,
  MAX(CASE WHEN label = 'GLUCOSE'     THEN valuenum END) AS GLUCOSE_max,
  MIN(CASE WHEN label = 'HEMATOCRIT'  THEN valuenum END) AS HEMATOCRIT_min,
  MAX(CASE WHEN label = 'HEMATOCRIT'  THEN valuenum END) AS HEMATOCRIT_max,
  MIN(CASE WHEN label = 'HEMOGLOBIN'  THEN valuenum END) AS HEMOGLOBIN_min,
  MAX(CASE WHEN label = 'HEMOGLOBIN'  THEN valuenum END) AS HEMOGLOBIN_max,
  MIN(CASE WHEN label = 'LACTATE'     THEN valuenum END) AS LACTATE_min,
  MAX(CASE WHEN label = 'LACTATE'     THEN valuenum END) AS LACTATE_max,
  MIN(CASE WHEN label = 'PLATELET'    THEN valuenum END) AS PLATELET_min,
  MAX(CASE WHEN label = 'PLATELET'    THEN valuenum END) AS PLATELET_max,
  MIN(CASE WHEN label = 'POTASSIUM'   THEN valuenum END) AS POTASSIUM_min,
  MAX(CASE WHEN label = 'POTASSIUM'   THEN valuenum END) AS POTASSIUM_max,
  MIN(CASE WHEN label = 'PTT'         THEN valuenum END) AS PTT_min,
  MAX(CASE WHEN label = 'PTT'         THEN valuenum END) AS PTT_max,
  MIN(CASE WHEN label = 'INR'         THEN valuenum END) AS INR_min,
  MAX(CASE WHEN label = 'INR'         THEN valuenum END) AS INR_max,
  MIN(CASE WHEN label = 'PT'          THEN valuenum END) AS PT_min,
  MAX(CASE WHEN label = 'PT'          THEN valuenum END) AS PT_max,
  MIN(CASE WHEN label = 'SODIUM'      THEN valuenum END) AS SODIUM_min,
  MAX(CASE WHEN label = 'SODIUM'      THEN valuenum END) AS SODIUM_max,
  MIN(CASE WHEN label = 'BUN'         THEN valuenum END) AS BUN_min,
  MAX(CASE WHEN label = 'BUN'         THEN valuenum END) AS BUN_max,
  MIN(CASE WHEN label = 'WBC'         THEN valuenum END) AS WBC_min,
  MAX(CASE WHEN label = 'WBC'         THEN valuenum END) AS WBC_max
FROM pvt
GROUP BY icustay_id;


-- ============================================================================
-- 1.3 gcs_si: worst GCS in the window (intubated = 15; previous value carried up to 6 h)
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.gcs_si` AS
WITH base AS (
  SELECT pvt.icustay_id, pvt.charttime,
     MAX(CASE WHEN pvt.itemid = 454 THEN pvt.valuenum END) AS GCSMotor,
     MAX(CASE WHEN pvt.itemid = 723 THEN pvt.valuenum END) AS GCSVerbal,
     MAX(CASE WHEN pvt.itemid = 184 THEN pvt.valuenum END) AS GCSEyes,
     CASE WHEN MAX(CASE WHEN pvt.itemid = 723 THEN pvt.valuenum END) = 0 THEN 1 ELSE 0 END AS EndoTrachFlag,
     ROW_NUMBER() OVER (PARTITION BY pvt.icustay_id ORDER BY pvt.charttime ASC) AS rn
  FROM (
     SELECT l.icustay_id,
        CASE WHEN l.itemid IN (723, 223900) THEN 723
             WHEN l.itemid IN (454, 223901) THEN 454
             WHEN l.itemid IN (184, 220739) THEN 184
             ELSE l.itemid END AS itemid,
        CASE WHEN l.itemid = 723    AND l.value = '1.0 ET/Trach'    THEN 0
             WHEN l.itemid = 223900 AND l.value = 'No Response-ETT' THEN 0
             ELSE l.valuenum END AS valuenum,
        CAST(l.charttime AS DATETIME) AS charttime
     FROM `physionet-data.mimiciii_clinical.chartevents` l
     INNER JOIN `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` ie ON l.icustay_id = ie.icustay_id
     WHERE l.itemid IN (184, 454, 723, 223900, 223901, 220739)
       AND COALESCE(l.error, 0) = 0
       AND CAST(l.charttime AS DATETIME) BETWEEN ie.si_starttime AND ie.si_endtime
  ) pvt
  GROUP BY pvt.icustay_id, pvt.charttime
),
gcs AS (
  SELECT b.*, b2.GCSVerbal AS GCSVerbalPrev, b2.GCSMotor AS GCSMotorPrev, b2.GCSEyes AS GCSEyesPrev,
     CASE
       WHEN b.GCSVerbal = 0 THEN 15
       WHEN b.GCSVerbal IS NULL AND b2.GCSVerbal = 0 THEN 15
       WHEN b2.GCSVerbal = 0 THEN COALESCE(b.GCSMotor, 6) + COALESCE(b.GCSVerbal, 5) + COALESCE(b.GCSEyes, 4)
       ELSE COALESCE(b.GCSMotor, COALESCE(b2.GCSMotor, 6))
          + COALESCE(b.GCSVerbal, COALESCE(b2.GCSVerbal, 5))
          + COALESCE(b.GCSEyes, COALESCE(b2.GCSEyes, 4))
     END AS GCS
  FROM base b
  LEFT JOIN base b2
    ON b.icustay_id = b2.icustay_id AND b.rn = b2.rn + 1
    AND b2.charttime > DATETIME_SUB(b.charttime, INTERVAL 6 HOUR)
),
gcs_final AS (
  SELECT gcs.*, ROW_NUMBER() OVER (PARTITION BY gcs.icustay_id ORDER BY gcs.GCS ASC) AS IsMinGCS
  FROM gcs
)
SELECT ie.icustay_id, gs.GCS AS MinGCS,
       COALESCE(gs.GCSMotor, gs.GCSMotorPrev)   AS GCSMotor,
       COALESCE(gs.GCSVerbal, gs.GCSVerbalPrev) AS GCSVerbal,
       COALESCE(gs.GCSEyes, gs.GCSEyesPrev)     AS GCSEyes,
       gs.EndoTrachFlag
FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` ie
LEFT JOIN gcs_final gs ON ie.icustay_id = gs.icustay_id AND gs.IsMinGCS = 1;


-- ============================================================================
-- 1.4 vitals_si: minimum MAP in the suspicion window
--   RECONSTRUCTED: the original query was not preserved. Same MAP itemids and
--   limits as _raw_vitals. Check against the existing table before rerunning.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.vitals_si` AS
SELECT s.icustay_id,
       MIN(ce.valuenum) AS MeanBP_Min,
       MAX(ce.valuenum) AS MeanBP_Max
FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
LEFT JOIN `physionet-data.mimiciii_clinical.chartevents` ce
  ON s.icustay_id = ce.icustay_id
  AND CAST(ce.charttime AS DATETIME) BETWEEN s.si_starttime AND s.si_endtime
  AND ce.itemid IN (456, 52, 6702, 443, 220052, 220181, 225312)
  AND ce.valuenum > 0 AND ce.valuenum < 300
  AND COALESCE(ce.error, 0) = 0
GROUP BY s.icustay_id;


-- ============================================================================
-- 1.5 bloodgasarterial_si: arterial blood gases, closest SpO2/FiO2 and PaO2/FiO2
--   Samples without specimen type are classified as arterial with the mimic-code
--   logistic model when the predicted probability is > 0.75.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.bloodgasarterial_si` AS
WITH cohort AS (
  SELECT s.icustay_id, s.si_starttime, s.si_endtime, ie.hadm_id
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
  INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON s.icustay_id = ie.icustay_id
),
bg AS (
  SELECT pvt.icustay_id, pvt.charttime,
    MAX(CASE WHEN label = 'SPECIMEN' THEN value END) AS SPECIMEN,
    MAX(CASE WHEN label = 'AADO2' THEN valuenum END) AS AADO2,
    MAX(CASE WHEN label = 'BASEEXCESS' THEN valuenum END) AS BASEEXCESS,
    MAX(CASE WHEN label = 'BICARBONATE' THEN valuenum END) AS BICARBONATE,
    MAX(CASE WHEN label = 'TOTALCO2' THEN valuenum END) AS TOTALCO2,
    MAX(CASE WHEN label = 'CARBOXYHEMOGLOBIN' THEN valuenum END) AS CARBOXYHEMOGLOBIN,
    MAX(CASE WHEN label = 'CHLORIDE' THEN valuenum END) AS CHLORIDE,
    MAX(CASE WHEN label = 'CALCIUM' THEN valuenum END) AS CALCIUM,
    MAX(CASE WHEN label = 'GLUCOSE' THEN valuenum END) AS GLUCOSE,
    MAX(CASE WHEN label = 'HEMATOCRIT' THEN valuenum END) AS HEMATOCRIT,
    MAX(CASE WHEN label = 'HEMOGLOBIN' THEN valuenum END) AS HEMOGLOBIN,
    MAX(CASE WHEN label = 'INTUBATED' THEN valuenum END) AS INTUBATED,
    MAX(CASE WHEN label = 'LACTATE' THEN valuenum END) AS LACTATE,
    MAX(CASE WHEN label = 'METHEMOGLOBIN' THEN valuenum END) AS METHEMOGLOBIN,
    MAX(CASE WHEN label = 'O2FLOW' THEN valuenum END) AS O2FLOW,
    MAX(CASE WHEN label = 'FIO2' THEN valuenum END) AS FIO2,
    MAX(CASE WHEN label = 'SO2' THEN valuenum END) AS SO2,
    MAX(CASE WHEN label = 'PCO2' THEN valuenum END) AS PCO2,
    MAX(CASE WHEN label = 'PEEP' THEN valuenum END) AS PEEP,
    MAX(CASE WHEN label = 'PH' THEN valuenum END) AS PH,
    MAX(CASE WHEN label = 'PO2' THEN valuenum END) AS PO2,
    MAX(CASE WHEN label = 'POTASSIUM' THEN valuenum END) AS POTASSIUM,
    MAX(CASE WHEN label = 'REQUIREDO2' THEN valuenum END) AS REQUIREDO2,
    MAX(CASE WHEN label = 'SODIUM' THEN valuenum END) AS SODIUM,
    MAX(CASE WHEN label = 'TEMPERATURE' THEN valuenum END) AS TEMPERATURE,
    MAX(CASE WHEN label = 'TIDALVOLUME' THEN valuenum END) AS TIDALVOLUME,
    MAX(CASE WHEN label = 'VENTILATIONRATE' THEN valuenum END) AS VENTILATIONRATE,
    MAX(CASE WHEN label = 'VENTILATOR' THEN valuenum END) AS VENTILATOR
  FROM (
    SELECT ie.icustay_id, CAST(le.charttime AS DATETIME) AS charttime, le.value,
      CASE
        WHEN itemid = 50800 THEN 'SPECIMEN'   WHEN itemid = 50801 THEN 'AADO2'
        WHEN itemid = 50802 THEN 'BASEEXCESS' WHEN itemid = 50803 THEN 'BICARBONATE'
        WHEN itemid = 50804 THEN 'TOTALCO2'   WHEN itemid = 50805 THEN 'CARBOXYHEMOGLOBIN'
        WHEN itemid = 50806 THEN 'CHLORIDE'   WHEN itemid = 50808 THEN 'CALCIUM'
        WHEN itemid = 50809 THEN 'GLUCOSE'    WHEN itemid = 50810 THEN 'HEMATOCRIT'
        WHEN itemid = 50811 THEN 'HEMOGLOBIN' WHEN itemid = 50812 THEN 'INTUBATED'
        WHEN itemid = 50813 THEN 'LACTATE'    WHEN itemid = 50814 THEN 'METHEMOGLOBIN'
        WHEN itemid = 50815 THEN 'O2FLOW'     WHEN itemid = 50816 THEN 'FIO2'
        WHEN itemid = 50817 THEN 'SO2'        WHEN itemid = 50818 THEN 'PCO2'
        WHEN itemid = 50819 THEN 'PEEP'       WHEN itemid = 50820 THEN 'PH'
        WHEN itemid = 50821 THEN 'PO2'        WHEN itemid = 50822 THEN 'POTASSIUM'
        WHEN itemid = 50823 THEN 'REQUIREDO2' WHEN itemid = 50824 THEN 'SODIUM'
        WHEN itemid = 50825 THEN 'TEMPERATURE' WHEN itemid = 50826 THEN 'TIDALVOLUME'
        WHEN itemid = 50827 THEN 'VENTILATIONRATE' WHEN itemid = 50828 THEN 'VENTILATOR'
        ELSE NULL END AS label,
      CASE
        WHEN valuenum <= 0 THEN NULL
        WHEN itemid = 50810 AND valuenum > 100 THEN NULL
        WHEN itemid = 50816 AND valuenum > 100 THEN NULL
        WHEN itemid = 50817 AND valuenum > 100 THEN NULL
        WHEN itemid = 50815 AND valuenum > 70 THEN NULL
        WHEN itemid = 50821 AND valuenum > 800 THEN NULL
        ELSE valuenum END AS valuenum
    FROM cohort ie
    LEFT JOIN `physionet-data.mimiciii_clinical.labevents` le
      ON le.hadm_id = ie.hadm_id
      AND CAST(le.charttime AS DATETIME) BETWEEN ie.si_starttime AND ie.si_endtime
      AND le.itemid IN (50800, 50801, 50802, 50803, 50804, 50805, 50806, 50807, 50808, 50809, 50810, 50811, 50812,
                        50813, 50814, 50815, 50816, 50817, 50818, 50819, 50820, 50821, 50822, 50823, 50824, 50825,
                        50826, 50827, 50828, 51545)
  ) pvt
  GROUP BY pvt.icustay_id, pvt.charttime
),
stg_spo2 AS (
  SELECT ie.icustay_id, CAST(ce.charttime AS DATETIME) AS charttime,
    MAX(CASE WHEN ce.valuenum <= 0 OR ce.valuenum > 100 THEN NULL ELSE ce.valuenum END) AS spo2
  FROM cohort ie
  INNER JOIN `physionet-data.mimiciii_clinical.chartevents` ce
    ON ie.icustay_id = ce.icustay_id
    AND CAST(ce.charttime AS DATETIME) BETWEEN DATETIME_SUB(ie.si_starttime, INTERVAL 2 HOUR) AND ie.si_endtime
  WHERE ce.itemid IN (646, 220277) AND COALESCE(ce.error, 0) = 0
  GROUP BY ie.icustay_id, CAST(ce.charttime AS DATETIME)
),
stg_fio2 AS (
  SELECT ie.icustay_id, CAST(ce.charttime AS DATETIME) AS charttime,
    MAX(CASE
      WHEN ce.itemid = 223835 THEN CASE
        WHEN ce.valuenum > 0 AND ce.valuenum <= 1 THEN ce.valuenum * 100
        WHEN ce.valuenum > 1 AND ce.valuenum < 21 THEN NULL
        WHEN ce.valuenum >= 21 AND ce.valuenum <= 100 THEN ce.valuenum ELSE NULL END
      WHEN ce.itemid IN (3420, 3422) THEN ce.valuenum
      WHEN ce.itemid = 190 AND ce.valuenum > 0.20 AND ce.valuenum < 1 THEN ce.valuenum * 100
      ELSE NULL END) AS fio2_chartevents
  FROM cohort ie
  INNER JOIN `physionet-data.mimiciii_clinical.chartevents` ce
    ON ie.icustay_id = ce.icustay_id
    AND CAST(ce.charttime AS DATETIME) BETWEEN DATETIME_SUB(ie.si_starttime, INTERVAL 4 HOUR) AND ie.si_endtime
  WHERE ce.itemid IN (3420, 190, 223835, 3422) AND COALESCE(ce.error, 0) = 0
  GROUP BY ie.icustay_id, CAST(ce.charttime AS DATETIME)
),
stg2 AS (
  SELECT bg.*, s1.spo2,
    ROW_NUMBER() OVER (PARTITION BY bg.icustay_id, bg.charttime ORDER BY s1.charttime DESC) AS lastRowSpO2
  FROM bg
  LEFT JOIN stg_spo2 s1
    ON bg.icustay_id = s1.icustay_id
    AND s1.charttime BETWEEN DATETIME_SUB(bg.charttime, INTERVAL 2 HOUR) AND bg.charttime
  WHERE bg.po2 IS NOT NULL
),
stg3 AS (
  SELECT bg.*, s2.fio2_chartevents,
    ROW_NUMBER() OVER (PARTITION BY bg.icustay_id, bg.charttime ORDER BY s2.charttime DESC) AS lastRowFiO2,
    1 / (1 + EXP(-(-0.02544 + 0.04598 * po2
      + COALESCE(-0.15356 * spo2, -0.15356 * 97.49420 + 0.13429)
      + COALESCE(0.00621 * fio2_chartevents, 0.00621 * 51.49550 - 0.24958)
      + COALESCE(0.10559 * hemoglobin, 0.10559 * 10.32307 + 0.05954)
      + COALESCE(0.13251 * so2, 0.13251 * 93.66539 - 0.23172)
      + COALESCE(-0.01511 * pco2, -0.01511 * 42.08866 - 0.01630)
      + COALESCE(0.01480 * fio2, 0.01480 * 63.97836 - 0.31142)
      + COALESCE(-0.00200 * aado2, -0.00200 * 442.21186 - 0.01328)
      + COALESCE(-0.03220 * bicarbonate, -0.03220 * 22.96894 - 0.06535)
      + COALESCE(0.05384 * totalco2, 0.05384 * 24.72632 - 0.01405)
      + COALESCE(0.08202 * lactate, 0.08202 * 3.06436 + 0.06038)
      + COALESCE(0.10956 * ph, 0.10956 * 7.36233 - 0.00617)
      + COALESCE(0.00848 * o2flow, 0.00848 * 7.59362 - 0.35803)
    ))) AS SPECIMEN_PROB
  FROM stg2 bg
  LEFT JOIN stg_fio2 s2
    ON bg.icustay_id = s2.icustay_id
    AND s2.charttime BETWEEN DATETIME_SUB(bg.charttime, INTERVAL 4 HOUR) AND bg.charttime
  WHERE bg.lastRowSpO2 = 1
)
SELECT icustay_id, charttime, SPECIMEN,
  CASE WHEN SPECIMEN IS NOT NULL THEN SPECIMEN WHEN SPECIMEN_PROB > 0.75 THEN 'ART' END AS SPECIMEN_PRED,
  SPECIMEN_PROB, SO2, spo2, PO2, PCO2, fio2_chartevents, FIO2, AADO2,
  CASE WHEN PO2 IS NOT NULL AND pco2 IS NOT NULL AND COALESCE(FIO2, fio2_chartevents) IS NOT NULL
       THEN (COALESCE(FIO2, fio2_chartevents) / 100) * (760 - 47) - (pco2 / 0.8) - po2 END AS AADO2_calc,
  CASE WHEN PO2 IS NOT NULL AND COALESCE(FIO2, fio2_chartevents) IS NOT NULL
       THEN 100 * PO2 / COALESCE(FIO2, fio2_chartevents) END AS PaO2FiO2,
  PH, BASEEXCESS, BICARBONATE, TOTALCO2, HEMATOCRIT, HEMOGLOBIN, CARBOXYHEMOGLOBIN, METHEMOGLOBIN,
  CHLORIDE, CALCIUM, TEMPERATURE, POTASSIUM, SODIUM, LACTATE, GLUCOSE,
  INTUBATED, TIDALVOLUME, VENTILATIONRATE, VENTILATOR, PEEP, O2Flow, REQUIREDO2
FROM stg3
WHERE lastRowFiO2 = 1
  AND (SPECIMEN = 'ART' OR SPECIMEN_PROB > 0.75);


-- ============================================================================
-- 1.6 sofa_si: SOFA in the suspicion window [-48 h, +24 h]
--   Renal component: creatinine only (no urine output table).
--   Cardiovascular: minimum MAP + maximum vasopressor rates in the window.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.sofa_si` AS
WITH weight AS (
  SELECT icustay_id, AVG(valuenum) AS weight
  FROM `physionet-data.mimiciii_clinical.chartevents`
  WHERE itemid IN (762, 763, 3723, 3580, 226512, 224639) AND valuenum > 0
  GROUP BY icustay_id
),
vaso_cv AS (
  SELECT ie.icustay_id,
    MAX(CASE WHEN cv.itemid = 30047 THEN cv.rate / COALESCE(wt.weight, 70.0) WHEN cv.itemid = 30120 THEN cv.rate END) AS rate_norepinephrine,
    MAX(CASE WHEN cv.itemid = 30044 THEN cv.rate / COALESCE(wt.weight, 70.0) WHEN cv.itemid IN (30119, 30309) THEN cv.rate END) AS rate_epinephrine,
    MAX(CASE WHEN cv.itemid IN (30043, 30307) THEN cv.rate END) AS rate_dopamine,
    MAX(CASE WHEN cv.itemid IN (30042, 30306) THEN cv.rate END) AS rate_dobutamine
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
  INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON s.icustay_id = ie.icustay_id
  INNER JOIN `physionet-data.mimiciii_clinical.inputevents_cv` cv
    ON s.icustay_id = cv.icustay_id AND CAST(cv.charttime AS DATETIME) BETWEEN s.si_starttime AND s.si_endtime
  LEFT JOIN weight wt ON ie.icustay_id = wt.icustay_id
  WHERE cv.itemid IN (30047, 30120, 30044, 30119, 30309, 30043, 30307, 30042, 30306) AND cv.rate IS NOT NULL
  GROUP BY ie.icustay_id
),
vaso_mv AS (
  SELECT ie.icustay_id,
    MAX(CASE WHEN mv.itemid = 221906 THEN mv.rate END) AS rate_norepinephrine,
    MAX(CASE WHEN mv.itemid = 221289 THEN mv.rate END) AS rate_epinephrine,
    MAX(CASE WHEN mv.itemid = 221662 THEN mv.rate END) AS rate_dopamine,
    MAX(CASE WHEN mv.itemid = 221653 THEN mv.rate END) AS rate_dobutamine
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
  INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON s.icustay_id = ie.icustay_id
  INNER JOIN `physionet-data.mimiciii_clinical.inputevents_mv` mv
    ON ie.icustay_id = mv.icustay_id AND CAST(mv.starttime AS DATETIME) BETWEEN s.si_starttime AND s.si_endtime
  WHERE mv.itemid IN (221906, 221289, 221662, 221653) AND mv.statusdescription != 'Rewritten'
  GROUP BY ie.icustay_id
),
pafi1 AS (
  SELECT bg.icustay_id, bg.charttime, bg.PaO2FiO2,
    CASE WHEN bg.PEEP IS NOT NULL OR bg.VENTILATOR IS NOT NULL OR bg.INTUBATED IS NOT NULL THEN 1 ELSE 0 END AS IsVent
  FROM `sepsis-mimic-def.sepsis3_workspace_def.bloodgasarterial_si` bg
),
pafi2 AS (
  SELECT icustay_id,
    MIN(CASE WHEN IsVent = 0 THEN PaO2FiO2 END) AS PaO2FiO2_novent_min,
    MIN(CASE WHEN IsVent = 1 THEN PaO2FiO2 END) AS PaO2FiO2_vent_min
  FROM pafi1 GROUP BY icustay_id
),
scorecomp AS (
  SELECT s.icustay_id,
    v.MeanBP_Min,
    COALESCE(cv.rate_norepinephrine, mv.rate_norepinephrine) AS rate_norepinephrine,
    COALESCE(cv.rate_epinephrine, mv.rate_epinephrine)       AS rate_epinephrine,
    COALESCE(cv.rate_dopamine, mv.rate_dopamine)             AS rate_dopamine,
    COALESCE(cv.rate_dobutamine, mv.rate_dobutamine)         AS rate_dobutamine,
    l.CREATININE_max AS Creatinine_Max,
    l.BILIRUBIN_max  AS Bilirubin_Max,
    l.PLATELET_min   AS Platelet_Min,
    pf.PaO2FiO2_novent_min,
    pf.PaO2FiO2_vent_min,
    gcs.MinGCS
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` s
  LEFT JOIN vaso_cv cv ON s.icustay_id = cv.icustay_id
  LEFT JOIN vaso_mv mv ON s.icustay_id = mv.icustay_id
  LEFT JOIN pafi2 pf   ON s.icustay_id = pf.icustay_id
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.vitals_si` v ON s.icustay_id = v.icustay_id
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.labs_si`   l ON s.icustay_id = l.icustay_id
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.gcs_si`  gcs ON s.icustay_id = gcs.icustay_id
),
scorecalc AS (
  SELECT icustay_id,
    CASE WHEN PaO2FiO2_vent_min < 100 THEN 4 WHEN PaO2FiO2_vent_min < 200 THEN 3
         WHEN PaO2FiO2_novent_min < 300 THEN 2 WHEN PaO2FiO2_novent_min < 400 THEN 1 ELSE 0 END AS respiration,
    CASE WHEN Platelet_Min < 20 THEN 4 WHEN Platelet_Min < 50 THEN 3
         WHEN Platelet_Min < 100 THEN 2 WHEN Platelet_Min < 150 THEN 1 ELSE 0 END AS coagulation,
    CASE WHEN Bilirubin_Max >= 12.0 THEN 4 WHEN Bilirubin_Max >= 6.0 THEN 3
         WHEN Bilirubin_Max >= 2.0 THEN 2 WHEN Bilirubin_Max >= 1.2 THEN 1 ELSE 0 END AS liver,
    CASE WHEN rate_dopamine > 15 OR rate_epinephrine > 0.1 OR rate_norepinephrine > 0.1 THEN 4
         WHEN rate_dopamine > 5 OR rate_epinephrine <= 0.1 OR rate_norepinephrine <= 0.1 THEN 3
         WHEN rate_dopamine > 0 OR rate_dobutamine > 0 THEN 2
         WHEN MeanBP_Min < 70 THEN 1 ELSE 0 END AS cardiovascular,
    CASE WHEN MinGCS BETWEEN 13 AND 14 THEN 1 WHEN MinGCS BETWEEN 10 AND 12 THEN 2
         WHEN MinGCS BETWEEN 6 AND 9 THEN 3 WHEN MinGCS < 6 THEN 4 ELSE 0 END AS cns,
    CASE WHEN Creatinine_Max >= 5.0 THEN 4
         WHEN Creatinine_Max BETWEEN 3.5 AND 4.9 THEN 3
         WHEN Creatinine_Max BETWEEN 2.0 AND 3.4 THEN 2
         WHEN Creatinine_Max BETWEEN 1.2 AND 1.9 THEN 1 ELSE 0 END AS renal
  FROM scorecomp
)
SELECT si.icustay_id,
  COALESCE(s.respiration, 0) + COALESCE(s.coagulation, 0) + COALESCE(s.liver, 0)
  + COALESCE(s.cardiovascular, 0) + COALESCE(s.cns, 0) + COALESCE(s.renal, 0) AS SOFA,
  s.respiration, s.coagulation, s.liver, s.cardiovascular, s.cns, s.renal
FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` si
LEFT JOIN scorecalc s ON si.icustay_id = s.icustay_id;


-- ============================================================================
-- 1.7 sepsis3_final_cohort: final cohort
--   SOFA >= 2, age >= 16 (> 89 recoded to 91.4), with chartevents, MetaVision,
--   first ICU stay, first service not CSURG, VSURG or TSURG.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.sepsis3_final_cohort` AS
WITH icu_stays_ranked AS (
  SELECT ie.subject_id, ie.hadm_id, ie.icustay_id, CAST(ie.intime AS DATETIME) AS intime, ie.dbsource,
         ROW_NUMBER() OVER (PARTITION BY ie.subject_id ORDER BY ie.intime) AS icu_seq
  FROM `physionet-data.mimiciii_clinical.icustays` ie
),
first_service_cte AS (
  SELECT hadm_id, curr_service AS first_service,
         ROW_NUMBER() OVER (PARTITION BY hadm_id ORDER BY transfertime) AS rn
  FROM `physionet-data.mimiciii_clinical.services`
),
cohort_base AS (
  SELECT sp.icustay_id, ie.subject_id, ie.hadm_id,
    sp.suspected_infection_time AS t0_sepsis,
    sofa.SOFA AS sofa_score,
    pat.gender, adm.ethnicity,
    CASE WHEN DATETIME_DIFF(ie.intime, CAST(pat.dob AS DATETIME), YEAR) > 89 THEN 91.4
         ELSE DATETIME_DIFF(ie.intime, CAST(pat.dob AS DATETIME), YEAR) END AS age,
    ie.dbsource, adm.has_chartevents_data, ie.icu_seq, srv.first_service
  FROM `sepsis-mimic-def.sepsis3_workspace_def.suspinfect_poe` sp
  INNER JOIN `sepsis-mimic-def.sepsis3_workspace_def.sofa_si` sofa ON sp.icustay_id = sofa.icustay_id
  INNER JOIN icu_stays_ranked ie ON sp.icustay_id = ie.icustay_id
  INNER JOIN `physionet-data.mimiciii_clinical.admissions` adm ON ie.hadm_id = adm.hadm_id
  INNER JOIN `physionet-data.mimiciii_clinical.patients`   pat ON ie.subject_id = pat.subject_id
  LEFT JOIN first_service_cte srv ON ie.hadm_id = srv.hadm_id AND srv.rn = 1
)
SELECT subject_id, hadm_id, icustay_id, t0_sepsis, sofa_score, age, gender, ethnicity, first_service
FROM cohort_base
WHERE sofa_score >= 2
  AND age >= 16
  AND has_chartevents_data = 1
  AND dbsource = 'metavision'
  AND icu_seq = 1
  AND first_service NOT IN ('CSURG', 'VSURG', 'TSURG')
ORDER BY subject_id;
-- sofa_score is NOT used as a predictor: it uses data up to +24 h after the
-- suspicion of infection (leakage).


-- ============================================================================
-- 1.8 hrv_analysis_mapping: time anchor for the HRV extraction
--   Anchor = suspicion time, or ICU admission + 1 h if the suspicion was before
--   admission or within its first hour. extraccion_hrv_v3.py derives the actual
--   analysis start (inicio_analisis_real) from this anchor and the ECG record.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.hrv_analysis_mapping` AS
SELECT cohorte.subject_id, cohorte.icustay_id, ie.intime, cohorte.t0_sepsis,
  CASE WHEN DATETIME_DIFF(cohorte.t0_sepsis, ie.intime, HOUR) < 1
       THEN DATETIME_ADD(ie.intime, INTERVAL 1 HOUR)
       ELSE cohorte.t0_sepsis END AS hrv_analysis_start_time
FROM `sepsis-mimic-def.sepsis3_workspace_def.sepsis3_final_cohort` cohorte
INNER JOIN `physionet-data.mimiciii_clinical.icustays` ie ON cohorte.icustay_id = ie.icustay_id;
