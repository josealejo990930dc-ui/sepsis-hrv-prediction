-- ============================================================================
-- 03_queries_vigentes_v4.sql: longitudinal tables (one row per prediction window)
--
-- PART A (independent of the window grid):
--   A1 _raw_vaso, A2 _raw_sedantes, A3 _raw_opioides: drug infusions as intervals.
-- PART B (run after generar_grilla_v3.py has uploaded hrv_prediction_points):
--   B0 grid check, B1 vitals, B2 drugs v4, B3 SOFA v3.3, B4 MAP and lactate,
--   B5 drug check, B6 outcome check (B0, B5 and B6 are read-only).
-- ============================================================================


-- ############################################################################
-- PART A: drug infusions as intervals
-- ############################################################################

-- ============================================================================
-- A1) _raw_vaso as intervals
--    MetaVision: [starttime, endtime) of each infusion segment with rate > 0.
--    CareVue: point records, assumed active for 1 h from charttime (residual:
--    the cohort is MetaVision).
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_vaso` AS
WITH weight AS (
  SELECT icustay_id, AVG(valuenum) AS weight
  FROM `physionet-data.mimiciii_clinical.chartevents`
  WHERE itemid IN (762, 763, 3723, 3580, 226512, 224639) AND valuenum > 0
  GROUP BY icustay_id
),
vaso_cv AS (
  SELECT
    cv.icustay_id,
    CAST(cv.charttime AS DATETIME) AS starttime,
    DATETIME_ADD(CAST(cv.charttime AS DATETIME), INTERVAL 1 HOUR) AS endtime,
    CASE WHEN cv.itemid = 30047 THEN cv.rate / COALESCE(wt.weight, 70.0)
         WHEN cv.itemid = 30120 THEN cv.rate END AS rate_norepi,
    CASE WHEN cv.itemid = 30044 THEN cv.rate / COALESCE(wt.weight, 70.0)
         WHEN cv.itemid IN (30119, 30309) THEN cv.rate END AS rate_epi,
    CASE WHEN cv.itemid IN (30043, 30307) THEN cv.rate END AS rate_dopa,
    CASE WHEN cv.itemid IN (30042, 30306) THEN cv.rate END AS rate_dobu
  FROM `physionet-data.mimiciii_clinical.inputevents_cv` cv
  LEFT JOIN weight wt ON cv.icustay_id = wt.icustay_id
  WHERE cv.itemid IN (30047, 30120, 30044, 30119, 30309, 30043, 30307, 30042, 30306)
    AND cv.rate IS NOT NULL AND cv.rate > 0
),
vaso_mv AS (
  SELECT
    mv.icustay_id,
    CAST(mv.starttime AS DATETIME) AS starttime,
    CAST(COALESCE(mv.endtime, mv.starttime) AS DATETIME) AS endtime,
    CASE WHEN mv.itemid = 221906 THEN mv.rate END AS rate_norepi,
    CASE WHEN mv.itemid = 221289 THEN mv.rate END AS rate_epi,
    CASE WHEN mv.itemid = 221662 THEN mv.rate END AS rate_dopa,
    CASE WHEN mv.itemid = 221653 THEN mv.rate END AS rate_dobu
  FROM `physionet-data.mimiciii_clinical.inputevents_mv` mv
  WHERE mv.itemid IN (221906, 221289, 221662, 221653)
    AND mv.statusdescription != 'Rewritten'
    AND mv.rate IS NOT NULL AND mv.rate > 0
)
SELECT * FROM vaso_cv
UNION ALL
SELECT * FROM vaso_mv;


-- ============================================================================
-- A2) _raw_sedantes as intervals
--    Propofol, midazolam, lorazepam and dexmedetomidine (itemid 225150; 222315 is
--    vasopressin). MetaVision: infusions (rate > 0) and boluses (amount > 0,
--    start = end). CareVue: point records.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_sedantes` AS
SELECT icustay_id,
       CAST(charttime AS DATETIME) AS starttime,
       CAST(charttime AS DATETIME) AS endtime
FROM `physionet-data.mimiciii_clinical.inputevents_cv`
WHERE itemid IN (30131, 30124, 30141)   -- propofol, midazolam, lorazepam
  AND rate IS NOT NULL AND rate > 0
UNION ALL
SELECT icustay_id,
       CAST(starttime AS DATETIME) AS starttime,
       CAST(COALESCE(endtime, starttime) AS DATETIME) AS endtime
FROM `physionet-data.mimiciii_clinical.inputevents_mv`
WHERE itemid IN (222168, 221668, 221385, 225150)   -- propofol, midazolam, lorazepam, dexmedetomidine
  AND statusdescription != 'Rewritten'
  AND ((rate IS NOT NULL AND rate > 0) OR (amount IS NOT NULL AND amount > 0));


-- ============================================================================
-- A3) _raw_opioides as intervals (fentanyl, morphine, hydromorphone)
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def._raw_opioides` AS
SELECT icustay_id,
       CAST(charttime AS DATETIME) AS starttime,
       CAST(charttime AS DATETIME) AS endtime
FROM `physionet-data.mimiciii_clinical.inputevents_cv`
WHERE itemid IN (30150, 30126, 30163)
  AND rate IS NOT NULL AND rate > 0
UNION ALL
SELECT icustay_id,
       CAST(starttime AS DATETIME) AS starttime,
       CAST(COALESCE(endtime, starttime) AS DATETIME) AS endtime
FROM `physionet-data.mimiciii_clinical.inputevents_mv`
WHERE itemid IN (221744, 225154, 221833)
  AND statusdescription != 'Rewritten'
  AND ((rate IS NOT NULL AND rate > 0) OR (amount IS NOT NULL AND amount > 0));




-- ############################################################################
-- PART B: run only after uploading the v3 window grid
-- ############################################################################

-- ============================================================================
-- B0) Check that hrv_prediction_points is the v3 grid (run before B1)
--     If the counts are those of the old grid (2981 windows / 499 patients),
--     the v3 grid has not been uploaded yet: stop.
-- ============================================================================
SELECT COUNT(*) AS n_ventanas, COUNT(DISTINCT icustay_id) AS n_pacientes,
       MIN(hour_index) AS min_h, MAX(hour_index) AS max_h
FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`;

-- ============================================================================
-- B1) vitals_longitudinal: values at T0 (MAP LOCF 1 h; RR, SpO2, temperature LOCF 6 h)
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.vitals_longitudinal` AS
WITH puntos AS (
  SELECT icustay_id, hour_index, t0_prediccion
  FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
),
locf AS (
  SELECT
    p.icustay_id, p.hour_index,
    ARRAY_AGG(v.rr_val   IGNORE NULLS ORDER BY v.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS rr_t0,
    ARRAY_AGG(v.spo2_val IGNORE NULLS ORDER BY v.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS spo2_t0,
    ARRAY_AGG(v.temp_val IGNORE NULLS ORDER BY v.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS temp_t0
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
    ON p.icustay_id = v.icustay_id
    AND v.charttime BETWEEN DATETIME_SUB(p.t0_prediccion, INTERVAL 6 HOUR) AND p.t0_prediccion
  GROUP BY p.icustay_id, p.hour_index
),
map_1h AS (
  SELECT
    p.icustay_id, p.hour_index,
    ARRAY_AGG(v.map_val IGNORE NULLS ORDER BY v.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS map_t0
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
    ON p.icustay_id = v.icustay_id
    AND v.charttime BETWEEN DATETIME_SUB(p.t0_prediccion, INTERVAL 1 HOUR) AND p.t0_prediccion
  GROUP BY p.icustay_id, p.hour_index
)
SELECT l.icustay_id, l.hour_index, m.map_t0, l.rr_t0, l.spo2_t0, l.temp_t0
FROM locf l
LEFT JOIN map_1h m USING (icustay_id, hour_index);


-- ============================================================================
-- B2) drugs_longitudinal v4
--    Rate of each vasopressor = rate of the infusion ACTIVE at each hourly point
--    (T0, T0-1h ... T0-11h); overlapping segments of one drug: highest rate.
--    NED and VIS at T0 = point k = 0; cumvp = sum of NED over the 12 points.
--    Sedative/opioid flag = any record overlapping (T0-12h, T0].
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.drugs_longitudinal` AS
WITH puntos AS (
  SELECT icustay_id, hour_index, t0_prediccion
  FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
),
horas AS (
  SELECT p.icustay_id, p.hour_index, k,
         DATETIME_SUB(p.t0_prediccion, INTERVAL k HOUR) AS t_h
  FROM puntos p, UNNEST(GENERATE_ARRAY(0, 11)) AS k
),
vaso_horario AS (
  SELECT
    h.icustay_id, h.hour_index, h.k,
    MAX(v.rate_norepi) AS rate_norepi,
    MAX(v.rate_epi)    AS rate_epi,
    MAX(v.rate_dopa)   AS rate_dopa,
    MAX(v.rate_dobu)   AS rate_dobu
  FROM horas h
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vaso` v
    ON h.icustay_id = v.icustay_id
    AND v.starttime <= h.t_h AND v.endtime > h.t_h
  GROUP BY h.icustay_id, h.hour_index, h.k
),
ned_horario AS (
  SELECT
    icustay_id, hour_index, k,
    rate_norepi, rate_epi, rate_dopa, rate_dobu,
    COALESCE(rate_norepi, 0) + COALESCE(rate_epi, 0) + COALESCE(rate_dopa, 0) / 100 AS ned,
    COALESCE(rate_dopa, 0) + COALESCE(rate_dobu, 0)
      + 100 * COALESCE(rate_epi, 0) + 100 * COALESCE(rate_norepi, 0) AS vis
  FROM vaso_horario
),
cum AS (
  SELECT icustay_id, hour_index, SUM(ned) AS cumvp
  FROM ned_horario
  GROUP BY icustay_id, hour_index
),
sed_flag AS (
  SELECT p.icustay_id, p.hour_index,
    MAX(CASE WHEN s.starttime IS NOT NULL THEN 1 ELSE 0 END) AS sedante_bin
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_sedantes` s
    ON p.icustay_id = s.icustay_id
    AND s.starttime <= p.t0_prediccion
    AND s.endtime >= DATETIME_SUB(p.t0_prediccion, INTERVAL 12 HOUR)
  GROUP BY p.icustay_id, p.hour_index
),
opi_flag AS (
  SELECT p.icustay_id, p.hour_index,
    MAX(CASE WHEN o.starttime IS NOT NULL THEN 1 ELSE 0 END) AS opioide_bin
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_opioides` o
    ON p.icustay_id = o.icustay_id
    AND o.starttime <= p.t0_prediccion
    AND o.endtime >= DATETIME_SUB(p.t0_prediccion, INTERVAL 12 HOUR)
  GROUP BY p.icustay_id, p.hour_index
)
SELECT
  n.icustay_id, n.hour_index,
  n.ned, n.vis, c.cumvp,
  n.rate_norepi, n.rate_epi, n.rate_dopa, n.rate_dobu,
  COALESCE(sf.sedante_bin, 0) AS sedante_bin,
  COALESCE(opi.opioide_bin, 0) AS opioide_bin
FROM ned_horario n
JOIN cum c USING (icustay_id, hour_index)
LEFT JOIN sed_flag sf USING (icustay_id, hour_index)
LEFT JOIN opi_flag opi USING (icustay_id, hour_index)
WHERE n.k = 0;

-- ============================================================================
-- B3) sofa_longitudinal v3.3
--    Evaluated at T0 and at T0 + 6 h. LOCF windows: respiratory 12 h, CNS 6 h,
--    renal and coagulation 24 h, liver 48 h; cardiovascular = MAP (LOCF 1 h) +
--    vasopressors active at the evaluation time.
--    Ties at the same charttime take the WORST value (deterministic).
--    *_valid_t0 / *_valid_6h show which component makes the SOFA non-evaluable.
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.sofa_longitudinal` AS

WITH puntos AS (
  SELECT icustay_id, hour_index, t0_prediccion
  FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
),

-- Two evaluation times per T0: T0 (offset 0) and T0 + 6 h (Mentasti horizon)
eval_points AS (
  SELECT icustay_id, hour_index, t0_prediccion AS eval_time, 0 AS eval_offset FROM puntos
  UNION ALL
  SELECT icustay_id, hour_index, DATETIME_ADD(t0_prediccion, INTERVAL 6 HOUR), 6 FROM puntos
),

-- ---------------- Respiratory (PaO2/FiO2, SpO2/FiO2 fallback), 12 h window ----------------
-- 12 h rather than 6 h: arterial gases are not hourly, and with 6 h the component was
-- available in only ~26 % of rows. 12 h equals the longest horizon considered in the
-- study and recovers ~70 % of the gain of a 24 h window.
resp_pf AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(STRUCT(bg.PaO2FiO2 AS pf,
                      CASE WHEN bg.PEEP IS NOT NULL OR bg.VENTILATOR IS NOT NULL
                                OR bg.INTUBATED IS NOT NULL THEN 1 ELSE 0 END AS is_vent)
              ORDER BY bg.charttime DESC, bg.PaO2FiO2 ASC LIMIT 1)[SAFE_OFFSET(0)] AS locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.bloodgasarterial_si` bg
    ON e.icustay_id = bg.icustay_id
    AND bg.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 12 HOUR) AND e.eval_time
    AND bg.PaO2FiO2 IS NOT NULL
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),
resp_sf_raw AS (
  -- SpO2/FiO2 is used only when PaO2/FiO2 is missing and SpO2 <= 97
  -- (Pandharipande et al. 2009 equivalence, applied regardless of ventilation status)
  SELECT
    v.icustay_id, v.charttime,
    100.0 * v.spo2_val / f.fio2_val AS sf_ratio
  FROM `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
  INNER JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_fio2` f
    ON v.icustay_id = f.icustay_id
    AND f.charttime BETWEEN DATETIME_SUB(v.charttime, INTERVAL 4 HOUR) AND v.charttime
  WHERE v.spo2_val IS NOT NULL AND v.spo2_val <= 97 AND f.fio2_val IS NOT NULL
),
resp_sf AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(sf.sf_ratio ORDER BY sf.charttime DESC, sf.sf_ratio ASC LIMIT 1)[SAFE_OFFSET(0)] AS sf_locf
  FROM eval_points e
  LEFT JOIN resp_sf_raw sf
    ON e.icustay_id = sf.icustay_id
    AND sf.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 12 HOUR) AND e.eval_time
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),
resp_score AS (
  SELECT
    p.icustay_id, p.hour_index, p.eval_offset,
    CASE
      WHEN p.locf.pf IS NOT NULL THEN
        CASE WHEN p.locf.is_vent = 1 AND p.locf.pf < 100 THEN 4
             WHEN p.locf.is_vent = 1 AND p.locf.pf < 200 THEN 3
             WHEN p.locf.pf < 300 THEN 2
             WHEN p.locf.pf < 400 THEN 1
             ELSE 0 END
      WHEN s.sf_locf IS NOT NULL THEN
        CASE WHEN s.sf_locf <= 66  THEN 4
             WHEN s.sf_locf <= 141 THEN 3
             WHEN s.sf_locf <= 220 THEN 2
             WHEN s.sf_locf <= 301 THEN 1
             ELSE 0 END
      ELSE NULL
    END AS resp,
    CASE WHEN p.locf.pf IS NOT NULL OR s.sf_locf IS NOT NULL THEN 1 ELSE 0 END AS resp_valid
  FROM resp_pf p
  LEFT JOIN resp_sf s USING (icustay_id, hour_index, eval_offset)
),

-- ---------------- CNS (GCS), 6 h window ----------------
neuro_score AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(g.gcs ORDER BY g.charttime DESC, g.gcs ASC LIMIT 1)[SAFE_OFFSET(0)] AS gcs_locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_gcs` g
    ON e.icustay_id = g.icustay_id
    AND g.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 6 HOUR) AND e.eval_time
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),

-- ---------------- Renal (creatinine) and coagulation (platelets), 24 h window ----------------
renal_coag_score AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(l.creatinine IGNORE NULLS ORDER BY l.charttime DESC, l.creatinine DESC LIMIT 1)[SAFE_OFFSET(0)] AS creatinine_locf,
    ARRAY_AGG(l.platelet   IGNORE NULLS ORDER BY l.charttime DESC, l.platelet ASC LIMIT 1)[SAFE_OFFSET(0)] AS platelet_locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_labs` l
    ON e.icustay_id = l.icustay_id
    AND l.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 24 HOUR) AND e.eval_time
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),

-- ---------------- Liver (bilirubin), 48 h window ----------------
hepatic_score AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(l.bilirubin IGNORE NULLS ORDER BY l.charttime DESC, l.bilirubin DESC LIMIT 1)[SAFE_OFFSET(0)] AS bilirubin_locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_labs` l
    ON e.icustay_id = l.icustay_id
    AND l.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 48 HOUR) AND e.eval_time
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),

-- ---------------- Cardiovascular: MAP (LOCF 1 h) + vasopressors ACTIVE at eval_time ----------------
-- A vasopressor counts when starttime <= t < endtime. MAP and vasopressors are
-- aggregated separately so their rows are not cross-joined.
map_ev AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    ARRAY_AGG(v.map_val IGNORE NULLS ORDER BY v.charttime DESC, v.map_val ASC LIMIT 1)[SAFE_OFFSET(0)] AS map_locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
    ON e.icustay_id = v.icustay_id
    AND v.charttime BETWEEN DATETIME_SUB(e.eval_time, INTERVAL 1 HOUR) AND e.eval_time
    AND v.map_val IS NOT NULL
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),
vaso_ev AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,
    MAX(va.rate_norepi) AS norepi_locf,
    MAX(va.rate_epi)    AS epi_locf,
    MAX(va.rate_dopa)   AS dopa_locf,
    MAX(va.rate_dobu)   AS dobu_locf
  FROM eval_points e
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vaso` va
    ON e.icustay_id = va.icustay_id
    AND va.starttime <= e.eval_time AND va.endtime > e.eval_time
  GROUP BY e.icustay_id, e.hour_index, e.eval_offset
),
cardio_score AS (
  SELECT m.icustay_id, m.hour_index, m.eval_offset, m.map_locf,
         v.norepi_locf, v.epi_locf, v.dopa_locf, v.dobu_locf
  FROM map_ev m
  JOIN vaso_ev v USING (icustay_id, hour_index, eval_offset)
),

-- ---------------- Consolidation and component scores ----------------
-- Start from eval_points and join every component on (icustay_id, hour_index,
-- eval_offset) so the two evaluation times are never mixed.
scored AS (
  SELECT
    e.icustay_id, e.hour_index, e.eval_offset,

    r.resp, r.resp_valid,

    CASE WHEN rc.platelet_locf < 20 THEN 4 WHEN rc.platelet_locf < 50 THEN 3
         WHEN rc.platelet_locf < 100 THEN 2 WHEN rc.platelet_locf < 150 THEN 1 ELSE 0 END AS coag,
    CASE WHEN rc.platelet_locf IS NOT NULL THEN 1 ELSE 0 END AS coag_valid,

    CASE WHEN he.bilirubin_locf >= 12.0 THEN 4 WHEN he.bilirubin_locf >= 6.0 THEN 3
         WHEN he.bilirubin_locf >= 2.0 THEN 2 WHEN he.bilirubin_locf >= 1.2 THEN 1 ELSE 0 END AS liver,
    CASE WHEN he.bilirubin_locf IS NOT NULL THEN 1 ELSE 0 END AS liver_valid,

    CASE WHEN ca.dopa_locf > 15 OR ca.epi_locf > 0.1 OR ca.norepi_locf > 0.1 THEN 4
         WHEN ca.dopa_locf > 5 OR ca.epi_locf > 0 OR ca.norepi_locf > 0 THEN 3
         WHEN ca.dopa_locf > 0 OR ca.dobu_locf > 0 THEN 2
         WHEN ca.map_locf < 70 THEN 1 ELSE 0 END AS cardio,
    CASE WHEN ca.map_locf IS NOT NULL THEN 1 ELSE 0 END AS cardio_valid,

    CASE WHEN n.gcs_locf BETWEEN 13 AND 14 THEN 1 WHEN n.gcs_locf BETWEEN 10 AND 12 THEN 2
         WHEN n.gcs_locf BETWEEN 6 AND 9 THEN 3 WHEN n.gcs_locf < 6 THEN 4 ELSE 0 END AS neuro,
    CASE WHEN n.gcs_locf IS NOT NULL THEN 1 ELSE 0 END AS neuro_valid,

    CASE WHEN rc.creatinine_locf >= 5.0 THEN 4
         WHEN rc.creatinine_locf BETWEEN 3.5 AND 4.9 THEN 3
         WHEN rc.creatinine_locf BETWEEN 2.0 AND 3.4 THEN 2
         WHEN rc.creatinine_locf BETWEEN 1.2 AND 1.9 THEN 1 ELSE 0 END AS renal,
    CASE WHEN rc.creatinine_locf IS NOT NULL THEN 1 ELSE 0 END AS renal_valid

  FROM eval_points e
  LEFT JOIN resp_score       r  USING (icustay_id, hour_index, eval_offset)
  LEFT JOIN renal_coag_score rc USING (icustay_id, hour_index, eval_offset)
  LEFT JOIN hepatic_score    he USING (icustay_id, hour_index, eval_offset)
  LEFT JOIN cardio_score     ca USING (icustay_id, hour_index, eval_offset)
  LEFT JOIN neuro_score      n  USING (icustay_id, hour_index, eval_offset)
),

totals AS (
  SELECT
    icustay_id, hour_index, eval_offset,
    resp, coag, liver, cardio, neuro, renal,
    resp_valid, coag_valid, liver_valid, cardio_valid, neuro_valid, renal_valid,
    (COALESCE(resp,0)+COALESCE(coag,0)+COALESCE(liver,0)
     +COALESCE(cardio,0)+COALESCE(neuro,0)+COALESCE(renal,0)) AS sofa_total,
    CASE WHEN resp_valid=1 AND coag_valid=1 AND liver_valid=1
          AND cardio_valid=1 AND neuro_valid=1 AND renal_valid=1
         THEN 1 ELSE 0 END AS sofa_valid
  FROM scored
)

SELECT
  t0.icustay_id, t0.hour_index,

  -- ---- Complete-case total SOFA (sensitivity analysis) ----
  CASE WHEN t0.sofa_valid = 1 THEN t0.sofa_total END AS sofa_t0,
  t0.sofa_valid AS sofa_t0_valid,
  CASE WHEN h6.sofa_valid = 1 THEN h6.sofa_total END AS sofa_t0_6h,
  h6.sofa_valid AS sofa_t0_6h_valid,

  -- ---- Components at T0 (NULL if not measured) ----
  CASE WHEN t0.resp_valid   = 1 THEN t0.resp   END AS sofa_resp_t0,
  CASE WHEN t0.coag_valid   = 1 THEN t0.coag   END AS sofa_coag_t0,
  CASE WHEN t0.liver_valid  = 1 THEN t0.liver  END AS sofa_liver_t0,
  CASE WHEN t0.cardio_valid = 1 THEN t0.cardio END AS sofa_cardio_t0,
  CASE WHEN t0.neuro_valid  = 1 THEN t0.neuro  END AS sofa_neuro_t0,
  CASE WHEN t0.renal_valid  = 1 THEN t0.renal  END AS sofa_renal_t0,

  -- ---- Model feature: partial SOFA at T0 = sum of available components ----
  (IF(t0.resp_valid=1, t0.resp, 0) + IF(t0.coag_valid=1, t0.coag, 0) + IF(t0.liver_valid=1, t0.liver, 0)
   + IF(t0.cardio_valid=1, t0.cardio, 0) + IF(t0.neuro_valid=1, t0.neuro, 0) + IF(t0.renal_valid=1, t0.renal, 0)) AS sofa_t0_parcial,
  (t0.resp_valid + t0.coag_valid + t0.liver_valid
   + t0.cardio_valid + t0.neuro_valid + t0.renal_valid) AS n_comp_validos_t0,

  -- ---- Paired delta (option C) ----
  -- Coagulation, cardiovascular, CNS, renal: required at T0 and T0 + 6 h.
  -- Respiratory and liver: included if measured at both times; contribute 0 if
  -- missing at both; window not evaluable if measured at only one.
  CASE WHEN t0.coag_valid*h6.coag_valid*t0.cardio_valid*h6.cardio_valid
            *t0.neuro_valid*h6.neuro_valid*t0.renal_valid*h6.renal_valid = 1
        AND t0.resp_valid  = h6.resp_valid
        AND t0.liver_valid = h6.liver_valid
       THEN 1 ELSE 0 END AS delta_pareado_valid,
  CASE WHEN t0.coag_valid*h6.coag_valid*t0.cardio_valid*h6.cardio_valid
            *t0.neuro_valid*h6.neuro_valid*t0.renal_valid*h6.renal_valid = 1
        AND t0.resp_valid  = h6.resp_valid
        AND t0.liver_valid = h6.liver_valid
       THEN (h6.coag - t0.coag) + (h6.cardio - t0.cardio)
          + (h6.neuro - t0.neuro) + (h6.renal - t0.renal)
          + IF(t0.resp_valid  = 1, h6.resp  - t0.resp,  0)
          + IF(t0.liver_valid = 1, h6.liver - t0.liver, 0)
  END AS delta_sofa_pareado,

  -- ---- Per-component validity at T0 and T0 + 6 h (diagnostic) ----
  t0.resp_valid AS resp_valid_t0,     h6.resp_valid AS resp_valid_6h,
  t0.coag_valid AS coag_valid_t0,     h6.coag_valid AS coag_valid_6h,
  t0.liver_valid AS liver_valid_t0,   h6.liver_valid AS liver_valid_6h,
  t0.cardio_valid AS cardio_valid_t0, h6.cardio_valid AS cardio_valid_6h,
  t0.neuro_valid AS neuro_valid_t0,   h6.neuro_valid AS neuro_valid_6h,
  t0.renal_valid AS renal_valid_t0,   h6.renal_valid AS renal_valid_6h
FROM totals t0
LEFT JOIN totals h6 ON t0.icustay_id=h6.icustay_id AND t0.hour_index=h6.hour_index AND h6.eval_offset=6
WHERE t0.eval_offset = 0;


-- ============================================================================
-- B4) cardio_lactato_longitudinal
--   Observation = (T0-12h, T0]; horizon = (T0, T0+6h].
--   Lactate ties: latest charttime, then highest value (deterministic).
-- ============================================================================
CREATE OR REPLACE TABLE `sepsis-mimic-def.sepsis3_workspace_def.cardio_lactato_longitudinal` AS
WITH puntos AS (
  SELECT icustay_id, hour_index, t0_prediccion
  FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
),
mbp_actual AS (
  SELECT p.icustay_id, p.hour_index,
    AVG(v.map_val)   AS mbp_mean_actual12h,
    COUNT(v.map_val) AS n_mbp_actual12h
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
    ON p.icustay_id = v.icustay_id AND v.map_val IS NOT NULL
    AND v.charttime > DATETIME_SUB(p.t0_prediccion, INTERVAL 12 HOUR)
    AND v.charttime <= p.t0_prediccion
  GROUP BY p.icustay_id, p.hour_index
),
mbp_siguiente AS (
  SELECT p.icustay_id, p.hour_index,
    AVG(v.map_val)          AS mbp_mean_next6h,
    COUNT(v.map_val)        AS n_mbp_next6h,
    COUNTIF(v.map_val > 65) AS n_mbp_next6h_gt65
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_vitals` v
    ON p.icustay_id = v.icustay_id AND v.map_val IS NOT NULL
    AND v.charttime > p.t0_prediccion
    AND v.charttime <= DATETIME_ADD(p.t0_prediccion, INTERVAL 6 HOUR)
  GROUP BY p.icustay_id, p.hour_index
),
lact_actual AS (
  SELECT p.icustay_id, p.hour_index,
    ARRAY_AGG(l.lactate ORDER BY l.charttime DESC, l.lactate DESC LIMIT 1)[SAFE_OFFSET(0)] AS lactate_last_actual12h
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_lactate` l
    ON p.icustay_id = l.icustay_id
    AND l.charttime > DATETIME_SUB(p.t0_prediccion, INTERVAL 12 HOUR)
    AND l.charttime <= p.t0_prediccion
  GROUP BY p.icustay_id, p.hour_index
),
lact_siguiente AS (
  SELECT p.icustay_id, p.hour_index,
    ARRAY_AGG(l.lactate ORDER BY l.charttime DESC, l.lactate DESC LIMIT 1)[SAFE_OFFSET(0)] AS lactate_last_next6h
  FROM puntos p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def._raw_lactate` l
    ON p.icustay_id = l.icustay_id
    AND l.charttime > p.t0_prediccion
    AND l.charttime <= DATETIME_ADD(p.t0_prediccion, INTERVAL 6 HOUR)
  GROUP BY p.icustay_id, p.hour_index
)
SELECT
  p.icustay_id, p.hour_index,
  ma.mbp_mean_actual12h, ma.n_mbp_actual12h,
  ms.mbp_mean_next6h, ms.n_mbp_next6h, ms.n_mbp_next6h_gt65,
  CASE WHEN ms.n_mbp_next6h > 0 THEN ms.n_mbp_next6h_gt65 / ms.n_mbp_next6h END AS mbp_pct_next6h_gt65,
  CASE WHEN ma.mbp_mean_actual12h > 0
       THEN (ms.mbp_mean_next6h - ma.mbp_mean_actual12h) / ma.mbp_mean_actual12h END AS mbp_incremento_pct,
  la.lactate_last_actual12h,
  ls.lactate_last_next6h
FROM puntos p
LEFT JOIN mbp_actual     ma USING (icustay_id, hour_index)
LEFT JOIN mbp_siguiente  ms USING (icustay_id, hour_index)
LEFT JOIN lact_actual    la USING (icustay_id, hour_index)
LEFT JOIN lact_siguiente ls USING (icustay_id, hour_index);

-- ============================================================================
-- B5) Check: drugs and cardiovascular component (analytic cohort, >= 6 h HRV)
-- ============================================================================
SELECT
  COUNT(*)                                             AS n_ventanas,
  ROUND(AVG(IF(d.ned > 0, 1, 0)) * 100, 1)             AS pct_ned_positivo,
  ROUND(AVG(d.sedante_bin) * 100, 1)                   AS pct_sedante,
  ROUND(AVG(d.opioide_bin) * 100, 1)                   AS pct_opioide,
  ROUND(AVG(IF(s.sofa_cardio_t0 >= 3, 1, 0)) * 100, 1) AS pct_cardio_3_4,
  APPROX_QUANTILES(IF(d.ned > 0, d.ned, NULL), 4)      AS cuartiles_ned_positivo
FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points` p
JOIN `sepsis-mimic-def.sepsis3_workspace_def.drugs_longitudinal` d USING (icustay_id, hour_index)
JOIN `sepsis-mimic-def.sepsis3_workspace_def.sofa_longitudinal`  s USING (icustay_id, hour_index)
WHERE p.n_horas_validas >= 6;

-- Duration of MetaVision infusion segments (endtime sanity check)
SELECT
  COUNT(*) AS n_tramos,
  COUNTIF(endtime <= starttime) AS n_sin_duracion,
  APPROX_QUANTILES(DATETIME_DIFF(endtime, starttime, MINUTE), 4) AS cuartiles_min
FROM `sepsis-mimic-def.sepsis3_workspace_def._raw_vaso`;

-- ============================================================================
-- B6) Check: outcomes with the paired delta (>= 6 h HRV)
--     Quick summary; the final labels (including death and discharge in the
--     horizon) are built by construir_ultimate_dataset_v3.py.
-- ============================================================================
WITH p AS (
  SELECT icustay_id, hour_index
  FROM `sepsis-mimic-def.sepsis3_workspace_def.hrv_prediction_points`
  WHERE n_horas_validas >= 6
),
lab AS (
  SELECT p.icustay_id,
    CASE WHEN s.delta_pareado_valid = 1 THEN IF(s.delta_sofa_pareado < 0, 1, 0)  END AS label_sofa_responder,
    CASE WHEN s.delta_pareado_valid = 1 THEN IF(s.delta_sofa_pareado >= 2, 1, 0) END AS label_sofa_deterioro_2pts,
    CASE WHEN c.n_mbp_actual12h > 0 AND c.n_mbp_next6h > 0
          AND c.lactate_last_actual12h IS NOT NULL AND c.lactate_last_next6h IS NOT NULL THEN
      IF( (COALESCE(c.mbp_pct_next6h_gt65, 0) > 0.5 OR COALESCE(c.mbp_incremento_pct, 0) >= 0.10)
          AND (   (c.lactate_last_actual12h < 2 AND c.lactate_last_next6h < 2)
               OR (c.lactate_last_actual12h BETWEEN 2 AND 4 AND c.lactate_last_next6h < 2)
               OR (c.lactate_last_actual12h > 4 AND c.lactate_last_next6h <= 4) ), 1, 0)
    END AS label_cardio
  FROM p
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.sofa_longitudinal` s USING (icustay_id, hour_index)
  LEFT JOIN `sepsis-mimic-def.sepsis3_workspace_def.cardio_lactato_longitudinal` c USING (icustay_id, hour_index)
)
SELECT
  COUNT(*) AS n_ventanas, COUNT(DISTINCT icustay_id) AS n_pacientes,
  COUNTIF(label_sofa_responder IS NOT NULL) AS n_eval_sofa,
  COUNT(DISTINCT IF(label_sofa_responder IS NOT NULL, icustay_id, NULL)) AS pac_eval_sofa,
  ROUND(AVG(label_sofa_responder) * 100, 1)      AS pct_responder_sofa,
  ROUND(AVG(label_sofa_deterioro_2pts) * 100, 1) AS pct_deterioro_2pts,
  COUNTIF(label_cardio IS NOT NULL) AS n_eval_cardio,
  COUNT(DISTINCT IF(label_cardio IS NOT NULL, icustay_id, NULL)) AS pac_eval_cardio,
  ROUND(AVG(label_cardio) * 100, 1) AS pct_responder_cardio,
  COUNTIF(label_sofa_responder IS NOT NULL AND label_cardio IS NOT NULL) AS n_eval_comun,
  COUNT(DISTINCT IF(label_sofa_responder IS NOT NULL AND label_cardio IS NOT NULL, icustay_id, NULL)) AS pac_eval_comun
FROM lab;
