# Pull Portugal QISAR2 data and process for GLOW analysis
library(RODBC)
library(tidyverse)
library(janitor)
library(rspiro)

# Notes ----

# To Do ----
# Visits - done
# Patient details - done
# Lifestyle - done
# Comorbidities
# Allergens
# HCRU
# Exacerbations
# Asthma control - done
# Spirometry - done
# Biomarkers
# Medications - probably better use to supplement table
# LTOCS
# Biologics - use biologic episodes data

# Pull Portugal data ----
## Connect to ISAR SQL tables
connection <- RODBC::odbcDriverConnect(
  'driver={SQL Server};
                              server=dsimport,4071;
                              database=QISAR2;
                              trusted_connection=true'
)

prt_data <- c(
  "RR_allergens",
  "RR_asthma_medications",
  "RR_biomarkers",
  "RR_exacerbation_history_detail",
  "RR_patient",
  "RR_patient_visits",
  "RR_safety",
  "RR_spirometry_predicted_values",
  "MT_Patient_Stats",
  "ST_biologic_episodes",
  "ST_antibiotic_courses",
  "data_dictionary",
  "data_dictionary_raw"
) |>
  set_names(nm = \(x) {
    str_remove(x, "_detail$") |>
      str_remove("^RR_") |>
      str_remove("^ST_UK_") |>
      str_remove("_predicted_values$") |>
      str_to_lower()
  }) |>
  map(\(table) {
    sqlQuery(
      connection,
      paste0(
        'SELECT * FROM QISAR2.prt.',
        table
      ),
      as.is = T
    )
  })

save(prt_data, file = "qisar2_prt_data.RData")

close(connection)

# Specifications ----
remove_empty_check_fields <- function(df) {
  # Find check columns where all values are NA
  all_empty_check <- df |>
    select(ends_with("_check")) |>
    select(where(\(x) all(is.na(x)))) |>
    names()

  # Build full set of columns to drop
  cols_to_drop <- c(
    all_empty_check
  )

  print("No check values in columns:")
  print(cols_to_drop)
  print("Columns Will be dropped.")

  # Drop columns
  df |> select(-any_of(cols_to_drop))
}

# Load data ----
load("qisar2_prt_data.RData")

#Split named list into objects
list2env(prt_data, envir = .GlobalEnv)


# Drop empty columns for sanity when reviewing DF
# patient_visits <- patient_visits |>
#   select(where(\(x) any(!is.na(x))))

# Dictionary ----
# variable_list <- gbr_data |>
#   imap(\(table, tablename) {
#     vars <- names(table)
#
#     vars |>
#       map(\(var) {
#         vals <- unique(table[[var]])
#         tibble(
#           qisar2_table = tablename,
#           qisar2_variable = var,
#           qisar2_values = paste(
#             if (length(vals) < 5) vals else vals[1:5],
#             collapse = ", "
#           )
#         )
#       }) |>
#       list_rbind()
#   }) |>
#   list_rbind() |>
#   mutate(
#     n = row_number(),
#     .by = qisar2_table
#   )
#
# write_xlsx(variable_list, "qisar2_gbr_data_dictionary.xlsx")

# Visits ----
visits <- patient |>
  # n = 377
  mutate(
    patient_id = structure(uidpat, label = "Patient ID"),
    date_of_birth = ymd(dob),
    # Provided age is age at baseline according to DD
    age_first_visit = structure(
      as.numeric(age),
      label = "Age at first visit, years"
    )
  ) |>
  select(patient_id, age_first_visit, date_of_birth) |>
  left_join(
    patient_visits |>
      select(
        patient_id = uidpat,
        visit_date = dov,
        uidpatvis
      ) |>
      distinct(),
    by = join_by(patient_id)
  ) |>
  mutate(
    visit_id = structure(
      case_when(
        str_detect(uidpatvis, "barm") ~ str_c(
          patient_id,
          "_baseline"
        ),
        str_detect(uidpatvis, "historic") ~ str_c(
          patient_id,
          "_historic"
        ),
        str_detect(uidpatvis, "f\\d{1,2}arm1") ~ str_c(
          patient_id,
          "_followup",
          str_replace(uidpatvis, ".*f(\\d{1,2})arm1.*", "\\1")
        ),
      ),
      label = "Visit ID"
    ),
    # Visit dates all plausible and no missing
    visit_date = structure(ymd(visit_date), label = "Visit date")
  ) |>
  distinct() |>
  mutate(
    num_of_visits = n(),
    first_visit_date = min(visit_date, na.rm = T),
    last_visit_date = max(visit_date, na.rm = T),
    age_at_visit = if_else(
      is.na(date_of_birth),
      age_first_visit +
        floor(interval(first_visit_date, visit_date) / years(1)),
      floor(interval(date_of_birth, visit_date) / years(1))
    ),
    .by = patient_id
  ) |>
  mutate(
    num_of_visits = structure(num_of_visits, label = "Number of visits"),
    first_visit_date = structure(first_visit_date, label = "First visit date"),
    last_visit_date = structure(last_visit_date, label = "Last visit date"),
    age_at_visit = structure(age_at_visit, label = "Age at visit, years")
  ) |>
  select(-age_first_visit, -date_of_birth) |>
  arrange(patient_id, visit_date, visit_id)
# Duplicate visit dates, keep first for now
# filter_out(
#   row_number() == 2,
#   .by = c(patient_id, visit_date)
# )

# Checks
visits_chks <- visits |>
  mutate(
    visit_order = case_when(
      str_detect(visit_id, "baseline") ~ 0,
      str_detect(visit_id, "historic") ~ -1,
      str_detect(visit_id, "followup") ~ as.numeric(str_extract(
        visit_id,
        "(?<=followup)\\d{1,2}$"
      ))
    )
  ) |>
  arrange(patient_id, visit_date) |>
  mutate(date_order = row_number() - 1, .by = patient_id)

# No patients with visits out of order
visits_chks |>
  filter(visit_order != date_order) |>
  distinct(patient_id) |>
  pull(patient_id)
visits_chks |>
  filter(str_detect(visit_id, "baseline") & date_order != 0) |>
  distinct(patient_id) |>
  pull(patient_id)

# Patient table ----
# Patient-level, i.e. non-varying details, no separate Demographic table
patient_table <- patient |>
  # Drop empty fields
  select(
    patient_id = uidpat,
    where(\(x) any(!is.na(x))),
    -contains("comstat"),
    -record_id
  ) |>
  # Add first and last visit
  left_join(
    visits |>
      distinct(
        patient_id,
        first_visit_date,
        last_visit_date,
        num_of_visits
      ),
    by = join_by(patient_id),
    relationship = "one-to-one"
  ) |>
  mutate(
    date_of_birth = structure(
      ymd(dob),
      label = "Date of birth, imputed as 15th"
    ),
    gender = structure(
      factor(gen, levels = c("Male", "Female")),
      label = "Gender"
    ),
    ethnicity = structure(eth, label = "Ethnicity"),
    asthma_onset = structure(
      as.numeric(asmage),
      label = "Age at asthma onset, years"
    ),
    country = structure("Portugal", label = "Country"),
    setting = structure("Portugal", label = "Setting"),
    site = structure(dagarm, label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system")
  ) |>
  select(
    patient_id,
    country,
    setting,
    site, # site is different to setting
    first_visit_date,
    last_visit_date,
    num_of_visits,
    data_capture_system,
    date_of_birth,
    gender,
    ethnicity,
    asthma_onset,
  ) |>
  arrange(first_visit_date)

# Notes:
# Aaron's SRD approach appears complicated because of need to get additional
# vars: index_date, max_beos_highest, highest_bec_imputed_flag, bt
# Better I think to process biologics, BEC, bt separately and pull in later

# GBR - gender only recorded once at baseline (non-varying in this case)

# Lifestyle ----
# VARS:
# weight - Weight is recorded at each visit and is used as recorded.
# height - Height is the visit's own, else the latest earlier visit's.
# BMI - BMI is recalculated from them, as the R did (QISAR's values are the recalculation), falling back to the recorded IMC; prt_bmi_recorded keeps the recorded figure.
# BSA - BSA is Du Bois, as the UK.
# smoking_status - carried forward from the latest visit that recorded it
# Errors:
# Only one patient appears to have varying height, PRT_57, and I think this is a mistake
# PRT_134 had swapped height/weight, used 70kg for weight as Dom did and use prt_height_cm = 1.62 instead of 0.70

lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    # Weight non-varying except six patients
    wt,
    bmi,
    smoke
  ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    patient |>
      select(
        patient_id = uidpat,
        # Height only recorded at baseline I think, non-varying except one patient
        ht
      ) |>
      distinct(),
    by = join_by(patient_id)
  ) |>
  mutate(
    weight = structure(as.numeric(wt), label = "Weight, kg"),
    height = structure(as.numeric(ht), label = "Height, m"),
    bmi = structure(as.numeric(bmi), label = "BMI, kg/m^2"),
    # PR_134 has error height/weight
    weight = case_when(
      uidpatvis == "PRT_134_barm1" ~ 70,
      .default = weight
    ),
    height = case_when(
      uidpatvis == "PRT_134_barm1" ~ 1.62,
      .default = height
    ),
    bmi = case_when(
      uidpatvis == "PRT_134_barm1" ~ round(weight / (height * height), 1),
      .default = bmi
    ),
    smoking_status = structure(
      factor(smoke, levels = c("Current Smoker", "Ex-Smoker", "Never Smoked")),
      label = "Smoking status, baseline"
    ),
    bmi_chk = round(weight / (height * height), 1)
  ) |>
  select(
    patient_id,
    visit_id,
    uidpatvis,
    visit_date,
    weight,
    height,
    bmi,
    bmi_chk,
    smoking_status
  )

# variable typing, check extreme values

# Comorbidity ----
# UK sends fewer comorbidities:
# allergic_rhinitis - recorded only at baseline, originally yes/no, coded as current/past for ISAR - Incorrect
# GS - don't use this variable as UK doesn't measure well? DP also likes it coded as ever/never

# nasal_polyps - recorded only at baseline, originally yes/no, but coded as current/never for ISAR - Incorrect
# eczema - recorded only at baseline, originally yes/no, but coded as current/never for ISAR - Incorrect
# anxiety - anxiety OR depression, single field, recorded only at baseline
# ocular - recorded only at baseline, UK only records cataracts
# ocular_type - recorded only at baseline, UK only records cataracts
# depression - identical to anxiety, except code 1/2, don't use
# type2_diabetes - UK field is diabetes, not specifically type 2, recorded only at baseline
# osteoporosis - recorded only at baseline, yes/no
# gerd - Gastro-oesophageal reflux, recorded only at baseline, yes/no
# atopic_disease - uk_atopic_disease, recorded only at baseline, yes/no, not calculated like REDCAP

# For other countries/redcap lets code to yes/no at visit rather than current, past, never

# SRD assumptions:
# Comorb-001 all comorbidities: If an “Ever” is recorded at a previous visit it is pulled forward to future visits.
# Comorb-002 [xx_stdate]: If the comorbidity start date is missing:
#   o	If the comorbidity is indicated as present:
#     §	Baseline visit, the visit date – 12 months will be used.
#     §	Follow-up visit, the visit date is used.
# Comorb-003 [xx_stdate]: If a partial date is provided, the middle of the month/year will be imputed.
#   o	If only the year is provided, the date will be imputed as YYYY-06-15
#     §	Unless the visit date is in the same year as the partial date provided and is before June 15th, then the visit date will be used
#   o	If the year and month are provided, the date will be imputed as YYYY-MM-15
# Comorb-004 [pneumonia_stdatex], [cataract_stdatex]: These are the only two comorbidities with multiple start dates.
#   o	Pneumonia is an acute disease and as such has multiple start dates for multiple instances, these are separate into pneumonia_stdate[1,2,3,4,5]
#   o	Cataracts, although a chronic disease, has multiple start dates, presumably for when the disease started in each eye, these are also separated.
# Comorb-005 [crs]: If nasal_polyps is “Ever” CRS will all be set to “Ever”.

comorbidities <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    allergic_rhinitis = allrhinitis,
    chronic_rhinosinusitis = chrrhiniosinusitis,
    nasal_polyps = npolyp,
    eczema,
    anxiety,
    anxiety_stdate = anxiety_dte,
    csd,
    csd_type,
    csd_type_oth,
    csd_dte,
    ocular,
    ocular_type,
    ocular_type_oth,
    ocular_dte,
    obstructive_sleep_apnoea = obsleepapnoea_yn,
    OSA_stdate = obsleepapnoea_stdt,
    renal_failure = renalf,
    renalf_stdate = renalf_dte,
    depression = depress,
    depression_stdate = depress_dte,
    type2_diabetes = t2diab,
    t2diab_stdate = t2diab_dte,
    peptic_ulcer = peptic,
    peptic_ulcer_stdate = peptic_dte,
    pneumonia = pneum,
    pneumonia_stdate = pneum_dte,
    osteoporosis = osteop,
    osteoporosis_stdate = osteop_dte,
    gerd,
    egpa,
    hypereosinophilic_syndrome = hyper,
    uk_atopic_disease
  ) |>
  select(
    where(\(x) any(!is.na(x))), # Keeps columns with any information
  ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  )
# Need to decide how to handle remaining SRD processing - is this necessary or helpful?

# Allergens ----
# Allergens are not provided for UK based on DD
allergens <- patient_visits |>
  select(patient_id = uidpat, uidpatvis, matches("sat|spt")) |>
  # Do we want to check the "no positive allergen" for SPT for SAT and vice-versa, excluded here but in Aaron's SRD script
  mutate(
    # Revised Aaron's code to update and make consistent - need to work on "other"?
    sat_none = sat_posallergen___8 == "No Positive Allergens",
    spt_none = spt_posallergen___12 == "No Positive Allergens",
    across(
      starts_with("sat_posallergen___"),
      \(x) {
        case_when(
          sat_yn == "No" | is.na(sat_yn) ~ NA_character_,
          x == "0" ~ "No",
          is.na(x) &
            sat_yn == "Yes" &
            sat_none ~ "No",
          !is.na(x) & x != "0" ~ "Yes",
          .default = x
        )
      }
      # .names = "{.col}_recode"
    ),
    across(
      starts_with("spt_posallergen___"),
      \(x) {
        case_when(
          spt_yn == "No" | is.na(spt_yn) ~ NA_character_,
          x == "0" ~ "No",
          is.na(x) &
            spt_yn == "Yes" &
            spt_none ~ "No",
          !is.na(x) & x != "0" ~ "Yes",
          .default = x
        )
      }
      # .names = "{.col}_recode"
    )
  ) |>
  select(-contains("_none")) |>
  rename_with(~ str_replace(., "_dte", "_date"), contains("dte")) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(contains("visit_"), .after = patient_id)

# HCRU ----
# exac_num - number of exacerbations in past 12 months
# emerg -	Number of A&E attendances for asthma in past 12 months
# hospadm - Number of hospital admissions for asthma in past 12 months
# invent - Number of episodes of invasive ventilation ever, recorded only at baseline
hcru <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    exac_yn,
    exac_num,
    emerg,
    hospadm,
    invent
  ) |>
  # UK does not provide any exacerbation details
  # left_join(
  #   exacerbation_history |> distinct(uidpatvis, exac_dte),
  #   by = "uidpatvis"
  # ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  mutate(
    exac_num = structure(
      as.numeric(exac_num),
      label = "Number of exacerbation in past 12 months"
    ),
    # All UK reporting is "In the past 12 months"
    exac_timeperiod = "In the past 12 months",
    hosp_adm_num = structure(
      as.numeric(hospadm),
      label = "Number of hospital admissions for asthma in past 12 months"
    ),
    hosp_timeperiod = "In the past 12 months",
    emerg_num = structure(
      as.numeric(emerg),
      label = "Number of A&E attendances for asthma in past 12 months"
    ),
    emerg_timeperiod = "In the past 12 months",
    inv_ventilation_num = structure(
      as.numeric(invent),
      label = "Number of episodes of invasive ventilation ever"
    ),
    inv_ventilation_timeperiod = "Ever, at baseline"
  ) |>
  select(
    patient_id,
    visit_id,
    visit_date,
    exac_num,
    exac_timeperiod,
    hosp_adm_num,
    hosp_timeperiod,
    emerg_num,
    emerg_timeperiod,
    inv_ventilation_num,
    inv_ventilation_timeperiod
  )

# Exacerbations ----
# UK does not provide any exacerbation details
# exacerbation_history |>
#   select(where(\(x) any(!is.na(x)))) |>
#   names()

# Asthma control ----
# Based on Aaron's SRD "Asthma_control_table.R"
# Portugal provides ACT and GINA, not sure how asthma_control derived in PRT data
# Recalculate to prioritise GINA > ACT > ACQ > AIRQ

# New asthma control function allows for not all variables to be present and calculates new asthma_control based on priority and identifies source
# Added to ISARanalysis
clean_asthma_control <- function(df = asthma_control_labelled) {
  gina_vars <- c(
    "gina_day_1",
    "gina_activity_2",
    "gina_night_3",
    "gina_reliever_4"
  )
  control_levels <- c("Uncontrolled", "Partly Controlled", "Well Controlled")
  source_labels <- c(gina_ac = "GINA", act_ac = "ACT", acq_ac = "ACQ")

  has_gina <- all(gina_vars %in% names(df))
  has_act <- "act_score" %in% names(df)
  has_acq <- "acq_score" %in% names(df)

  if (!any(has_gina, has_act, has_acq)) {
    warning(
      "Asthma control variables not found (need all four GINA items, act_score, or acq_score). Returning data unchanged.",
      call. = FALSE
    )
    return(df)
  }
  if (!"country" %in% names(df)) {
    warning(
      "Country variable not found. Returning data unchanged.",
      call. = FALSE
    )
    return(df)
  }

  # GINA: requires all four items
  if (has_gina) {
    df <- df |>
      mutate(
        # Recode gina control questions as yes = 1 and no = 0
        across(
          all_of(gina_vars),
          ~ case_when(
            . == "Yes" ~ 1,
            . == "1" ~ 1,
            . == "No" ~ 0,
            . == "2" ~ 0,
            .default = NA
          )
        ),
        # Calculate gina control score
        gina_score = structure(
          rowSums(pick(all_of(gina_vars)), na.rm = FALSE),
          label = "GINA symptom control score"
        ),
        # Range check
        gina_score_rangefail = gina_score < 0 | gina_score > 4,
        # Classify gina control score
        gina_ac = structure(
          case_when(
            gina_score_rangefail ~ NA,
            gina_score == 0 ~ "Well Controlled",
            gina_score <= 2 ~ "Partly Controlled",
            gina_score <= 4 ~ "Uncontrolled"
          ) |>
            factor(levels = control_levels),
          label = "Asthma control, GINA"
        )
      )
  } else if (any(gina_vars %in% names(df))) {
    message(
      "GINA control not derived: missing ",
      paste(setdiff(gina_vars, names(df)), collapse = ", ")
    )
  }

  # ACT
  if (has_act) {
    df <- df |>
      mutate(
        act_score = structure(
          as.numeric(act_score),
          label = "ACT score"
        ),
        # Range check
        act_score_rangefail = act_score < 5 | act_score > 25,
        act_ac = structure(
          case_when(
            act_score_rangefail ~ NA,
            act_score < 16 ~ "Uncontrolled",
            act_score < 20 ~ "Partly Controlled",
            act_score >= 20 ~ "Well Controlled",
            .default = NA
          ) |>
            factor(levels = control_levels),
          label = "Asthma control, ACT"
        )
      )
  }

  # ACQ (version of test used is unavailable)
  if (has_acq) {
    df <- df |>
      mutate(
        acq_score = structure(
          as.numeric(acq_score),
          label = "ACQ score"
        ),
        # Range check, exclude ACQ scores outside range, 0 - 6, ACQ is mean of 7 questions with 7-point scale (0 = no impairment; 6 = maximum impairment)
        acq_score_rangefail = acq_score < 0 | acq_score > 6,
        acq_ac = structure(
          case_when(
            acq_score_rangefail ~ NA,
            acq_score <= 0.75 ~ "Well Controlled",
            acq_score < 1.5 ~ "Partly Controlled",
            acq_score >= 1.5 ~ "Uncontrolled",
            .default = NA
          ) |>
            factor(levels = control_levels),
          label = "Asthma control, ACQ"
        ),
        acq_version = structure(
          case_when(
            is.na(acq_score) ~ NA,
            country == "UK" & str_detect(visit_id, "baseline") ~ "ACQ7",
            country == "UK" & str_detect(visit_id, "followup") ~ "ACQ6",
            country %in% c("Poland", "Mexico", "Brazil", "Italy") ~ "ACQ7",
            country %in% c("Belgium", "Denmark") ~ "ACQ6",
            !is.na(acq_score) ~ as.character("Unspecified ACQ")
          ),
          label = "ACQ version"
        ),
      )
  }

  # Combined control: gina_ac > act_ac > acq_ac (first non-missing, per row)
  priority <- intersect(names(source_labels), names(df))

  df <- df |>
    mutate(
      # reduce same as coalesce(coalesce(gina_ac, act_ac), acq_ac) but won't fail if column not present
      asthma_control = structure(
        reduce(pick(all_of(priority)), coalesce) |>
          factor(levels = control_levels),
        label = "Asthma control"
      ),
      # identify source of control measure selected
      asthma_control_source = pick(all_of(priority)) |>
        imap(\(x, nm) {
          if_else(is.na(x), NA_character_, source_labels[[nm]])
        }) |>
        reduce(coalesce)
    )
  source_levels <- unname(source_labels[priority])

  if (has_acq) {
    df <- df |>
      mutate(
        asthma_control_source = if_else(
          asthma_control_source == "ACQ",
          acq_version,
          asthma_control_source
        )
      )
    source_levels <- c(
      setdiff(source_levels, "ACQ"),
      c("ACQ7", "ACQ6", "Unspecified ACQ")
    )
  }

  df |>
    mutate(
      asthma_control_source = structure(
        factor(asthma_control_source, levels = source_levels),
        label = "Asthma control source"
      )
    ) |>
    select(
      where(\(x) any(!is.na(x))),
      -country
    ) |>
    relocate(
      patient_id,
      uidpatvis,
      visit_id,
      visit_date,
      any_of(
        c(
          "gina_day_1",
          "gina_activity_2",
          "gina_night_3",
          "gina_reliever_4",
          "gina_lf_5",
          "gina_score",
          "gina_score_rangefail",
          "gina_ac",
          "act_score",
          "act_score_rangefail",
          "act_ac",
          "acq_score",
          "acq_score_rangefail",
          "acq_version",
          "acq_ac",
          "asthma_control",
          "asthma_control_source"
        )
      ),
      everything()
    )
}

asthma_control <- patient_table |>
  select(patient_id, country) |>
  left_join(
    patient_visits |>
      select(
        patient_id = uidpat,
        uidpatvis,
        contains("gina"),
        any_of(
          c(
            "act_score",
            "acq_score"
          )
        ),
        contains("asthma_control")
      ),
    by = join_by(patient_id),
    relationship = "one-to-many"
  ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = patient_id) |>
  rename(
    gina_activity_2 = gina_activity,
    prt_asthma_control = asthma_control
  ) |>
  clean_asthma_control()

# Spirometry ----
# Modification of clean_spirometry function from ISARanalysis
clean_spirometry <- function(
  df = spirometry_labelled,
  fev1_range_low = 0.15,
  fev1_range_high = 8,
  fvc_range_low = 0.15,
  fvc_range_high = 8
) {
  df |>
    rename_with(
      ~ str_remove(.x, "observation_|visit_|spirometry_"),
      .cols = contains("_date")
    ) |>
    # Drop rows with missing measurement date
    filter_out(is.na(date)) |>
    # If post-bronchodilator measures are not available at all time points for a patient, use pre-bronchodilator only to be consistent over time
    mutate(
      exclude_post_fev1_measures = !all(!is.na(post_fev1)),
      exclude_post_fvc_measures = !all(!is.na(post_fvc)),
      .by = patient_id
    ) |>
    # Range checks against values specified above
    mutate(
      across(c(pre_fev1, post_fev1, pre_fvc, post_fvc), \(var) as.numeric(var)),
      across(
        c(pre_fev1, post_fev1),
        \(var) {
          !is.na(var) &
            (var < fev1_range_low | var > fev1_range_high)
        },
        .names = "{.col}_rangefail"
      ),
      across(
        c(pre_fvc, post_fvc),
        \(var) {
          !is.na(var) &
            (var < fvc_range_low | var > fvc_range_high)
        },
        .names = "{.col}_rangefail"
      ),
      # FEV1 should be less than FVC, but could be the same
      pre_fev1_fvc_implausible = !is.na(pre_fev1) &
        !is.na(pre_fvc) &
        pre_fev1 > pre_fvc,
      post_fev1_fvc_implausible = !is.na(post_fev1) &
        !is.na(post_fvc) &
        post_fev1 > post_fvc,
      # Flag records to exclude based on above checks
      pre_fev1_exclude = pre_fev1_rangefail | pre_fev1_fvc_implausible,
      pre_fvc_exclude = pre_fvc_rangefail | pre_fev1_fvc_implausible,
      post_fev1_exclude = post_fev1_rangefail |
        post_fev1_fvc_implausible |
        exclude_post_fev1_measures,
      post_fvc_exclude = post_fvc_rangefail |
        post_fev1_fvc_implausible |
        exclude_post_fvc_measures
    ) |>
    mutate(
      # Create single fev1 measure for baseline and longitudinal description
      fev1 = structure(
        case_when(
          pre_fev1_exclude & post_fev1_exclude ~ NA,
          !is.na(post_fev1) & !post_fev1_exclude ~ post_fev1,
          !is.na(post_fev1) & post_fev1_exclude & !pre_fev1_exclude ~ pre_fev1,
          !is.na(pre_fev1) & !pre_fev1_exclude ~ pre_fev1,
          .default = NA
        ) |>
          as.numeric(),
        label = "FEV1 (L)"
      ),
      # Flag value used for overall measure
      fev1_flag = case_when(
        pre_fev1_exclude & post_fev1_exclude ~ "none",
        !is.na(post_fev1) & !post_fev1_exclude ~ "post",
        !is.na(post_fev1) & post_fev1_exclude & !pre_fev1_exclude ~ "pre",
        !is.na(pre_fev1) & !pre_fev1_exclude ~ "pre",
        .default = NA_character_ # Empty if no FEV1 measure
      ),
      # Create single fvc measure for baseline and longitudinal description
      fvc = structure(
        case_when(
          pre_fvc_exclude & post_fvc_exclude ~ NA,
          !is.na(post_fvc) & !post_fvc_exclude ~ post_fvc,
          !is.na(post_fvc) & post_fvc_exclude & !pre_fvc_exclude ~ pre_fvc,
          !is.na(pre_fvc) & !pre_fvc_exclude ~ pre_fvc,
          .default = NA
        ) |>
          as.numeric(),
        label = "FVC (L)"
      ),
      fvc_flag = case_when(
        pre_fvc_exclude & post_fvc_exclude ~ "none",
        !is.na(post_fvc) & !post_fvc_exclude ~ "post",
        !is.na(post_fvc) & post_fvc_exclude & !pre_fvc_exclude ~ "pre",
        !is.na(pre_fvc) & !pre_fvc_exclude ~ "pre",
        .default = NA_character_ # Empty if no FVC measure
      ),
      # Calculate ratio measure for baseline and longitudinal description
      fev1_fvc_ratio = structure(
        fev1 / fvc,
        label = "FEV1/FVC ratio"
      ),
      # Create categorical ratio measure for baseline and longitudinal description
      fev1_fvc_ratio_cat = structure(
        case_when(
          fev1_fvc_ratio < 0.70 ~ "<0.70",
          fev1_fvc_ratio >= 0.70 ~ "≥0.70",
          .default = NA
        ) |>
          factor(levels = c("<0.70", "≥0.70")),
        label = "FEV1/FVC ratio"
      )
    )
}

# Only one spirometry record provided per visit (i.e. no duplicated uidpatvis)
# sp_date provided for some, use this or visit date if missing
spirometry_cln <- spirometry |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date, age_at_visit),
    by = join_by(uidpatvis)
  ) |>
  select(
    patient_id = uidpat,
    visit_id,
    uidpatvis,
    visit_date,
    sp_date,
    age_at_visit,
    age_lf,
    ht_lf,
    pre_fev1 = bdfev1_pre,
    pre_fvc = bdfvc_pre,
    post_fev1 = bdfev1_post,
    post_fvc = bdfvc_post
  ) |>
  left_join(
    patient_table |>
      select(patient_id, gender, date_of_birth),
    by = join_by(patient_id)
  ) |>
  left_join(
    # Add height at visit
    lifestyle |>
      select(uidpatvis, height),
    by = join_by(uidpatvis)
  ) |>
  mutate(
    # some differences between my calculated age at visit and age_lf, see prt_age_source, recalc from completed sp_date (with visit_date)
    err_age = age_at_visit != as.numeric(age_lf),
    # height looks okay, except digits
    err_ht = height != as.numeric(ht_lf)
  ) |>
  mutate(
    # Use visit date if test not available
    spirometry_date = if_else(is.na(sp_date), visit_date, ymd(sp_date)),
    # Use height at visit (non-varying, except one patient)
    ht_lf = height,
    # Calculate age at test/visit
    age_lf = floor(interval(date_of_birth, spirometry_date) / years(1)),
    across(c(pre_fev1, pre_fvc, post_fev1, post_fvc), as.numeric)
  ) |>
  select(
    patient_id,
    visit_id,
    spirometry_date,
    age_lf,
    ht_lf,
    gender,
    pre_fev1,
    pre_fvc,
    post_fev1,
    post_fvc
  ) |>
  arrange(patient_id, spirometry_date) |>
  clean_spirometry() |>
  # Add derived spirometry measures using Rspiro package
  mutate(
    fev1_GLIgl_lln = rspiro::LLN_GLIgl(age_lf, ht_lf, gender, param = "fev1"),
    fev1_GLIgl_pctpred = rspiro::pctpred_GLIgl(
      age_lf,
      ht_lf,
      gender,
      FEV1 = fev1
    ),
    fev1_GLIgl_zscore = rspiro::zscore_GLIgl(
      age_lf,
      ht_lf,
      gender,
      FEV1 = fev1
    ),
    fvc_GLIgl_lln = rspiro::LLN_GLIgl(age_lf, ht_lf, gender, param = "fvc"),
    fvc_GLIgl_pctpred = rspiro::pctpred_GLIgl(age_lf, ht_lf, gender, FVC = fvc),
    fvc_GLIgl_zscore = rspiro::zscore_GLIgl(age_lf, ht_lf, gender, FVC = fvc),
    fev1_GLIgl_pctpred_cat = structure(
      case_when(
        fev1_GLIgl_pctpred < 80 ~ "<80%",
        fev1_GLIgl_pctpred >= 80 ~ "≥80%",
        .default = NA
      ) |>
        factor(levels = c("<80%", "≥80%")),
      label = "percent predicted FEV1 (%)"
    ),
    fev1_GLIgl_zscore_cat = structure(
      case_when(
        fev1_GLIgl_zscore < -1.64 ~ "< -1.64",
        fev1_GLIgl_zscore >= -1.64 ~ "≥ -1.64",
        .default = NA
      ) |>
        factor(levels = c("< -1.64", "≥ -1.64")),
      label = "FEV1 Z-score"
    )
  )

# None with FVC < FEV1
spirometry_cln |>
  filter(if_any(contains("implausible")) | if_any(contains("rangefail"))) |>
  pull(patient_id) |>
  noquote()

# Biomarkers ----
# PRT doesn't provide test-specific dates for IgE, use visit date
biomarkers_cln <- biomarkers |>
  select(
    patient_id = uidpat,
    uidpatvis,
    where(\(x) any(!is.na(x))), # Keeps columns with any information
    -uiddb,
    -record_id,
    -contains("redcap"),
    -dagarm,
    -beos_yn,
    beos_result = beos_cnt,
    beos_date = beos_dte,
    ige_result = ige_cnt,
    feno_result = feno_rlt,
    -contains("prt_")
  ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  pivot_longer(
    cols = !all_of(c("patient_id", "uidpatvis", "visit_id", "visit_date")),
    names_to = c("test", ".value"),
    names_sep = "_"
  ) |>
  filter_out(is.na(result)) |>
  mutate(
    result = as.numeric(result),
    result = case_when(
      # Truncate IgE and FeNO counts (same as NJH)
      test == "ige" & result <= 2 ~ 2,
      test == "ige" & result >= 5000 ~ 5000,
      test == "feno" & result <= 5 ~ 5,
      test == "feno" & result >= 300 ~ 300,
      .default = result
    ),
    unit = case_when(
      test == "ige" ~ "IU/L", # Confirm this
      test == "beos" ~ "x10^9/L",
      test == "feno" ~ "ppb"
    )
  )

# left_join(visits |> select(visit_id, visit_date), by = join_by(visit_id)) |>
# relocate(visit_date, .after = visit_id)

# result = case_when(
#   # Convert to cells/uL (multiply K/uL or 10*9/L by 1000)
#   unit == "K/uL" ~ eos_value * 1000,
#   unit == "10*9/L" ~ eos_value * 1000,
#   unit == "Unknown" ~ eos_value * 1000, # 6 unknown values, all other tests for these patients in K/uL and results look similar
#   .default = eos_value
# ),
# orig_unit = unit,
# unit = "cells/uL"
# )

# Medications ----
medications <- asthma_medications |>
  # Keep non-empty columns
  select(
    patient_id = uidpat,
    uidpatvis,
    where(\(x) any(!is.na(x))),
    -uiddb,
    -record_id,
    -dagarm,
    -contains("redcap"),
  ) |>
  # Bare drug columns (no underscore) -> "<drug>_name"
  rename_with(
    ~ paste0(.x, "_name"),
    .cols = !all_of(c(
      "patient_id",
      "uidpatvis",
      "MedicationCategory"
    )) &
      !contains("_")
  ) |>
  # All "oth" columns (LAMA, LTRA, MABX) contain only "Unknown"
  mutate(
    # No further information under lama_oth, all "Unknown"
    lama_name = case_when(
      lama_oth == "Unknown" ~ "Other",
      .default = lama_name
    )
    # No further information under ltra_oth and mabx_oth, all "Unknown"
  ) |>
  select(
    -contains("_oth"),
    -contains("uk_"),
    -icslabalama_explicit,
    -contains("reason"),
    -MedicationCategory #Not sure how helpful this variable as it doesn't align with all biologics
  ) |>
  pivot_longer(
    cols = !all_of(c(
      "patient_id",
      "uidpatvis"
    )),
    names_to = c("drug", ".value"),
    names_sep = "_"
  ) |>
  # All (name, ldose, stdate, ongoing, endate) missing where yn == NA & yn == "No"
  # filter(yn == "Yes") |>
  # Range checks
  mutate(
    drug = case_when(
      drug %in% c("antiige", "antiil4", "antiil5") ~ "Biologic",
      drug == "mos" ~ "LTOCS",
      drug == "ics" ~ "ICS only",
      drug == "mabx" ~ "Macrolides",
      .default = str_to_upper(drug)
    ),
    across(contains("date"), ymd)
  ) |>
  left_join(
    visits |> select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  distinct() |>
  filter_out(drug == "Biologic") |>
  arrange(patient_id, visit_date, drug)
# Earliest medication start date is always a year before baseline visit
# Ongoing flag is derived (so not from UK), not sure how helpful this is.
# chk <- medications |> mutate(date_chk = min(stdate, na.rm = T) == (min(visit_date) %m-% months(12)), .by = patient_id, drug)

# Biologics ----
# Use ST biologic episodes file as provided by PRT rather than patient_visits
bx_cln <- st_biologic_episodes |>
  select(
    patient_id = uidpat,
    slot,
    drug_name,
    # Raw start/stop dates checked for parsing errors - none currently
    bx_start_date = start_date,
    bx_end_date = stop_date,
    ongoing,
    contains("stop"),
    contains("switched"),
    bx_entry_type = biologic_entry_type
  ) |>
  # Flag and drop rows with duplicated start dates and drug names, keep record with end date
  # mutate(
  #   dup_start_name = n() > 1,
  #   .by = c(patient_id, start_date, drug_name)
  # ) |>
  # Exclude some duplicated rows
  filter_out(
    when_any(
      patient_id == "GBR01_1116" & is.na(bx_end_date),
      patient_id == "GBR01_1131" & is.na(bx_end_date),
      patient_id == "GBR01_413" & slot == 6,
      patient_id == "GBR01_415" & is.na(bx_end_date),
      patient_id == "GBR03_1108" &
        drug_name == "Mepolizumab" &
        is.na(bx_end_date),
      patient_id == "GBR04_1087" & is.na(bx_end_date),
      # GBR03_1108 has a random benra prescription with no end date
      # GBR01_82 has a lot of blank rows and then random end date = last_visit_date that should be moved up?
      patient_id == "GBR01_413" &
        drug_name == "Reslizumab" &
        bx_end_date == ymd("2020-08-17") # Duplicated episode, keep earlier start
    )
  ) |>
  # mutate(
  #   dup_start_name2 = n() > 1,
  #   .by = c(patient_id, start_date)
  # )
  mutate(
    bx_start_date = ymd(bx_start_date),
    bx_end_date = ymd(bx_end_date),
    drug_name = case_when(
      str_detect(drug_name, "Omalizumab") ~ "Omalizumab",
      str_detect(drug_name, "99") ~ NA,
      str_detect(drug_name, "Other") ~ NA,
      .default = drug_name
    ),
    bx_name_short = structure(
      factor(
        drug_name,
        levels = c(
          "Benralizumab",
          "Dupilumab",
          "Mepolizumab",
          "Omalizumab",
          "Reslizumab",
          "Tezepelumab"
        )
      ),
      label = "Biologic name"
    ),
    bx_type = structure(
      case_when(
        bx_name_short %in% c("Mepolizumab", "Reslizumab") ~ "Anti-IL5",
        bx_name_short %in% c("Benralizumab") ~ "Anti-IL5R",
        bx_name_short == "Omalizumab" ~ "Anti-IgE",
        bx_name_short == "Dupilumab" ~ "Anti-IL4Ralpha",
        bx_name_short == "Tezepelumab" ~ "Anti-TSLP"
      ) |>
        factor(
          levels = c(
            "Anti-IgE",
            "Anti-IL4Ralpha",
            "Anti-IL5",
            "Anti-IL5R",
            "Anti-TSLP"
          ),
          labels = c(
            "Anti-IgE",
            "Anti-IL4Ralpha",
            "Anti-IL5",
            "Anti-IL5R",
            "Anti-TSLP"
          )
        ),
      label = "Biologic class"
    ),
    bx_ongoing = structure(
      factor(ongoing, levels = c("Yes", "No")),
      label = "Biologic ongoing"
    ),
    bx_switchreason = structure(
      factor(
        stop_detail,
        levels = c(
          "Suboptimal clinical response",
          "Side-effects",
          "Patient preference",
          "Other"
        )
      ),
      label = "Stop or switch reason"
    ),
    .by = patient_id
  ) |>
  select(
    patient_id,
    bx_start_date,
    bx_end_date,
    bx_ongoing,
    bx_name_short,
    bx_type,
    bx_switchreason
    # switched - IGNORE FOR NOW
  ) |>
  left_join(
    patient_table |>
      select(
        patient_id,
        data_capture_system,
        country,
        setting,
        site,
        first_visit_date,
        last_visit_date,
        num_of_visits
      ),
    by = join_by(patient_id),
    relationship = "many-to-one"
  ) |>
  relocate(
    data_capture_system,
    country,
    setting,
    site,
    contains("visit"),
    .after = patient_id
  )

# Check prior_bx_name
bx_chk1 <- bx_cln |>
  # Add and check status of biologics reported prior to UK registry entry
  right_join(
    patient |>
      filter(!is.na(uk_prior_biologic_name)) |>
      select(
        patient_id = uidpat,
        # uk_prior_biologic, # 0/1 flag doesn't align with availability of uk_prior_biologic_name
        prior_bx_name = uk_prior_biologic_name
      ) |>
      mutate(prior_bx_name = str_to_title(prior_bx_name)),
    by = join_by(patient_id),
    relationship = "many-to-one"
  ) |>
  relocate(prior_bx_name, .after = bx_name_short) |>
  mutate(
    bx_list = paste(unique(bx_name_short), collapse = ", "),
    # See "uksar prior biologic patients.xlsx" for manual review (aligns with code below)
    prior_biologic_chk = case_when(
      is.na(prior_bx_name) ~ "No prior biologic listed",
      # Exclude - Prior biologic flagged and named but not listed in episodes
      !str_detect(bx_list, prior_bx_name) ~ "prior biologic not listed",
      # Include - Prior biologic flagged and named and likely listed in episodes because same drug start date prior to first visit
      first(bx_name_short) == prior_bx_name &
        first(bx_start_date) <
          first_visit_date ~ "prior start before first visit",
      # Maybe - Prior biologic flagged and named and likely listed in episodes because same drug start date prior on or after first visit
      first(bx_name_short) == prior_bx_name &
        first(bx_start_date) >=
          first_visit_date ~ "prior start after first visit",
      # Exclude - Missing start date
      any(is.na(bx_start_date)) ~ "missing start date"
    ),
    .by = patient_id
  )

bx_cln_exclusions <- patient |>
  distinct(patient_id = uidpat) |>
  left_join(
    bx_cln |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  left_join(
    bx_chk1 |>
      distinct(patient_id, prior_biologic_chk),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date) |>
  mutate(
    # Look at overlaps (1 day probably okay)
    overlap = interval(
      dplyr::lag(bx_end_date, default = first(bx_start_date)),
      bx_start_date
    ) /
      days(1) <
      -1,
    first_biologic = bx_start_date == min(bx_start_date),
    follow_up_mos = interval(bx_start_date, last_visit_date) / months(1),
    exclusions = case_when(
      all(is.na(bx_data)) ~ "No biologic information",
      any(is.na(bx_name_short)) ~ "Missing biologic name",
      any(is.na(bx_start_date)) ~ "Missing start date",
      prior_biologic_chk == "prior biologic not listed" ~ "Non-biologic naive",
      n_distinct(bx_start_date) < n() ~ "Duplicate start dates",
      any(bx_end_date < bx_start_date) ~ "End date before start date",
      all(bx_start_date > last_visit_date) ~ "All future treatments",
      any(bx_start_date > last_visit_date) ~ "Any future treatments",
      any(is.na(bx_end_date)) ~ "Missing end date",
      any(overlap) ~ "Overlapping episodes (>1 day)",
      any(
        first_biologic & year(bx_start_date) < 2016
      ) ~ "Started biologic before 2016",
      any(
        first_biologic & follow_up_mos < 3
      ) ~ "<3 months post-initiation follow-up",
      .default = NA
    ) |>
      factor(
        levels = c(
          "No biologic information",
          "Missing biologic name",
          "Missing start date",
          "Non-biologic naive",
          "Duplicate start dates",
          "End date before start date", #Look like imputations
          "All future treatments",
          "Any future treatments",
          "Missing end date",
          "Overlapping episodes (>1 day)",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up"
        )
      ),
    exclusions_table = case_when(
      exclusions == "No biologic information" ~ "No biologic information",
      exclusions == "Non-biologic naive" ~ "Non-biologic naive",
      exclusions ==
        "Started biologic before 2016" ~ "Started biologic before 2016",
      exclusions ==
        "<3 months post-initiation follow-up" ~ "<3 months post-initiation follow-up",
      # exclusions == "<11 months post-initiation follow-up" ~ "Eligible",
      is.na(exclusions) ~ "Eligible",
      .default = "Data queries"
    ) |>
      factor(
        levels = c(
          "No biologic information",
          "Non-biologic naive",
          "Data queries",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up",
          "Eligible"
        )
      ),
    .by = patient_id
  )

# See "20261014_uk_feasibility.R"
bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions) |>
  adorn_totals() |>
  arrange(!exclusions == "Total")

bx_cln_exclusions |>
  distinct(patient_id, exclusions_table) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  arrange(!exclusions_table == "Total")

# 20261027 sent queries to Dominic P

# 10 patients with no biologic episode information have info under medications
# 20261026 DP - These patients are "intent only" visit listed bx info but none in episodes, see uk_biologic_unplaced
no_bx_info <- bx_cln_exclusions |>
  filter(exclusions == "No biologic information") |>
  distinct(patient_id) |>
  pull(patient_id)
forreview <- asthma_medications |>
  filter(uidpat %in% no_bx_info) |>
  select(uidpat, contains("anti"), -uk_anti_ige_treat) |>
  filter(if_any(contains("anti"), \(x) !is.na(x))) |>
  distinct(uidpat) |>
  pull(uidpat) |>
  noquote()
forreview

# GBR01_1028 – single episode of benra with missing end date?
# GBR01_82 – rows with missing and last row missing drug_name, end date in wrong row?
# GBR01_1049 – missing end date
# GBR01_1053 – missing start and end dates
# GBR01_413 – overlapping records for two resli spells, different starts but same end
# GBR01_556 – 9 day overlap, probably okay
# GBR01_634 – Xolair overlaps with later mepo and benra intervals, probably bad end date
# GBR03_683 – 3 month overlap mepo > dupi
# GBR04_1149 – looks like end date for Benra should be 2024-12-02 not 2025
# GBR03_1108 - started Benra on same day as Mepo?

## Analysis table ----
bx_analysis <- bx_cln_exclusions |>
  filter(is.na(exclusions)) |>
  select(-bx_data, -contains("exclusions"), -prior_biologic_chk, -overlap) |>
  # UK outcome =
  # UK continuers: 1 biologic episode, flagged as ongoing, or stop indicated more than a month after last visit
  # UK stoppers: 1 biologic episode, flagged as stopped before or within a month of last visit; >1 biologic and gap between 1st and 2nd >6 months
  # UK switchers: >1 biologic episode, <=6 months between 1 and 2nd
  mutate(
    index_date = min(bx_start_date),
    follow_up_mos = interval(index_date, last_visit_date) / months(1),
    bx_ct = n(),
    bx_nbr = row_number(),
    prev_end_date = dplyr::lag(bx_end_date),
    gap_months = case_when(
      bx_nbr > 1 ~ interval(
        prev_end_date,
        bx_start_date
      ) /
        months(1),
      .default = NA
    ),
    gt6mos_gap = bx_nbr == 2 & gap_months > 6,
    outcome = case_when(
      all(bx_ct == 1 & bx_ongoing == "Yes") ~ "Continue", # all end dates same as last visit date
      all(
        bx_ct == 1 &
          bx_ongoing == "No" &
          bx_end_date > last_visit_date %m+% months(1)
      ) ~ "Continue", # Censored not stop
      # any(bx_ct == 1 & bx_ongoing == "No" & bx_end_date <= last_visit_date) ~ "Stop",
      all(
        bx_ct == 1 &
          bx_ongoing == "No" &
          (bx_end_date < last_visit_date |
            bx_end_date <= last_visit_date %m+% months(1))
      ) ~ "Stop", # Recorded stop before or within a month of last visit
      # any(bx_ct > 1 & gt6mos_gap) ~ "Stop", # Drop this for ATS abstract, need to define gap
      any(bx_ct > 1) ~ "Switch"
    ),
    # Calculate episode length in months
    episode_months = case_when(
      outcome == "Continue" ~ interval(bx_start_date, last_visit_date) /
        months(1),
      outcome == "Stop" ~ interval(bx_start_date, bx_end_date) / months(1),
      outcome == "Switch" ~ interval(bx_start_date, bx_end_date) / months(1)
    ),
    .by = patient_id
  )

save(
  bx_cln,
  bx_cln_exclusions,
  bx_analysis,
  file = "uk_biologic_data.Rdata"
)

# Save SRD output ----
# save(
#   visits,
#   patient_table,
#   lifestyle,
#   comorbidities,
#   allergens,
#   hcru,
#   asthma_control,
#   spirometry_cln,
#   biomarkers,
#   medications,
#   bx_cln,
#   bx_cln_exclusions,
#   bx_analysis,
#   file = "srd_uk_data.Rdata"
# )
