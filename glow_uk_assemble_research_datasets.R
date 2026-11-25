# Pull UKSAR QISAR2 data and process for ISAR analysis
library(RODBC)
library(tidyverse)
library(janitor)
library(ISARanalysis)
library(kableExtra)
library(rspiro)

# Notes ----
# Draws from AB's SRD scripts but dropped Imputation flagging, etc
# 20261027 Incorporated eligibility steps from UK feasibility calculations and sent data queries to DP
# 20261027 "switched" variable in biologic_episodes is under review, not sure it's useful
# 20261124 Some small changes to data and counts

# To Do ----
# Visits - done
# Patients - done
# Lifestyle - done
# Comorbidities - focus on T2 and OCS-related for now, needs processing
# Allergens - allergens are not provided for UK
# HCRU - check timeframes
# Exacerbations - exacerbation details not provided for UK, can only use exac_num
# Asthma control - done
# Spirometry - done
# Biomarkers - done
# Medications - done
# Biologics - done

# Pull UK data ----
## Connect to ISAR SQL tables
connection <- RODBC::odbcDriverConnect(
  'driver={SQL Server};
                              server=dsimport,4071;
                              database=QISAR2;
                              trusted_connection=true'
)

uk_data <- c(
  "RR_asthma_medications",
  "RR_biomarkers",
  "RR_exacerbation_history_detail",
  "RR_patient",
  "RR_patient_visits",
  "RR_spirometry_predicted_values",
  "MT_Patient_Stats",
  "st_biologic_episodes", # use for biologic analyses
  "data_dictionary",
  "data_dictionary_raw"
) |>
  set_names(nm = \(x) {
    str_remove(x, "_detail$") |>
      str_remove("^RR_") |>
      str_remove("^st_") |>
      str_remove("_predicted_values$") |>
      str_to_lower()
  }) |>
  map(\(table) {
    sqlQuery(
      connection,
      paste0(
        'SELECT * FROM QISAR2.uk.',
        table
      ),
      as.is = T
    )
  })

save(uk_data, file = here::here("datasets", "qisar2_uk_data.RData"))

close(connection)

# Specifications ----

# Load data ----
load(here::here("datasets", "qisar2_uk_data.RData"))

#Split named list into objects
list2env(uk_data, envir = .GlobalEnv)

# Visits ----
uk_visits <- patient |>
  # n = 912
  mutate(
    patient_id = structure(uidpat, label = "Patient ID"),
    date_of_birth = NA,
    age_first_visit = as.numeric(age)
  ) |>
  select(patient_id, age_first_visit, date_of_birth) |>
  left_join(
    patient_visits |>
      distinct(
        patient_id = uidpat,
        visit_date = dov,
        uidpatvis
      ),
    by = join_by(patient_id),
    relationship = "one-to-many"
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
    visit_date = structure(dmy(visit_date), label = "Visit date")
  ) |>
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

# Checks
visits_chks <- uk_visits |>
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
rm(visits_chks)

# Patient table ----
# Patient-level, i.e. non-varying details, no separate Demographic table
uk_patients <- patient |>
  # Drop empty fields
  select(
    patient_id = uidpat,
    where(\(x) any(!is.na(x))),
    -contains("comstat"),
    -record_id
  ) |>
  # Add first and last visit
  left_join(
    uk_visits |>
      distinct(patient_id, first_visit_date, last_visit_date, num_of_visits),
    by = join_by(patient_id),
    relationship = "one-to-one"
  ) |>
  mutate(
    # date_of_birth = structure(
    #   ymd(str_replace(dob, "UNK", "15")),
    #   label = "Date of birth, imputed as 15th"
    # ),
    gender = structure(
      factor(gen, levels = c("Male", "Female")),
      label = "Gender"
    ),
    age_first_visit = structure(
      as.numeric(age),
      label = "Age at first visit, years"
    ),
    ethnicity = structure(eth, label = "Ethnicity"),
    asthma_onset = structure(
      as.numeric(asmage),
      label = "Age at asthma onset, years"
    ),
    asthma_onset_cat = structure(
      case_when(
        asthma_onset < 12 ~ "<12",
        asthma_onset <= 40 ~ "12-40",
        asthma_onset > 40 ~ ">40",
        .default = NA
      ) |>
        factor(levels = c("<12", "12-40", ">40")),
      label = "Age at asthma onset, years"
    ),
    country = structure("UK", label = "Country"),
    setting = structure("UK", label = "Setting"),
    site = structure(dagarm, label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system"),
    # ht = structure(as.numeric(ht), label = "Height, metres"),
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
    gender,
    age_first_visit,
    ethnicity,
    asthma_onset,
    asthma_onset_cat
  ) |>
  arrange(first_visit_date)

# Notes:
# AB's SRD approach appears complicated because of need to get additional
# vars: index_date, max_beos_highest, highest_bec_imputed_flag, bt
# Better I think to process biologics, BEC, bt separately and pull in later

# Lifestyle ----
# VARS:
# height - only recorded at baseline
# weight - only recorded at baseline
# BMI - recorded at each followup, can disagree with carried-forward weight
# smoking_status - only recorded at baseline
uk_lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    wt,
    bmi,
    smoke
  ) |>
  left_join(
    uk_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    patient |>
      select(
        patient_id = uidpat,
        # Height only recorded at baseline
        ht
      ) |>
      distinct(),
    by = join_by(patient_id)
  ) |>
  # Range check height and weight and calculate BMI, obesity
  ISARanalysis::clean_lifestyle() |>
  # Function cleans smoking status to plausible values over time, UK only recorded at baseline
  ISARanalysis::clean_smoking() |>
  mutate(
    smoking_status = structure(
      factor(
        smoking_status_cln,
        levels = c("Current Smoker", "Ex-Smoker", "Never Smoked")
      ),
      label = "Smoking status"
    )
  ) |>
  # count(smoke, smoking_status_cln)
  select(
    patient_id,
    uidpatvis,
    visit_id,
    visit_date,
    weight,
    height,
    bmi,
    obese,
    smoking_status
  )

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

# Focus on T2 and OCS related conditions for now
uk_comorbidities <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    allergic_rhinitis = allrhinitis,
    crs = chrrhiniosinusitis,
    nasal_polyps = npolyp,
    eczema,
    osteoporosis = osteop,
    osteoporosis_stdate = osteop_dte,
    osteop_dte_imput,
    depress,
    depress_dte,
    depress_dte_imput,
    anxiety,
    anxiety_dte,
    anxiety_dte_imput,
    # t2diab,
    # t2diab_dte,
    # t2diab_dte_imput,
    # hfailure,
    # hfailure_stdate,
    # mi,
    # mi_stdate,
    # stroke,
    # stroke_stdate,
    # embolism,
    # csd = circulatory system disease?
    # csd,
    # csd_type, #"1-5?"
    # csd_type_oth,
    # csd_stdate = csd_dte,
    # csd_dte_imput,
    # glaucoma,
    # glaucoma_stdate
    # ocular contains glaucoma and cataracts?
    # ocular,
    # ocular_type,
    # ocular_type_oth,
    # ocular_dte,
    # ocular_dte_imput,
    # osa = obsleepapnoea_yn,
    # osa_stdate = obsleepapnoea_stdt,
    # osa_stdate_imputed = obsleepapnoea_stdt_imput,
    # renalf,
    # renalf_dte,
    # renalf_dte_imput,
    # peptic_ulcer = peptic,
    # peptic_ulcer_stdate = peptic_dte,
    # peptic_dte_imput,
    # pneumonia = pneum,
    # pneumonia_dte = pneum_dte,
    # pneumonia_stdate_imputed = pneum_dte_imput,
    # gerd,
    # gerd_stdate,
    # copd,
    # egpa,
    # egpa_stdate,
    # hypereos = hyper,
    # hypereos_stdate
    # uk_atopic_disease
  ) |>
  rename_with(~ str_replace(., "imput$", "imputed"), contains("imput")) |>
  rename_with(~ str_replace(., "dte", "stdate"), contains("dte")) |>
  # select(
  #   where(\(x) any(!is.na(x))), # Keeps columns with any information
  # ) |>
  left_join(
    uk_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  mutate(
    # For allergic_rhinitis, crs, eczema Italy provides current, past, never, which we can code as "Ever", "Never" for SRD
    across(c(allergic_rhinitis, crs, eczema), \(x) {
      replace_values(x, "1" ~ "Current", "2" ~ "Past", "3" ~ "Never") |>
        factor(levels = c("Current", "Past", "Never"))
    }),
    # For nasal polyps, csd, anxiety, depression, ocular, osa, renalfailure,
    # t2diabetes, peptic ulcers, pneumonia, osteoporosis all Yes/No
    # For SRD, yes/no mapped as current/never, directly but not sure if it's built longitudinally
    #   hfailure = case_when(
    #     csd == "Yes" & csd_type == "1" ~ "Ever",
    #     csd == "Yes" ~ "Never",
    #     csd == "No" ~ "Never",
    #     .default = NA
    #   ) |>
    #     factor(levels = c("Ever", "Never")),
    #   mi = case_when(
    #     csd == "Yes" & csd_type == "2" ~ "Ever",
    #     csd == "Yes" ~ "Never",
    #     csd == "No" ~ "Never",
    #     .default = NA
    #   ) |>
    #     factor(levels = c("Ever", "Never")),
    #   stroke = case_when(
    #     csd == "Yes" & csd_type == "3" ~ "Ever",
    #     csd == "Yes" ~ "Never",
    #     csd == "No" ~ "Never",
    #     .default = NA
    #   ) |>
    #     factor(levels = c("Ever", "Never")),
    #   embolism = case_when(
    #     csd == "Yes" & csd_type == "4" ~ "Ever",
    #     csd == "Yes" ~ "Never",
    #     csd == "No" ~ "Never",
    #     .default = NA
    #   ) |>
    #     factor(levels = c("Ever", "Never")),
    #   csd_other = case_when(
    #     csd == "Yes" & csd_type == "5" ~ "Ever",
    #     csd == "Yes" ~ "Never",
    #     csd == "No" ~ "Never",
    #     .default = NA
    #   ) |>
    #     factor(levels = c("Ever", "Never"))
  )

# Allergens ----
# Allergens are not provided for UK based on DD

# HCRU ----
# exac_num - number of exacerbations in past 12 months
# emerg -	Number of A&E attendances for asthma in past 12 months
# hospadm - Number of hospital admissions for asthma in past 12 months
# invent - Number of episodes of invasive ventilation ever, recorded only at baseline
uk_hcru <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    exac_yn,
    exac_num,
    exac_timeframe,
    emerg,
    emerg_timeframe,
    hospadm,
    hospadm_timeframe,
    invent,
    invent_timeframe
  ) |>
  # UK does not provide any exacerbation details
  # left_join(
  #   exacerbation_history |> distinct(uidpatvis, exac_dte),
  #   by = "uidpatvis"
  # ) |>
  left_join(
    uk_visits |>
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
uk_exacerbations <- exacerbation_history

# Asthma control ----
# UK uses ACQ7 at baseline, ACQ6 at follow-up
uk_asthma_control <- uk_patients |>
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
            "acq_score",
            "airq_score"
          )
        ),
        contains("asthma_control"),
        -uk_gina_classification # overall baseline GINA classification
      ),
    by = join_by(patient_id),
    relationship = "one-to-many"
  ) |>
  left_join(
    uk_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = country) |>
  rename(
    gina_activity_2 = gina_activity
  ) |>
  ISARanalysis::clean_asthma_control()

# Spirometry ----
uk_spirometry <- spirometry |>
  rename(patient_id = uidpat) |>
  left_join(
    uk_visits |>
      select(uidpatvis, visit_id, visit_date, age_at_visit),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    uk_patients |>
      select(patient_id, gender),
    by = join_by(patient_id)
  ) |>
  left_join(
    uk_lifestyle |>
      select(uidpatvis, height),
    by = join_by(uidpatvis)
  ) |>
  select(
    patient_id,
    visit_id,
    uidpatvis,
    visit_date,
    gender,
    age_lf = age_at_visit,
    ht_lf = height,
    sp_date,
    sp_date_imput,
    pre_fev1 = bdfev1_pre,
    pre_fvc = bdfvc_pre,
    post_fev1 = bdfev1_post,
    post_fvc = bdfvc_post
  ) |>
  arrange(patient_id, visit_date) |>
  ISARanalysis::clean_spirometry() |>
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

# Still a few with FVC < FEV1
uk_spirometry |>
  filter(if_any(contains("implausible")) | if_any(contains("rangefail"))) |>
  pull(patient_id) |>
  noquote()

# Biomarkers ----
# UK doesn't provide test-specific dates only visit measures
uk_biomarkers <- biomarkers |>
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
    ige_result = ige_cnt,
    feno_result = feno_rlt,
    -uk_feno_suppression_cat
  ) |>
  left_join(
    uk_visits |>
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
    # Use visit date if test date is missing
    date = visit_date,
    result = case_when(
      # Convert BEC to cells/mcL
      test == "beos" ~ as.numeric(result) * 1000,
      # IgE kU/L is same as IU/mL
      .default = as.numeric(result)
    ),
    unit = case_when(
      test == "beos" ~ "cells/μl",
      test == "ige" ~ "IU/mL",
      test == "feno" ~ "ppb"
    )
  ) |>
  # Ensure units converted before apply range checks
  clean_biomarkers(
    drop_duplicate = T,
    drop_rangefail = T,
    reduce_function = "max",
    flag_only = F,
    keep_flags = F
  ) |>
  mutate(
    result = case_when(
      # Truncate IgE and FeNO counts (same as NJH)
      test == "ige" & result <= 2 ~ 2,
      test == "ige" & result >= 5000 ~ 5000,
      test == "feno" & result <= 5 ~ 5,
      test == "feno" & result >= 300 ~ 300,
      .default = result
    )
  )

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
uk_medications <- asthma_medications |>
  # Keep non-empty columns
  select(
    patient_id = uidpat,
    uidpatvis,
    where(\(x) any(!is.na(x))),
    -uiddb,
    -record_id,
    -dagarm,
    -contains("redcap"),
    -contains("imput")
  ) |>
  # Bare drug columns (no underscore) -> "<drug>_name"
  rename_with(
    ~ paste0(.x, "_name"),
    .cols = !all_of(c(
      "patient_id",
      "uidpatvis"
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
  ) |>
  pivot_longer(
    cols = !all_of(c(
      "patient_id",
      "uidpatvis"
    )),
    names_to = c("therapy", ".value"),
    names_sep = "_"
  ) |>
  # All (name, ldose, stdate, ongoing, endate) missing where yn == NA & yn == "No"
  # filter(yn == "Yes") |>
  # Range checks
  filter(!is.na(name)) |>
  mutate(
    therapy = case_when(
      therapy %in% c("antiige", "antiil4", "antiil5") ~ "Biologic",
      therapy == "mos" ~ "LTOCS",
      therapy == "ics" ~ "ICS only",
      therapy == "mabx" ~ "Macrolides",
      .default = str_to_upper(therapy)
    ),
    across(contains("date"), dmy)
  ) |>
  left_join(
    uk_visits |> select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  distinct() |>
  filter_out(therapy == "Biologic") |>
  arrange(patient_id, visit_date, therapy)

# Earliest medication start date is always a year before baseline visit
# Ongoing flag is derived (so not from UK), not sure how helpful this is.
# chk <- medications |> mutate(date_chk = min(stdate, na.rm = T) == (min(visit_date) %m-% months(12)), .by = patient_id, drug)

# Biologics ----
# Use ST biologic episodes file as provided by UKSAR rather than patient_visits
# Note - some fields in patient visits appear populated with bx info but not all, incomplete?
uk_biologics <- biologic_episodes |>
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
    bx_start_date = dmy(bx_start_date),
    bx_end_date = dmy(bx_end_date),
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
    uk_patients |>
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
bx_chk1 <- uk_biologics |>
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

# Add prior biologic conclusion
uk_biologics <- uk_biologics |>
  left_join(
    bx_chk1 |>
      distinct(patient_id, prior_biologic_chk),
    by = join_by(patient_id)
  )

# Save datasets ----
save(
  uk_visits,
  uk_patients,
  uk_lifestyle,
  uk_comorbidities,
  # uk_allergens,
  uk_hcru,
  uk_exacerbations,
  uk_asthma_control,
  uk_spirometry,
  uk_biomarkers,
  uk_medications,
  uk_biologics,
  file = here::here("datasets", "qisar2_uk_srd.Rdata")
)
