# Pull Portugal QISAR2 data and process for GLOW analysis
library(RODBC)
library(tidyverse)
library(janitor)
library(ISARanalysis)
library(kableExtra)
library(rspiro)

# Notes ----
# Pending queries around HCRU (exac) timeframes
# PRT_57 has problematic height and weight and other records - exclude for now
# Some visits with duplicate dates
# PRT_167 has LABA Indacaterol recorded as LAMA
# Looks to be an issue with LABA processing as a lot of "No" records with named drugs

# To Do ----
# Visits - done
# Patients - done
# Lifestyle - done
# Comorbidities - focus on T2 and OCS-related for now, needs processing
# Allergens - needs processing to get final variable
# HCRU - done for now
# Exacerbations - done for now
# Asthma control - done
# Spirometry - done
# Biomarkers - done
# Medications - done, needs further cleaning of "other" values and ICS LABA LAMA harmonisation
# Biologics - done, with queries to send, 100+ patients with missing start date

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

# Load data ----
load(here::here("datasets", "qisar2_prt_data.RData"))

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
prt_visits <- patient |>
  # n = 377
  mutate(
    patient_id = structure(uidpat, label = "Patient ID"),
    date_of_birth = ymd(dob),
    # Provided age is age at baseline according to DD
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
visits_chks <- prt_visits |>
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
  pull(patient_id) |>
  noquote() |>
  kable() |>
  kable_styling()

# Find patients where baseline is not first visit
visits_chks |>
  filter(str_detect(visit_id, "baseline") & date_order != 0) |>
  distinct(patient_id) |>
  pull(patient_id) |>
  noquote() |>
  kable() |>
  kable_styling()

# Patient table ----
# Patient-level, i.e. non-varying details, no separate Demographic table
prt_patients <- patient |>
  # Drop empty fields
  select(
    patient_id = uidpat,
    where(\(x) any(!is.na(x))),
    -contains("comstat"),
    -record_id
  ) |>
  # Add first and last visit
  left_join(
    prt_visits |>
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
    country = structure("Portugal", label = "Country"),
    setting = structure("Portugal", label = "Setting"),
    site = structure(dagarm, label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system")
  ) |>
  select(
    patient_id,
    country,
    setting,
    site,
    first_visit_date,
    last_visit_date,
    num_of_visits,
    data_capture_system,
    gender,
    date_of_birth,
    ethnicity,
    asthma_onset,
    asthma_onset_cat
  ) |>
  arrange(first_visit_date)

# Notes:
# Aaron's SRD approach appears complicated because of need to get additional
# vars: index_date, max_beos_highest, highest_bec_imputed_flag, bt
# Better I think to process biologics, BEC, bt separately and pull in later

# Lifestyle ----
# Visit level information
# VARS:
# weight - Weight is recorded at each visit and is used as recorded.
# height - Height is the visit's own, else the latest earlier visit's.
# BMI - BMI is recalculated from them, as the R did (QISAR's values are the recalculation), falling back to the recorded IMC; prt_bmi_recorded keeps the recorded figure.
# BSA - BSA is Du Bois, as the UK.
# smoking_status - carried forward from the latest visit that recorded it
# Errors:
# Only one patient appears to have varying height, PRT_57, and I think this is a mistake
# PRT_134 had swapped height/weight, used 70kg for weight as Dom did and use prt_height_cm = 1.62 instead of 0.70

prt_lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    # Weight non-varying except six patients
    wt,
    # Height non-varying except one patient
    ht = prt_height_cm,
    bmi,
    smoke
  ) |>
  mutate(ht = as.numeric(ht) / 100) |>
  left_join(
    prt_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  # left_join(
  #   patient |>
  #     select(
  #       patient_id = uidpat,
  #       # Height only recorded at baseline I think, non-varying except one patient
  #       ht
  #     ) |>
  #     distinct(),
  #   by = join_by(patient_id)
  # ) |>
  # Range check height and weight and calculate BMI, obesity
  ISARanalysis::clean_lifestyle() |>
  # Function cleans smoking status to plausible values over time
  ISARanalysis::clean_smoking() |>
  mutate(
    smoking_status = structure(
      factor(
        smoking_status_cln,
        levels = c("Current smoker", "Ex-smoker", "Never smoked")
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

# Check smoking status
# lifestyle |>
#   filter(any(smoking_status!=smoking_status_cln),.by= patient_id) |>
#   select(patient_id, uidpatvis, visit_date, contains("smoking_")) |>
#   # filter(any(n_distinct(visit_date) < n()),.by= patient_id) |>
#   view()

# Comorbidity ----
# Focus on T2 and OCS related conditions for now
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

prt_comorbidities <- patient_visits |>
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
  ) |>
  rename_with(~ str_replace(., "imput$", "imputed"), contains("imput")) |>
  rename_with(~ str_replace(., "dte", "stdate"), contains("dte")) |>
  # select(
  #   where(\(x) any(!is.na(x))), # Keeps columns with any information
  # ) |>
  left_join(
    prt_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  mutate(
    # For allergic_rhinitis, crs, eczema Italy provides current, past, never, which we can code as "Ever", "Never" for SRD
    # across(c(allergic_rhinitis, crs, eczema), \(x) {
    #   replace_values(x, "1" ~ "Current", "2" ~ "Past", "3" ~ "Never") |>
    #     factor(levels = c("Current", "Past", "Never"))
    # }),
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
# Use AB's SRD approach for processing for now, though I'm not sure the QISAR structure makes senses
# 20261111 AS/AT - sat_posallergen___0 and spt_posallergen___0 correspond to "No data" - all NA for Italy at least
# Not sure if "No data" and "No positive allergen" make sense.
# AB's approach probably most sensible but unclear whether to use "positive" or
# "clinically relevant" - this variable was dropped from SRD and is only available for serum tests.
prt_allergens <- allergens |>
  select(
    everything(),
    patient_id = uidpat,
    -uiddb,
    -dagarm
  ) |>
  # Do we want to check the "no positive allergen" for SPT for SAT and vice-versa, excluded here but in Aaron's SRD script
  mutate(
    # Revised Aaron's code to update and make consistent - need to work on "other"?
    sat_none = sat_posallergen___8 == "No Positive Allergens",
    spt_none = spt_posallergen___12 == "No Positive Allergens",
    # Recode columns as Yes/No
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
      },
      .names = "{.col}_yn"
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
      },
      .names = "{.col}_yn"
    )
  ) |>
  select(-contains("_none"))
# rename(
#   # Using _pos instead of _yn as in SRD because this is the flag for a positive test
#   sat_pos_dust_result = sat_pos_dust,
#   sat_pos_dust = sat_pos_dust_yn,
#   sat_pos_grass_result = sat_pos_grass,
#   sat_pos_dust = sat_pos_grass_yn,
#
#   sat_cat_pos = sat_posallergen___3,
#   sat_mould_pos = sat_posallergen___4,
#   sat_dog_pos = sat_posallergen___5,
#   sat_aspergillu_pos = sat_posallergen___6,
#   sat_other_pos = sat_posallergen___7,
#   sat_other_specify_1 = sat_pos_othspecify,
#   sat_other_result_1 = sat_pos_oth,
#   # sat_other_specify_2 = sat_pos_othspecify_2,
#   # sat_other_result_2 = sat_pos_oth_2,
#   # sat_other_specify_3 = sat_pos_othspecify_3,
#   # sat_other_result_3 = sat_pos_oth_3,
#   # sat_other_specify_4 = sat_pos_othspecify_4,
#   # sat_other_result_4 = sat_pos_oth_4
#   spt_dust_yn = spt_posallergen___7,
#   spt_grass_yn = spt_posallergen___1,
#   spt_cat_yn = spt_posallergen___9,
#   spt_mould_yn = spt_posallergen___5,
#   spt_dog_yn = spt_posallergen___10,
#   spt_aspergillus_yn = spt_posallergen___4,
#   spt_tree_yn = spt_posallergen___2,
# ) |>
#   rename_with(~ str_replace(., "_dte", "_date"), contains("dte")) |>
#   left_join(
#     prt_visits |>
#       select(uidpatvis, visit_id, visit_date),
#     by = join_by(uidpatvis)
#   ) |>
#   relocate(contains("visit_"), .after = patient_id) |>
#   # Add steps to create allergens variable used for analysis, modifying SRD - as in "glow_baseline.R"
#   mutate(
#     any_positive = if_any(matches("^(sat|spt)_.+_yn$"), ~ .x %in% "Yes"),
#     any_tested = if_any(matches("^(sat|spt)_.+_yn$"), ~ !is.na(.x)),
#     allergy = case_when(
#       any_positive ~ "Yes",
#       any_tested & !any_positive ~ "No",
#       .default = NA
#     )
#   )
# filter_out(is.na(allergy)) |>
# GS previously used if any positive allergen results then the patient is positive
# summarise(
#   allergy = if_else(any(allergy == "Yes"), "Yes", "No"),
#   .by = patient_id
# ) |>
#   mutate(
#     allergy = structure(
#       factor(allergy, levels = c("Yes", "No")),
#       label = "Any Positive Allergen Test"
#     )
#   )

# HCRU ----
# exac_num - number of exacerbations in past 12 months, as in UK
# emerg -	Number of A&E attendances for asthma in past 12 months
# hospadm - Number of hospital admissions for asthma in past 12 months
# invent - Number of episodes of invasive ventilation ever, recorded only at baseline
prt_hcru <- patient_visits |>
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
  left_join(
    prt_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  mutate(
    exac_num = structure(
      as.numeric(exac_num),
      label = "Number of exacerbation in past 12 months"
    ),
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
    uidpatvis,
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
# Use number for now
prt_exacerbations <- exacerbation_history

# Asthma control ----
# Portugal provides ACT and GINA, not sure how asthma_control derived in PRT data
# Recalculate to prioritise GINA > ACT > ACQ > AIRQ
prt_asthma_control <- prt_patients |>
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
    prt_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = country) |>
  rename(
    gina_activity_2 = gina_activity,
    prt_asthma_control = asthma_control
  ) |>
  ISARanalysis::clean_asthma_control()

# Spirometry ----
# Only one spirometry record provided per visit (i.e. no duplicated uidpatvis)
# sp_date provided for some, use this or visit date if missing
prt_spirometry <- spirometry |>
  rename(patient_id = uidpat) |>
  left_join(
    prt_visits |>
      select(uidpatvis, visit_id, visit_date, age_at_visit),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    prt_patients |>
      select(patient_id, gender),
    by = join_by(patient_id)
  ) |>
  left_join(
    prt_lifestyle |>
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

# None with FVC < FEV1
prt_spirometry |>
  filter(if_any(contains("implausible")) | if_any(contains("rangefail"))) |>
  pull(patient_id) |>
  noquote()

# Biomarkers ----
# Harmonising data sources:
# IgE all results in kU/L
# BEOS convert to cells/uL (multiply K/uL or 10*9/L by 1000)
# FeNO all results in ppb
# PRT doesn't provide test-specific dates for IgE, use visit date if test date missing
prt_biomarkers <- biomarkers |>
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
    prt_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  mutate(
    ige_date = as.character(NA)
  ) |>
  pivot_longer(
    cols = !all_of(c("patient_id", "uidpatvis", "visit_id", "visit_date")),
    names_to = c("test", ".value"),
    names_sep = "_"
  ) |>
  filter_out(is.na(result)) |>
  mutate(
    # Use visit date if test date is missing
    date = if_else(!is.na(date), ymd(date), visit_date), # Almost all test dates, where available for FeNO, are same as visit date, except ~20 earlier ones
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

# Medications ----
prt_medications <- asthma_medications |>
  # Keep non-empty columns
  select(
    patient_id = uidpat,
    uidpatvis,
    where(\(x) any(!is.na(x))),
    -uiddb,
    -record_id,
    -dagarm,
    -contains("redcap"),
    -contains("imput"), # 1 if day imputed as 15th or day/month as 15/6
    -icslabalama_explicit
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
  # Need to clean "oth" columns: mos, laba, icslaba, icslabalama, lama
  mutate(
    mos_name = case_when(
      mos_oth == "prednisolona" ~ "Prednisolone",
      .default = mos_name
    ),
    # PRT_178, PRT_286 using mometasone as ICS or triple therapy?
    # ics_yn = case_when(
    #   uidpatvis == "" & str_to_snake(icslaba_oth) == "mometasona" ~ "Yes",
    #   .default = ics_yn
    # ),
    # ics_name = case_when(
    #   uidpatvis == "" & str_to_snake(icslaba_oth) == "mometasona" ~ "Yes",
    #
    # ),
    icslaba_name = case_when(
      str_to_snake(icslaba_oth) %in%
        c(
          "indacaterol_e_furoato_de_mometasona",
          "indacaterol_e_mometasona",
          "indacaterol_furoato_de_mometasona",
          "mometasona_e_indacaterol",
          "mometasona_indacaterol"
        ) ~ "Mometasone furoate + Indacaterol",
      .default = icslaba_name
    ),
    # PRT_167 has LABA Indacaterol recorded as LAMA
    laba_yn = case_when(
      uidpatvis == "PRT_167_f2arm1" & lama_oth == "indacaterol" ~ "Yes",
      .default = laba_yn
    ),
    laba_name = case_when(
      laba_oth == "Indacaterol" ~ "Indacaterol",
      laba_oth == "Vilanterol" ~ "Vilanterol",
      uidpatvis == "PRT_167_f2arm1" & lama_oth == "indacaterol" ~ "Indacaterol",
      .default = lama_name
    ),
    laba_ongoing = case_when(
      uidpatvis == "PRT_167_f2arm1" & lama_oth == "indacaterol" ~ 1,
      .default = laba_ongoing
    ),
    lama_yn = case_when(
      uidpatvis == "PRT_167_f2arm1" & lama_oth == "indacaterol" ~ "No",
      .default = lama_yn
    ),
    lama_name = case_when(
      lama_oth == "umeclidinio" ~ "Umeclidinium",
      lama_oth == "indacaterol" ~ NA_character_,
      .default = lama_name
    ),
    lama_ongoing = case_when(
      uidpatvis == "PRT_167_f2arm1" & lama_oth == "indacaterol" ~ NA,
      .default = lama_ongoing
    )
  ) |>
  select(
    -contains("_oth"),
    -contains("prt_"),
    -contains("switch")
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
  # filter(is.na(yn)) |> select(where(\(x) any(!is.na(x)))) |>
  # filter(yn == "No") |> select(where(\(x) any(!is.na(x)))) |>
  # Looks to be an issue with LABA processing as a lot of "No" records with named drugs
  filter(yn == "Yes") |>
  mutate(
    therapy = case_when(
      therapy %in% c("antiige", "antiil4", "antiil5") ~ "Biologic",
      therapy == "mos" ~ "LTOCS",
      therapy == "ics" ~ "ICS only",
      therapy == "mabx" ~ "Macrolides",
      .default = str_to_upper(therapy)
    ),
    across(contains("date"), ymd)
  ) |>
  left_join(
    prt_visits |> select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
  filter_out(therapy == "Biologic") |>
  arrange(patient_id, visit_date, therapy)

# Biologics ----
# Use ST biologic episodes file as provided by PRT rather than patient_visits
prt_biologics <- st_biologic_episodes |>
  select(
    patient_id = uidpat,
    where(\(x) any(!is.na(x))),
    # # Raw start/stop dates checked for parsing errors - none currently
    bx_start_date = start_date,
    bx_end_date = stop_date,
  ) |>
  # Flag and drop rows with duplicated start dates and drug names, keep record with end date
  # mutate(
  #   dup_start_name = n() > 1,
  #   .by = c(patient_id, start_date, drug_name)
  # ) |>
  mutate(
    recorded_visit_date = ymd(recorded_visit_date)
  ) |>
  # Some duplicated visit dates (many-to-many), but looks to work to add visit ID
  left_join(
    prt_visits |>
      select(
        patient_id,
        uidpatvis,
        visit_id,
        visit_date
      ),
    by = join_by(patient_id, recorded_visit_date == visit_date)
  ) |>
  mutate(
    bx_start_date = case_when(
      start_date_precision == "Y" ~ NA,
      .default = ymd(bx_start_date)
    ),
    bx_end_date = case_when(
      bx_end_date == "Y" ~ NA,
      .default = ymd(bx_end_date)
    ),
    drug_name = case_when(
      str_detect(drug_name, "Omalizumab") ~ "Omalizumab",
      str_detect(drug_name, "99") ~ NA,
      str_detect(drug_name, "Other") ~ NA,
      .default = str_to_title(drug_name)
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
      factor(ongoing, levels = c("1", "0"), labels = c("Yes", "No")),
      label = "Biologic ongoing"
    ),
    bx_switchreason = structure(
      factor(
        switch_reason,
        levels = c(
          "Lack of efficacy",
          "Side effects",
          "Access restriction",
          "Patient preference",
          "Other"
        ),
        labels = c(
          "Suboptimal clinical response",
          "Side-effects",
          "Access restriction",
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
    uidpatvis,
    visit_id,
    visit_date = recorded_visit_date,
    bx_start_date,
    bx_end_date,
    bx_ongoing,
    bx_name_short,
    bx_type,
    bx_switchreason
  ) |>
  # A lot of duplicated records
  distinct() |>
  left_join(
    prt_patients |>
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
    uidpatvis,
    .after = patient_id
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date, bx_switchreason) |>
  # mutate(
  #   dup = n(),
  #   .by = c(patient_id, bx_start_date, bx_name_short)
  # ) |>
  # Exclude some duplicated rows
  filter_out(
    when_any(
      patient_id == "PRT_103" & is.na(bx_switchreason),
      patient_id == "PRT_106" & is.na(bx_switchreason),
      patient_id == "PRT_107" & is.na(bx_end_date),
      patient_id == "PRT_111" & is.na(bx_end_date),
      patient_id == "PRT_117" & is.na(bx_end_date),
      patient_id == "PRT_122" & is.na(bx_switchreason),
      patient_id == "PRT_139" & is.na(bx_end_date),
      patient_id == "PRT_140" & is.na(bx_switchreason),
      patient_id == "PRT_147" & is.na(bx_end_date),
      patient_id == "PRT_164" & bx_switchreason == "Patient preference", #Prioritise efficacy
      patient_id == "PRT_79" & is.na(bx_switchreason),
      patient_id == "PRT_89" &
        is.na(bx_end_date) &
        bx_name_short == "Omalizumab",
      patient_id == "PRT_89" &
        is.na(bx_switchreason) &
        bx_name_short == "Mepolizumab",
      patient_id == "PRT_89" &
        is.na(bx_end_date) &
        bx_name_short == "Benralizumab",
      patient_id == "PRT_90" & bx_switchreason == "Patient preference", #Prioritise efficacy
      patient_id == "PRT_159" & is.na(bx_switchreason),
      patient_id == "PRT_159" & is.na(bx_end_date),
      patient_id == "PRT_184" & is.na(bx_end_date) # Two rows with same drug with different end dates remain
    )
  )

# Save datasets ----
save(
  prt_visits,
  prt_patients,
  prt_lifestyle,
  prt_comorbidities,
  prt_allergens,
  prt_hcru,
  prt_exacerbations,
  prt_asthma_control,
  prt_spirometry,
  prt_biomarkers,
  prt_medications,
  prt_biologics,
  file = here::here("datasets", "qisar2_prt_srd.Rdata")
)
