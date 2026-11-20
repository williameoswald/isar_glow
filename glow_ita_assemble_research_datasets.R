# Pull Italy QISAR2 data and process for ISAR analysis
library(RODBC)
library(tidyverse)
library(janitor)
library(ISARanalysis)
library(kableExtra)
library(rspiro)

# Notes ----

# To Do ----
# Visits - done
# Patients - done
# Lifestyle - done
# Comorbidities - focus on T2 and OCS-related for now, needs processing
# Allergens - needs processing to get final variable
# HCRU - check timeframes
# Exacerbations - very few dates, compare numbers against reported counts
# Asthma control - done
# Spirometry - done
# Biomarkers - done, dates across visits look funny and some beos values may already have been converted
# Medications - done
# Biologics - done, queries sent, exclude patients for now

# Pull Italy data ----
## Connect to ISAR SQL tables
connection <- RODBC::odbcDriverConnect(
  'driver={SQL Server};
                              server=dsimport,4071;
                              database=QISAR2;
                              trusted_connection=true'
)

ita_data <- c(
  "RR_allergens",
  "RR_asthma_control",
  "RR_asthma_medications",
  "RR_biomarkers",
  "RR_comorbidities",
  "RR_exacerbation_history_detail",
  "RR_patient",
  "RR_patient_visits",
  "RR_spirometry_predicted_values",
  "ST_asthma_medications",
  "ST_exits_reason",
  "data_dictionary",
  "data_dictionary_raw"
) |>
  set_names(nm = \(x) {
    str_remove(x, "_detail$") |>
      str_remove("^RR_") |>
      str_replace("ST_asthma_", "ST_") |>
      str_remove("^ST_") |>
      str_remove("_predicted_values$") |>
      str_to_lower()
  }) |>
  map(\(table) {
    sqlQuery(
      connection,
      paste0(
        'SELECT * FROM QISAR2.ita.',
        table
      ),
      as.is = T
    )
  })

save(ita_data, file = "qisar2_ita_data.RData")

close(connection)

# Specifications ----

# Load data ----
load("qisar2_ita_data.RData")

#Split named list into objects
list2env(ita_data, envir = .GlobalEnv)

# Visits ----
ita_visits <- patient |>
  # n = 4138
  # Italy provides DOB so use this instead of age (not sure where this comes from), with day imputed as 15th
  mutate(
    patient_id = structure(uidpat, label = "Patient ID"),
    date_of_birth = ymd(str_replace(dob, "UNK", "15")),
    age_first_visit = as.numeric(NA)
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

# Check visit ID versus date order
visits_chks <- ita_visits |>
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

# Find patients with visits out of order
visits_chks |>
  filter(any(visit_order != date_order), .by = patient_id) |>
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
rm(visits_chks)

# Patients ----
# Patient-level, i.e. non-varying details, no separate Demographic table
ita_patients <- patient |>
  # Drop empty fields
  select(
    patient_id = uidpat,
    where(\(x) any(!is.na(x))),
    -contains("comstat"),
    -record_id
  ) |>
  # Add first and last visit
  left_join(
    ita_visits |>
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
      ymd(str_replace(dob, "UNK", "15")),
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
    country = structure("Italy", label = "Country"),
    setting = structure("Italy", label = "Setting"),
    site = structure(dagarm, label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system"),
    # Height recorded at each visit for Italy, ita_height_at_visit
    # ht = structure(as.numeric(ht), label = "Height, metres"),
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
# AB's SRD approach appears complicated because of need to get additional
# vars: index_date, max_beos_highest, highest_bec_imputed_flag, bt
# Better I think to process biologics, BEC, bt separately and pull in later

# Lifestyle ----
# Visit level information
# VARS:
# weight - Weight is recorded at each visit and is used as recorded.
# height - Height is recorded at each visit and is used as recorded.
# BMI - Recalculate
# smoking_status - carried forward from the latest visit that recorded it, with logic rules

# Looks like smoking only recorded once, at baseline presumably
ita_lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    wt,
    ht = ita_height_at_visit, # Italian height variable at visit
    smoke # Varies by visit, a lot of Never smokers convert to Ex-smokers
  ) |>
  left_join(
    ita_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  # Fill height (not weight) missings across visits, last measure forward and also back
  fill(ht, .direction = "downup", .by = patient_id) |>
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

# Comorbidities ----
# Focus on T2 and OCS related conditions for now
ita_comorbidities <- comorbidities |>
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
    ita_visits |>
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
# Use AB's SRD approach for processing for now, though I'm not sure the QISAR structure makes senses
# 20261111 AS/AT - sat_posallergen___0 and spt_posallergen___0 correspond to "No data" - all NA for Italy at least
# Not sure if "No data" and "No positive allergen" make sense.
# AB's approach probably most sensible but unclear whether to use "positive" or
# "clinically relevant" - this variable was dropped from SRD and is only available for serum tests.
ita_allergens <- allergens |>
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
rename(
  # Using _pos instead of _yn as in SRD because this is the flag for a positive test
  sat_pos_dust_result = sat_pos_dust,
  sat_pos_dust = sat_pos_dust_yn,
  sat_pos_grass_result = sat_pos_grass,
  sat_pos_dust = sat_pos_grass_yn,

  sat_cat_pos = sat_posallergen___3,
  sat_mould_pos = sat_posallergen___4,
  sat_dog_pos = sat_posallergen___5,
  sat_aspergillu_pos = sat_posallergen___6,
  sat_other_pos = sat_posallergen___7,
  sat_other_specify_1 = sat_pos_othspecify,
  sat_other_result_1 = sat_pos_oth,
  # sat_other_specify_2 = sat_pos_othspecify_2,
  # sat_other_result_2 = sat_pos_oth_2,
  # sat_other_specify_3 = sat_pos_othspecify_3,
  # sat_other_result_3 = sat_pos_oth_3,
  # sat_other_specify_4 = sat_pos_othspecify_4,
  # sat_other_result_4 = sat_pos_oth_4
  spt_dust_yn = spt_posallergen___7,
  spt_grass_yn = spt_posallergen___1,
  spt_cat_yn = spt_posallergen___9,
  spt_mould_yn = spt_posallergen___5,
  spt_dog_yn = spt_posallergen___10,
  spt_aspergillus_yn = spt_posallergen___4,
  spt_tree_yn = spt_posallergen___2,
) |>
  rename_with(~ str_replace(., "_dte", "_date"), contains("dte")) |>
  left_join(
    ita_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(contains("visit_"), .after = patient_id) |>
  # Add steps to create allergens variable used for analysis, modifying SRD - as in "glow_baseline.R"
  mutate(
    any_positive = if_any(matches("^(sat|spt)_.+_yn$"), ~ .x %in% "Yes"),
    any_tested = if_any(matches("^(sat|spt)_.+_yn$"), ~ !is.na(.x)),
    allergy = case_when(
      any_positive ~ "Yes",
      any_tested & !any_positive ~ "No",
      .default = NA
    )
  )
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
# Need timeframes
ita_hcru <- patient_visits |>
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
    ita_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = uidpatvis) |>
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
ita_exacerbations <- exacerbation_history

# Asthma control ----
# Italy doesn't include GINA question 1 on activity, based on ACT > ACQ
ita_asthma_control <- ita_patients |>
  select(patient_id, country) |>
  left_join(
    asthma_control |>
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
        contains("asthma_control")
      ),
    by = join_by(patient_id),
    relationship = "one-to-many"
  ) |>
  left_join(
    ita_visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = country) |>
  ISARanalysis::clean_asthma_control()

# Spirometry ----
ita_spirometry <- spirometry |>
  rename(patient_id = uidpat) |>
  # mutate(
  #   p = str_match(uidpat, "^([^_]+)_([^_]+)_(\\d+)$"),
  #   patient_id = str_c(
  #     p[, 2],
  #     p[, 3],
  #     "_",
  #     str_pad(p[, 4], width = 5, pad = "0")
  #   ),
  #   q = str_match(uidpatvis, "^([^_]+)_([^_]+)_(\\d+)_(.*)$"),
  #   uidpatvis = str_c(
  #     q[, 2],
  #     q[, 3],
  #     "_",
  #     str_pad(q[, 4], width = 5, pad = "0"),
  #     "_",
  #     q[, 5]
  #   )
  # ) |>
  # select(-p, -q, -uidpat) |>
  left_join(
    ita_visits |>
      select(uidpatvis, visit_id, visit_date, age_at_visit),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    ita_patients |>
      select(patient_id, gender),
    by = join_by(patient_id)
  ) |>
  left_join(
    ita_lifestyle |>
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
ita_spirometry |>
  filter(if_any(contains("implausible")) | if_any(contains("rangefail"))) |>
  pull(patient_id) |>
  noquote()

# Biomarkers ----
# Question - impute missing dates with visit date, then check for duplicates?
# Italy doesn't provide IgE dates
# Dates across visits look funny and some beos values may already have been converted
ita_biomarkers <- biomarkers |>
  select(
    patient_id = uidpat,
    everything(),
    -uiddb,
    -dagarm,
    -contains("_yn"),
    -contains("ita_"),
    -contains("_imput"),
    -contains("_timeframe"),
    -import_log_id
  ) |>
  rename_with(~ str_replace(., "_cnt", "_result"), everything()) |>
  rename_with(~ str_replace(., "_dte", "_date"), everything()) |>
  mutate(
    feno_unit = "ppb"
  ) |>
  left_join(
    ita_visits |>
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
    date = if_else(!is.na(date), dmy(date), visit_date), # Almost all test dates, where available for FeNO, are same as visit date, except ~20 earlier ones
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
  )

# Medications ----
# GS - therapy data provided by Italy is not linked to visits. The file lists
# start dates (and potentially en dates if treatment has been discontinued) for
# all asthma-related drugs. It should be a lot easier to use this data without
# taking visits into consideration. We need to be careful not assuming « No »
# when the drug is not listed: eg, for LTOCS, I am assuming missing when only
# bx is listed there (it’s not plausible that a patient would initiate bx with
# no other asthma-related drug). I am also looking at the time window covered by
# the data, ie data can be entered in the CRF but the therapy data not entered.

# Use supplement table rather than qisar-style table for querying and creating biologic table below
ita_medications <- medications |>
  select(
    -episode_id,
    -medication_class, # Eventually same as derived as "therapy" below but latter more specific
    -contains("imput"),
    -contains("import"),
    start_date = therapy_start,
    end_date = therapy_end,
  ) |>
  mutate(
    # One UNK day in date to impute as 15th
    start_date = case_when(
      str_detect(start_date, "^UNK-") ~ str_replace(start_date, "^UNK-", "15-"),
      .default = start_date
    )
  ) |>
  # Bare therapy columns (no underscore) -> "<therapy>_name"
  rename_with(
    ~ str_remove_all(., "_") |>
      str_replace("labeldose", "_labeldose") |>
      str_replace("freq", "_freq"),
    matches("mos|ics|anti|rhin|bron")
  ) |>
  rename_with(
    ~ paste0(.x, "_name"),
    .cols = !all_of(c("patient_id", "start_date", "end_date")) &
      !matches("dose|freq|biologic")
  ) |>
  pivot_longer(
    cols = !all_of(c(
      "patient_id",
      "start_date",
      "end_date",
      "biologic_switch_reason",
      "biologic_switch_adverse_event_grade"
    )),
    names_to = c("therapy", ".value"),
    names_sep = "_"
  ) |>
  filter(!is.na(name)) |>
  mutate(
    therapy = case_when(
      therapy %in% c("antiige", "antiil4", "antiil5", "antitslp") ~ "Biologic",
      therapy == "mos" ~ "LTOCS",
      therapy == "ics" ~ "ICS only",
      therapy == "icslaba" ~ "ICS LABA",
      therapy == "icslabalama" ~ "ICS LABA LAMA",
      therapy == "macrolideantibiotic" ~ "Macrolides",
      therapy == "rhinitistherapy" ~ "Rhinitis",
      therapy == "bronchialthermoplasty" ~ "Thermoplasty",
      .default = str_to_upper(therapy)
    ),
    across(contains("date"), dmy)
  ) |>
  relocate(therapy, name, labeldose, freq, .after = patient_id) |>
  arrange(patient_id, start_date, therapy)

# Biologics ----
# No rows with multiple biologics listed on same row
# medications |> filter(rowSums(!is.na(pick(starts_with("anti")))) > 1)

# Use supplement therapy file as provided rather than asthma_medications table
ita_biologics <- ita_medications |>
  filter(therapy == "Biologic") |>
  select(
    patient_id,
    bx_start_date = start_date,
    bx_end_date = end_date,
    bx_name = name,
    bx_switchreason = biologic_switch_reason,
    bx_switch_adverse_event_grade = biologic_switch_adverse_event_grade
  ) |>
  mutate(
    bx_name_short = structure(
      factor(
        str_trim(str_remove(bx_name, " \\(.*\\)$")),
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
      if_else(is.na(bx_end_date), "Yes", "No") |>
        factor(levels = c("Yes", "No")),
      label = "Biologic ongoing"
    ),
    # reason_no_clinical_efficacy = str_detect(switchreason, "No clinical efficacy"),
    # reason_patients_choice = str_detect(switchreason, "Patient's choice"),
    # reason_side_effects = str_detect(switchreason, "Side effects"),
    # reason_access_restriction = str_detect(switchreason, "Restrictions on access to biological drugs"),
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
  ) |>
  left_join(
    ita_patients |>
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
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date)

# Save datasets ----
save(
  ita_visits,
  ita_patients,
  ita_lifestyle,
  ita_comorbidities,
  ita_allergens,
  ita_hcru,
  ita_exacerbations,
  ita_asthma_control,
  ita_spirometry,
  ita_biomarkers,
  ita_medications,
  ita_biologics,
  file = here::here("datasets", "ita_srd.Rdata")
)
