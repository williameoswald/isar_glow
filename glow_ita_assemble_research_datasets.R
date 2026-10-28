# Pull Italy QISAR2 data and process for GLOW analysis
library(RODBC)
library(tidyverse)
library(janitor)

# Notes ----
# 20261009 AT has been asked to revisit rebuild to ensure she pulls in 
# all information and uses all transferred data so things may change. 
# Worth finishing processing here, so we can make baseline and compare with G's
# 20261015 Italy provides DOB as YYYY-MM-UNK -> Impute as 15th and calculate


# !!!!! - Exclude 419 patients with 2="No" for both inclusion criteria?

# To Do ----
# Medications - probably better use supplement table
# Exacerbations
# Biologics
# HCRU
# Comorbidities
# LTOCS - pending

# Pull NJH data ----
## Connect to ISAR SQL tables
connection <- RODBC::odbcDriverConnect(
  'driver={SQL Server};
                              server=dsimport,4071;
                              database=QISAR2;
                              trusted_connection=true'
)

# Italy ----
ita_data <- c(
  "RR_asthma_medications",
  "RR_blood_and_sputum_eosinophils_detail",
  "RR_exacerbation_history_detail",
  "RR_patient",
  "RR_patient_visits",
  "RR_spirometry_predicted_values",
  "ST_medications", # Medications as shared by Italy pre-ISAR processing
  "ST_exitreason" # Lost to follow up flags (Italy only)
  ) |>
  set_names(nm = \(x) {
    str_remove(x, "_detail$") |>
      str_remove("^RR_") |>
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
load("qisar2_ita_data.RData")

#Split named list into objects
list2env(ita_data, envir = .GlobalEnv)

# Dictionary ----
# variable_list <- ita_data |>
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
# write_xlsx(variable_list, "qisar2_ita_data_dictionary.xlsx")

# Visits ----
visits <- patient |> 
  # n = 4120, 3 duplicated patient records
  # Italy provides DOB so use this instead of age (not sure where this comes from), with day imputed as 15th
  distinct(patient_id = uidpat) |>
  # Right join because some patients in visits are missing from patient
  right_join(
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
    .by = patient_id
  ) |>
  mutate(
    num_of_visits = structure(num_of_visits, label = "Number of visits"),
    first_visit_date = structure(first_visit_date, label = "First visit date"),
    last_visit_date = structure(last_visit_date, label = "Last visit date"),
    # age_at_visit = structure(age_at_visit, label = "Age at visit, years")
  ) |>
  # select(-age_first_visit) |>
  arrange(patient_id, visit_date)

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

# 69 patients with visits out of order
visits_chks |>
  filter(visit_order != date_order) |>
  distinct(patient_id) |>
  pull(patient_id)
# ITA_SITE13_13 baseline after followup
# ITA_SITE156_2 no follow up 3, and 4 is six months after 2
# ITA_SITE16_128 followup 1 out of order
# ITA_SITE16_9 followup 8 and 9 out of order
# ITA_SITE17_2 followup 6 and 7 out of order
# ITA_SITE20_27 followup 6 and 7 out of order
# ITA_SITE21_14 followup 8 and 9 out of order
# ITA_SITE22_11 no visit 4
# ITA_SITE22_24 followup 4 and 5 out of order
# ITA_SITE22_39 followup 3 and 4 out of order
# ITA_SITE23_16 no visit 7
# 13 patients where baseline is not first visit
visits_chks |>
  filter(str_detect(visit_id, "baseline") & date_order != 0) |>
  distinct(patient_id) |>
  pull(patient_id) |>
  noquote() |>
  kable() |>
  kable_styling()

# Patient table ----
# Patient-level, i.e. non-varying details, no separate Demographic table
patient_table <- patient |> 
  # n = 4120, 3 duplicated patient records
  distinct(patient_id = uidpat, site = dagarm, .keep_all = TRUE) |>
  # Drop empty fields
  select(
    patient_id,
    site,
    where(\(x) any(!is.na(x))),
    -contains("comstat"),
    -record_id
  ) |> 
  # Add first and last visit
  left_join(
    visits |> 
      distinct(patient_id,
             first_visit_date,
             last_visit_date,
             num_of_visits),
    by = join_by(patient_id),
    relationship = "one-to-one"
  ) |>
  mutate(
    gender = structure(
      factor(gen, levels = c("Male", "Female")),
      label = "Gender"
    ),
    date_of_birth = structure(ymd(str_replace(dob, "UNK", "15")),
                              label = "Date of birth, day imputed"
    ),
    ethnicity = structure(eth, label = "Ethnicity"),
    asthma_onset = structure(as.numeric(asmage), label = "Age at asthma onset, years"),
    # Convert height to metres
    ht = structure(as.numeric(ht) / 100, label = "Height, metres"),
    country = structure("Italy", label = "Country"),
    site = structure(dagarm, label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system"),
    first_visit_date = structure(first_visit_date, label = "First visit date"),
    last_visit_date = structure(last_visit_date, label = "Last visit date"),
    num_of_visits = structure(num_of_visits, label = "Number of visits")
  ) |> 
  select(
    patient_id,
    country,
    site, # site is different to setting
    first_visit_date,
    last_visit_date,
    num_of_visits,
    data_capture_system,
    gender,
    date_of_birth,
    ethnicity,
    asthma_onset,
    ht
  ) |> 
  arrange(first_visit_date)
# Impute DOB for Italy (add 15th for day) and calculate age at index_date and at spirometry test

# Age check
# agechk <- patient |>
#   select(patient_id = uidpat, dob, age) |>
#   right_join(patient_visits |>
#                select(patient_id = uidpat, visit_id = uidpatvis, dov),
#              by = join_by(patient_id)) |>
#   mutate(
#     year_of_birth = str_sub(dob, 1, 4),
#     dov = ymd(dov),
#     date_of_birth = ymd(str_replace(dob, "UNK", "15")),
#     age_at_visit = floor(interval(date_of_birth, dov)/years(1))
#   ) |>
#   arrange(patient_id, dov) |> 
#   filter(any(age != min(age_at_visit, na.rm = T)), .by = patient_id)
# "Age" appears to be age at baseline visit
# With this filter provided age is usually one year higher than calculated so we know DOB is that month AND prior to visit date

# Notes:
# Aaron's SRD approach appears complicated because of need to get additional 
# vars: index_date, max_beos_highest, highest_bec_imputed_flag, bt
# Better I think to process biologics, BEC, bt separately and pull in later

# Lifestyle ----
# Q - Need to know WHEN height, weight, BMI, smoking_status are recorded, asked AT 20261015
lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    wt,
    bmi,
    smoke
    ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  # Keep earliest visit for smoke, wt, bmi - this aligns with UK
  slice_min(order_by = visit_date, by = patient_id, with_ties = T) |> 
  # get_dupes(patient_id) |> view()
  # ITA_SITE22_42, ITS_SITE368_208, ITA_SITE368_97 have baseline and followup visits on the same day
  filter_out(
    patient_id %in% c("ITA_SITE22_42", "ITA_SITE368_208", "ITA_SITE368_97") & str_detect(visit_id, "followup")
  ) |> 
  left_join(
    patient_table |> 
      select(patient_id, 
             ht) |> 
      distinct(),
    by = join_by(patient_id)
    ) |> 
  mutate(
    weight = structure(as.numeric(wt), label = "Baseline weight, kg"),
    height = structure(as.numeric(ht), label = "Baseline height, m"),
    bmi = structure(weight/ (height^2), label = "Baseline BMI, kg/m^2"),
    smoking_status = structure(
      factor(smoke, levels = c("Current Smoker", "Ex-Smoker", "Never Smoked")),
      label = "Smoking status, baseline"
    )
  ) |>
  select(
    patient_id,
    weight,
    height,
    bmi,
    smoking_status
  )

# Comorbidity ----
# 20261015 STOPPED HERE!!!  Trying to align Italy and UK approaches but AT delays
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

# Allergens ----
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
# Use number for now

# Asthma control ----
# Based on Aaron's SRD "Asthma_control_table.R"
# Italy doesn't include GINA question 1 on activity, based on ACT > ACQ
asthma_control <- total_patients |>
  left_join(
    patient_visits |>
      select(
        patient_id = uidpat,
        visit_id = uidpatvis,
        visit_date = dov,
        asthma_control,
        acq_score,
        act_score,
        gina_day_1,
        gina_activity_2 = gina_activity,
        gina_night_3,
        gina_reliever_4,
        gina_lf_5
      ),
    by = join_by(patient_id)
  ) |>
  mutate(
    visit_id = case_when(
      str_detect(visit_id, "barm") ~ as.character(paste0(
        patient_id,
        "_baseline"
      )),
      str_detect(visit_id, "historic") ~ as.character(paste0(
        patient_id,
        "_historic"
      )),
      str_detect(visit_id, "f\\d{1,2}arm1") ~ as.character(paste0(
        patient_id,
        "_followup",
        sub(".*f(\\d{1,2})arm1.*", "\\1", visit_id)
      )),
      str_detect(visit_id, "adeseguimiearm") ~ str_replace(
        visit_id,
        "_adeseguimiearm1([a-z]?)$",
        ~ {
          letter <- str_remove(.x, "_adeseguimiearm1")
          idx <- if_else(letter == "", 1L, match(letter, letters))
          paste0("_followup", idx)
        }
      )
    ),
    # Visit dates all plausible and no missing
    visit_date = ymd(visit_date),
    across(contains("score"), as.numeric),
    acq_version = case_when(
      !is.na(acq_score) &
        str_detect(patient_id, "POL|MEX|BRZ|GBR|ITA") ~ as.character("ACQ7"),
      !is.na(acq_score) & str_detect(patient_id, "BEL|DEN") ~ as.character(
        "ACQ6"
      ),
      !is.na(acq_score) ~ as.character("Unspecified ACQ")
    )
  ) |>
  # Exclude ACT scores outside range, 5 <= ACT <= 25
  filter_out(act_score > 25 | act_score < 5) |>
  # Exclude ACQ scores outside range, 0 - 6, ACQ is mean of 7 questions with 7-point scale (0 = no impairment; 6 = maximum impairment)
  filter_out(acq_score > 6 | acq_score < 0) |>
  mutate(
    # Classify ACQ score (version of test used is unavailable)
    acq_ac = case_when(
      acq_score <= 0.75 ~ "Well Controlled",
      acq_score < 1.5 ~ "Partly Controlled",
      acq_score >= 1.5 ~ "Uncontrolled",
      .default = NA
    ) |>
      factor(
        levels = c("Uncontrolled", "Partly Controlled", "Well Controlled")
      ),
    # Classify ACT score
    act_ac = case_when(
      act_score < 16 ~ "Uncontrolled",
      act_score < 20 ~ "Partly Controlled",
      act_score >= 20 ~ "Well Controlled",
      .default = NA
    ) |>
      factor(
        levels = c("Uncontrolled", "Partly Controlled", "Well Controlled")
      ),
    # Results are usually prioritised GINA > ACT > ACQ > AIRQ.
    # Italy has missing responses for GINA questions 1 and 4, so prioritise ACT.
    asthma_control = case_when(
      !is.na(act_ac) ~ act_ac,
      !is.na(acq_ac) ~ acq_ac,
      .default = NA
    )
  ) |>
  filter_out(is.na(asthma_control)) |>
  select(
    patient_id,
    visit_id,
    visit_date,
    asthma_control,
    act_score,
    act_ac,
    acq_score,
    acq_version,
    acq_ac
  ) |>
  # Remove duplicates
  distinct()


group_by(patient_id, visit_id) |>
  tidyr::fill(everything(), .direction = "updown") |>
  distinct |>
  ungroup |>
  mutate(
    visit_date_imputed_flag = case_when(
      visit_date < as.Date("2000-01-01") ~ as.character(
        "Visit date nulled as it was prior to 2000-01-01, likely a birth date"
      ),
      visit_date > Sys.Date() ~ as.character(
        "Visit date nulled as it was in the future"
      ),
      T ~ NA
    ),
    visit_date = ifelse(
      visit_date < as.Date("2000-01-01") | visit_date > Sys.Date(),
      NA,
      as.character(visit_date)
    ),
    visit_date_imputed = ifelse(!is.na(visit_date_imputed_flag), 1, 0),

    acq_version = case_when(
      !is.na(acq_score) &
        grepl("POL|MEX|BRZ|GBR|ITA", patient_id) ~ as.character("ACQ7"),
      !is.na(acq_score) & grepl("BEL|DEN", patient_id) ~ as.character("ACQ6"),
      !is.na(acq_score) ~ as.character("Unspecified ACQ")
    ),
  )
select(
  patient_id,
  visit_id,
  visit_id_imputed,
  visit_id_imputed_flag,
  visit_date,
  visit_date_imputed,
  visit_date_imputed_flag,
  asthma_control,
  asthma_control_imputed,
  asthma_control_imputed_flag,
  acq_score,
  acq_score_imputed,
  acq_score_imputed_flag,
  acq_version,
  acq_version_imputed,
  acq_version_imputed_flag,
  act_score,
  act_score_imputed,
  act_score_imputed_flag,
  gina_day_1,
  gina_day_1_imputed,
  gina_day_1_imputed_flag,
  gina_activity_2,
  gina_activity_2_imputed,
  gina_activity_2_imputed_flag,
  gina_night_3,
  gina_night_3_imputed,
  gina_night_3_imputed_flag,
  gina_reliever_4,
  gina_reliever_4_imputed,
  gina_reliever_4_imputed_flag,
  gina_lf_5,
  gina_lf_5_imputed,
  gina_lf_5_imputed_flag
) |>
  distinct |>
  mutate(
    visit_id = case_when(
      grepl("barm", visit_id) ~ as.character(paste0(patient_id, "_baseline")),
      grepl("historic", visit_id) ~ as.character(paste0(
        patient_id,
        "_historic"
      )),
      grepl("f\\d{1,2}arm1", visit_id) ~ as.character(paste0(
        patient_id,
        "_followup",
        sub(".*f(\\d{1,2})arm1.*", "\\1", visit_id)
      )),
      grepl("adeseguimiearm", visit_id) ~ str_replace(
        visit_id,
        "_adeseguimiearm1([a-z]?)$",
        ~ {
          letter <- str_remove(.x, "_adeseguimiearm1")
          idx <- if_else(letter == "", 1L, match(letter, letters))
          paste0("_followup", idx)
        }
      )
    ),
    visit_date = as.Date(visit_date),
    acq_score = as.numeric(acq_score),
    act_score = as.numeric(act_score),
    gina_day_1 = as.numeric(gina_day_1),
    gina_activity_2 = as.numeric(gina_activity_2),
    gina_night_3 = as.numeric(gina_night_3),
    gina_reliever_4 = as.numeric(gina_reliever_4),
    gina_lf_5 = as.numeric(gina_lf_5),
    asthma_control = case_when(
      asthma_control == "Partially Controlled" ~ as.character(
        "Partly Controlled"
      ),
      T ~ as.character(asthma_control)
    )
  ) |>
  inner_join(Visit |> distinct(visit_id), by = "visit_id") |>
  mutate(
    across(
      matches("gina_"),
      ~ case_when(
        . == "0" ~ NA,
        . == "1" ~ "Yes",
        . == "2" ~ "No"
      )
    )
  )

# Spirometry ----

# Biomarkers ----
# Based on Aaron's SRD "Biomarker_table.R"
# Italy doesn't provide IgE dates
biomarkers <- total_patients |>
  left_join(
    blood_and_sputum_eosinophils |>
      select(
        patient_id = uidpat,
        visit_id = uidpatvis,
        where(\(x) any(!is.na(x))),
        -uiddb,
        -record_id,
        -dagarm,
        -contains("redcap"),
        -beos_yn,
        beos_result = beos_cnt,
        beos_date = beos_dte,
        ige_result = ige_cnt,
        feno_result = feno_rlt
      ),
    by = join_by(patient_id)
  ) |>
  mutate(
    ige_date = NA
  ) |>
  pivot_longer(
    cols = !all_of(c("patient_id", "visit_id")),
    names_to = c("test", ".value"),
    names_sep = "_"
  ) |>
  filter_out(is.na(result), is.na(date)) |>
  mutate(
    result = as.numeric(result),
    date = ymd(date),
    result = 
      case_when(
        # Truncate IgE and FeNO counts (same as NJH)
        test == "ige" & result <= 2 ~ 2,
        test == "ige" & result >= 5000 ~ 5000,
        test == "feno" & result <= 5 ~ 5,
        test == "feno" & result <= 300 ~ 300,
        test == "beos"
        .default = result
      ),
    unit = 
      case_when(
        test == "ige" ~ 
      )
  ) |> 
  
  
# left_join(visits |> select(visit_id, visit_date), by = join_by(visit_id)) |>
# relocate(visit_date, .after = visit_id)

# Need range checks on IGE and FeNO and units for BEC
# Some values recorded as >5000 or <2, truncate all the same way
result = case_when(
  ige_value <= 2 ~ 2,
  ige_value >= 5000 ~ 5000,
  .default = round(ige_value, 0)
)
# Some values recorded as >300 or <5, truncate all the same way
result = case_when(
  feno_value <= 5 ~ 5,
  feno_value >= 300 ~ 300,
  .default = round(feno_value, 0)
)
# Exclude ACT scores outside range, 5 <= ACT <= 25
filter_out(act_score > 25 | act_score < 5) |>
  # Exclude ACQ scores outside range, 0 - 6, ACQ is mean of 7 questions with 7-point scale (0 = no impairment; 6 = maximum impairment)
  filter_out(acq_score > 6 | acq_score < 0) |>
  # Some values recorded as >300 or <5, truncate all the same way
  result = case_when(
    feno_value <= 5 ~ 5,
    feno_value >= 300 ~ 300,
    .default = round(feno_value, 0)
  )
result = case_when(
  ige_value <= 2 ~ 2,
  ige_value >= 5000 ~ 5000,
  .default = round(ige_value, 0)
)
result = case_when(
  # Convert to cells/uL (multiply K/uL or 10*9/L by 1000)
  unit == "K/uL" ~ eos_value * 1000,
  unit == "10*9/L" ~ eos_value * 1000,
  unit == "Unknown" ~ eos_value * 1000, # 6 unknown values, all other tests for these patients in K/uL and results look similar
  .default = eos_value
),
orig_unit = unit,
unit = "cells/uL"
) |>


# Medications ----
# GS - therapy data provided by Italy is not linked to visits. The file lists 
# start dates (and potentially en dates if treatment has been discontinued) for 
# all asthma-related drugs. It should be a lot easier to use this data without 
# taking visits into consideration. We need to be careful not assuming « No » 
# when the drug is not listed: eg, for LTOCS, I am assuming missing when only 
# bx is listed there (it’s not plausible that a patient would initiate bx with 
# no other asthma-related drug). I am also looking at the time window covered by
# the data, ie data can be entered in the CRF but the therapy data not entered.

# Table with pre-processed medications - this will be a supplement table
# SELECT * FROM QISAR2.ita.example_medications

# Based on Aaron's SRD "Medication_table.R"
medications <- asthma_medications |>
  # Keep non-empty columns
  select(
    patient_id = uidpat,
    visit_id = uidpatvis,
    where(\(x) any(!is.na(x))),
    -uiddb,
    -record_id,
    -dagarm,
    -contains("redcap"),
    # All "other" tezepelumab recorded under antiil5
    -antiil5_oth
  ) |>
  # Bare drug columns (no underscore) -> "<drug>_name"
  rename_with(
    ~ paste0(.x, "_name"),
    .cols = !all_of(c("patient_id", "visit_id")) & !contains("_")
  ) |>
  mutate(
    # Recode "1" to "Yes"
    across(c(lama_yn, ltra_yn, theo_yn, mabx_yn), \(x) {
      replace_values(x, "1" ~ "Yes")
    }),
    visit_id = case_when(
      str_detect(visit_id, "barm") ~ as.character(paste0(
        patient_id,
        "_baseline"
      )),
      str_detect(visit_id, "historic") ~ as.character(paste0(
        patient_id,
        "_historic"
      )),
      str_detect(visit_id, "f\\d{1,2}arm1") ~ as.character(paste0(
        patient_id,
        "_followup",
        sub(".*f(\\d{1,2})arm1.*", "\\1", visit_id)
      )),
      str_detect(visit_id, "adeseguimiearm") ~ str_replace(
        visit_id,
        "_adeseguimiearm1([a-z]?)$",
        ~ {
          letter <- str_remove(.x, "_adeseguimiearm1")
          idx <- if_else(letter == "", 1L, match(letter, letters))
          paste0("_followup", idx)
        }
      )
    ),
    # Recode "other" macrolides: Ciprofloxacin; Levofloxacin
    mabx_name = case_when(
      mabx_name == "Other" & !is.na(mabx_oth) ~ mabx_oth,
      .default = mabx_name
    ),
    # Recode "other" LABA: Terbutaline
    laba_name = case_when(
      laba_name == "Other" & !is.na(laba_oth) ~ laba_oth,
      .default = laba_name
    )
  ) |>
  select(-contains("_oth")) |>
  pivot_longer(
    cols = !all_of(c("patient_id", "visit_id")),
    names_to = c("drug", ".value"),
    names_sep = "_"
  ) |>
  # All missing where yn == NA & yn == "No"
  filter(yn == "Yes") |>
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
  left_join(visits |> select(visit_id, visit_date), by = join_by(visit_id)) |>
  relocate(visit_date, .after = visit_id) |>
  arrange(patient_id, drug, visit_date)
# Need to get first instance of biologics and build spells
# Can we assume they're on a drug if it's listed against the last visit?
# Other drugs similarly duplicated

# Medications Supplement ----
s_medications <- supplement_medications |> 
  rename(start_date = therapy_start,
         end_date = therapy_end,
         name = drug) |> 
  mutate(
    drug_type = case_when(
      drug_type %in% c("antiige", "antiil4", "antiil5", "antitslp") ~ "Biologic",
      drug_type == "mos" ~ "LTOCS",
      drug_type == "ics" ~ "ICS only",
      drug_type == "macrolide" ~ "Macrolides",
      .default = str_to_upper(drug_type)
    ),
    across(contains("date"), ymd)
  ) |> 
  filter_out(
    if_all(c(drug_type, name, ldose, freq, switchreason), is.na)
    )

# Biologics ----
bx_cln <- medications |> 
  filter(drug_type == "Biologic") |> 
  arrange(patient_id, start_date) |> 
  mutate(bx_ct = n(),
         bx_nbr = row_number(),
         bx_name = structure(
           factor(name,
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
             bx_name %in% c("Mepolizumab", "Reslizumab") ~ "Anti-IL5",
             bx_name %in% c("Benralizumab") ~ "Anti-IL5R",
             bx_name == "Omalizumab" ~ "Anti-IgE",
             bx_name == "Dupilumab" ~ "Anti-IL4Ralpha",
             bx_name == "Tezepelumab" ~ "Anti-TSLP"
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
           if_else(is.na(end_date), "Yes", "No") |> 
             factor(levels = c("Yes", "No")),
           label = "Biologic ongoing"
         ),
         reason_no_clinical_efficacy = str_detect(switchreason, "No clinical efficacy"),
         reason_patients_choice = str_detect(switchreason, "Patient's choice"),
         reason_side_effects = str_detect(switchreason, "Side effects"),
         reason_access_restriction = str_detect(switchreason, "Restrictions on access to biological drugs"),
         .by = patient_id) |> 
  select(patient_id, contains("date"), bx_ongoing, bx_ct, bx_nbr, bx_name, bx_type, contains("reason_"))

single <- bx_cln |> 
  filter(bx_ct == 1)

# Of patients with end dates, 50% have switch reasons
# Assume they remain on treatment, unless there's an end date? Or assume they stop only where there's a switchreason?
single |> filter(!is.na(end_date)) |> tabyl(switchreason)

multiple <- bx_cln |> 
  filter(bx_ct > 1) |> 
  mutate(
    benra_reliever = !is.na(end_date) & bx_name %in% "Benralizumab" & start_date==end_date
  ) |> 
  group_by(patient_id) 
keys  <- multiple |> group_keys() |> pull(patient_id)

results <- multiple |> 
  group_split() |>
  set_names(nm = keys) |>
  map(\(patient) process_patient_v3(patient), .progress = T)

spells <- results |>
  map("spells") |>
  list_rbind() |> 
  select(-effective_end) # not needed



ITA_SITE13_4 two ongoing?
ITA_SITE13_19 two months then stopped for several month - look within first year?
ITA_SITE13_11 4 month gap between oma
ITA_SITE10_45 4 month gap between mepo
ITA_SITE366_23 Benra injection 1day for exacerbation?


  pivot_wider(id_cols = c(patient_id, bx_ct),
              names_from = bx_nbr,
              values_from = -c(patient_id, bx_ct),
              names_vary = "slowest") |> 
  mutate(overlaps =
           case_when(
             end_date_1 > start_date_2 ~ TRUE,
             end_date_2 > start_date_3 ~ TRUE,
             end_date_3 > start_date_4 ~ TRUE,
             end_date_4 > start_date_5 ~ TRUE,
             
           ))
# 1) pivot wide
# 2) overlaps?
# 3) same drug - minimal gap < 3 months?


process_patient_v3 <- function(patient_data, month_gap_allowed = 3) {
  patient_id <- patient_data |> distinct(patient_id) |> pull()
  
  spells <- patient_data |>
    filter_out(is.na(start_date)) |>
    arrange(start_date, end_date) |>
    # ID instances of biologic use in chronological order
    mutate(
      change_med = bx_name != dplyr::lag(bx_name, default = first(bx_name)),
      change = if_else(
        bx_name != dplyr::lag(bx_name, default = first(bx_name)),
        1,
        0
      ),
      bx_id = cumsum(as.integer(change)) + 1
    ) |>
    # Within instance flag rows to create spells
    mutate(
      running_max = purrr::accumulate(end_date, max),
      prev_max = dplyr::lag(running_max, default = first(end_date)),
      new_spell = start_date > prev_max %m+% months(month_gap_allowed),
      .by = bx_id
    ) |>
    mutate(spell_id = cumsum(new_spell), .by = bx_id) |>
    # Within instance summarise rows to create spells
    summarise(
      spell_start = min(start_date),
      spell_end = max(end_date),
      .by = c(bx_id, bx_name, spell_id)
    ) |>
    mutate(
      patient_id = patient_id,
      .before = 1
    ) |>
    arrange(bx_id, spell_start) |>
    mutate(
      next_spell_start = lead(spell_start),
      effective_end = if_else(is.na(spell_end), next_spell_start, spell_end),
      .by = bx_id
    ) |>
    select(
      patient_id,
      bx_name,
      bx_id,
      spell_start,
      spell_end,
      effective_end
    ) |>
    arrange(spell_start)
  
  # Rows with no start_date or end_date: confirm  drug is present in listed spells.
  dropped_rows_missing_all_dates <- patient_data |>
    filter(when_all(is.na(start_date), is.na(end_date))) |>
    distinct(bx_name, .keep_all = T) |>
    rowwise() |>
    mutate(
      biologic_listed = bx_name %in% spells$bx_name
    ) |>
    ungroup()
  
  # Rows with no start_date: check whether end_date inside spell for same drug,
  # using effective_end (start of next spell for that drug) for intervals with no end date.
  # TRUE = already covered, safe to drop; FALSE = not covered by any known
  # treatment window, review manually.
  dropped_rows_missing_start_dates <- patient_data |>
    filter(when_all(is.na(start_date), !is.na(end_date))) |>
    rowwise() |>
    mutate(
      end_date_in_known_spell = any(
        spells$bx_name == bx_name &
          spells$spell_start <= end_date & end_date <= spells$effective_end,
        na.rm = TRUE
      )
    ) |>
    ungroup()
  
  list(
    spells = spells,
    dropped_rows_missing_all_dates = dropped_rows_missing_all_dates,
    dropped_rows_missing_start_dates = dropped_rows_missing_start_dates
  )
}
