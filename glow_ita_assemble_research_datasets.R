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
# Patient details - done
# Lifestyle - done
# Comorbidities - focus on T2 and OCS-related for now
# Allergens - needs processing to get final variable
# HCRU - check timeframes*
# Exacerbations - very few dates, compare numbers against reported counts
# Asthma control - done
# Spirometry - done
# Biomarkers - finish function
# Medications - done
# Biologics - done, queries sent to NE for site, exclude for now

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

# Visits ----
visits <- patient |>
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

clean_smoking <- function(df) {
  # Check where duplicated visit dates
  # Check: df |> count(smoke, smoking_status_cln)

  # Fill NA with previous,
  # Ex-smoker > Current = ok, keep current
  # Ex-smoker > Never = ignore never, bring ex forward
  # Ex-smoker > Ex-smoker = ok, keep ex-smoker
  # Current > Current = ok
  # Current > Ex_smoker = ok
  # Current > Never = ignore never, bring current forward
  # Never smoker > Current = ok
  # Never smoker > Ex-smoker = never changes to ex-smoker, except never>current>ex, current left alone
  # Never smoker > Never = ok

  df |>
    mutate(smoke = as.character(smoke)) |>
    arrange(patient_id, visit_date) |>
    mutate(
      prior_smoked = dplyr::lag(
        cumany(smoke %in% c("Ex-smoker", "Current smoker")),
        default = FALSE
      ),
      # Blank out rows to be replaced (NA, or an impossible never)
      smoking_status_cln = if_else(
        smoke == "Never smoked" & prior_smoked,
        NA_character_,
        smoke
      ),
      .by = patient_id
    ) |>
    fill(smoking_status_cln, .direction = "down", .by = patient_id) |>
    mutate(
      # Find earliest ex/current record
      first_smk = first(smoke[smoke %in% c("Ex-smoker", "Current smoker")]),
      # If earliest ex/current record is "Ex-smoker" then replace "Never" with "Ex-smoker"
      # Leaves all "Never" and any "Current" alone
      smoking_status_cln = if_else(
        smoking_status_cln == "Never smoked" & first_smk %in% "Ex-smoker",
        "Ex-smoker",
        smoking_status_cln
      ),
      .by = patient_id
    ) |>
    select(-prior_smoked, -first_smk)
}
clean_lifestyle <- function(df = lifestyle_labelled) {
  df |>
    mutate(
      wt = as.numeric(wt),
      ht = as.numeric(ht),
      # Range checks
      across(
        c(wt),
        \(var) {
          var < 20 | var > 635
        },
        .names = "{.col}_rangefail"
      ),
      across(
        c(ht),
        \(var) {
          var < 1.3 | var > 2.2
        },
        .names = "{.col}_rangefail"
      ),
      weight = structure(if_else(!wt_rangefail, wt, NA), label = "Weight, kg"),
      height = structure(if_else(!ht_rangefail, ht, NA), label = "Height, m"),
      bmi = structure(
        case_when(
          ht_rangefail | wt_rangefail ~ NA,
          is.na(weight) | is.na(height) ~ NA,
          .default = weight / (height^2)
        ),
        label = "BMI, kg/m^2"
      ),
      obese = structure(
        case_when(
          bmi < 30 ~ "No",
          bmi >= 30 ~ "Yes",
          .default = NA
        ) |>
          factor(levels = c("Yes", "No")),
        label = "Obese (BMI ≥30 kg/m2)"
      )
    )
}

# Looks like smoking only recorded once, at baseline presumably
lifestyle <- patient_visits |>
  select(
    patient_id = uidpat,
    uidpatvis,
    wt,
    ht = ita_height_at_visit, # Italian height variable at visit
    smoke # Varies by visit, a lot of Never smokers convert to Ex-smokers
  ) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  # Fill height (not weight) missings across visits, last measure forward and also back
  fill(ht, .direction = "downup", .by = patient_id) |>
  # Range check height and weight and calculate BMI, obesity
  clean_lifestyle() |>
  # Function cleans smoking status to plausible values over time
  clean_smoking() |>
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
comorbidities_cln <- comorbidities |>
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
    t2diab,
    t2diab_dte,
    t2diab_dte_imput,
    # hfailure,
    # hfailure_stdate,
    # mi,
    # mi_stdate,
    # stroke,
    # stroke_stdate,
    # embolism,
    # csd = circulatory system disease?
    csd,
    csd_type, #"1-5?"
    csd_type_oth,
    csd_stdate = csd_dte,
    csd_dte_imput,
    # glaucoma,
    # glaucoma_stdate
    # ocular contains glaucoma and cataracts?
    ocular,
    ocular_type,
    ocular_type_oth,
    ocular_dte,
    ocular_dte_imput,
    depress,
    depress_dte,
    depress_dte_imput,
    anxiety,
    anxiety_dte,
    anxiety_dte_imput,
    osa = obsleepapnoea_yn,
    osa_stdate = obsleepapnoea_stdt,
    osa_stdate_imputed = obsleepapnoea_stdt_imput,
    renalf,
    renalf_dte,
    renalf_dte_imput,
    peptic_ulcer = peptic,
    peptic_ulcer_stdate = peptic_dte,
    peptic_dte_imput,
    pneumonia = pneum,
    pneumonia_dte = pneum_dte,
    pneumonia_stdate_imputed = pneum_dte_imput,
    gerd,
    # gerd_stdate,
    copd,
    egpa,
    # egpa_stdate,
    hypereos = hyper,
    # hypereos_stdate
  ) |>
  rename_with(~ str_replace(., "imput$", "imputed"), contains("imput")) |>
  rename_with(~ str_replace(., "dte", "stdate"), contains("dte")) |>
  # select(
  #   where(\(x) any(!is.na(x))), # Keeps columns with any information
  # ) |>
  left_join(
    visits |>
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
    hfailure = case_when(
      csd == "Yes" & csd_type == "1" ~ "Ever",
      csd == "Yes" ~ "Never",
      csd == "No" ~ "Never",
      .default = NA
    ) |>
      factor(levels = c("Ever", "Never")),
    mi = case_when(
      csd == "Yes" & csd_type == "2" ~ "Ever",
      csd == "Yes" ~ "Never",
      csd == "No" ~ "Never",
      .default = NA
    ) |>
      factor(levels = c("Ever", "Never")),
    stroke = case_when(
      csd == "Yes" & csd_type == "3" ~ "Ever",
      csd == "Yes" ~ "Never",
      csd == "No" ~ "Never",
      .default = NA
    ) |>
      factor(levels = c("Ever", "Never")),
    embolism = case_when(
      csd == "Yes" & csd_type == "4" ~ "Ever",
      csd == "Yes" ~ "Never",
      csd == "No" ~ "Never",
      .default = NA
    ) |>
      factor(levels = c("Ever", "Never")),
    csd_other = case_when(
      csd == "Yes" & csd_type == "5" ~ "Ever",
      csd == "Yes" ~ "Never",
      csd == "No" ~ "Never",
      .default = NA
    ) |>
      factor(levels = c("Ever", "Never"))
  )

# Allergens ----
# Use AB's SRD approach for processing for now, though I'm not sure the QISAR structure makes senses
# 20261111 AS/AT - sat_posallergen___0 and spt_posallergen___0 correspond to "No data" - all NA for Italy at least
# Not sure if "No data" and "No positive allergen" make sense.
# AB's approach probably most sensible but unclear whether to use "positive" or
# "clinically relevant" - this variable was dropped from SRD and is only available for serum tests.
allergens_cln <- allergens |>
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
  select(-contains("_none")) |>

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
    visits |>
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
summarise(
  allergy = if_else(any(allergy == "Yes"), "Yes", "No"),
  .by = patient_id
) |>
  mutate(
    allergy = structure(
      factor(allergy, levels = c("Yes", "No")),
      label = "Any Positive Allergen Test"
    )
  )

# HCRU ----
# Need timeframes
hcru <- patient_visits |>
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
    visits |>
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

# Asthma control ----
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
          "airq_score",
          "asthma_control",
          "asthma_control_source"
        )
      ),
      everything(),
      -country
    )
}
# Italy doesn't include GINA question 1 on activity, based on ACT > ACQ
asthma_control_cln <- patient_table |>
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
    visits |>
      select(uidpatvis, visit_id, visit_date),
    by = join_by(uidpatvis)
  ) |>
  relocate(visit_id, visit_date, .after = country) |>
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
    # Drop rows with missing measurement date is missing
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

spirometry_cln <- spirometry |>
  mutate(
    p = str_match(uidpat, "^([^_]+)_([^_]+)_(\\d+)$"),
    patient_id = str_c(
      p[, 2],
      p[, 3],
      "_",
      str_pad(p[, 4], width = 5, pad = "0")
    ),
    q = str_match(uidpatvis, "^([^_]+)_([^_]+)_(\\d+)_(.*)$"),
    uidpatvis = str_c(
      q[, 2],
      q[, 3],
      "_",
      str_pad(q[, 4], width = 5, pad = "0"),
      "_",
      q[, 5]
    )
  ) |>
  select(-p, -q, -uidpat) |>
  left_join(
    visits |>
      select(uidpatvis, visit_id, visit_date, age_at_visit),
    by = join_by(uidpatvis)
  ) |>
  left_join(
    patient_table |>
      select(patient_id, gender),
    by = join_by(patient_id)
  ) |>
  left_join(
    lifestyle |>
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

# Still a few with FVC < FEV1
spirometry_cln |>
  filter(if_any(contains("implausible")) | if_any(contains("rangefail"))) |>
  pull(patient_id) |>
  noquote()

# Biomarkers ----
# Question - impute missing dates with visit date, then check for duplicates?
# Italy doesn't provide IgE dates

clean_biomarkers <- function(
  df = biomarkers_labelled,
  beos_range_low = 0,
  beos_range_high = 5000,
  feno_range_low = 0,
  feno_range_high = 300,
  ige_range_low = 0,
  ige_range_high = 100000,
  reduce_function = "max",
  drop_duplicate = TRUE,
  drop_rangefail = FALSE
) {
  if (!reduce_function %in% c("min", "max", "mean")) {
    stop('De-duplicating function must be "min", "mean", or "max"')
  }

  df <- df |>
    # remove_empty_imputation_fields() |>
    # rename_with(~ str_replace(.x, "_rlt", "_cnt"), contains("feno")) |>
    # pivot_longer(
    #   cols = -c(
    #     patient_id,
    #     visit_id,
    #     visit_id_imputed,
    #     visit_id_imputed_flag,
    #     result_date,
    #     result_date_imputed,
    #     result_date_imputed_flag
    #   ),
    #   names_to = c("test", ".value"),
    #   names_pattern = "^(.+?)_(cnt|cnt_imputed_flag|cnt_imputed)$"
    # ) |>
    # filter_out(is.na(cnt)) |>
    # rename_with(~ str_replace(.x, "cnt", "result"), contains("cnt")) |>
    # mutate(result_date = ymd(result_date)) |>
    # Range checks against values specified above
    mutate(
      exclude_rangefail = case_when(
        test == "beos" &
          (result < beos_range_low | result > beos_range_high) ~ TRUE,
        test == "feno" &
          (result < feno_range_low | result > feno_range_high) ~ TRUE,
        test == "ige" &
          (result < ige_range_low | result > ige_range_high) ~ TRUE,
        .default = FALSE
      ),
    ) |>
    mutate(
      duplicate_result = case_when(
        n() > 1 ~ TRUE,
        .default = FALSE
      ),
      keep_result = case_when(
        duplicate_result &
          reduce_function == "min" &
          result != min(result, na.rm = TRUE) ~ FALSE,
        duplicate_result &
          reduce_function == "mean" &
          result != mean(result, na.rm = TRUE) ~ FALSE,
        duplicate_result &
          reduce_function == "max" &
          result != max(result, na.rm = TRUE) ~ FALSE,
        .default = TRUE
      ),
      .by = c(patient_id, date, test)
    )

  if (drop_duplicate) {
    df <- df |>
      filter(keep_result) |>
      # Some patients still have identical rows - these can be dropped, keeps first
      distinct(patient_id, date, test, .keep_all = TRUE) |>
      select(-keep_result, -duplicate_result)
  }

  if (drop_rangefail) {
    df <- df |>
      filter_out(exclude_rangefail) |>
      select(-exclude_rangefail)
  }
}

biomarkers_cln <- biomarkers |>
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
    date = if_else(!is.na(date), dmy(date), visit_date), # Almost all test dates, where available for FeNO, are same as visit date, except ~20 earlier ones
    result = as.numeric(result)
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

# Use supplement table for querying and creating biologic table below
medications_cln <- medications |>
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
bx_cln <- medications_cln |>
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
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date)

# Italy-specific queries, sent to NE 20261120:
# Non-biologic naive - None found - can we confirm?
# Duplicate start dates -	Patient has more than one record with same start date, regardless of whether therapies different or not
# End date before start date - Row-specific, end and start dates swapped?
# All future treatments - All of patient's start dates are after their last visit date
# Any future treatments - Any of patient's start dates are after their last visit date
# Missing end date - Did not include this query for Italy - can we confirm that these are "ongoing" therapies?
# Overlapping episodes (>1 day) - same	Prior end date overlaps with next start date by more than 1 day (on the same therapy) - can we combine these episodes?
# Overlapping episodes (>1 day) - different	Prior end date overlaps with next start date by more than 1 day (on a different therapy) - can we confirm if this is plausible, potentially "bridging" therapies? What would be the maximum overlap allowable for a check?
# Gap on same therapy (>0 months) -	Apparent gap in treatment based on prior end date and next start date on the same therapy - can we confirm if this is plausible and why? Is there a length of time we can use to combine episodes or to identify stops?

## Exclusions ----
# Based on those used for UK
bx_cln_exclusions <- patient_table |>
  select(patient_id, date_of_birth) |>
  left_join(
    bx_cln |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    index_age = interval(date_of_birth, min(bx_start_date)) / years(1),
    age_at_baseline = interval(
      date_of_birth,
      min(bx_start_date) %m-% months(12)
    ) /
      years(1),
    # Look at overlaps (1 day probably okay)
    overlap_same = bx_name_short ==
      dplyr::lag(bx_name_short, default = first(bx_name_short)) &
      interval(
        dplyr::lag(bx_end_date, default = first(bx_start_date)),
        bx_start_date
      ) /
        days(1) <
        -1,
    overlap_diff = bx_name_short !=
      dplyr::lag(bx_name_short, default = first(bx_name_short)) &
      interval(
        dplyr::lag(bx_end_date, default = first(bx_start_date)),
        bx_start_date
      ) /
        days(1) <
        -1,
    tx_gap = bx_name_short ==
      dplyr::lag(bx_name_short, , default = first(bx_name_short)) &
      interval(
        dplyr::lag(bx_end_date, default = first(bx_start_date)),
        bx_start_date
      ) /
        months(1) >
        0,
    first_biologic = bx_start_date == min(bx_start_date),
    follow_up_mos = interval(bx_start_date, last_visit_date) / months(1),
    exclusions = case_when(
      all(is.na(bx_data)) ~ "No biologic information",
      any(is.na(bx_name_short)) ~ "Missing biologic name",
      any(is.na(bx_start_date)) ~ "Missing start date",
      # prior_biologic_chk == "prior biologic not listed" ~ "Non-biologic naive",
      n_distinct(bx_start_date) < n() ~ "Duplicate start dates",
      any(bx_end_date < bx_start_date) ~ "End date before start date",
      all(bx_start_date > last_visit_date) ~ "All future treatments",
      any(bx_start_date > last_visit_date) ~ "Any future treatments",
      # any(is.na(bx_end_date)) ~ "Missing end date", # Missing end date currently ongoing tx...
      any(overlap_same) ~ "Overlapping episodes (>1 day) - same",
      any(overlap_diff) ~ "Overlapping episodes (>1 day) - different",
      any(tx_gap) ~ "Gap on same therapy (>0 months)",
      age_at_baseline <
        18 ~ "Under-18 at baseline, 12 months prior to initiation",
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
          "Overlapping episodes (>1 day) - same",
          "Overlapping episodes (>1 day) - different",
          "Gap on same therapy (>0 months)",
          "Under-18 at baseline, 12 months prior to initiation",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up"
        )
      ),
    exclusions_table = case_when(
      exclusions == "No biologic information" ~ "No biologic information",
      exclusions ==
        "Under-18 at baseline, 12 months prior to initiation" ~ "Under 18 years-old",
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
          "Under 18 years-old",
          "Non-biologic naive",
          "Data queries",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up",
          "Eligible"
        )
      ),
    .by = patient_id
  )

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

writexl::write_xlsx(
  bx_cln_exclusions |>
    filter(exclusions_table == "Data queries") |>
    distinct(patient_id, exclusions),
  here::here("data_cleaning", "ita_bx_data_queries_20261120.xlsx")
)

## Analysis table ----
bx_analysis <- bx_cln_exclusions |>
  filter(is.na(exclusions)) |>
  select(
    patient_id,
    data_capture_system,
    country,
    setting,
    site,
    contains("visit"),
    bx_start_date,
    bx_end_date,
    bx_ongoing,
    bx_name_short,
    bx_type,
    bx_switchreason
  ) |>
  # ITA outcome =
  # ITA continuers: 1 biologic episode, flagged as ongoing, or stop indicated more than a month after last visit
  # ITA stoppers: 1 biologic episode, flagged as stopped before or within a month of last visit; >1 biologic and gap between 1st and 2nd >6 months
  # ITA switchers: >1 biologic episode, <=6 months between 1 and 2nd
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
  file = "ita_biologic_data.Rdata"
)
