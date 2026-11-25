# Pull NJH QISAR2 data and process for ISAR analysis
library(RODBC)
library(tidyverse)
library(janitor)
library(ISARanalysis)
library(kableExtra)
library(rspiro)


# To Do ----
# Visits -
# Patients -
# Lifestyle - done, need to confirm smoking status values 1,2,3,4
# Comorbidities - focus on T2 and OCS-related for now, needs processing
# Allergens -
# HCRU -
# Exacerbations -
# Asthma control -
# Spirometry -
# Biomarkers -
# Medications -
# Biologics - done, queries sent, exclude patients for now

# 20261109 - need to convert patients to bx_cln_exclusions to match other sites?

# Notes ----
# 20261109 - realised that needed to comment out "any missing start dates" in
# exclusions steps to get 20261015 feasibility numbers.

# Pull NJH data ----
## Connect to ISAR SQL tables
connection <- RODBC::odbcDriverConnect(
  'driver={SQL Server};
                              server=dsimport,4071;
                              database=QISAR2;
                              trusted_connection=true'
)

njh_data <- c(
  "AIRQ",
  "Allergens",
  "AsthmaControlTest",
  "Biologics",
  "Comorbidities",
  "Demographics",
  "Diagnostics",
  "Eosinophils",
  "ExacerbationDetails",
  "FeNO",
  "HeightWeight",
  "IgE",
  "Medication",
  "OCS",
  "Spirometry"
) |>
  set_names(nm = \(x) str_to_lower(x)) |>
  map(\(table) {
    sqlQuery(
      connection,
      paste0(
        'SELECT * FROM QISAR2.usa_njh.NJH_',
        table
      ),
      as.is = T
    )
  })
# str_remove(x, "^Table\\d+_")

# Pending: bronchial thermoplasty; LTOCS

save(njh_data, file = "qisar2_njh_data.RData")

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
load(here::here("datasets", "qisar2_njh_data.RData"))

#Split named list into objects
list2env(njh_data, envir = .GlobalEnv)

# NJH specific ----
# All patients and whether in respective tables
patient_tables <- njh_data |>
  map(\(df) distinct(df, patient_id)) |>
  list_rbind(names_to = "table") |>
  mutate(present = TRUE) |>
  pivot_wider(names_from = table, values_from = present, values_fill = FALSE)

# patient_tables |>
#   filter(when_all(
#     demographics,
#     heightweight,
#     eosinophils,
#     ige,
#     feno,
#     biologics,
#     medication,
#     airq,
#     asthmacontroltest,
#     diagnostics,
#     spirometry,
#     exacerbationdetails
#   ))

# Visits ----
all_observation_dates <- njh_data |>
  map(\(df) select(df, any_of(c("patient_id", "observation_date")))) |>
  list_rbind(names_to = "table") |>
  mutate(observation_date = ymd(observation_date)) |>
  select(patient_id, table, observation_date) |>
  arrange(patient_id, observation_date)

visit_summary <- all_observation_dates |>
  summarise(
    num_of_visits = n(),
    first_visit_date = min(observation_date, na.rm = TRUE),
    last_visit_date = max(observation_date, na.rm = TRUE),
    .by = patient_id
  )

njh_visits <- all_observation_dates |>
  left_join(
    visit_summary,
    by = join_by(patient_id),
    relationship = "many-to-one"
  ) |>
  relocate(contains("visit"), .after = observation_date) |>
  left_join(
    demographics |>
      select(patient_id, observation_date, year_of_birth, age_at_visit) |>
      mutate(observation_date = ymd(observation_date)),
    by = join_by(patient_id, observation_date)
  )

### Numbers of records ----
# nbr_observations <- all_observation_dates |>
#   summarise(
#     # events_by_date = paste0(table, collapse = ", "), When not grouping by table
#     n_obs = n(),
#     .by = c(patient_id, observation_date, table)
#   ) |>
#   # complete(patient_id, observation_date, table, fill = list(n_obs = 0))  #
#   pivot_wider(
#     id_cols = c(patient_id, observation_date),
#     names_from = table,
#     values_from = n_obs,
#     values_fill = 0
#   )

# Look at patient data timeline
# examine <- all_observation_dates |>
#   distinct(patient_id) |>
#   slice_sample(n = 1) |>
#   pull(patient_id)
#
# all_observation_dates |>
#   filter(patient_id %in% examine) |>
#   summarise(
#     # events_by_date = paste0(table, collapse = ", "), When not grouping by table
#     n_obs = n(),
#     .by = c(patient_id, observation_date, table)
#   ) |>
#   ggplot() +
#   geom_point(aes(x = observation_date, y = table, colour = table, size = n_obs))
#
# all_examine_data <- njh_data |>
#   map(\(df) filter(df, patient_id %in% examine))

# Patients ----
# "current age" provided by NJH in earlier transfers, understood to be age at extraction
# "year of birth" = year of extraction - current age = "Inferred" under year_of_birth_source
# "year of birth" provided in later transfers = "Actual" or "Actual_Supplement" under year_of_birth_source
# "age at visit" = visit date - year of birth
# Ignore ethnicity, occupation, and smoking details for now
# "asthma_onset" age not available for most NJH patients
njh_patients <- demographics |>
  remove_empty_check_fields() |>
  # right_join(index_dates, by = join_by(patient_id)) |>
  # relocate(index_date, .after = last_visit_date) |>
  mutate(
    observation_date = ymd(observation_date),
    date_of_birth = make_date(year = year_of_birth, month = 7, day = 2),
    gender = structure(
      factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
      label = "Gender"
    ),
    smoking_status = structure(
      factor(
        smoke_stat,
        levels = c(1, 2, 3),
        labels = c("Current", "Ex-Smoker", "Never")
      ),
      label = "Smoking Status"
    ),
    country = structure("USA", label = "Country"),
    setting = structure("USA - NJH", label = "Setting"),
    site = structure("USA-NJH", label = "Site"),
    data_capture_system = structure("Bespoke", label = "Data capture system"),
  ) |>
  arrange(patient_id, observation_date) |>
  # Exclude patients for now until NJH confirms
  filter_out(any(!is.na(multiplegender_check)), .by = patient_id) |>
  select(
    -multiplegender_check,
    -contains("ethnicity"),
    -contains("smoke_"),
    -contains("occupation")
  ) |>
  distinct(patient_id, gender, date_of_birth, smoking_status)

get_dupes(demographics_cln, patient_id) |> view()

# Multiple gender check:
# US_NJH_11230 duplicated records with 1 and 2
# US_NJH_11771 duplicated records with 1 and 2, then a block of 1s, last records are 2s
# US_NJH_15221 all 1s except a single visit records a 2
# US_NJH_24302 almost all 2 except 4 records with 1
# US_NJH_24324 all 2 except most recent 1 after 5 years
# US_NJH_26887 all 2 except a couple of more recent visits as 1
# US_NJH_29342 duplicated records with 1 and 2

# Look for change in fixed variables across visits
# tochk <- demographics |>
#   mutate(
#     across(
#       all_of(c(
#         "gender"
#       )),
#       \(x) x != dplyr::lag(x, default = first(x)),
#       .names = "{.col}_changed"
#     ),
#     any_changed = if_any(ends_with("_changed")),
#     .by = patient_id
#   )

# Lifestyle ----
# Flag outlying values
is_outlier <- function(
  x,
  method = c("iqr", "3sd"),
  low = TRUE,
  high = TRUE
) {
  method <- match.arg(method)
  # 3sd assumes normality, is influenced by outliers, and has small sample problems
  # iqr robust to outliers and normality and works with small samples

  if (!low && !high) {
    stop("At least one of `low` or `high` must be TRUE.")
  }

  if (method == "iqr") {
    q <- quantile(x, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    spread <- 1.5 * (q[2] - q[1])
    lower <- q[1] - spread
    upper <- q[2] + spread
  } else {
    centre <- mean(x, na.rm = TRUE)
    spread <- 3 * sd(x, na.rm = TRUE)
    lower <- centre - spread
    upper <- centre + spread
  }

  is_low <- if (low) x < lower else FALSE
  is_high <- if (high) x > upper else FALSE

  is_low | is_high
}

njh_lifestyle <- heightweight |>
  select(patient_id, observation_date, wt = weight_kg, ht = height_m) |>
  # Range checks on height and weight, as in other sites
  clean_lifestyle() |> # Converts wildly out of range values to missing
  # Remove outliers within patients (may be mis-entered) using IQR (small-sample)
  mutate(
    weight_clean = if_else(
      is_outlier(weight, method = "iqr") %in% TRUE,
      NA_real_,
      weight
    ),
    height_clean = if_else(
      is_outlier(height, method = "iqr") %in% TRUE,
      NA_real_,
      height
    ),
    .by = patient_id
  ) |>
  # Reduce multiple "clean" measurements on same day as average
  summarise(
    weight = mean(weight_clean, na.rm = TRUE),
    height = mean(height_clean, na.rm = TRUE),
    # n_readings = n(),
    .by = c(patient_id, observation_date)
  ) |>
  # Get rid of "not-a-number", convert to missing
  mutate(across(c(weight, height), \(x) if_else(is.nan(x), NA_real_, x))) |>
  mutate(
    weight = structure(weight, label = "Weight, kg"),
    height = structure(height, label = "Height, m"),
    bmi = structure(
      weight / (height^2),
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
  ) |>
  full_join(
    demographics |>
      select(patient_id, visit_date = observation_date, smoke_stat) |>
      mutate(
        smoke = factor(
          smoke_stat,
          levels = c("1", "2", "3"),
          labels = c("Current smoker", "Ex-smoker", "Never smoked")
        )
      ) |>
      clean_smoking() |>
      rename(observation_date = visit_date),
    by = join_by(patient_id, observation_date)
  ) |>
  mutate(
    observation_date = ymd(observation_date),
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
    observation_date,
    weight,
    height,
    bmi,
    obese,
    smoking_status
  )

# Comorbidities ----
comorbidities_cln <- comorbidities |>
  rename_with(~ str_replace(., "comorbidity_", ""), everything()) |>
  mutate(
    across(contains("date"), ymd),
    name = case_when(
      str_detect(name, "Allergic Rhinitis") ~ "Allergic rhinitis",
      str_detect(name, "Chronic Rhinosinusitis") ~ "Chronic rhinosinusitis",
      str_detect(name, "Diabetes mellitus") ~ "Diabetes",
      str_detect(name, "Nasal Polyps") ~ "Nasal polyps",
      .default = name
    ),
    status = factor(
      status,
      levels = c(
        "Active",
        "Auto Resolved",
        "Current",
        "Deleted",
        "Denied",
        "Never",
        "Past",
        "Resolved",
        "Transitioned"
      ),
      labels = c(
        "Current",
        "Past",
        "Current",
        "Past",
        "Past",
        "Never",
        "Past",
        "Past",
        "Past"
      )
    )
  ) |>
  arrange(patient_id, name, observation_date) |>
  mutate(
    change = any(status != dplyr::lead(status, default = last(status))),
    change_error = case_when(
      n() == 1 ~ "one record",
      row_number() == 1 ~ "ok",
      dplyr::lag(status) == "Never" & status == "Never" ~ "ok",
      dplyr::lag(status) %in%
        c("Never", "Past", "Current") &
        status == "Current" ~ "ok",
      dplyr::lag(status) %in% c("Past", "Current") & status == "Past" ~ "ok",
      dplyr::lag(status) %in%
        c("Past", "Current") &
        status == "Never" ~ "error",
      dplyr::lag(status) == "Never" & status == "Past" ~ "error"
    ),
    .by = c(patient_id, name)
  )

# Look at records with change
comorbidities_cln |> filter(change) |> view()
# TO DO
# Determine direction of change, later visit --> earlier visit:
# current --> never = error
# past --> never = error

# never --> past = error

# current --> past = ok
# past --> past = ok

# current --> current = ok
# past --> current = ok
# never --> current = ok,

# never --> never = ok
# Could collapse to vector of status and check patterns to flag patients with errors

# Allergens ----

# HCRU ----

# Exacerbations ----
# For now assume all listed are exacerbation - SCORING IN PROGRESS
exac_cln <- exacerbationdetails |>
  mutate(
    exac_date = ymd(exac_steroid_startdate),
  ) |>
  # Exclude missing start dates - Exclude at patient level for now
  filter_out(any(is.na(exac_date)), .by = patient_id) |>
  # Add visit dates
  left_join(
    visit_dates |> select(patient_id, first_visit_date),
    by = join_by(patient_id)
  ) |>
  relocate(contains("visit"), .after = patient_id) |>
  # Duplicated exacerbation dates
  select(patient_id, contains("visit"), exac_date) |>
  distinct() |>
  arrange(patient_id, exac_date)

### Baseline ----
bl_exac_num <- exac_cln |>
  # Add query date and select eligible patients
  inner_join(
    index_dates |>
      # Keep eligible patient start dates
      filter(is.na(exclusions)) |>
      select(patient_id, index_date),
    by = join_by(patient_id)
  ) |>
  mutate(
    followup_to_date = interval(first_visit_date, index_date),
    followup_to_date_mos = followup_to_date / months(1),
    # Only include exacerbations if followup window at least 11 months
    followup_ok = followup_to_date_mos >= 11,
    # Search looking back 12 months
    search_start = index_date %m-% months(12),
    search_window = interval(search_start, index_date),
    # Count exacerbations if enough followup and occurred within 12 months of index
    include_exac = followup_ok & exac_date %within% search_window
  ) |>
  # Only count among patients with sufficient followup
  filter(followup_ok) |>
  summarise(bl_exac_num = sum(include_exac), .by = patient_id)


# Asthma control ----
asthma_control <- asthmacontroltest |>
  clean_names() |>
  # Exclude invalid values for available ACT total score <5 or >25
  filter_out(!is.na(act_overallscore_check)) |>
  mutate(
    chk = rowSums(across(c(
      act1_activitylimitation,
      act2_shortnessofbreath,
      act3_nighttimesymptoms,
      act4_rescueinhaleruse,
      act5_overallasthmacontrolrating
    ))),
    act_status = factor(
      act_status,
      levels = c("Good Control", "Not Well Controlled", "Poorly Controlled"),
      labels = c("Well Controlled", "Partly Controlled", "Uncontrolled")
    )
  ) |>
  relocate(chk, .after = act_overall_score) |>
  select(
    patient_id,
    ac_date = observation_date,
    asthma_control = act_status,
    act_score = act_overall_score
  )

# airq2 <- airq |>
#   clean_names() |>
#   # Exclude invalid values
#   filter_out(!is.na(airq_overallscore_check)) |>
#   mutate(
#     chk = rowSums(across(matches("^airq\\d+_")), na.rm = T),
#     chk2 = if_any(matches("^airq\\d+_"), is.na),
#     airq_status = factor(
#       airq_status,
#       levels = c(
#         "Well Controlled",
#         "Not Well Controlled",
#         "Very Poorly Controlled"
#       ),
#       labels = c("Well Controlled", "Partly Controlled", "Uncontrolled")
#     )
#   ) |>
#   relocate(chk, chk2, .after = airq_overall_score)
# select(
#   patient_id,
#   ac_date = observation_date,
#   asthma_control = airq_status,
#   airq_score = airq_overall_score
# )

# Spirometry ----
# TO DO - figure out how/where to summarise average for multiple records on same day
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
clean_spirometry <- function(
  df = spirometry_labelled,
  fev1_range_low = 0.15,
  fev1_range_high = 8,
  fvc_range_low = 0.15,
  fvc_range_high = 8
) {
  df |>
    # Drop rows with missing measurement date is missing
    filter_out(is.na(observation_date)) |>
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

njh_spirometry <- spirometry |>
  rename_with(~ str_remove(., "reported_"), contains("reported")) |>
  # remove_empty_check_fields() |>
  # Drop perc_pred fields, as we'll recalculate these
  select(-contains("percpred")) |>
  ISARanalysis::clean_spirometry(fvc_range_high = 10) |> #Two records with FVC >8 but <10
  select(
    patient_id,
    date,
    fev1,
    fev1_flag,
    fvc,
    fvc_flag,
    fev1_fvc_ratio,
    fev1_fvc_ratio_cat
  ) |>
  arrange(patient_id, observation_date) |>
  # Exclude duplicates
  distinct()

# Checks
# Merged age and height values are correct
# Range, low = 0.15, high = 8
# FEV1 < FVC

# Questions
# Are the percent predicted calculated and reported or have we calculated? --> All reported...

### Baseline ----
bl_spiro <- spirometry_cln |>
  mutate(spirometry_date = ymd(observation_date)) |>
  rename_with(
    ~ str_replace(.x, "visit_|spirometry_", ""),
    .cols = any_of("spirometry_date")
  ) |>
  inner_join(
    index_dates |>
      distinct() |>
      mutate(
        # https://lubridate.tidyverse.org/reference/mplus.html
        search_start = index_date %m-% months(12),
        search_end = index_date, # query date (inclusive)
        search_interval = interval(search_start, search_end)
      ) |>
      select(patient_id, contains("search_")),
    by = join_by(patient_id),
    relationship = "many-to-one"
  ) |>
  # Test results within search interval
  filter(
    date %within% search_interval
  ) |>
  arrange(patient_id, desc(date)) |>
  mutate(
    "spirometry_date" := date,
    "spirometry_gap" := interval(date, search_end) / months(1)
  ) |>
  slice_head(n = 1, by = patient_id)

# Biomarkers ----
# Process to get in same
## IgE ----
# All results in kU/L
ige_cln <- ige |>
  # Following excludes source_column to de-duplicate values
  distinct(
    patient_id,
    observation_date,
    result_operator,
    ige_value,
    unit,
    .keep_all = T
  ) |>
  mutate(
    test = "ige",
    result_date = ymd(observation_date),
    ige_value = as.numeric(ige_value),
    # Some values recorded as >5000 or <2, truncate all the same way
    result = case_when(
      ige_value <= 2 ~ 2,
      ige_value >= 5000 ~ 5000,
      .default = round(ige_value, 0)
    )
  ) |>
  # Excludes 2 duplicates and selects key variables
  distinct(patient_id, test, result_date, result, unit) |>
  arrange(patient_id, result_date)

## FeNO ----
# All results in ppb
feno_cln <- feno |>
  # Drop missing results
  filter_out(is.na(feno_value)) |>
  # Only single source "Exhaled Nitric Oxide"
  select(-source_column) |>
  mutate(
    test = "feno",
    result_date = ymd(observation_date),
    feno_value = as.numeric(feno_value),
    # Some values recorded as >300 or <5, truncate all the same way
    result = case_when(
      feno_value <= 5 ~ 5,
      feno_value >= 300 ~ 300,
      .default = round(feno_value, 0)
    )
  ) |>
  # Excludes any duplicates and selects key variables
  distinct(patient_id, test, result_date, result, unit) |>
  arrange(patient_id, result_date)

# feno_cln |>
#   mutate(
#     n_obs = n(),
#     lag = interval(result_date, dplyr::lag(result_date)) / months(1),
#     mean_lag = mean(lag, na.rm = T),
#     .by = patient_id
#   ) |>
#   summarise(
#     n_obs = mean(n_obs),
#     mean_lag = mean(mean_lag, na.rm = T)
#   )

## EOS ----
# DF 24/09/26 - confirmed okay to combine eos source values as "all are bec"
beos_cln <- eosinophils |>
  # Exclude sputum samples, not used for ISAR
  filter_out(sample_type == "Sputum") |>
  mutate(
    test = "beos",
    result_date = ymd(observation_date),
    eos_value = as.numeric(eos_value),
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
  # Excludes any duplicates and selects key variables
  distinct(patient_id, test, result_date, result, unit) |>
  arrange(patient_id, result_date)

# Combine biomarkers for querying by date and test
njh_biomarkers <- bind_rows(
  beos_cln,
  ige_cln,
  feno_cln
) |>
  arrange(patient_id, result_date)

rm(list = c("eosinophils", "beos_cln", "feno", "feno_cln", "ige", "ige_cln"))

# Medications ----
# Do we want to include nasal steroids?
njh_medications <- medication |>
  rename_with(~ str_replace(., "medication_", ""), everything()) |>
  rename(therapy = category, name = class) |>
  mutate(
    across(contains("date"), ymd),
    # Start date is the same as observation date for 90% of records, use if start_date is missing
    start_date = if_else(
      is.na(start_date),
      observation_date,
      start_date
    ),
    # Replace implausible date
    end_date = if_else(
      end_date == "1900-01-01",
      NA,
      end_date
    )
  ) |>
  arrange(patient_id, start_date) |>
  # Exclude non-asthma uses of medications
  filter_out(reason_nonasthma == "1") |>
  # Exclude if start and end date are missing
  filter_out(when_all(
    is.na(start_date),
    is.na(end_date)
  )) |>
  # End date is before start date
  filter_out(start_date > end_date) |>
  mutate(
    tx_months = interval(start_date, end_date) /
      months(1)
  )

# Biologics ----
# Discontinuation date is between 1 and 39 days after the end_date, can probably ignore
njh_biologics <- biologics |>
  left_join(njh_visits, by = join_by(patient_id)) |>
  relocate(contains("visit"), .after = patient_id) |>
  rename_with(~ str_remove(.x, "medication_"), contains("medication")) |>
  rename_with(~ str_replace(.x, "biologic_", "bx_"), contains("biologic")) |>
  mutate(
    across(contains("date"), ymd),
    bx_name_short = factor(
      bx_class,
      levels = c(
        "Depemokimab",
        "Tezepelumab",
        "Benralizumab",
        "Dupilumab",
        "Reslizumab",
        "Mepolizumab",
        "Omalizumab"
      )
    ),
    bx_type = case_when(
      bx_name_short == "Benralizumab" ~ "Anti-IL5R",
      bx_name_short == "Dupilumab" ~ "Anti-IL4Ralpha",
      .default = bx_category
    ) |>
      factor(
        levels = c(
          "Anti-IgE",
          "Anti-IL4Ralpha",
          "Anti-IL5",
          "Anti-IL5R",
          "Anti-TSLP"
        )
      ),
    start_date = case_when(
      # 1900-01-01
      start_date == ymd("1900-01-01") ~ NA,
      # if start date missing and is a prescription: assume prescription date
      is.na(start_date) & historicalorreported == 0 ~ observation_date,
      # set to missing if medication reported by patient and start date same as observation date: e.g. GS - patient comes to NJH, reports taking bx, and somehow they enter start date as "today's date"
      # historicalorreported == 1 & start_date == observation_date ~ NA,
      .default = start_date
    ),
    end_date = case_when(
      # 1900-01-01
      end_date == ymd("1900-01-01") ~ NA,
      # Future dates
      # end_date > last_visit_date ~ NA, # Keep future dates for now
      .default = end_date
    ),
    # tx_months = interval(start_date, end_date) / months(1)
  ) |>
  relocate(bx_name_short, bx_type, .after = end_date) |>
  # relocate(tx_months, .after = end_date) |>
  relocate(contains("visit_date"), .after = patient_id) |>
  arrange(patient_id, start_date) |>
  rename(
    bx_start_date = start_date,
    bx_end_date = end_date
  ) |>
  select(-bx_category, -bx_class)

# Save datasets ----
save(
  njh_visits,
  njh_patients,
  njh_lifestyle,
  njh_comorbidities,
  njh_allergens,
  njh_hcru,
  njh_exacerbations,
  njh_asthma_control,
  njh_spirometry,
  njh_biomarkers,
  njh_medications,
  njh_biologics,
  file = here::here("datasets", "qisar2_njh_srd.Rdata")
)
