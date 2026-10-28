library(RODBC)
library(tidyverse)
library(janitor)


# To Do ----
# Medications
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

njh_data <- c(
  "AIRQ",
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

# Pending: allergens; bronchial thermoplasty; LTOCS

save(njh_data, file = "qisar2_njh_data.RData")

close(connection)

# Functions ----
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
load("qisar2_njh_data.RData")

#Split named list into objects
list2env(njh_data, envir = .GlobalEnv)

# Dictionary ----
# variable_list <- njh_data |>
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
# write_xlsx(variable_list, "qisar2_njh_data_dictionary.xlsx")

# Process data
# Patients ----
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

## Observation and visit dates ----
all_observation_dates <- njh_data |>
  map(\(df) select(df, any_of(c("patient_id", "observation_date")))) |>
  list_rbind(names_to = "table") |>
  mutate(observation_date = ymd(observation_date)) |>
  select(patient_id, table, observation_date) |>
  arrange(patient_id, observation_date)

### Visit dates ----
visit_dates <- all_observation_dates |>
  summarise(
    first_visit_date = min(observation_date, na.rm = TRUE),
    last_visit_date = max(observation_date, na.rm = TRUE),
    .by = patient_id
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
examine <- all_observation_dates |>
  distinct(patient_id) |>
  slice_sample(n = 1) |>
  pull(patient_id)

all_observation_dates |>
  filter(patient_id %in% examine) |>
  summarise(
    # events_by_date = paste0(table, collapse = ", "), When not grouping by table
    n_obs = n(),
    .by = c(patient_id, observation_date, table)
  ) |>
  ggplot() +
  geom_point(aes(x = observation_date, y = table, colour = table, size = n_obs))

all_examine_data <- njh_data |>
  map(\(df) filter(df, patient_id %in% examine))

# Biologics ----
# Discontinuation date is between 1 and 39 days after the end_date, can probably ignore
bx_cln <- biologics |>
  left_join(visit_dates, by = join_by(patient_id)) |>
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
  relocate(bx_name_short, bx_type, .before = bx_type) |>
  # relocate(tx_months, .after = end_date) |>
  relocate(contains("visit_date"), .after = patient_id) |>
  mutate(
    exclusions = case_when(
      any(
        year(observation_date) < 2016
      ) ~ "Started biologic before 2016",
      all(
        reason_nonasthma == 1
      ) ~ "Prescribed a biologic only for non-asthma reason",
      any(
        reason_nonasthma == 1
      ) ~ "Prescribed a biologic for non-asthma reason",
      # all(
      #   historicalorreported == 1
      # ) ~ "All patient reported dates",
      # No dates provided
      all(
        is.na(start_date) & is.na(end_date)
      ) ~ "No start or end dates provided",
      all(is.na(start_date)) ~ "No start dates provided",
      # Exclude if medication reported by patient and start date same as observation date
      any(
        start_date == observation_date & historicalorreported == 1
      ) ~ "Unreliable start date - Reported", # Exclude for now
      # January first, very unlikely and probably imputed by NJH
      any(
        month(start_date) == 1 & day(start_date) == 1
      ) ~ "Unreliable start date - Imputed",
      # any(is.na(start_date)) ~ "Missing a biologic initiation date",
      .default = NA_character_
    ) |>
      factor(
        levels = c(
          "Started biologic before 2016",
          "Prescribed a biologic only for non-asthma reason",
          "Prescribed a biologic for non-asthma reason",
          # "All patient reported dates", # catch patients here or with missing start dates
          "No start or end dates provided",
          "No start dates provided",
          "Unreliable start date - Reported",
          "Unreliable start date - Imputed"
          # "Missing a biologic initiation date"
        )
      ),
    .by = patient_id
  ) |>
  arrange(patient_id, start_date) |>
  rename(
    bx_start_date = start_date,
    bx_end_date = end_date
  ) |>
  select(-bx_category, -bx_class)

# Initial exclusions
bx_cln |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions)

## Exclusions and initiation dates ----
# Use first biologic initiation date (n=3663 patients)
# 20261015 - refined to use 3 months' follow-up as minimum for Obj1

# Defines all eligible biologic patients for further analysis
patients <- bx_cln |>
  # Keep one row per excluded patient
  filter_out(
    when_all(!is.na(exclusions), row_number() > 1),
    .by = patient_id
  ) |>
  # Keep earliest start date, keeps ties, ie started multiple on same day
  slice_min(bx_start_date, by = patient_id) |>
  distinct(
    patient_id,
    first_visit_date,
    last_visit_date,
    exclusions,
    index_date = bx_start_date,
    bx_name_short,
  ) |>
  mutate(
    country = "USA",
    setting = "NJH",
    data_capture_system = "Bespoke",
    exclusions = case_when(
      is.na(exclusions) &
        index_date <
          first_visit_date ~ "Started first biologic before first observation",
      is.na(exclusions) & n() > 1 ~ "Duplicated start date - multiple biologic",
      is.na(exclusions) &
        interval(index_date, last_visit_date) / months(1) <
          3 ~ "<3 months post-initiation follow-up",
      .default = exclusions
    ) |>
      factor(
        levels = c(
          "Started biologic before 2016",
          "Prescribed a biologic only for non-asthma reason",
          "Prescribed a biologic for non-asthma reason",
          "No start or end dates provided",
          "No start dates provided",
          "Unreliable start date - Reported",
          "Unreliable start date - Imputed",
          "Started first biologic before first observation",
          "Duplicated start date - multiple biologic",
          "<3 months post-initiation follow-up"
        )
      ),
    exclusions_table = case_when(
      exclusions ==
        "Prescribed a biologic only for non-asthma reason" ~ "Non-biologic naive",
      exclusions ==
        "Prescribed a biologic for non-asthma reason" ~ "Non-biologic naive",
      exclusions ==
        "Started first biologic before first observation" ~ "Non-biologic naive",
      exclusions ==
        "Started biologic before 2016" ~ "Started biologic before 2016",
      exclusions ==
        "<3 months post-initiation follow-up" ~ "<3 months post-initiation follow-up",
      is.na(exclusions) ~ "Eligible",
      .default = "Data queries pending"
    ) |>
      factor(
        levels = c(
          "Started biologic before 2016",
          "Non-biologic naive",
          "Data queries pending",
          "<3 months post-initiation follow-up",
          "Eligible"
        )
      ),
    .by = patient_id
  ) |>
  # Keep one row per excluded patient
  filter_out(
    when_all(!is.na(exclusions), row_number() > 1),
    .by = patient_id
  )
# Set values to missing for patients with an exclusion
# mutate(
#   across(c(bx_name, bx_type, index_date), \(x) {
#     replace(x, which(!is.na(exclusions) & exclusions != "E"), NA)
#   })
# )

patients |>
  tabyl(exclusions) |>
  adorn_totals()

patients |>
  tabyl(exclusions) |>
  adorn_pct_formatting() |>
  kable() |>
  kable_styling()

eligible <- patients |>
  filter(is.na(exclusions)) |>
  pull(patient_id)

index_dates <- patients |>
  # Restrict to eligible for obj1
  filter(is.na(exclusions)) |>
  select(patient_id, index_date)

## Switch dates ----
# 20261006 Use minimum start date per biologic - no stopping information
# switch_dates <- bx_cln |>
#   # Selects eligible patients
#   right_join(index_dates, by = join_by(patient_id)) |>
#   relocate(index_date, .after = last_visit_date) |>
#   arrange(patient_id, start_date, end_date) |>
#   # Exclude start dates 24 + 3 months after initiation
#   # filter_out(start_date > index_date %m+% months(28)) |>
#   slice_min(
#     start_date,
#     by = c(patient_id, bx_name),
#     with_ties = T,
#     na_rm = T
#   ) |>
#   distinct(
#     patient_id,
#     first_visit_date,
#     last_visit_date,
#     index_date,
#     bx_name,
#     bx_type,
#     start_date
#   ) |>
#   mutate(n_obs = n(), .by = patient_id)

# 20261015 Use spells approach rather than first start
bx_episodes <- bx_cln |>
  # Selects eligible patients
  right_join(index_dates, by = join_by(patient_id)) |>
  relocate(
    first_visit_date,
    last_visit_date,
    index_date,
    .after = last_visit_date
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    change_med = bx_name_short !=
      dplyr::lag(bx_name_short, default = first(bx_name_short)),
    change = if_else(
      bx_name_short !=
        dplyr::lag(bx_name_short, default = first(bx_name_short)),
      1,
      0
    ),
    bx_nbr = cumsum(as.integer(change)) + 1,
    # Base change on biologic names
    any_change = any(change == 1),
    .by = c(patient_id)
  ) |>
  summarise(
    any_change = first(any_change),
    bx_start_date = min(bx_start_date, na.rm = T),
    bx_name_short = first(bx_name_short),
    bx_type = first(bx_type),
    .by = c(patient_id, bx_nbr)
  ) |>
  mutate(
    bx_ct = max(bx_nbr),
    # Outcome based on number of spells
    outcome = case_when(
      bx_ct == 1 ~ "Continuer",
      bx_ct > 1 ~ "Switch"
    ),
    .by = patient_id
  ) |>
  left_join(
    patients |>
      select(
        patient_id,
        data_capture_system,
        country,
        setting,
        first_visit_date,
        last_visit_date,
        index_date
      ),
    by = join_by(patient_id)
  )

bx_episodes |>
  distinct(patient_id, .keep_all = T) |>
  tabyl(outcome) |>
  adorn_totals()

## Switch dates ----
# Keep second record for switchers
switch_dates <- spells |>
  filter(
    spell_ct > 1,
    spell_id == 2
  ) |>
  left_join(index_dates, by = join_by(patient_id)) |>
  mutate(
    followup_mos = interval(switch_date, last_visit_date) / months(1),
    eligible_fup = followup_mos >= 11,
    post_switch_date = switch_date %m+% months(12)
  ) |>
  select(
    patient_id,
    first_visit_date,
    last_visit_date,
    index_date,
    switch_date,
    post_switch_date,
    followup_mos,
    eligible_fup
  ) |>
  filter(eligible_fup)


# US_NJH_14691 exclude for now - looks like they had a historic reported use of Oma

# 615 patients with only 1 medicine
# NOTE - May have gaps in treatment on same med
# Need to define stops
# Restrict to first 2 years
single_med <- switch_dates |>
  filter(all(change == 0), .by = patient_id) |>
  arrange(start_date, observation_date, end_date) |>
  mutate(
    # Ignore starts 24 + 3 months after first initiation
    ignore_record = !is.na(start_date) &
      (start_date > index_date %m+% months(27))
  )
summarise(
  n_obs = n(),
  first_visit_date = first(first_visit_date),
  last_visit_date = last(last_visit_date),
  start_date = min(start_date, na.rm = T),
  end_date1 = last(end_date), #Don't exclude missing
  end_date2 = max(end_date, na.rm = F), #Don't exclude missing
  end_date3 = max(end_date, na.rm = T), #Exclude missing
  bx_name = paste(unique(na.omit(bx_name)), collapse = ", "),
  bx_type = paste(unique(na.omit(bx_type)), collapse = ", "),
  .by = patient_id
) |>
  mutate(
    stop1 = !is.na(end_date1),
    stop2 = !is.na(end_date2),
    stop3 = !is.na(end_date3),
    censor_date = start_date %m+% months(12)
  )

# Single med spells
single_med_results <- switch_dates |>
  filter(all(change == 0), .by = patient_id) |>
  arrange(start_date, observation_date, end_date) |>
  # mutate(
  #   # Ignore starts 24 + 3 months after first initiation
  #   ignore_record = !is.na(start_date) &
  #     (start_date > index_date %m+% months(27))
  # ) |>
  group_by(patient_id) |>
  group_split() |>
  map(\(patient) process_patient_v2(patient), .progress = T)

single_med_spells <- single_med_results |>
  map("spells") |>
  list_rbind()
write_xlsx(single_med_spells, "single_med_spells_20261001.xlsx")
dropped_rows_all <- single_med_results |>
  map("dropped_rows_missing_all_dates") |>
  list_rbind()
dropped_rows_start <- single_med_results |>
  map("dropped_rows_missing_start_dates") |>
  list_rbind()


# Discontinuation date only different for 49 patients (1.5%)
# stop1 = stop2 = stop3 for 110 patients, meaning
# end_date approach doesn't matter if only 1 row
# 75 (52%) of patients with 1 row have missing end_date
# 68 of patients with 1 row have all the same end_date summarised results

# TO DO - define stoppers as end date before observation date?

switchers <- switch_dates |>
  filter(any(change != 0), .by = patient_id) |>
  summarise(
    start_date = min(start_date, na.rm = T),
    bx_type = first(bx_type),
    .by = c(patient_id, bx_name)
  )

## Bx spells ----
# Row missing all dates (among patients with some dates) - check that biologic is already listed for patient
missing_all_dates <- bx_cln |>
  filter(is.na(exclusions)) |>
  filter(when_all(is.na(start_date), is.na(end_date))) |>
  distinct(patient_id, bx_name, bx_type, .keep_all = T)
# Note - all of these records are reported by patient

# Check treatment rows with only end dates based on missing start dates are already in spells
missing_start_dates <- bx_cln |>
  filter(is.na(exclusions)) |>
  filter(when_all(is.na(start_date), !is.na(end_date))) |>
  select(patient_id, bx_name, bx_type, start_date, end_date)

# Treatment spells defined excluding rows with missing start dates
# NB - Need to decide how to approach concurrent treatments
spells <- bx_cln |>
  filter(is.na(exclusions)) |>
  # Exclude rows with missing start dates
  filter_out(is.na(start_date)) |>
  arrange(patient_id, start_date) |>
  mutate(
    bx_interval = interval(start_date, end_date),
    change_med = bx_name != dplyr::lag(bx_name, default = first(bx_name)),
    date_overlap = start_date %within% dplyr::lag(bx_interval),
    change = if_else(
      bx_name != dplyr::lag(bx_name, default = first(bx_name)),
      1,
      0
    ),
    spell_id = cumsum(as.integer(change)) + 1,
    .by = c(patient_id)
  ) |>
  relocate(
    bx_interval,
    change_med,
    date_overlap,
    change,
    spell_id,
    .after = end_date
  ) |>
  mutate(
    spell_start_date = min(start_date, na.rm = TRUE),
    spell_end_date = max(end_date, na.rm = TRUE),
    .by = c(patient_id, spell_id)
  ) |>
  mutate(
    spell_ct = max(spell_id),
    .by = patient_id
  ) |>
  # summarise(start_date = min(start_date, na.rm = TRUE),
  #           end_date = max(end_date, na.rm = TRUE),
  #           .by = c(patient_id, first_visit_date, last_visit_date, bx_name, bx_type, spell_id, exclusions)) |>
  arrange(patient_id, start_date) |>
  relocate(contains("spell"), .after = end_date) |>
  distinct(
    patient_id,
    first_visit_date,
    last_visit_date,
    bx_name,
    bx_type,
    spell_id,
    spell_ct,
    spell_start_date,
    spell_end_date
  )

# Check 1 - row only lists biologic, no dates, so see if biologic already listed, exclude patient if not
# Patient reported drug already in list
ok_to_include_patient <- spells |>
  left_join(
    missing_all_dates |>
      select(patient_id, bx_name, start_date, end_date, historicalorreported) |>
      mutate(ok = TRUE),
    by = join_by(patient_id, bx_name),
    relationship = "many-to-one"
  ) |>
  filter(any(ok), .by = patient_id)
# Patient reported drug not in list
# Note - These look to be cases where either a different drug is reported by patient
# or the same drug is already listed but without a start date so it's been excluded
# from the spells dataset
patient_to_exclude <- missing_all_dates |>
  anti_join(
    spells,
    by = join_by(patient_id, bx_name)
  ) |>
  distinct(patient_id) |>
  pull(patient_id)
checking <- bx_cln |>
  filter(patient_id %in% to_exclude)

# Check 2 - biologic listed with only an end date. See if biologic is listed
# and date overlaps with existing interval.
# Test rows with missing start dates against existing spells to see if they're accounted for already
ok_to_include_patient2 <- missing_start_dates |>
  inner_join(
    spells,
    by = join_by(
      patient_id,
      bx_name,
      between(end_date, spell_start_date, spell_end_date)
    )
  )
patient_to_exclude2 <- missing_start_dates |>
  anti_join(ok_to_include_patient2, by = join_by(patient_id, bx_name)) |>
  pull(patient_id)

## Clean spells ----
spells_cln <- spells |>
  filter_out(patient_id %in% patient_to_exclude) |>
  filter_out(patient_id %in% patient_to_exclude2) |>
  mutate(
    spell_end_date = if_else(
      is.infinite(spell_end_date),
      last_visit_date,
      spell_end_date
    ),
    tx_months = round(interval(spell_start_date, spell_end_date) / months(1), 2)
  ) |>
  mutate(
    spell_overlap = int_overlaps(
      interval(spell_start_date, spell_end_date),
      dplyr::lag(interval(spell_start_date, spell_end_date))
    ),
    spell_overlap_mos = if_else(
      spell_overlap,
      interval(spell_start_date, dplyr::lag(spell_end_date)) / months(1),
      NA
    ),
    spell_overlap_same_bx = if_else(
      spell_overlap & spell_overlap_mos > 0,
      bx_name == dplyr::lag(bx_name),
      NA
    ),
    .by = patient_id
  )
# Note - need to think about overlaps - are these switches or should we look at concurrent treatments?

#### Plot spells ----
# Define biologic order
bx_levels <- c(
  "Depemokimab",
  "Tezepelumab",
  "Benralizumab",
  "Dupilumab",
  "Reslizumab",
  "Mepolizumab",
  "Omalizumab"
)

# Named colours keyed to those levels
bx_colours <- setNames(
  RColorBrewer::brewer.pal(7, "Set2"),
  bx_levels
)

totest <- spells_cln |>
  filter(spell_ct > 1) |>
  distinct(patient_id) |> #n=507 patients with multiple spells
  slice_sample(n = 10) |>
  pull(patient_id)
spells_cln |>
  filter(patient_id %in% totest) |>
  mutate(
    spell_id = as.factor(spell_id),
    spell_end_date_imputed = if_else(spell_end_date == "-Inf", "Yes", "No"),
    spell_end_date = case_when(
      spell_end_date == "-Inf" ~ last_visit_date,
      .default = spell_end_date
    )
  ) |>
  ggplot() +
  geom_segment(
    aes(
      y = spell_id,
      yend = spell_id,
      x = spell_start_date,
      xend = spell_end_date,
      colour = bx_name
    ),
    size = 6,
    position = position_nudge(y = 0.12)
  ) +
  geom_point(
    aes(
      y = spell_id,
      x = spell_end_date,
      fill = bx_name,
      colour = bx_name,
      shape = spell_end_date_imputed
    ),
    size = 6,
    position = position_nudge(y = 0.12),
    show.legend = c(shape = TRUE, fill = FALSE, colour = FALSE)
  ) +
  scale_colour_manual(values = bx_colours) +
  scale_fill_manual(values = bx_colours) +
  scale_y_discrete() +
  scale_shape_manual(values = c(NA, 4)) +
  scale_x_date(
    date_labels = "%Y",
    date_breaks = "year",
    limits = as.Date(c("2004-01-01", "2026-12-31"))
  ) +
  facet_wrap(vars(patient_id), ncol = 1, scales = "free_y")

# Look at start dates by FDA approval + six months
# Check patients on multiple biologics manually

# Demographics ----
# "current age" provided by NJH in earlier transfers, understood to be age at extraction
# "year of birth" = year of extraction - current age = "Inferred" under year_of_birth_source
# "year of birth" provided in later transfers = "Actual" or "Actual_Supplement" under year_of_birth_source
# "age at visit" = visit date - year of birth
# Ignore ethnicity, occupation, and smoking details for now
demographics_cln <- demographics |>
  remove_empty_check_fields() |>
  right_join(index_dates, by = join_by(patient_id)) |>
  relocate(index_date, .after = last_visit_date) |>
  mutate(
    observation_date = ymd(observation_date),
    gender = structure(
      factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
      label = "Gender"
    ),
    date_of_birth = make_date(year = year_of_birth, month = 7, day = 2),
    index_age = floor(interval(date_of_birth, start_date) / years(1)),
    smoking_status = structure(
      factor(
        smoke_stat,
        levels = c(1, 2, 3),
        labels = c("Current", "Ex-Smoker", "Never")
      ),
      label = "Smoking Status"
    )
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
  distinct(patient_id, gender, year_of_birth, smoking_status)

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

### Baseline ----
bl_ac <- asthma_control |>
  mutate(ac_date = ymd(ac_date)) |>
  rename_with(
    ~ str_replace(.x, "ac_|visit_|spirometry_", ""),
    .cols = any_of("ac_date")
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
    "ac_date" := date,
    "ac_gap" := interval(date, search_end) / months(1)
  ) |>
  slice_head(n = 1, by = patient_id)

### Switch date ----
# Closest to switch
### Post-switch ----
# Closest to 12 months after switch

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

### Switch date ----
### Post-switch ----

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
# Medications ----
# Do we want to include nasal steroids?
medication_cln <- medication |>
  rename_with(~ str_replace(., "medication_", ""), everything()) |>
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


# Biomarkers ----
# Process to get in same
## IGE ----
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

## FENO ----
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
biomarkers <- bind_rows(
  beos_cln,
  ige_cln,
  feno_cln
) |>
  arrange(patient_id, result_date)

rm(list = c("eosinophils", "beos_cln", "feno", "feno_cln", "ige", "ige_cln"))

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

spirometry_cln <- spirometry |>
  rename_with(~ str_remove(., "reported_"), contains("reported")) |>
  remove_empty_check_fields() |>
  # Drop perc_pred fields, as we'll recalculate these
  select(-contains("percpred")) |>
  clean_spirometry(fvc_range_high = 10) |> #Two records with FVC >8 but <10
  select(
    patient_id,
    observation_date,
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

### Switch date ----
# Closest to switch
### Post-switch ----
# Closest to 12 months after switch
