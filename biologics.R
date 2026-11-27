# 20260812
# Assemble biologics dataset, SRD no longer separating biologics,
# so reproduce steps here using the full medications table.
# 20261015 - revised to use 3 months' follow-up, not 12
# 20261125 - revised to match objects to bespoke processing

isar_biologics <- isar$medications_labelled |>
  remove_empty_imputation_fields() |>
  filter(med_class == "Biologic" & med_yn == "Yes") |>
  # filter_out(is.na(med_name)) |>
  rename_with(
    .cols = contains("med_"),
    ~ str_replace(., "med_", "bx_")
  ) |>
  rename_with(
    .cols = contains("bx_st"),
    ~ str_replace(., "bx_st", "bx_start_")
  ) |>
  rename_with(
    .cols = contains("bx_en"),
    ~ str_replace(., "bx_en", "bx_end_")
  ) |>
  select(
    patient_id,
    bx_name,
    contains("bx_start_"),
    contains("bx_end_"),
    contains("bx_ongoing"),
    bx_switchreason = bio_switchreason
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
    #Add steps from Biologics_table.R SRD script
    bx_ongoing_imputed_flag = case_when(
      row_number() == n() &
        is.na(bx_end_date) &
        is.na(bx_ongoing) ~ as.character(
        "Last biologic with missing end date, assume ongoing"
      )
    ),
    bx_ongoing = case_when(
      row_number() == n() &
        is.na(bx_end_date) &
        is.na(bx_ongoing) ~ as.character("Ongoing - Maintain"),
      .default = bx_ongoing
    ),
    bx_ongoing_imputed = ifelse(!is.na(bx_ongoing_imputed_flag), 1, 0)
  ) |>
  relocate(bx_name_short, bx_type, .after = bx_name) |>
  # Added distinct 20260812, this is the final step in the biologics_table.R script
  distinct()

# index_dates <- visit_labelled |>
#   slice_min(visit_date, by = patient_id, with_ties = F) |>
#   select(patient_id, first_visit_date = visit_date) |>
#   left_join(
#     meds |>
#       slice_min(bx_start_date, by = patient_id, with_ties = F) |>
#       select(patient_id, bx_start_date),
#     by = join_by(patient_id)
#   ) |>
#   mutate(
#     index_date = coalesce(bx_start_date, first_visit_date)
#   )

## Exclusions ----
isar_bx_cln_exclusions <- isar$patient_labelled |>
  filter_out(country %in% c("Denmark", "UK", "Italy", "Portugal")) |>
  filter_out(str_detect(patient_id, "NJH")) |>
  filter_out(str_detect(patient_id, "SITETEZ")) |>
  remove_empty_imputation_fields() |>
  select(-data_capture_system) |>
  #   Demographics - add age and index_date, confirm the latter
  left_join(
    isar$demographics_labelled |>
      # Can't use clean function for now because no index_age
      select(
        patient_id,
        date_of_birth
      ),
    by = join_by(patient_id),
    relationship = "one-to-one"
  ) |>
  #   Left join to keep patients without biologic data for consort
  left_join(
    isar_biologics |>
      mutate(
        #   Exclude imputed months AND days
        exclude_imp_date = case_when(
          bx_start_date_imputed_flag == "MonthsDays" ~ TRUE,
          bx_end_date_imputed_flag == "MonthsDays" ~ TRUE
        ),
        bx_data = TRUE
      ),
    by = "patient_id",
    relationship = "one-to-many"
  ) |>
  group_by(patient_id) |>
  #   Flag issues
  arrange(bx_start_date) |>
  mutate(
    # Missing bx name but not dates - empty rows dropped above
    err1 = case_when(
      any(bx_data) &
        is.na(bx_name) ~ "missing_bx_name_not_dates"
    ),
    # Missing start date
    err2 = case_when(
      any(bx_data) &
        is.na(bx_start_date) ~ "missing_bx_start_date"
    ),
    # err3 = case_when(
    #   any(bx_flag) &
    #     n() > 1 &
    #     all(is.na(bx_end_date)) ~ "no_end_dates"
    # ),
    # Later biologic starts before end of previous
    err4 = case_when(
      any(bx_data) &
        bx_start_date <
          dplyr::lag(bx_end_date) ~ "later_bx_start_before_prior_end"
    ),
    # Biologic ends before start
    err5 = case_when(
      any(bx_data) &
        bx_end_date < bx_start_date ~ "end_before_start"
    ),
    err6 = case_when(
      any(bx_data) &
        any(
          duplicated(
            pick(bx_start_date) |>
              filter(!is.na(bx_start_date))
          )
        ) ~ "duplicate_bx_start_date"
    ),
    # Duplicate start dates within patient and biologic
    err7 = case_when(
      any(bx_data) &
        any(
          duplicated(
            pick(bx_name, bx_start_date) |>
              filter(!is.na(bx_name) & !is.na(bx_start_date))
          )
        ) ~ "duplicate_bx_start_date_and_name"
    ),
    # Multiple recordings of same biologic (ie doses), ignores possible gaps
    err8 = case_when(
      any(bx_data) &
        n() > 1 &
        n_distinct(bx_name, na.rm = TRUE) == 1 ~ "multiple_same_bx"
    ),
    err9 = case_when(
      any(bx_data) &
        n() > 1 &
        all(is.na(bx_end_date)) ~ ">1 record, no end dates"
    ),
    err10 = case_when(
      any(bx_data) &
        n() > 1 &
        !last(is.na(bx_end_date)) &
        any(is.na(bx_end_date)) ~ ">1 record, missing not-last end date"
    ),
    err11 = case_when(
      last(is.na(bx_end_date)) &
        last(bx_start_date) >
          last_visit_date ~ "missing end date and started last biologic after last visit"
    )
  ) |>
  fill(contains("err"), .direction = "downup") |>
  unite("errors", contains("err"), sep = ", ", na.rm = T, remove = F) |>
  mutate(
    bx_initiation_date = if (all(is.na(bx_start_date))) {
      as.Date(NA)
    } else {
      min(bx_start_date, na.rm = T)
    },
    bx_initiation_year = year(bx_initiation_date),
    bx_initiation_age = interval(ymd(date_of_birth), bx_initiation_date) /
      years(1),
    age_at_baseline = interval(
      ymd(date_of_birth),
      bx_initiation_date %m-% months(12)
    ) /
      years(1),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
    exclusions = case_when(
      all(is.na(bx_data)) ~ "No biologic information",
      all(num_of_visits == 0) ~ "No visits recorded", # means no visit dates for imputation
      country %in% c("Ecuador", "Qatar") ~ "Fewer than two biologics available",
      any(is.na(bx_name)) ~ "Missing biologic name",
      any(is.na(bx_start_date)) ~ "Missing start date",
      any(exclude_imp_date) ~ "Imputed start/end MonthsDays",
      any(errors != "") ~ "Flagged errors",
      age_at_baseline <
        18 ~ "Under-18 at baseline, 12 months prior to initiation",
      any(bx_initiation_year < 2016) ~ "Started biologic before 2016",
      # max(1:n()) > 6 ~ "More than 6 listed treatments",
      follow_up_mos < 3 ~ "<3 months post-initiation follow-up",
      # Add patient-specific data queries here by setting as NA before flagging as eligible:
      patient_id == "CAN_S11_14" ~ NA, # 2017 Teze patient
      patient_id == "CAN_S2_20" ~ NA, # Bad DOB
      patient_id == "BEL_SITE1_363" ~ NA, # Ended initial and switched to next biologic long after last visit
      patient_id == "MEX_SITE13_12282" ~ NA, # DOB in 1899
      .default = "Eligible"
    ) |>
      factor(
        levels = c(
          "No biologic information",
          "No visits recorded",
          "Fewer than two biologics available",
          "Missing biologic name",
          "Missing start date",
          "Imputed start/end MonthsDays",
          # "More than 6 listed treatments",
          "Flagged errors",
          "Under-18 at baseline, 12 months prior to initiation",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up",
          "Eligible"
        )
      ),
    exclusions_table = case_when(
      exclusions == "No biologic information" ~ "No biologic information",
      exclusions ==
        "Under-18 at baseline, 12 months prior to initiation" ~ "Under-18",
      exclusions == "Non-biologic naive" ~ "Non-biologic naive",
      exclusions ==
        "Started biologic before 2016" ~ "Started biologic before 2016",
      exclusions ==
        "<3 months post-initiation follow-up" ~ "<3 months post-initiation follow-up",
      exclusions == "Eligible" ~ "Eligible",
      .default = "Data queries"
    ) |>
      factor(
        levels = c(
          "No biologic information",
          "Under-18",
          "Non-biologic naive",
          "Data queries",
          "Started biologic before 2016",
          "<3 months post-initiation follow-up",
          "Eligible"
        )
      )
    # country = if_else(
    #   is.na(country),
    #   str_extract(patient_id, "^[^_]+"),
    #   country
    # ),
  ) |>
  ungroup()

# Check that index date aligns with first biologic date
bx_cln |>
  mutate(
    index_year = year(index_date),
    index_year_chk = first(bx_start_year),
    .by = patient_id
  ) |>
  filter(index_year != index_year_chk) # empty set

error_summary <- isar_bx_cln_exclusions |>
  distinct(patient_id, .keep_all = T) |>
  filter(errors != "") |>
  select(-errors) |>
  pivot_longer(
    cols = contains("err"),
    names_to = "error_column",
    values_to = "error_flag"
  ) |>
  filter(!is.na(error_flag))
error_summary |>
  count(error_flag, sort = TRUE) |>
  arrange(error_flag) |>
  kableExtra::kable() |>
  kableExtra::kable_styling()
error_summary |>
  count(country, error_flag, sort = TRUE) |>
  arrange(country, error_flag)
error_summary |> view()

### Table ----
# See "20261015_srd_feasibility.R"
# Excluding only NJH, Denmark, and UK matches earlier feasibility numbers
isar_bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions) |>
  adorn_totals() |>
  arrange(!exclusions == "Total")

isar_bx_cln_exclusions |>
  distinct(patient_id, exclusions_table) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  arrange(!exclusions_table == "Total")

## Analysis ----
isar_bx_analysis <- isar_bx_cln_exclusions |>
  # Exclude patients with errors
  filter(
    when_any(
      exclusions == "Eligible"
      # exclusions == "Fewer than two biologics available",
      # exclusions == "Biologic start date before 2016",
      # exclusions == "Less than 3 months of follow-up"
    )
  ) |>
  relocate(data_capture_system = version, .after = patient_id) |>
  select(
    -bx_name,
    -bx_data,
    -contains("err"),
    # -contains("exclu")
  ) |>
  arrange(patient_id, bx_start_date) |>
  # Assign earlier end date as subsequent start date in these cases
  mutate(
    bx_initiation_date = min(bx_start_date),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
    bx_nbr = row_number(),
    bx_ct = n(),
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
      all(bx_ct == 1 & is.na(bx_end_date)) ~ "Continue",
      all(
        bx_ct == 1 &
          bx_end_date > last_visit_date %m+% months(1)
      ) ~ "Continue", # Censored not stop
      # any(bx_ct == 1 & bx_ongoing == "No" & bx_end_date <= last_visit_date) ~ "Stop",
      all(
        bx_ct == 1 &
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
  ) |>
  arrange(patient_id, bx_start_date) |>
  rename(site = setting) |>
  mutate(
    setting = case_when(
      str_detect(site, "UNC") ~ "USA - UNC",
      str_detect(site, "Michigan") ~ "USA - Michigan",
      str_detect(site, "NJH") ~ "USA - NJH",
      str_detect(site, "Tallaght") ~ "IRL - Tallaght",
      str_detect(site, "Beaumont") ~ "IRL - Beaumont",
      .default = country
    )
  )

save(
  isar_biologics,
  isar_bx_cln_exclusions,
  isar_bx_analysis,
  file = here::here("datasets", "biologic_data_isar.RData")
)
