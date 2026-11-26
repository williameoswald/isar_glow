# ELIGIBILITY ----
# Revises cleaned biologic data in respective environment
# Adds exclusions to create "*_bx_cln_exclusions"
# Filters to eligible to creae "*_bx_analysis"
# Uses bx_analysis datasets to get full list of eligible patients with dates

# UK ----
## Exclusions ----
uk_biologics <- uk$uk_biologics
uk_bx_cln_exclusions <- uk$uk_patients |>
  select(patient_id, age_first_visit) |>
  left_join(
    uk_biologics |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    bx_initiation_age = age_first_visit +
      floor(interval(first_visit_date, min(bx_start_date)) / years(1)),
    age_at_baseline = age_first_visit +
      floor(
        interval(first_visit_date, min(bx_start_date) %m-% months(12)) /
          years(1)
      ),
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
      dplyr::lag(bx_name_short, default = first(bx_name_short)) &
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
      prior_biologic_chk == "prior biologic not listed" ~ "Non-biologic naive",
      n_distinct(bx_start_date) < n() ~ "Duplicate start dates",
      any(bx_end_date < bx_start_date) ~ "End date before start date",
      all(bx_start_date > last_visit_date) ~ "All future treatments",
      any(bx_start_date > last_visit_date) ~ "Any future treatments",
      any(is.na(bx_end_date)) ~ "Missing end date",
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

### Table ----
# See "20261014_uk_feasibility.R"
uk_bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions) |>
  adorn_totals() |>
  arrange(!exclusions == "Total")

uk_bx_cln_exclusions |>
  distinct(patient_id, exclusions_table) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  arrange(!exclusions_table == "Total")

# 20261027 sent queries to Dominic P

# 10 patients with no biologic episode information have info under medications
# 20261026 DP - These patients are "intent only" visit listed bx info but none in episodes, see uk_biologic_unplaced
no_bx_info <- uk_bx_cln_exclusions |>
  filter(exclusions == "No biologic information") |>
  distinct(patient_id) |>
  pull(patient_id)
# forreview <- asthma_medications |>
#   filter(uidpat %in% no_bx_info) |>
#   select(uidpat, contains("anti"), -uk_anti_ige_treat) |>
#   filter(if_any(contains("anti"), \(x) !is.na(x))) |>
#   distinct(uidpat) |>
#   pull(uidpat) |>
#   noquote()
# forreview

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

## Analysis ----
uk_bx_analysis <- uk_bx_cln_exclusions |>
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
  # UK outcome =
  # UK continuers: 1 biologic episode, flagged as ongoing, or stop indicated more than a month after last visit
  # UK stoppers: 1 biologic episode, flagged as stopped before or within a month of last visit; >1 biologic and gap between 1st and 2nd >6 months
  # UK switchers: >1 biologic episode, <=6 months between 1 and 2nd
  mutate(
    bx_initiation_date = min(bx_start_date),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
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

### Save ----
save(
  uk_biologics,
  uk_bx_cln_exclusions,
  uk_bx_analysis,
  file = here::here("datasets", "biologic_data_uk.Rdata")
)
# rm(list = ls()[!str_detect(ls(), "analysis|exclusions")])

# Italy ----

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
ita_biologics <- ita$ita_biologics
ita_bx_cln_exclusions <- ita$ita_patients |>
  select(patient_id, date_of_birth) |>
  left_join(
    ita_biologics |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    bx_initiation_age = interval(date_of_birth, min(bx_start_date)) / years(1),
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
      dplyr::lag(bx_name_short, default = first(bx_name_short)) &
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

### Table ----
ita_bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions) |>
  adorn_totals() |>
  arrange(!exclusions == "Total")

ita_bx_cln_exclusions |>
  distinct(patient_id, exclusions_table) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  arrange(!exclusions_table == "Total")

# Queries sent to NE 20261120
# writexl::write_xlsx(
#   bx_cln_exclusions |>
#     filter(exclusions_table == "Data queries") |>
#     distinct(patient_id, exclusions),
#   here::here("data_cleaning", "ita_bx_data_queries_20261120.xlsx")
# )

## Analysis ----
ita_bx_analysis <- ita_bx_cln_exclusions |>
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
    bx_initiation_date = min(bx_start_date),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
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

### Save ----
save(
  ita_biologics,
  ita_bx_cln_exclusions,
  ita_bx_analysis,
  file = here::here("datasets", "biologic_data_ita.Rdata")
)
# rm(list = ls()[!str_detect(ls(), "analysis|exclusions")])

# See "20261120_ita_feasibility.R"
ita_bx_analysis |>
  distinct(patient_id, .keep_all = T) |>
  tabyl(outcome) |>
  adorn_totals() |>
  adorn_pct_formatting(1)

# Portugal ----

# Portugal-specific queries, sent to CG 20261123:

## Exclusions ----
prt_biologics <- prt$prt_biologics
prt_bx_cln_exclusions <- prt$prt_patients |>
  select(patient_id, date_of_birth) |>
  left_join(
    prt_biologics |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    bx_initiation_age = interval(date_of_birth, min(bx_start_date)) / years(1),
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
      dplyr::lag(bx_name_short, default = first(bx_name_short)) &
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
      n_distinct(bx_start_date) < n() ~ "Duplicate start dates",
      any(bx_end_date < bx_start_date) ~ "End date before start date",
      all(bx_start_date > last_visit_date) ~ "All future treatments",
      any(bx_start_date > last_visit_date) ~ "Any future treatments",
      any(is.na(bx_end_date) & is.na(bx_ongoing)) ~ "Missing end date",
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

### Table ----
prt_bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions) |>
  adorn_totals() |>
  arrange(!exclusions == "Total")

prt_bx_cln_exclusions |>
  distinct(patient_id, exclusions_table) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  arrange(!exclusions_table == "Total")

# Queries sent to NE 20261120
# writexl::write_xlsx(
#   prt_bx_cln_exclusions |>
#     filter(exclusions_table == "Data queries") |>
#     distinct(patient_id, exclusions),
#   here::here("data_cleaning", "prt_bx_data_queries_20261123.xlsx")
# )

## Analysis ----
prt_bx_analysis <- prt_bx_cln_exclusions |>
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
  # PRT outcome =
  # PRT continuers: 1 biologic episode, flagged as ongoing, or stop indicated more than a month after last visit
  # PRT stoppers: 1 biologic episode, flagged as stopped before or within a month of last visit; >1 biologic and gap between 1st and 2nd >6 months
  # PRT switchers: >1 biologic episode, <=6 months between 1 and 2nd
  mutate(
    bx_initiation_date = min(bx_start_date),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
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
          is.na(bx_ongoing) &
          bx_end_date > last_visit_date %m+% months(1)
      ) ~ "Continue", # Censored not stop
      # any(bx_ct == 1 & bx_ongoing == "No" & bx_end_date <= last_visit_date) ~ "Stop",
      all(
        bx_ct == 1 &
          is.na(bx_ongoing) &
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

### Save ----
save(
  prt_biologics,
  prt_bx_cln_exclusions,
  prt_bx_analysis,
  file = here::here("datasets", "biologic_data_prt.Rdata")
)
# rm(list = ls()[!str_detect(ls(), "analysis|exclusions")])

# See "20261120_prt_feasibility.R"
prt_bx_analysis |>
  distinct(patient_id, .keep_all = T) |>
  tabyl(outcome) |>
  adorn_totals() |>
  adorn_pct_formatting(1)

# NJH ----
# No revised SRD dataset yet
load(here::here("datasets", "qisar2_srd_njh.Rdata"))

# bx_cln now called njh_biologics, removed exclusions steps to make study agnostic and match others
# njh_bx_cln_exclusions now contains same exclusions
## Exclusions ----
njh_bx_cln_exclusions <- njh_patients |>
  select(patient_id, date_of_birth) |>
  right_join(
    njh_biologics |>
      mutate(bx_data = TRUE),
    by = join_by(patient_id)
  ) |>
  arrange(patient_id, bx_start_date, bx_end_date) |>
  mutate(
    bx_initiation_age = interval(date_of_birth, min(bx_start_date)) / years(1),
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
      age_at_baseline <
        18 ~ "Under-18 at baseline, 12 months prior to initiation",
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
        is.na(bx_start_date) & is.na(bx_end_date)
      ) ~ "No start or end dates provided",
      all(is.na(bx_start_date)) ~ "No start dates provided",
      # Exclude if medication reported by patient and start date same as observation date
      any(
        bx_start_date == observation_date & historicalorreported == 1
      ) ~ "Unreliable start date - Reported", # Exclude for now
      # January first, very unlikely and probably imputed by NJH
      any(
        month(bx_start_date) == 1 & day(bx_start_date) == 1
      ) ~ "Unreliable start date - Imputed",
      # any(is.na(start_date)) ~ "Missing a biologic initiation date", #Comment out this line and in factor below to get 20261015 feasibility number
      .default = NA_character_
    ) |>
      factor(
        levels = c(
          "Under-18 at baseline, 12 months prior to initiation",
          "Started biologic before 2016",
          "Prescribed a biologic only for non-asthma reason",
          "Prescribed a biologic for non-asthma reason",
          # "All patient reported dates", # catch patients here or with missing start dates
          "No start or end dates provided",
          "No start dates provided",
          "Unreliable start date - Reported",
          "Unreliable start date - Imputed"
          # "Missing a biologic initiation date" #Comment out this line and in factor below to get 20261015 feasibility number
        )
      ),
    .by = patient_id
  )

# Initial exclusions
njh_bx_cln_exclusions |>
  distinct(patient_id, exclusions) |>
  tabyl(exclusions)

# Use first biologic initiation date (n=3663 patients)
# 20261015 - refined to use 3 months' follow-up as minimum for Obj1

# Defines all eligible biologic patients for further analysis
patients <- njh_bx_cln_exclusions |>
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
    setting = "USA - NJH",
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
          "Under-18 at baseline, 12 months prior to initiation",
          "Started biologic before 2016",
          "Prescribed a biologic only for non-asthma reason",
          "Prescribed a biologic for non-asthma reason",
          "No start or end dates provided",
          "No start dates provided",
          "Unreliable start date - Reported",
          "Unreliable start date - Imputed",
          # "Missing a biologic initiation date",
          "Started first biologic before first observation",
          "Duplicated start date - multiple biologic",
          "<3 months post-initiation follow-up"
        )
      ),
    exclusions_table = case_when(
      exclusions ==
        "Under-18 at baseline, 12 months prior to initiation" ~ "Under 18 years-old",
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
          "Under 18 years-old",
          "Non-biologic naive",
          "Data queries pending",
          "Started biologic before 2016",
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

### Table ----
# Aligns with earlier number
patients |>
  tabyl(exclusions) |>
  adorn_totals()
patients |>
  tabyl(exclusions_table) |>
  adorn_totals()

patients |>
  tabyl(exclusions) |>
  adorn_pct_formatting() |>
  kable() |>
  kable_styling()

index_dates <- patients |>
  # Restrict to eligible for obj1
  filter(is.na(exclusions)) |>
  select(patient_id, bx_initiation_date = index_date)

## Analysis ----
# 20261015 Use spells approach rather than first start
njh_bx_analysis <- njh_bx_cln_exclusions |>
  # Selects eligible patients
  right_join(index_dates, by = join_by(patient_id)) |>
  # relocate(
  #   first_visit_date,
  #   last_visit_date,
  #   bx_initiation_date,
  #   .after = last_visit_date
  # ) |>
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
    .by = c(patient_id)
  ) |>
  summarise(
    bx_start_date = min(bx_start_date, na.rm = T),
    bx_name_short = first(bx_name_short),
    bx_type = first(bx_type),
    .by = c(patient_id, bx_initiation_date, bx_nbr)
  ) |>
  left_join(
    njh_patients |>
      distinct(
        patient_id,
        data_capture_system,
        country,
        setting,
        site,
        first_visit_date,
        last_visit_date,
        num_of_visits
      ),
    by = join_by(patient_id)
  ) |>
  relocate(
    data_capture_system,
    country,
    setting,
    site,
    first_visit_date,
    last_visit_date,
    num_of_visits,
    bx_initiation_date,
    .after = patient_id
  ) |>
  mutate(
    bx_end_date = as.Date(NA),
    follow_up_mos = interval(bx_initiation_date, last_visit_date) / months(1),
    bx_ct = max(bx_nbr),
    # Outcome based on number of spells
    outcome = case_when(
      bx_ct == 1 ~ "Continue",
      bx_ct > 1 ~ "Switch"
    ),
    .by = patient_id
  ) |>
  relocate(bx_nbr, .before = bx_ct)

njh_bx_analysis |>
  distinct(patient_id, .keep_all = T) |>
  tabyl(outcome) |>
  adorn_totals()

### Save ----
save(
  njh_bx_cln_exclusions,
  patients,
  njh_bx_analysis,
  file = here::here("datasets", "biologic_data_njh.Rdata")
)

## Switch dates ----
# Keep second record for switchers
# 20261110 - this matches 20261015 feasibility numbers, n=299
switch_dates <- njh_bx_analysis |>
  filter(
    bx_ct > 1,
    bx_nbr == 2
  ) |>
  # left_join(index_dates, by = join_by(patient_id)) |>
  mutate(
    followup_mos = interval(bx_start_date, last_visit_date) / months(1),
    eligible_fup = followup_mos >= 11,
    post_switch_date = bx_start_date %m+% months(12)
  ) |>
  select(
    patient_id,
    first_visit_date,
    last_visit_date,
    bx_initiation_date,
    switch_date = bx_start_date,
    post_switch_date,
    followup_mos,
    eligible_fup
  ) |>
  filter(eligible_fup)

# Eligible ----
# File for querying baseline and pre/post analyses
# Actually couldn't this be the analysis dataset for obj 1 and used for querying baseline?
# patient_id, initiation_date, switchstop_date (pre), index_date + 12 months

get_eligible <- function(df) {
  return(
    df |>
      # All biologic patients
      select(
        any_of(c(
          "patient_id",
          "data_capture_system",
          "country",
          "setting",
          "site",
          "first_visit_date",
          "last_visit_date",
          "num_of_visits",
          "bx_initiation_date",
          "follow_up_mos",
          "bx_nbr",
          "bx_ct",
          "outcome",
          "bx_start_date",
          "bx_end_date",
          "bx_ongoing",
          "bx_name_short",
          "bx_type",
          "bx_switchreason"
        ))
      ) |>
      pivot_wider(
        names_from = bx_nbr,
        values_from = any_of(c(
          "bx_start_date",
          "bx_end_date",
          "bx_ongoing",
          "bx_name_short",
          "bx_type",
          "bx_switchreason"
        )),
        names_vary = "slowest"
      ) |>
      mutate(
        # Date of biologic change
        switchstop_date = case_when(
          outcome == "Switch" ~ bx_start_date_2,
          outcome == "Stop" ~ bx_end_date_1,
          .default = NA
        ),
        # Months from biologic initiation to biologic change
        switchstop_mos = interval(bx_initiation_date, switchstop_date) /
          months(1),
        # Months of follow-up from biologic change to last visit
        post_switchstop_followup_mos = interval(
          switchstop_date,
          last_visit_date
        ) /
          months(1),
        # >= 11 months of follow-up from biologic change to last visit
        eligible_fup = post_switchstop_followup_mos >= 11,
        # Query date for post-change interval
        post_switchstop_date = switchstop_date %m+% months(12)
      )
  )
}

eligible <- bind_rows(
  get_eligible(uk_bx_analysis),
  get_eligible(ita_bx_analysis),
  get_eligible(prt_bx_analysis),
  get_eligible(njh_bx_analysis),
  get_eligible(isar_bx_analysis)
) |>
  # Relabel after binding
  sjlabelled::copy_labels(df_origin = uk_bx_analysis)

# Eligible numbers pretty close to previous ones with some new exclusions (<18 for NJH & ISAR)
eligible |>
  tabyl(eligible_fup, outcome, data_capture_system) |>
  adorn_totals("row")

eligible |>
  filter(data_capture_system == "Bespoke") |>
  tabyl(eligible_fup, outcome, country) |>
  adorn_totals("row")
