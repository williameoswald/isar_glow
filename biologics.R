# 20260401
# Clean using new biologics dataset

bx_cln <- biologics_labelled |>
  rename_with(
    .cols = contains("bio_"),
    ~ str_replace(., "bio_", "bx_")
  ) |>
  mutate(
    bx_name_short = structure(
      factor(
        str_trim(str_remove(bx_name, " \\(.*\\)$")),
        levels = c(
          "Benralizumab",
          "Mepolizumab",
          "Reslizumab",
          "Dupilumab",
          "Omalizumab",
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
            "Anti-IL5",
            "Anti-IL5R",
            "Anti-IL4Ralpha",
            "Anti-IgE",
            "Anti-TSLP"
          ),
          labels = c(
            "Anti-IL5",
            "Anti-IL5R",
            "Anti-IL4Ralpha",
            "Anti-IgE",
            "Anti-TSLP"
          )
        ),
      label = "Biologic class"
    )
  ) |>
  relocate(bx_name_short, bx_type, .after = bx_name) |>
  #   Reverse imputation for now
  mutate(
    # bx_start_date = structure(
    #   case_when(
    #     bx_stdate_imputed_flag == "Date missing, visit date imputed" ~ NA,
    #     TRUE ~ bx_stdate
    #   ),
    #   label = "Biologic start date"
    # ),
    # bx_end_date = structure(
    #   case_when(
    #     bx_endate_imputed_flag ==
    #       "New prescription at next visit, next prescription start date imputed" ~ NA,
    #     TRUE ~ bx_endate
    #   ),
    #   label = "Biologic end date"
    # ),
    #   Exclude imputed months AND days
    exclude_imp_date = case_when(
      bx_stdate_imputed_flag == "MonthsDays" ~ TRUE,
      bx_endate_imputed_flag == "MonthsDays" ~ TRUE
    ),
    bx_flag = TRUE
  ) |>
  group_by(patient_id) |>
  # fill(c("exclude_imp_date")) |>
  ungroup() |>
  rename(bx_start_date = bx_stdate, bx_end_date = bx_endate)
#   Imputed dates reverted to missing, only accepting imputed days within month for dates
#   Drop imputed flags and keep distinct rows
# select(-contains("imputed"))  #n=51747

## Analysis ----
# 23480 patients without biologic data
patient_labelled |>
  anti_join(bx_cln, by = "patient_id") |>
  nrow()
# 0 patients in biologic data without patient data
bx_cln |>
  anti_join(patient_labelled, by = "patient_id") |>
  distinct(patient_id) |>
  nrow()

bx_err <- patient_labelled |>
  remove_empty_imputation_fields() |>
  select(-data_capture_system) |>
  #   Demographics - add age and index_date, confirm the latter
  left_join(
    clean_demographics() |>
      select(
        patient_id,
        index_date,
        index_age
      ),
    by = join_by(patient_id),
    relationship = "one-to-one"
  ) |>
  #   Left join to keep patients without biologic data for consort
  left_join(
    biologics_labelled |>
      clean_biologics() |>
      remove_empty_imputation_fields() |>
      #   Reverse imputation for now
      mutate(
        # bx_start_date = structure(
        #   case_when(
        #     bx_stdate_imputed_flag == "Date missing, visit date imputed" ~ NA,
        #     TRUE ~ bx_stdate
        #   ),
        #   label = "Biologic start date"
        # ),
        # bx_end_date = structure(
        #   case_when(
        #     bx_endate_imputed_flag ==
        #       "New prescription at next visit, next prescription start date imputed" ~ NA,
        #     TRUE ~ bx_endate
        #   ),
        #   label = "Biologic end date"
        # ),
        #   Exclude imputed months AND days
        exclude_imp_date = case_when(
          bx_stdate_imputed_flag == "MonthsDays" ~ TRUE,
          bx_endate_imputed_flag == "MonthsDays" ~ TRUE
        ),
        bx_flag = TRUE
      ),
    by = "patient_id",
    relationship = "one-to-many"
  ) |>
  # filter(country %in% c("UK", "Denmark", "Spain", "Italy", "USA")) |>
  group_by(patient_id) |>
  #   Flag issues
  arrange(bx_start_date) |>
  mutate(
    # Missing bx name but not dates - empty rows dropped above
    err1 = case_when(
      any(bx_flag) &
        is.na(bx_name) ~ "missing_bx_name_not_dates"
    ),
    # Missing start date
    err2 = case_when(
      any(bx_flag) &
        is.na(bx_start_date) ~ "missing_bx_start_date"
    ),
    # err3 = case_when(
    #   any(bx_flag) &
    #     n() > 1 &
    #     all(is.na(bx_end_date)) ~ "no_end_dates"
    # ),
    # Later biologic starts before end of previous
    err4 = case_when(
      any(bx_flag) &
        bx_start_date <
          dplyr::lag(bx_end_date) ~ "later_bx_start_before_prior_end"
    ),
    # Biologic ends before start
    err5 = case_when(
      any(bx_flag) &
        bx_end_date < bx_start_date ~ "end_before_start"
    ),
    err6 = case_when(
      any(bx_flag) &
        any(
          duplicated(
            pick(bx_start_date) |>
              filter(!is.na(bx_start_date))
          )
        ) ~ "duplicate_bx_start_date"
    ),
    # Duplicate start dates within patient and biologic
    err7 = case_when(
      any(bx_flag) &
        any(
          duplicated(
            pick(bx_name, bx_start_date) |>
              filter(!is.na(bx_name) & !is.na(bx_start_date))
          )
        ) ~ "duplicate_bx_start_date_and_name"
    ),
    # Multiple recordings of same biologic (ie doses), ignores possible gaps
    err8 = case_when(
      any(bx_flag) &
        n() > 1 &
        n_distinct(bx_name, na.rm = TRUE) == 1 ~ "multiple_same_bx"
    ),
    err9 = case_when(
      any(bx_flag) &
        n() > 1 &
        all(is.na(bx_end_date)) ~ ">1 record, no end dates"
    ),
    err10 = case_when(
      any(bx_flag) &
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
    bx_start_year = year(bx_start_date),
    exclusions = case_when(
      all(is.na(bx_flag)) ~ "Not in biologic data",
      all(num_of_visits == 0) ~ "No visits recorded", # means no visit dates for imputation
      country %in% c("Ecuador", "Qatar") ~ "Fewer than two biologics available",
      any(bx_start_year < 2016) ~ "Biologic start date before 2016",
      index_age < 18 ~ "Under-18",
      any(is.na(bx_name)) ~ "Biologic name not listed",
      any(is.na(bx_start_date)) ~ "Biologic start date not listed",
      any(exclude_imp_date) ~ "Imputed start/end MonthsDays",
      # max(1:n()) > 6 ~ "More than 6 listed treatments",
      any(errors != "") ~ "Flagged errors",
      .default = "Eligible"
    ) |>
      factor(
        levels = c(
          "Not in biologic data",
          "No visits recorded",
          "Fewer than two biologics available",
          "Biologic start date before 2016",
          "Under-18",
          "Biologic name not listed",
          "Biologic start date not listed",
          "Imputed start/end MonthsDays",
          # "More than 6 listed treatments",
          "Flagged errors",
          "Eligible"
        )
      ),
    # country = if_else(
    #   is.na(country),
    #   str_extract(patient_id, "^[^_]+"),
    #   country
    # ),
  ) |>
  ungroup()

# Check that index date aligns with first biologic date
bx_err |>
  mutate(
    index_year = year(index_date),
    index_year_chk = first(bx_start_year),
    .by = patient_id
  ) |>
  filter(index_year != index_year_chk) # empty set

error_summary <- bx_err |>
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
  arrange(error_flag)
error_summary |>
  count(country, error_flag, sort = TRUE) |>
  arrange(country, error_flag)
error_summary |> view()

# Biologic history ----
bx_history <- bx_err |>
  # Drop patients with errors
  filter(
    exclusions == "Eligible" |
      exclusions == "Fewer than two biologics available" |
      exclusions == "Biologic start date before 2016"
  ) |>
  select(
    -bx_flag,
    -contains("err"),
    # -contains("exclu")
  ) |>
  arrange(patient_id, bx_start_date) |>
  group_by(patient_id) |>
  # Assign earlier end date as subsequent start date in these cases
  mutate(
    bx_nbr = row_number(),
    bx_ct = n(),
    # Use last visit date where valid for missing end dates
    bx_end_date2 = case_when(
      # If 1 bx listed and no end date available, assume on treatment until last visit, unless last visit is before start date
      bx_ct == 1 &
        is.na(bx_end_date) &
        last_visit_date >= bx_start_date ~ last_visit_date, #switched to allow same last visit and start
      # If last bx and no end date available, assume end date is last visit
      row_number() == bx_ct &
        last(is.na(bx_end_date)) &
        last_visit_date >= bx_start_date ~ last_visit_date, #switched to allow same last visit and start
      TRUE ~ bx_end_date
    ),
    # Avoid 0 length intervals by adding 1 day to end
    bx_end_date2 = if_else(
      bx_start_date == bx_end_date2,
      bx_end_date2 + days(1),
      bx_end_date2
    ),
    # Calculate treatment months
    tx_months = interval(bx_start_date, bx_end_date2) / months(1),
    # Add gap in months to allow visualisation
    # Calculate gap between end of current treatment row and start of next
    next_bx_start_date = if_else(
      bx_nbr < bx_ct,
      dplyr::lead(bx_start_date),
      NA
    ),
    tx_gap = interval(bx_end_date2, next_bx_start_date),
    tx_gap_months = tx_gap / months(1),
    # Flag if gap row needs to be added
    tx_gap_row = coalesce(tx_gap > 0, FALSE),
  ) |>
  relocate(c(bx_start_date, bx_end_date2, tx_months), .after = bx_end_date) |>
  # Duplicate row where treatment gap indicated
  tidyr::uncount(
    weights = if_else(tx_gap_row, 2L, 1L),
    .id = "gap", # specify new column name to indicate "gap" rows
    .remove = TRUE # remove "weights" column if present, not needed
  ) |>
  # Insert gaps into treatment months
  mutate(
    tx_months = if_else(gap == 2, tx_gap_months, tx_months),
    bx_start_date = if_else(
      gap == 2,
      dplyr::lag(bx_end_date2),
      bx_start_date
    ),
    bx_end_date2 = if_else(gap == 2, dplyr::lead(bx_start_date), bx_end_date2),
    bx_type = if_else(gap == 2, "No treatment", bx_type),
    bx_name = if_else(gap == 2, "Treatment gap", bx_name),
    bx_name_short = if_else(gap == 2, "Treatment gap", bx_name_short),
    bx_nbr = if_else(gap == 2, NA, bx_nbr),
    # Calculate months since index date per treatment and gaps for visualisation
    start_months = dplyr::lag(cumsum(tx_months), default = 0),
    end_months = cumsum(tx_months)
  ) |>
  ungroup() |>
  arrange(patient_id, bx_start_date)

save(bx_cln, bx_err, bx_history, file = "biologic_data.RData")
