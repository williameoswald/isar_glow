# Based on "analysis_ats_exploratory.R"

load("biologic_data.RData")

# Specifications ----
library(patchwork)

# Define biologic order
bx_levels <- c(
  "Tezepelumab",
  "Benralizumab",
  "Dupilumab",
  "Reslizumab",
  "Mepolizumab",
  "Omalizumab"
)

# Named colours keyed to those levels
bx_colours <- setNames(
  RColorBrewer::brewer.pal(6, "Set2"),
  bx_levels
)

ftsize = 14
ftscale = 0.8

# Patient flow ----
bx_cln |>
  filter_out(exclusions == "Not in biologic data") |>
  droplevels() |>
  distinct(patient_id, .keep_all = TRUE) |>
  mutate(
    exclusions_table = case_when(
      exclusions == "Biologic name not listed" ~ "Incomplete or incorrect data",
      exclusions ==
        "Biologic start date not listed" ~ "Incomplete or incorrect data",
      exclusions ==
        "Imputed start/end MonthsDays" ~ "Incomplete or incorrect data",
      exclusions == "Flagged errors" ~ "Incomplete or incorrect data",
      exclusions ==
        "Less than 3 months of follow-up" ~ "Less than 3 months of follow-up",
      .default = exclusions
    ) |>
      factor(
        levels = c(
          "No visits recorded",
          "Fewer than two biologics available",
          "Biologic start date before 2016",
          "Under-18",
          "Incomplete or incorrect data",
          "Less than 3 months of follow-up",
          "Eligible"
        )
      )
  ) |>
  tabyl(exclusions_table) |>
  adorn_totals() |>
  adorn_pct_formatting(1) |>
  arrange(exclusions_table != "Total")

# Overall followup ----
fup <- bx_cln |>
  # filter(exclusions != "Not in biologic data") |>
  filter(exclusions == "Eligible") |>
  slice_min(bx_start_date, by = patient_id, with_ties = F) |>
  mutate(
    follow_up = follow_up_mos / 12
  ) |>
  select(
    patient_id,
    first_visit_date,
    last_visit_date,
    bx_start_date,
    follow_up
  )

# Distribution of followup
c("follow_up") |>
  map_dfr(\(x) distfunc(data = fup, var = x))
distfunc(
  bx_cln |>
    distinct(patient_id, .keep_all = T) |>
    filter(exclusions %in% c("Eligible", "Less than 12 months of follow-up")),
  var = "follow_up_mos"
)
fu <- fup |>
  # filter_out(outcome == "Continue") |>
  ggplot(aes(
    x = follow_up,
  )) +
  geom_boxplot(alpha = 0.5) +
  scale_y_discrete(limits = rev) +
  # scale_fill_brewer(palette = "Set2", name = "Outcome") +
  # scale_colour_brewer(palette = "Set2", name = "Outcome") +
  labs(
    x = "Months from biologic end date (or start date for continuers) to date of last recorded visit",
    fill = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )
fu
# preview_or_save(fu, "ers_fup_by_outcome", w = 18, h = 18, dpi = 300)

fup_plot <- fup |>
  ggplot(aes(
    x = follow_up,
    y = 1
  )) +
  # geom_point(
  #   # aes(
  #   #   fill = gap_date_imputed,
  #   #   colour = gap_date_imputed
  #   # ),
  #   position = position_jitterdodge(dodge.width = 0.60),
  #   alpha = 0.4,
  #   size = 3,
  #   stroke = 0,
  #   na.rm = TRUE
  # ) +
  # https://stackoverflow.com/questions/52217852/how-to-add-means-to-a-ggplot-geom-point-plot
  stat_summary(
    fun = "mean",
    geom = "point",
    col = "grey30",
    size = 3,
    shape = 24,
    fill = "darkred"
  ) +
  geom_boxplot(
    width = 0.3,
    alpha = 0.2,
    na.rm = TRUE,
    outlier.shape = NA,
    colour = "grey30"
  ) +
  geom_violin(
    width = 0.5,
    alpha = 0.2,
    na.rm = TRUE,
    colour = "grey30",
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 25, by = 1),
    limits = c(0, 25)
  ) +
  # scale_y_discrete(limits = rev) +
  # scale_fill_brewer(palette = "Set2", name = "Outcome") +
  # scale_colour_brewer(palette = "Set2", name = "Outcome") +
  labs(
    x = "Years from first biologic start to last visit",
    y = NULL,
    fill = NULL,
    colour = NULL
  ) +
  # theme_minimal() +
  theme(
    # aspect.ratio = 1,
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.position.inside = c(0.85, 0.85),
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    # plot.title = element_text(size = ftsize, colour = "grey30", hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = ftsize, colour = "grey30"),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = ftsize, colour = "grey30")
  )
fup_plot

# Country representation ----
bx_cln |>
  distinct(patient_id, .keep_all = T) |>
  filter(exclusions %in% c("Eligible", "Less than 12 months of follow-up")) |>
  droplevels() |>
  count(country) |>
  adorn_totals() |>
  adorn_percentages("col") |>
  adorn_pct_formatting(1) |>
  adorn_ns("front") |>
  arrange(country != "Total")

# Lollipop chart of eligible patient counts by country
country_counts <- bx_cln |>
  distinct(patient_id, .keep_all = T) |>
  filter(exclusions %in% c("Eligible", "Less than 12 months of follow-up")) |>
  distinct(patient_id, .keep_all = TRUE) |>
  mutate(
    country = case_when(
      country == "UK" ~ "United Kingdom",
      country == "USA" ~ "United States of America",
      country == "UAE" ~ "United Arab Emirates",
      TRUE ~ country
    )
  ) |>
  count(country) |>
  mutate(p = n / sum(n)) |>
  ggplot() +
  geom_segment(
    aes(
      y = country,
      yend = country,
      x = 0,
      xend = n
    ),
    color = "grey"
  ) +
  geom_point(
    aes(
      y = country,
      x = n
    ),
    size = 6,
    colour = "#073763"
  ) +
  scale_x_continuous(
    limits = c(0, 3000)
    # expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  # scale_colour_gradient(
  #   low = "#07376350",        # lightest shade
  #   high = "#073763"   # your chosen colour
  # ) +
  # scale_colour_viridis_d(limits = rev) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.text = element_markdown(size = 24, colour = "grey30"),
    axis.title.x = element_markdown(size = 24, colour = "grey30")
  ) +
  labs(x = "Eligible Participants", y = NULL)
country_counts

preview_and_save(country_counts, "country_counts", h = 30, w = 30, dpi = 300)

# Approval dates ----
datefile <- "C:/Users/William.Oswald/Optimum Patient Care Global Limited/ISAR Team - Documents/1. Research Projects/1. Active/GLOW/3. Data/2. Dataset/approval_data/biologics_asthma_approval_reimbursement.xlsx"
excel_sheets(datefile)

dates <- c("Europe", "Americas", "Asia") |>
  set_names() |>
  map(\(sheet) read_excel(datefile, sheet = sheet)) |>
  list_rbind(names_to = "region") |>
  clean_names() |>
  unite(note, notes) |>
  mutate(
    approval_date2 = excel_numeric_to_date(as.numeric(approval_date)),
    reimbursement_date2 = excel_numeric_to_date(as.numeric(reimbursement_year)),
    bx_name_short = case_when(
      str_detect(biologic, "Oma") ~ "Omalizumab",
      str_detect(biologic, "Mep") ~ "Mepolizumab",
      str_detect(biologic, "Res") ~ "Reslizumab",
      str_detect(biologic, "Ben") ~ "Benralizumab",
      str_detect(biologic, "Dup") ~ "Dupilumab",
      str_detect(biologic, "Tez") ~ "Tezepelumab",
      str_detect(biologic, "Dep") ~ "Depemokimab"
    )
  ) |>
  select(
    region,
    country,
    biologic,
    bx_name_short,
    approval_date2,
    reimbursement_date2
  ) |>
  left_join(
    bx_history |>
      mutate(index_year = year(bx_start_date)) |>
      filter_out(is.na(bx_name_short) | is.na(index_year)) |>
      filter_out(bx_name_short == "Treatment gap") |>
      slice_min(
        n = 1,
        order_by = bx_start_date,
        by = c(country, bx_name_short),
        with_ties = F
      ) |>
      select(country, bx_name_short, bx_start_date),
    by = join_by(country, bx_name_short)
  ) |>
  mutate(
    bx_name_short = factor(
      bx_name_short,
      levels = c(
        "Omalizumab",
        "Mepolizumab",
        "Reslizumab",
        "Benralizumab",
        "Dupilumab",
        "Tezepelumab",
        "Depemokimab"
      )
    ),
    across(contains("date"), ymd)
  ) |>
  pivot_longer(
    cols = contains("date"),
    names_to = "type",
    values_to = "date"
  ) |>
  mutate(
    type = replace_values(
      type,
      "approval_date2" ~ "approval",
      "reimbursement_date2" ~ "reimbursement",
      "bx_start_date" ~ "isar_date"
    ) |>
      factor(
        levels = c("approval", "reimbursement", "isar_date"),
        labels = c("Approval", "Reimbursement", "First ISAR")
      )
  ) |>
  filter_out(
    when_any(is.na(date), bx_name_short == "Depemokimab", year(date) < 2000)
  ) |>
  droplevels()

date_plot <- dates |>
  filter(country %in% c("Spain", "Italy", "UK", "USA")) |>
  ggplot() +
  geom_point(
    aes(y = bx_name_short, x = date, colour = type),
    size = 5,
    position = position_jitter(height = 0.05)
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_date(
    # expand = expansion(add = c(0, 0)),
    # limits = c(2000,)
  ) +
  facet_wrap(vars(country), scales = "fixed") +
  labs(
    x = NULL,
    y = NULL,
    colour = "Date"
  ) +
  theme(
    # panel.background = element_rect(fill="transparent", colour=NA_character_), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(1.2, "cm"),
    legend.key = element_rect(fill = NA, colour = NULL),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = ftsize, colour = "grey30"),
    axis.text.y = element_text(size = ftsize, colour = "grey30"),
    # axis.title = element_text(size = ftsize, colour = "grey30"),
    strip.text = element_text(size = ftsize, colour = "grey30"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    # strip.text.y.bottom = element_text()
  )
date_plot
preview_and_save(
  date_plot,
  "ers_approval_dates",
  h = 19,
  w = 33,
  dpi = 300
)

## ISAR date by biologic and country ----
isar_dates <- bx_history |>
  filter_out(is.na(bx_name_short) | is.na(bx_start_date)) |>
  filter_out(bx_name_short == "Treatment gap") |>
  filter(country %in% c("Italy", "Spain", "UK", "USA")) |>
  select(country, bx_name_short, bx_start_date) |>
  mutate(bx_name_short = factor(bx_name_short, levels = bx_levels))

isar_dates_plot <- isar_dates |>
  filter(country == "Spain") |>
  arrange(country, bx_name_short) |>
  ggplot(aes(
    x = bx_start_date,
    y = bx_name_short,
    fill = bx_name_short,
    colour = bx_name_short
  )) +
  geom_boxplot(alpha = 0.5) +
  scale_y_discrete(limits = rev) +
  scale_fill_brewer(palette = "Set2", name = "Outcome") +
  scale_colour_brewer(palette = "Set2", name = "Outcome") +
  labs(
    x = "Months from biologic end date (or start date for continuers) to date of last recorded visit",
    fill = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )
isar_dates_plot

date_plot <- function(countryname) {
  isar_sub <- isar_dates |>
    filter(country == countryname) |>
    arrange(bx_name_short)

  dates <- dates |>
    filter(country == countryname, type %in% c("Approval", "Reimbursement")) |>
    select(bx_name_short, type, date) |>
    filter(bx_name_short %in% unique(isar_sub$bx_name_short))

  isar_sub |>
    ggplot(aes(
      x = bx_start_date,
      y = bx_name_short
    )) +
    geom_vline(
      data = dates,
      aes(xintercept = date, colour = bx_name_short, linetype = type),
      linewidth = 1.1,
      alpha = 0.8,
      key_glyph = "path"
    ) +
    geom_point(
      aes(
        fill = bx_name_short,
        colour = bx_name_short
      ),
      position = position_jitterdodge(dodge.width = 0.60),
      alpha = 0.2,
      size = 3,
      stroke = 0,
      na.rm = TRUE
    ) +
    # https://stackoverflow.com/questions/52217852/how-to-add-means-to-a-ggplot-geom-point-plot
    stat_summary(
      fun = "mean",
      geom = "point",
      col = "grey30",
      size = 3,
      shape = 24,
      fill = "darkred"
    ) +
    geom_boxplot(
      width = 0.3,
      alpha = 0.2,
      na.rm = TRUE,
      outlier.shape = NA,
      colour = "grey30"
    ) +
    geom_violin(
      width = 0.5,
      alpha = 0.2,
      na.rm = TRUE,
      colour = "grey30",
      drop = FALSE
    ) +
    scale_x_date(
      date_breaks = "2 year",
      date_labels = "%Y",
      limits = as.Date(c("2000-07-02", "2027-07-02")),
      expand = c(0, 0)
      # guide = guide_axis(n.dodge = 2)
    ) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = bx_colours, guide = NULL) +
    scale_colour_manual(values = bx_colours, guide = NULL) +
    scale_linetype_manual(values = c(1, 2)) +
    labs(
      x = "Biologic start dates",
      y = NULL,
      fill = NULL,
      colour = NULL,
      linetype = NULL
    ) +
    facet_wrap(vars(bx_name_short), scales = "free_y", ncol = 1) +
    guides(
      linetype = guide_legend(
        nrow = 2,
        reverse = F
        # override.aes = list(size = 5)
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_blank(),
      legend.position = "right",
      # legend.position.inside = c(0.85, 0.85),
      legend.text = element_text(size = ftsize, colour = "grey30"),
      legend.title = element_text(
        size = ftsize,
        colour = "grey30",
        face = "bold"
      ),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(1.2, "cm"),
      # plot.title = element_text(size = ftsize, colour = "grey30", hjust = 0.5, face = "bold"),
      # axis.text.x = element_text(size = ftsize, colour = "grey30", angle = 45, vjust = 1, hjust = 1),
      axis.text = element_text(size = ftsize, colour = "grey30"),
      # axis.text.y = element_blank(),
      axis.title = element_text(size = ftsize, colour = "grey30")
    )
}

preview_and_save(
  date_plot("Spain"),
  "ers_isar_dates_spain",
  w = 33,
  h = 19,
  dpi = 300
)
preview_and_save(
  date_plot("Italy"),
  "ers_isar_dates_italy",
  w = 33,
  h = 19,
  dpi = 300
)
preview_and_save(
  date_plot("UK"),
  "ers_isar_dates_uk",
  w = 33,
  h = 19,
  dpi = 300
)
preview_and_save(
  date_plot("USA"),
  "ers_isar_dates_usa",
  w = 33,
  h = 19,
  dpi = 300
)


# Prescriptions ----
## By country ----
# Based on 20260422
year_labels <- bind_rows(
  tibble(
    bx_name_short = rev(bx_levels),
    year = c(2003, 2015, 2016, 2017, 2017, 2021),
    label = "FDA"
  ),
  tibble(
    bx_name_short = rev(bx_levels),
    year = c(2013, 2017, 2017, 2021, 2019, 2023),
    label = "NICE"
  ),
  tibble(
    bx_name_short = rev(bx_levels),
    year = c(2005, 2015, 2016, 2017, 2018, 2022),
    label = "EMA"
  )
) |>
  arrange(label, year)

# Identify number of biologics per country and year, includes all years
df_countries <- bx_history |>
  mutate(index_year = year(bx_start_date)) |>
  filter_out(is.na(bx_name_short) | is.na(index_year)) |>
  filter_out(bx_name_short == "Treatment gap") |>
  group_by(country, bx_name_short) |>
  mutate(
    outlier = if_else(is_outlier(index_year, low = T, high = F), "Yes", "No") |>
      factor(levels = c("Yes", "No")),
    z = (index_year - mean(index_year, na.rm = TRUE)) /
      sd(index_year, na.rm = TRUE),
    outlier2 = if_else(z < -2, "Yes", "No") |>
      factor(levels = c("Yes", "No"))
  ) |>
  left_join(
    year_labels |>
      pivot_wider(
        id_cols = bx_name_short,
        values_from = year,
        names_from = label
      ) |>
      clean_names(),
    by = "bx_name_short"
  ) |>
  mutate(
    outlier3 = if_else(index_year < fda, "Yes", "No") |>
      factor(levels = c("Yes", "No"))
  ) |>
  #   Filter out statistical and patients prior to FDA (earliest)
  filter_out(outlier == "Yes" | outlier2 == "Yes" | outlier3 == "Yes") |>
  summarise(index_year = min(index_year), .groups = "drop") |>
  mutate(
    bx_name_short = factor(
      bx_name_short,
      levels = c(
        "Omalizumab",
        "Mepolizumab",
        "Reslizumab",
        "Benralizumab",
        "Dupilumab",
        "Tezepelumab"
      )
    ),
    region = case_when(
      country %in%
        c(
          "Argentina",
          "Brazil",
          "Canada",
          "Colombia",
          "Ecuador",
          "Mexico",
          "USA"
        ) ~ "Americas",
      country %in%
        c(
          "India",
          "Japan",
          "Kuwait",
          "Qatar",
          "Saudi Arabia",
          "Singapore",
          "South Korea",
          "Taiwan",
          "UAE/Dubai"
        ) ~ "Asia",
      .default = "Europe"
    )
  )

filler <- crossing(
  country = df_countries |> distinct(country) |> pull(),
  bx_name_short = df_countries |> distinct(bx_name_short) |> pull(),
  all_years = 2000:2026,
)

# Gantt style
df_countries2 <- df_countries |>
  left_join(filler, by = c("country", "bx_name_short")) |>
  filter_out(all_years < index_year) |>
  arrange(country, bx_name_short, all_years)

# Add gap lines manually on gantt style (geom_tile)
lines <- tibble(x = seq(1999.5, 2025.5, by = 1))

fda_band <- year_labels |>
  group_by(label, year) |>
  mutate(
    total = n(),
    y = (row_number() - 0.5) / total * 0.5 # scales to fit within 0-0.5 regardless of n
  ) |>
  ungroup() |>
  ggplot(aes(y = year, x = y, fill = bx_name_short)) +
  geom_tile(aes(width = 0.5 / total, height = 1)) +
  scale_y_reverse(
    limits = c(1999.5, 2026.5),
    breaks = seq(2000, 2026, by = 2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = bx_colours) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text.x.bottom = element_text(
      angle = 60,
      hjust = 0.5,
      vjust = 0.60,
      size = ftsize * .4,
      colour = "grey30"
    )
  ) +
  facet_grid(cols = vars(label), switch = "x")

fda_band

stagger_labels <- function(labels) {
  # Find odd number labels using modulo, which returns remainder
  list(map_chr(seq_along(labels[[1]]), \(i) {
    if_else(i %% 2 == 0, str_c("\n\n", labels[[1]][i]), labels[[1]][i])
  }))
}

main_plot <- df_countries2 |>
  # filter(region == "Europe") |>
  mutate(
    bx_name_short = factor(bx_name_short, levels = bx_levels),
    country = case_when(
      country == "Saudi Arabia" ~ "Saudi\nArabia",
      country == "South Korea" ~ "South\nKorea",
      country == "UAE/Dubai" ~ "UAE\nDubai",
      .default = country
    )
  ) |>
  ggplot(aes(y = all_years, x = bx_name_short, fill = bx_name_short)) +
  geom_tile() +
  geom_hline(data = lines, aes(yintercept = x), colour = "grey95") +
  scale_x_discrete(limits = rev) +
  scale_y_reverse(
    breaks = seq(2000, 2026, by = 2),
    limits = c(1999.5, 2026.5),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = bx_colours) +
  labs(y = NULL, x = NULL, fill = NULL) +
  guides(
    fill = guide_legend(
      nrow = 1,
      reverse = TRUE
      # override.aes = list(size = 5)
    )
  ) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey90", colour = "grey90"),
    axis.text.y = element_text(
      angle = 0,
      size = ftsize,
      colour = "grey30"
    ),
    axis.text.x = element_blank(),
    # strip.text.x.bottom = element_text(angle = 60, hjust = 0.5, vjust = 0.65, size = ftsize*.5),
    strip.text.x.bottom = element_text(colour = "grey30", size = ftsize * 0.25),
    legend.position = "bottom",
    legend.title = element_text(colour = "grey30", size = ftsize),
    legend.text = element_text(colour = "grey30", size = ftsize * 0.8)
    # legend.key.width = unit(1.2, "cm"),
    # legend.key.height = unit(1.2, "cm")
  ) +
  facet_grid(
    cols = vars(country),
    switch = "x",
    scales = "free_x",
    labeller = stagger_labels
  )
main_plot
full <- fda_band +
  main_plot +
  plot_layout(widths = c(3, 55)) &
  theme(
    # panel.background = element_rect(fill="transparent", colour=NA_character_), #transparent panel bg
    # plot.background = element_rect(fill="transparent", colour=NA_character_) #transparent plot bg
  )

preview_and_save(full, w = 34, h = 17, dpi = 300)
preview_and_save(
  full,
  "ers_bx_bycountry_h_no_outliers",
  w = 34,
  h = 17,
  dpi = 300
)

## By region ----
# Vertical plot
fda_band <- year_labels |>
  group_by(label, year) |>
  mutate(
    total = n(),
    y = (row_number() - 0.5) / total * 0.5 # scales to fit within 0-0.5 regardless of n
  ) |>
  ungroup() |>
  ggplot(aes(x = year, y = y, fill = bx_name_short)) +
  geom_tile(aes(height = 0.5 / total, , width = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    limits = c(1999.5, 2026.5),
    breaks = seq(2000, 2026, by = 2)
  ) +
  scale_fill_manual(values = bx_colours) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_text(angle = 0, size = ftsize * 0.45)
  ) +
  facet_grid(rows = vars(label), switch = "y")

# Vertical plot by region
c("Asia", "Americas", "Europe") |>
  map(\(reg) {
    main_plot <- df_countries2 |>
      filter(region == reg) |>
      mutate(
        bx_name_short = factor(bx_name_short, levels = bx_levels)
      ) |>
      ggplot(aes(x = all_years, y = bx_name_short, fill = bx_name_short)) +
      geom_tile() +
      geom_vline(data = lines, aes(xintercept = x), colour = "grey95") +
      scale_y_discrete(limits = rev) +
      scale_x_continuous(
        breaks = seq(2000, 2026, by = 2),
        limits = c(1999.5, 2026.5),
        expand = c(0, 0)
      ) +
      scale_fill_manual(values = bx_colours) +
      labs(x = NULL, y = NULL, fill = NULL) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey90", colour = "grey90"),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = ftsize
        ),
        axis.text.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = ftsize * 0.75),
        legend.title = element_text(size = ftsize),
        legend.text = element_text(size = ftsize)
      ) +
      facet_grid(rows = vars(country), switch = "y", scales = "free_y")

    full <- main_plot / fda_band + plot_layout(heights = c(40, 3))
    # theme(
    #   # panel.background = element_rect(fill="transparent", colour=NA_character_), #transparent panel bg
    #   plot.background = element_rect(fill="transparent", colour=NA_character_), #transparent plot bg
    #   panel.background = element_rect(fill="transparent", colour=NA_character_)
    # )

    preview_and_save(
      full,
      paste0("ers_bx_bycountry_", reg),
      w = 34,
      h = 17,
      dpi = 300
    )
  })

# First prescriptions by year ----
bx_ct_byyear <- bx_history |>
  filter_out(
    exclusions == "Biologic start date before 2016" |
      exclusions == "Fewer than two biologics available"
  ) |>
  mutate(index_year = year(bx_start_date)) |>
  filter_out(is.na(bx_name) | is.na(index_year)) |>
  # Biologic initiation, ie first
  filter(bx_nbr == 1) |>
  count(index_year, bx_name_short) |>
  mutate(
    bx_name_short = factor(bx_name_short, levels = bx_levels)
  )

a <- bx_ct_byyear |>
  summarise(total = sum(n), .by = index_year) |>
  ggplot(aes(x = index_year, y = total)) +
  geom_bar(stat = "identity", colour = "grey30") +
  geom_text(
    aes(label = total),
    size = ftsize * ftscale / .pt,
    colour = "grey30",
    vjust = -0.5
  ) +
  scale_y_continuous(
    limits = c(0, 1500),
  ) +
  scale_x_continuous(
    breaks = seq(2016, 2026, by = 1),
    expand = c(0, 0)
  ) +
  labs(
    y = NULL,
    x = NULL,
    title = "Number of first biologic prescriptions by year"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = ftsize, colour = "grey30"),
    axis.text.y = element_text(size = ftsize, colour = "grey30"),
    axis.text.x = element_blank()
  )

b <- bx_ct_byyear |>
  mutate(pct = number_helper(n / sum(n) * 100, digits = 1), .by = index_year) |>
  ggplot(aes(x = index_year, y = n, fill = bx_name_short)) +
  geom_bar(position = "fill", stat = "identity") +
  # geom_text(
  #   aes(label =  paste0(pct,"%"), y = n),
  #   size = ftsize * ftscale / .pt,
  #   colour = "grey30",
  #   stat = "identity",
  #   position = position_fill(vjust = 0.5)
  # ) +
  scale_y_continuous(
    # limits = c(0, 1),
    labels = scales::label_percent(),
    expand = expansion(mult = c(0.025, 0))
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  ) +
  scale_x_continuous(
    breaks = seq(2016, 2026, by = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = bx_colours) +
  labs(
    y = NULL,
    x = NULL,
    fill = NULL,
    title = "First biologic prescription by year"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = ftsize, colour = "grey30"),
    axis.text = element_text(size = ftsize, colour = "grey30")
  )

bx_initiation <- a /
  b +
  plot_layout(heights = c(1, 3)) &
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(fill = "transparent", colour = NA_character_) #transparent plot bg
  )

bx_initiation
ggsave(
  here("figures", paste0("ers_bx_initiation.png")),
  bx_initiation,
  width = 33,
  height = 19,
  units = "cm",
  dpi = 300
)

# First prescriptions by country ----
bx_ct_bycountry <- bx_history |>
  filter_out(
    exclusions == "Biologic start date before 2016" |
      exclusions == "Fewer than two biologics available"
  ) |>
  filter_out(is.na(bx_name)) |>
  # Biologic initiation, ie first
  filter(bx_nbr == 1) |>
  count(country, bx_name_short) |>
  # add_count(country, wt = n, name = "N")
  mutate(prop = n / sum(n), .by = country) |>
  mutate(
    bx_name_short = factor(bx_name_short, levels = bx_levels)
  )

country_order <- bx_ct_bycountry |>
  select(country, bx_name_short, prop) |>
  pivot_wider(
    names_from = bx_name_short,
    values_from = prop,
    values_fill = 0
  ) |>
  arrange(
    desc(Omalizumab),
    desc(Mepolizumab),
    desc(Reslizumab),
    desc(Dupilumab),
    desc(Benralizumab),
    desc(Tezepelumab)
  ) |>
  pull(country)

a <- bx_ct_bycountry |>
  mutate(country = factor(country, levels = country_order)) |>
  summarise(total = sum(n), .by = country) |>
  ggplot(aes(x = country, y = total)) +
  geom_bar(stat = "identity", colour = "grey30") +
  geom_text(
    aes(label = total),
    size = ftsize * ftscale / .pt,
    colour = "grey30",
    vjust = -0.5
  ) +
  scale_y_continuous(
    limits = c(0, 3000),
    expand = expansion(add = c(0, 300))
  ) +
  # scale_x_continuous(
  #   breaks = seq(2016, 2026, by = 1),
  #   expand = c(0, 0)
  # ) +
  labs(
    y = NULL,
    x = NULL,
    title = "Number of first biologic prescriptions by country"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = ftsize, colour = "grey30"),
    axis.text.y = element_text(size = ftsize, colour = "grey30"),
    axis.text.x = element_blank()
  )
a
b <- bx_ct_bycountry |>
  mutate(country = factor(country, levels = country_order)) |>
  # mutate(pct = number_helper(n / sum(n) * 100, digits = 1), .by = country) |>
  ggplot(aes(x = country, y = n, fill = bx_name_short)) +
  geom_bar(position = "fill", stat = "identity") +
  # geom_text(
  #   aes(label =  paste0(pct,"%"), y = n),
  #   size = ftsize * ftscale / .pt,
  #   colour = "grey30",
  #   stat = "identity",
  #   position = position_fill(vjust = 0.5)
  # ) +
  scale_y_continuous(
    # limits = c(0, 1),
    labels = scales::label_percent(),
    expand = expansion(mult = c(0.025, 0))
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  ) +
  # scale_x_continuous(
  #   breaks = seq(2016, 2026, by = 1),
  #   expand = c(0, 0)
  # ) +
  scale_fill_manual(values = bx_colours) +
  labs(
    y = NULL,
    x = NULL,
    fill = NULL,
    title = "First biologic prescription by country"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.position.inside = c(0.85, 0.85),
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize, colour = "grey30"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    axis.text = element_text(
      size = ftsize,
      colour = "grey30",
      angle = 45,
      vjust = 0,
      hjust = 1
    ),
    axis.title = element_text(size = ftsize, colour = "grey30")
  )
b
bx_initiation <- a /
  b +
  plot_layout(heights = c(1, 3)) &
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(fill = "transparent", colour = NA_character_) #transparent plot bg
  )

bx_initiation
ggsave(
  here("figures", paste0("ers_bx_initiation_bycountry.png")),
  bx_initiation,
  width = 33,
  height = 19,
  units = "cm",
  dpi = 300
)

# Time to cessation of treatment ----
# GS - Don't show Ns at risk or events for ERS
df_survival <- bx_history |>
  mutate(index_year = year(bx_start_date)) |>
  filter_out(is.na(bx_name_short) | is.na(index_year)) |>
  filter_out(bx_name_short == "Treatment gap") |>
  filter(
    when_any(
      exclusions == "Eligible",
      exclusions == "Less than 12 months of follow-up"
    )
  ) |>
  # 20260902 exclude patient with 2017 Teze start date
  filter_out(patient_id == "CAN_S11_14") |>
  # 20260902 exclude NJH for now
  filter_out(str_detect(patient_id, "NJH")) |>
  group_by(patient_id) |>
  # Treatment gaps
  # filter_out(any(tx_gap_row)) |>
  # Keep only first two switches for everybody, including those with multiple switches
  filter(bx_nbr <= 2) |>
  arrange(patient_id, bx_nbr) |>
  # censoring for switch includes stopping as well as loss to followup, censoring for stop includes switch as well as loss to follow-up
  mutate(
    # Flag switchers - ok to use logical with Surv()
    switch = any(bx_nbr == 2),
    # Assume if bx_ongoing is missing or maintain then they continue
    stop = case_when(
      bx_nbr == 1 & bx_ongoing == "Stop" & !switch ~ TRUE,
      .default = FALSE
    ),
    switchstop = case_when(
      switch | stop ~ TRUE,
      .default = FALSE
    ),
    outcome = case_when(
      switch ~ "Switch",
      stop ~ "Stop",
      .default = "Censor"
    ) |>
      factor(levels = c("Censor", "Switch", "Stop")),
    index_year_cat = case_when(
      index_year <= 2019 ~ "2016-2019",
      index_year <= 2021 ~ "2020-2021",
      .default = "2022-2026"
    ) |>
      factor(levels = c("2016-2019", "2020-2021", "2022-2026"))
  ) |>
  filter_out(bx_nbr == 2) |>
  ungroup() |>
  select(
    country,
    patient_id,
    index_year,
    index_year_cat,
    bx_ct,
    bx_nbr,
    bx_switchreason,
    bx_name_short,
    bx_type,
    contains("date"),
    tx_months,
    bx_ongoing,
    bx_ongoing_imputed,
    bx_ongoing_imputed_flag,
    switch,
    stop,
    switchstop,
    outcome
  )

library(survival)

## Composite outcome KM ----
biologic_switchstop <- survfit2(
  Surv(tx_months, switchstop) ~ bx_name_short,
  data = df_survival
)
#Get median times
biologic_switchstop |>
  tidy() |>
  filter(estimate < 0.5) |>
  # Returns smallest time below 0.5
  slice_min(time, n = 1, with_ties = FALSE, by = strata) |>
  mutate(
    strata = str_split_i(strata, "=", 2),
    time = number_helper(time, digits = 1)
  ) |>
  select(strata, time) |>
  kable() |>
  kable_styling()

year_switchstop <- survfit2(
  Surv(tx_months, switchstop) ~ index_year_cat,
  data = df_survival
)
#Get median times
year_switchstop |>
  tidy() |>
  filter(estimate < 0.5) |>
  # Returns smallest time below 0.5
  slice_min(time, n = 1, with_ties = FALSE, by = strata) |>
  mutate(
    strata = str_split_i(strata, "=", 2),
    time = number_helper(time, digits = 1)
  ) |>
  select(strata, time) |>
  kable() |>
  kable_styling()

biologic_switchstop_19 <- survfit2(
  Surv(tx_months, switchstop) ~ bx_name_short,
  data = df_survival |> filter(index_year <= 2019)
)

biologic_switchstop_21 <- survfit2(
  Surv(tx_months, switchstop) ~ bx_name_short,
  data = df_survival |> filter(index_year > 2019 & index_year <= 2021)
)
biologic_switchstop_26 <- survfit2(
  Surv(tx_months, switchstop) ~ bx_name_short,
  data = df_survival |> filter(index_year > 2021)
)

km_plot <- function(fit, bx = TRUE) {
  # Pull the last observed point of each curve for label placement
  label_df <- broom::tidy(fit) |>
    # Truncate for labeling
    filter(when_all(time <= 60, !is.na(estimate))) |>
    slice_max(time, n = 1, by = strata) |>
    mutate(labeltxt = str_split_i(strata, "=", 2), risk = 1 - estimate)

  plot <- fit |>
    ggsurvfit("survival") +
    # add_confidence_interval() +
    # add_risktable(
    #   risktable_stats = c(
    #     "n.risk",
    #     # "cum.censor",
    #     "cum.event"
    #   )
    # ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = time, y = estimate, label = labeltxt, colour = labeltxt),
      # xlim = c(max(label_df$time) * 1.02, NA),  # push labels to the right of the plot
      fontface = "bold",
      size = ftsize / .pt,
      direction = "y",
      # xlim = c(5.5, NA),
      # hjust = 1,
      segment.size = .7,
      segment.alpha = .5,
      segment.linetype = "dotted",
      box.padding = .4,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20,
      max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
      show.legend = F
    ) +
    theme_ggsurvfit_KMunicate() +
    # Restrict to first five years
    coord_cartesian(xlim = c(0, 60), ylim = c(0, 1)) +
    scale_x_continuous(
      # expand = expansion(mult = c(0.02, 0)),
      limits = c(0, 60),
      breaks = seq(0, 60, by = 12),
      # labels = 1:5
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 1),
      breaks = seq(0, 1, by = .1),
      labels = scales::label_percent()
    ) +
    theme(
      # aspect.ratio = 1,
      legend.position = "none",
      axis.text = element_text(size = ftsize, colour = "grey30"),
      axis.title = element_text(size = ftsize, colour = "grey30")
    ) +
    labs(
      x = "Months since biologic initiation",
      y = "Continuing initial treatment (%)"
    )

  if (bx) {
    plot <- plot +
      scale_colour_manual(values = bx_colours)

    plot
  } else {
    plot
  }
}

preview_and_save(
  km_plot(biologic_switchstop),
  "ers_km_bx_switchstop_nonjh",
  h = 19,
  w = 17,
  dpi = 300
)

preview_and_save(
  km_plot(biologic_switchstop_19),
  "ers_km_bx_switchstop_nonjh_19",
  h = 19,
  w = 17,
  dpi = 300
)
preview_and_save(
  km_plot(biologic_switchstop_21),
  "ers_km_bx_switchstop_nonjh_21",
  h = 19,
  w = 17,
  dpi = 300
)
preview_and_save(
  km_plot(biologic_switchstop_26),
  "ers_km_bx_switchstop_nonjh_26",
  h = 19,
  w = 17,
  dpi = 300
)

preview_and_save(
  km_plot(year_switchstop, bx = F),
  "ers_km_yr_switchstop_nonjh",
  h = 19,
  w = 17,
  dpi = 300
)

# Baseline characteristics ----
# Use SHINE dataset for now
datapath <- 'C:/Users/William.Oswald/Optimum Patient Care Global Limited/ISAR Team - Documents/1. Research Projects/Active/SHINE/3. Data/2. Dataset'

# Load list of analysis datasets for each definition
load(here(datapath, "shine_dataset_list.Rdata"))

table_vars <- c(
  "setting",
  "index_date",
  "bx_type",
  "index_age",
  "sex",
  "ethnicity",
  "smoking_status",
  "bmi",
  "obese",
  "exac_bl",
  "ltocs_yn_bl",
  "cont_bl",
  "fev1_percpred_bl",
  "fev1_percpred_bl_cat",
  "fev1_fvc_ratio_bl",
  "fev1_fvc_ratio_bl_cat",
  "age_asthma_onset",
  "age_asthma_onset_cat",
  "asthma_duration",
  "asthma_duration_cat",
  "bec_bl",
  "ige_bl",
  "feno_bl",
  "allergy",
  "eosinophilic_grade",
  "ml_ics_bl",
  "ml_ics_bl_cat",
  "lama_bl",
  "ltra_bl",
  "theophylline_bl",
  "macrolide_bl",
  "ar",
  "crs",
  "np",
  "ad",
  "osteoporosis_bl",
  "anxdep_bl"
)

bl_df <- readRDS(here(datapath, "SHINE_clean.rds")) |>
  select(patient_id, any_of(table_vars), -bx_type) |>
  distinct()

bl_analysis <- df_survival |>
  left_join(bl_df, by = join_by(patient_id)) |>
  mutate(
    outcome = case_when(
      !switch & !stop ~ "Continue",
      switch ~ "Switch",
      stop ~ "Stop"
    ) |>
      factor(levels = c("Continue", "Switch", "Stop"))
  )

# Steps to tidy up summary tables before converting to flextable
format_table <- \(df) {
  df |>
    filter(
      value != "No",
      value != "Negative",
      value != "Never",
      value != "Male",
      value != "≥80%",
      value != "≥0.7"
    ) |>
    mutate(
      value = case_when(
        !str_detect(value, "median") ~ paste0(value, ": n (%)"),
        TRUE ~ value
      ),
      label = str_replace_all(
        label,
        c(
          "Bmi" = "BMI",
          "Kg/M2" = "kg/m^2",
          "Fev1" = "FEV1",
          "Fvc" = "FVC",
          "Gina" = "GINA",
          "Ltocs" = "LTOCS",
          "Feno" = "FENO",
          "Ige" = "IgE",
          "Ics" = "ICS",
          "Lama" = "LAMA",
          "Ltra" = "LTRA",
          "Ppb" = "ppb",
          "Cells/Mcl" = "cells/mcL",
          "Iu/Ml" = "IU/mL"
        )
      ),
      label = str_replace(label, " By Biologic Initiation Date", ""),
      label = str_replace(label, " Up To Biologic Initiation", ""),
      # label = str_replace(label, " In Past Year", ""),
      label = str_replace(label, ", Continuous", ""),
      label = str_replace(
        label,
        " In The Year Preceding Biologic Initiation",
        ""
      ),
      label = str_replace(label, "History of ", ""),
      label = str_replace(label, " Closest To Biologic Initiation", ""),
      label = case_when(
        varname == "ml_ics_bl" ~ "ICS dose (BDP-equivalent mcg/day)",
        varname == "ml_ics_bl_cat" ~ "ICS dose (BDP-equivalent mcg/day)",
        TRUE ~ label
      )
      #     label == "setting_short_list" ~ "Geographic setting*",
      #     label == "Age At Biologic Initiation" ~ "Age At Biologic Initiation (Years)",
      #     label == "exac_bl" ~ "Asthma Exacerbations In Past Year",
      #     label == "obese" ~ "Obese (BMI ≥30 kg/m^2)",
      #     label == "age_asthma_onset_cat" ~ "Age At Asthma Onset (Years)",
      #     label == "asthma_duration_cat" ~ "Asthma Duration (Years)",
      #     label == "Chronic Rhinosinusitis" ~ "CRS (+/-NP)",
      #     label == "Nasal Polyposis" ~ "CRSwNP",
      #     label == "Chronic Rhinosinusitis" ~ "CRS (+/-NP)",
      #     label == "fev1_percpred_bl_cat" ~ "Most Recent Percent Predicted FEV1 (%)",
      #     label == "fev1_fvc_ratio_bl_cat" ~ "Most Recent FEV1/FVC Ratio",
      #     label == "exac_bl" ~ "Number Of Asthma Exacerbations",
      #     TRUE ~ label
      # )
    ) |>
    group_by(label) |>
    mutate(
      label = case_when(
        1:n() != 1 ~ "",
        TRUE ~ label
      ),
      across(
        contains("n"),
        ~ case_when(
          1:n() != 1 ~ NA,
          TRUE ~ .
        )
      )
    ) |>
    ungroup() |>
    select(-varname)
}

bl_table <- table_vars |>
  map_df(summary_table2, group_var = "outcome", data = bl_analysis)
format_table() |>
  flextable() |>
  align(align = "right", part = "body") |>
  align(align = "center", part = "header") |>
  # merge_v(j = "varname", part = "body") |>
  valign(j = "label", valign = "top", part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Aptos", part = "all") |>
  bold(bold = TRUE, part = "header") |>
  bold(i = NULL, j = "label", bold = TRUE, part = "body") |>
  align(i = NULL, j = "label", align = "left", part = "body") |>
  set_table_properties(layout = "autofit") |>
  height(height = 0.7, part = "body", unit = "cm") |>
  set_header_labels(
    "label" = "Characteristic",
    "n" = "N",
    "value" = "Measure/category",
    "result" = "Result"
  ) |>
  # Remove spacing before and after paragraph
  padding(padding.top = 0, padding.bottom = 0, part = "all") |>
  border_remove() |>
  hline_top(part = "body") |>
  hline_bottom(part = "body")
bl_table

# Patterns ----
# Based on 20260415
# Shows switch from first biologic to second
pattern_analysis <- bx_history |>
  filter_out(
    exclusions == "Biologic start date before 2016" |
      exclusions == "Fewer than two biologics available"
  ) |>
  filter(!is.na(bx_nbr), bx_ct > 1) |>
  mutate(
    # index_year = min(year(bx_start_date)),
    index_bx_name = first(bx_name_short),
    index_bx_type = first(bx_type),
    .by = patient_id
  ) |>
  #   By class
  # count(index_bx_type, bx_nbr, bx_type) |>
  #   By product
  count(index_bx_name, bx_nbr, bx_name_short) |>
  filter_out(bx_nbr == 1) |>
  mutate(
    index_bx_name = factor(index_bx_name, levels = bx_levels),
    bx_name_short = factor(bx_name_short, levels = bx_levels)
  )

# switch1 <- pattern_analysis |>
#   # First switch
#   filter(bx_nbr == 2) |>
#   add_count(index_bx_type, wt = n, name = "total") |>
#   mutate(group_label = glue::glue("{index_bx_type} N={total}")) |>
#   ggplot(aes(y = bx_type, x = n, fill = bx_type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_brewer(palette = "Set2") +
#   scale_y_discrete(limits = rev) +
#   scale_x_continuous(expand = c(0,0)) +
#   labs(y = "Next biologic class", x = "Number of patients", fill = NULL) +
#   facet_wrap(vars(group_label), nrow = 5) +
#   theme_minimal() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     legend.position = "none",
#     axis.text = element_text(size = 12, colour = "grey30"),
#     axis.title = element_text(size = 12, colour = "grey30"),
#     strip.text = element_text(size = 12, colour = "grey30")
#   )
# switch1

switch1_h <- pattern_analysis |>
  # First switch
  filter(bx_nbr == 2) |>
  arrange(desc(index_bx_name), desc(bx_name_short)) |>
  add_count(index_bx_name, wt = n, name = "total") |>
  mutate(
    group_label = glue::glue("{index_bx_name} N={total}"),
    group_label = factor(group_label, levels = unique(group_label)),
    prop = n / total
  ) |>
  ggplot(aes(x = bx_name_short, y = prop, fill = bx_name_short)) +
  # geom_bar(aes(x = 3.5, y = total, fill = index_bx_name),
  #          stat = "identity",
  #          width = 6,
  #          alpha = 0.7) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bx_colours) +
  # scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(limits = rev) +
  # scale_y_continuous(
  #   expand = expansion(add = c(0, 50))
  # ) +
  scale_y_continuous(
    expand = expansion(add = c(0, 0)),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::label_percent()
  ) +
  labs(
    x = "First biologic",
    y = "Percentage of patients",
    fill = "Next biologic",
    title = NULL
  ) +
  facet_grid(cols = vars(group_label), switch = "x") +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = TRUE
      # override.aes = list(size = 5)
    )
  ) +
  # theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize * 1.2, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1, size = ftsize, colour = "grey30"),
    axis.text.y = element_text(size = ftsize, colour = "grey30"),
    axis.title = element_text(size = ftsize, colour = "grey30"),
    strip.text = element_text(size = ftsize, colour = "grey30"),
    # legend.key.width = unit(1, "cm"),
    # legend.key.height = unit(1, "cm"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    panel.spacing = unit(2, "lines")
    # strip.text.y.bottom = element_text()
  )
switch1_h

preview_and_save(
  switch1_h,
  "ers_bx_switch1_h",
  w = 40,
  h = 19,
  dpi = 300
)

# Shows switch from second biologic to third, among double switchers
pattern_analysis <- bx_history |>
  filter_out(
    exclusions == "Biologic start date before 2016" |
      exclusions == "Fewer than two biologics available"
  ) |>
  filter(!is.na(bx_nbr), bx_ct > 2) |>
  filter_out(bx_nbr == 1) |>
  mutate(
    # index_year = min(year(bx_start_date)),
    index_bx_name = first(bx_name_short),
    index_bx_type = first(bx_type),
    .by = patient_id
  ) |>
  #   By class
  # count(index_bx_type, bx_nbr, bx_type) |>
  #   By product
  count(index_bx_name, bx_nbr, bx_name_short) |>
  mutate(
    index_bx_name = factor(index_bx_name, levels = bx_levels),
    bx_name_short = factor(bx_name_short, levels = bx_levels)
  )

switch2_h <- pattern_analysis |>
  # First switch
  filter(bx_nbr == 3) |>
  arrange(desc(index_bx_name), desc(bx_name_short)) |>
  add_count(index_bx_name, wt = n, name = "total") |>
  mutate(
    group_label = glue::glue("{index_bx_name} N={total}"),
    group_label = factor(group_label, levels = unique(group_label)),
    prop = n / total
  ) |>
  ggplot(aes(x = bx_name_short, y = prop, fill = bx_name_short)) +
  # geom_bar(aes(x = 3.5, y = total, fill = index_bx_name),
  #          stat = "identity",
  #          width = 6,
  #          alpha = 0.7) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bx_colours) +
  # scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(limits = rev) +
  # scale_y_continuous(
  #   expand = expansion(add = c(0, 50))
  # ) +
  scale_y_continuous(
    expand = expansion(add = c(0, 0)),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::label_percent()
  ) +
  labs(
    x = "Second biologic",
    y = "Percentage of patients",
    fill = "Next biologic",
    title = NULL
  ) +
  facet_grid(cols = vars(group_label), switch = "x") +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = TRUE
      # override.aes = list(size = 5)
    )
  ) +
  # theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent panel bg
    plot.background = element_rect(
      fill = "transparent",
      colour = NA_character_
    ), #transparent plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = ftsize, colour = "grey30"),
    legend.title = element_text(size = ftsize * 1.2, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1, size = ftsize, colour = "grey30"),
    axis.text.y = element_text(size = ftsize, colour = "grey30"),
    axis.title = element_text(size = ftsize, colour = "grey30"),
    strip.text = element_text(size = ftsize, colour = "grey30"),
    # legend.key.width = unit(1.2, "cm"),
    # legend.key.height = unit(1.2, "cm"),
    plot.title = element_text(
      size = ftsize,
      colour = "grey30",
      hjust = 0.5,
      face = "bold"
    ),
    panel.spacing = unit(2, "lines")
    # strip.text.y.bottom = element_text()
  )
switch2_h

preview_and_save(
  switch2_h,
  "ers_bx_switch2_h",
  w = 40,
  h = 19,
  dpi = 300
)

# Alluvial plot ----
# Shows all transitions, maybe start with this?
source("r_sankey4.R")

df_alluvial <- bx_history |>
  filter_out(
    exclusions == "Biologic start date before 2016" |
      exclusions == "Fewer than two biologics available"
  ) |>
  filter(!is.na(bx_nbr), bx_ct > 1) |>
  mutate(index_year = min(year(bx_start_date)), .by = patient_id) |>
  mutate(
    bx_type2 = case_when(
      str_detect(bx_type, "IL5") ~ "Anti-IL5/5R",
      .default = bx_type
    )
  )

alluvial_bx_name <- sankeytime(
  df_alluvial,
  "bx_name_short",
  time_var = "bx_nbr",
  show_prevalence = F,
  show_flow_counts = F
)

preview_and_save(
  alluvial_bx_name,
  "ers_alluvial_biologic_px",
  w = 33,
  h = 19,
  dpi = 300
)

## By class ----
alluvial_bx_type <- sankeytime(
  df_alluvial,
  "bx_type2",
  time_var = "bx_nbr",
  show_prevalence = F,
  show_flow_counts = F
)
alluvial_bx_type
preview_and_save(
  alluvial_bx_type,
  "ers_alluvial_bx_type",
  w = 33,
  h = 19,
  dpi = 300
)

## By class, 2022+ ----
alluvial_bx_type_22 <- sankeytime(
  df_alluvial |> filter(index_year > 2021),
  "bx_type2",
  time_var = "bx_nbr",
  show_prevalence = F,
  show_flow_counts = F
)
alluvial_bx_type_22
preview_and_save(
  alluvial_bx_type_22,
  "ers_alluvial_bx_type_22",
  w = 33,
  h = 19,
  dpi = 300
)

# Treatment gaps ----
# See 20260814 tx gaps by switch for original
# Present distribution of treatment gaps between switches
switches <- bx_cln |>
  # Drop patients with errors for now
  filter(
    when_any(
      exclusions == "Eligible",
      exclusions == "Less than 12 months of follow-up"
    )
  ) |>
  select(
    -bx_flag,
    -contains("err"),
    # -contains("exclu")
  ) |>
  arrange(patient_id, bx_start_date) |>
  mutate(
    bx_row_nbr = row_number(),
    bx_ct = n(),
    .by = patient_id
  ) |>
  # Restrict to patients with >1 biologic
  filter(bx_ct > 1) |>
  select(-exclusions) |>
  arrange(patient_id, bx_row_nbr, .by = patient_id) |>
  pivot_wider(
    id_cols = 1:9,
    names_from = bx_row_nbr,
    values_from = contains("bx_"),
    names_vary = "slowest"
  ) |>
  mutate(
    switch_gap_1 = interval(bx_end_date_1, bx_start_date_2) / months(1),
    switch_gap_2 = interval(bx_end_date_2, bx_start_date_3) / months(1),
    switch_gap_3 = interval(bx_end_date_3, bx_start_date_4) / months(1),
    switch_gap_4 = interval(bx_end_date_4, bx_start_date_5) / months(1),
    switch_gap_5 = interval(bx_end_date_5, bx_start_date_6) / months(1)
  ) |>
  select(contains("switch_gap"), bx_start_year_1, country) |>
  pivot_longer(
    cols = -c(country, bx_start_year_1),
    names_to = "switch",
    names_prefix = "switch_gap_",
    values_to = "gap_months"
  ) |>
  filter_out(is.na(gap_months)) |>
  mutate(
    switch = factor(switch, levels = rev(c("1", "2", "3", "4", "5"))),
    index_year = factor(bx_start_year_1)
  )

distfunc(data = switches, var = "gap_months", groupvar = "switch")

gap_plot <- function(varname, y_axis_name) {
  switches |>
    ggplot(aes(
      x = gap_months,
      y = {{ varname }}
    )) +
    geom_point(
      aes(
        fill = {{ varname }},
        colour = {{ varname }}
      ),
      position = position_jitterdodge(dodge.width = 0.60),
      alpha = 0.2,
      size = 3,
      stroke = 0,
      na.rm = TRUE
    ) +
    # https://stackoverflow.com/questions/52217852/how-to-add-means-to-a-ggplot-geom-point-plot
    stat_summary(
      fun = "mean",
      geom = "point",
      col = "grey30",
      size = 3,
      shape = 24,
      fill = "darkred"
    ) +
    geom_boxplot(
      width = 0.3,
      alpha = 0.2,
      na.rm = TRUE,
      outlier.shape = NA,
      colour = "grey30"
    ) +
    geom_violin(
      width = 0.5,
      alpha = 0.2,
      na.rm = TRUE,
      colour = "grey30",
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq(0, 96, by = 12)
    ) +
    scale_y_discrete(limits = rev) +
    # scale_fill_brewer(palette = "Set2") +
    # scale_colour_brewer(palette = "Set2") +
    labs(
      x = "Months from previous biologic end date to next biologic start date",
      y = y_axis_name,
      fill = NULL,
      colour = NULL
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.position.inside = c(0.85, 0.85),
      legend.text = element_text(size = ftsize, colour = "grey30"),
      legend.title = element_text(size = ftsize, colour = "grey30"),
      # plot.title = element_text(size = ftsize, colour = "grey30", hjust = 0.5, face = "bold"),
      axis.text = element_text(size = ftsize, colour = "grey30"),
      # axis.text.y = element_blank(),
      axis.title = element_text(size = ftsize, colour = "grey30")
    )
}

preview_and_save(
  gap_plot(switch, y_axis_name = "Switch number"),
  "ers_tg_by_switch",
  w = 18,
  h = 18,
  dpi = 300
)
preview_and_save(
  gap_plot(index_year, y_axis_name = "First biologic initiation year"),
  "ers_tg_by_year",
  w = 18,
  h = 18,
  dpi = 300
)
preview_and_save(
  gap_plot(country, y_axis_name = "Country"),
  "ers_tg_by_country",
  w = 18,
  h = 18,
  dpi = 300
)


# Follow-up time ----
# 20260501
# How much follow-up is available after a stop or a switch?

df_followup <- bx_history |>
  mutate(index_year = year(bx_start_date)) |>
  filter_out(is.na(bx_name_short) | is.na(index_year)) |>
  filter_out(bx_name_short == "Treatment gap") |>
  filter(exclusions == "Eligible") |>
  group_by(patient_id) |>
  # Treatment gaps
  # filter_out(any(tx_gap_row)) |>
  # Keep only first two biologics for everybody, including those with multiple switches
  filter(bx_nbr <= 2) |>
  arrange(patient_id, bx_nbr) |>
  mutate(
    # Flag switchers - ok to use logical with Surv()
    switch = any(bx_nbr == 2),
    # Assume if bx_ongoing is missing or maintain then they continue
    stop = case_when(
      bx_nbr == 1 & bx_ongoing == "Stop" & !switch ~ TRUE,
      .default = FALSE
    ),
    switchstop = case_when(
      switch | stop ~ TRUE,
      .default = FALSE
    ),
    outcome = case_when(
      switch ~ "Switch",
      stop ~ "Stop",
      .default = "Continue"
    ) |>
      factor(levels = c("Continue", "Switch", "Stop")),
    follow_up_mos = case_when(
      outcome == "Continue" ~ interval(
        start = bx_start_date,
        end = last_visit_date
      ) /
        months(1),
      .default = interval(start = bx_end_date, end = last_visit_date) /
        months(1)
    ),
    follow_up_mos2 = case_when(
      outcome == "Continue" ~ interval(
        start = bx_start_date,
        end = last_visit_date
      ) /
        months(1),
      .default = interval(start = bx_end_date2, end = last_visit_date) /
        months(1)
    )
  ) |>
  filter_out(bx_nbr == 2) |>
  ungroup() |>
  select(-contains("tx_gap")) |>
  # 20260903 Recode negative followup to missing for ERS slides, mostly NJH
  mutate(across(contains("follow_up"), ~ if_else(.x < 0, NA, .x)))

# Distribution of followup
c("follow_up_mos") |>
  map_dfr(\(x) distfunc(data = df_followup, var = x, groupvar = "outcome"))
# kable() |>
# kableExtra::kable_styling()

# fu <- df_followup |>
#   # filter_out(outcome == "Continue") |>
#   ggplot(aes(x = follow_up_mos, fill = outcome, colour = outcome)) +
#   geom_density(alpha = 0.3) +
#   scale_fill_brewer(palette = "Set2", name = "Outcome") +
#   scale_colour_brewer(palette = "Set2", name = "Outcome") +
#   labs(
#     x = "Months from biologic end date to date of last recorded visit",
#     fill = NULL,
#     colour = NULL
#   ) +
#   theme_minimal() +
#   theme(
#     aspect.ratio = 1
#   )
# fu

fu <- df_followup |>
  # filter_out(outcome == "Continue") |>
  arrange(outcome) |>
  ggplot(aes(
    x = follow_up_mos,
    y = outcome,
    fill = outcome,
    colour = outcome
  )) +
  geom_boxplot(alpha = 0.5) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = seq(0, 136, by = 12)) +
  scale_fill_brewer(palette = "Set2", name = "Outcome") +
  scale_colour_brewer(palette = "Set2", name = "Outcome") +
  labs(
    x = "Months from biologic end date (or start date for continuers) to date of last recorded visit",
    fill = NULL,
    colour = NULL
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )
fu
preview_and_save(fu, "ers_fup_by_outcome", w = 25, h = 18, dpi = 300)

# Base cohort ----
# For basecohort origins see 20260512
basecohort <- bx_cln |>
  # Drop patients with errors for now
  filter(
    exclusions == "Eligible"
  ) |>
  select(
    -bx_flag,
    -contains("err"),
    # -contains("exclu")
  ) |>
  arrange(patient_id, bx_start_date) |>
  mutate(
    bx_row_nbr = row_number(),
    bx_ct = n(),
    .by = patient_id
  ) |>
  # Filter out later biologics, >2
  filter_out(bx_row_nbr > 2) |>
  select(-exclusions) |>
  group_by(patient_id) |>
  arrange(patient_id, bx_row_nbr) |>
  pivot_wider(
    id_cols = 1:9,
    names_from = bx_row_nbr,
    values_from = contains("bx_"),
    names_vary = "slowest"
  ) |>
  ungroup() |>
  rename_with(~ paste0("first_", str_remove(.x, "_1")), ends_with("_1")) |>
  rename_with(~ paste0("switch_", str_remove(.x, "_2")), ends_with("_2")) |>
  rename(bx_ct = first_bx_ct) |>
  mutate(
    # Flag switchers (all dose records w/in 3 months should have been combined), >2 biologic treatments
    switch = case_when(
      switch_bx_row_nbr == 2 ~ TRUE,
      .default = FALSE
    ),
    # It is OK to assume if bx_ongoing is missing then they continue, Aaron will update this in dataset
    stop = case_when(
      first_bx_ongoing == "Stop" & !switch ~ TRUE,
      .default = FALSE
    ),
    switch_date = if_else(switch, switch_bx_start_date, NA),
    stop_date = if_else(stop, first_bx_end_date, NA),
    outcome = case_when(
      !switch & !stop ~ "Continue",
      stop ~ "Stop",
      switch ~ "Switch"
    ) |>
      factor(levels = c("Continue", "Switch", "Stop")),
    treatment_mos = case_when(
      outcome == "Continue" ~ interval(
        start = first_bx_start_date,
        end = last_visit_date
      ) /
        months(1),
      .default = interval(
        start = first_bx_start_date,
        end = first_bx_end_date
      ) /
        months(1)
    ),
    treatment_gap_mos = case_when(
      outcome == "Switch" ~ interval(
        start = first_bx_end_date,
        end = switch_date
      ) /
        months(1),
      .default = NA
    ),
    # Ignores treatment gaps for switchers
    treatment_mos_alt = case_when(
      outcome == "Switch" ~ interval(
        start = first_bx_start_date,
        end = switch_date
      ) /
        months(1),
      .default = treatment_mos
    ),
    follow_up_mos = case_when(
      outcome == "Continue" ~ interval(
        start = first_bx_start_date,
        end = last_visit_date
      ) /
        months(1),
      outcome == "Switch" ~ interval(
        start = switch_date,
        end = last_visit_date
      ) /
        months(1),
      outcome == "Stop" ~ interval(
        start = stop_date,
        end = last_visit_date
      ) /
        months(1)
    ),
    gap_date_imputed = case_when(
      first_bx_end_date_imputed == 1 |
        switch_bx_start_date_imputed == 1 ~ "Imputed date",
      .default = "No dates imputed"
    ) |>
      factor(levels = c("Imputed date", "No dates imputed"))
  ) |>
  select(
    -contains("row_nbr"),
    -contains("class"),
    -switch_bx_name,
    -switch_bx_ct,
    -first_bx_name,
    -contains("name_imputed"),
    -contains("ongoing_imputed"),
    -contains("switchreason_imputed")
  )
