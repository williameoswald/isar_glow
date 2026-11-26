# GLOW ISAR Switchers & Stoppers assembly and analysis ----
# William Oswald
# 20260528

# Modifications/Decisions:
# 20260311 - missing last end dates will be imputed to be last visit date

# To do:

# SET UP ----
rm(list = ls())

options(scipen = 999)

# Change directory to source folder
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
if (!require("pacman")) {
  install.packages("pacman")
}
pkgs <- c(
  "tidyverse",
  "haven",
  "labelled",
  "ggtext",
  "ggsurvfit",
  "patchwork",
  "readxl",
  "writexl",
  "here",
  "conflicted",
  "flextable",
  "kableExtra",
  "officer",
  "janitor",
  "arsenal",
  "survey",
  "binom",
  "sandwich",
  "lme4",
  "lmtest",
  "broom",
  "broom.mixed",
  "broom.helpers",
  "marginaleffects",
  "remotes",
  "rvg",
  "quarto"
) |>
  pacman::p_load(pkgs, character.only = T)
rm(pkgs)
conflicted::conflict_prefer("filter", "dplyr")

pak::pak("williameoswald/ISARanalysis")
library(ISARanalysis)

# SPECIFICATIONS ----
# Location of data folder on Sharepoint
# https://stackoverflow.com/questions/17605563/efficiently-convert-backslash-to-forward-slash-in-r
# Type pp + tab to paste copied filepath using snippet
datapath <- "C:/Users/william.oswald/Optimum Patient Care Global Limited/ISAR Team - Documents/Dataset production"
# The latest folder with a date name is the same as "latest", use dated folder for a useful version tag in dataset
versiondate <- "20260805"

# FUNCTIONS ----

# Load helper functions
source("functions.R")

source(
  "C:/Users/William.Oswald/OneDrive - Optimum Patient Care Global Limited/Documents/Github/useful/r_summary_table_v2.R"
)

# LOAD DATA ----
## SRD v1 ----
# Use this for countries with pending rebuilds
dictionary <- readxl::read_xlsx(
  here::here(datapath, "documentation", "ISAR SRD Data Dictionary v1.0.xlsx")
)

# Create list of dataset filenames (without full path)
filenames <- list.files(
  here::here(datapath, "datasets", versiondate),
  pattern = "\\.csv$"
)

# Map over all listed files or select specific datasets
# filenames <- filenames[filenames != "7. HCRU.csv"] or filenames[c(1, 4)]
# Will create new objects in "isar" environment for each dataset
isar <- new.env()
filenames |>
  purrr::set_names() |>
  purrr::map(
    \(x) import_datasets(x, keep_raw = F, output_env = isar),
    .progress = T
  )

### SRD v1 biologic data ----
source("biologics.R")

## QISAR2 ----
# See "isar_srd" scripts
# Load QISAR2/SRD data into respective environments and process
# for GLOW eligibility, objective 1 analysis, and baseline
uk <- new.env()
# See "isar_srd/assemble_research_datasets_uk.R"
load(here::here("datasets", "qisar2_srd_uk.Rdata"), envir = uk)
ita <- new.env()
# See "isar_srd/assemble_research_datasets_ita.R"
load(here::here("datasets", "qisar2_srd_ita.Rdata"), envir = ita)
prt <- new.env()
# See "isar_srd/assemble_research_datasets_prt.R"
load(here::here("datasets", "qisar2_srd_prt.Rdata"), envir = prt)
njh <- new.env()
# See "isar_srd/assemble_research_datasets_njh.R"
load(here::here("datasets", "qisar2_srd_njh.Rdata"))

### QISAR2 biologic data ----
# 1) Creates biologics analysis dataset and flag eligible patients
# 2) Assemble eligible initiation date and index (switch/stop) dates for querying SRD
source("glow_eligibility.R")

## BASELINE ASSEMBLY ----
# Query datasets for baseline measures, using bx_initiation_date
source("glow_baseline.R")

# OBJ1 ANALYSIS ----
source("glow_analysis_obj1.R")

# SWITCH/STOP ASSEMBLY ----
# Query datasets for pre/post measures, using switchstop_date
# source("glow_prepost.R")

# OBJ2 ANALYSIS ----
# source("glow_analysis_obj2.R")

# REPORT ----

rmarkdown::render(
  "GLOW_internal_results_tables_and_figures.rmd",
  encoding = "UTF-8",
  output_file = paste0(
    "GLOW_internal_results_tables_and_figures_",
    format(Sys.Date(), "%Y%m%d"),
    ".docx"
  ),
  output_dir = "report"
)
