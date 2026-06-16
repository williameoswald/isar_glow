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
versiondate <- "20260603"

# FUNCTIONS ----

# Load helper functions
source("functions.R")

source(
  "C:/Users/William.Oswald/OneDrive - Optimum Patient Care Global Limited/Documents/Github/useful/r_summary_table_v2.R"
)

# ASSEMBLE DATA ----

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
# Will create new objects in global environment for each dataset
filenames |>
  purrr::set_names() |>
  purrr::map(\(x) import_datasets(x, keep_raw = F), .progress = T)

# BASELINE ASSEMBLY ----
# source("1_baseline.R")
