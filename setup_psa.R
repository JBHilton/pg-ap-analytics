# In this script we load parameters and assemble all the parameter objects we
# will need for the analysis.

# Flag to determine whether to convert rates in the data to probabilities
CONVERSIONS_REQUIRED <- FALSE

# Set to true to directly calculate probabilities for the no-LOF subpopulation,
# or false to use values from other subpopulations
UNIQUE_NO_LOF_PROBS <- FALSE

# Step correction for death in decision tree and Markov model, i.e. average time
# to death given death occurs within a given year.
AVE_TIME_TO_EVENT <- 0.5

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

# Set time horizon
time_hor <- 40.
time_step <- 1.

#### Functions for calculating model probabilities ####

# Function for calculating probabilities from odds ratios and baseline odds
prob_from_or <- function(or, base_odds){
  return (or * base_odds / (1 + or * base_odds))
} %>% Vectorize()

odds_from_prob <- function(prob){
  return (prob / (1 - prob))
} %>% Vectorize()


#### Start by setting up names ####
# There are four subpopulations based on genotype and drug, and six health
# states in the Markov cohort model.

# Define subpopulations

# Keeping this for now in case we revert and want a quick reference point for
# exact formatting
old_subpop_names <- c("clo_lof",
                      "clo_no_lof",
                      "tic",
                      "pra")

subpop_names <- c("ac_lof",
                  "ac_no_lof",
                  "at",
                  "ap")

n_subpops <- length(subpop_names)

# Get list of states for Markov cohort model
markov_states <- c("no_event",
                   "stroke",
                   "post_stroke",
                   "mi",
                   "post_mi",
                   "death")
n_states <- length(markov_states)

# The following pastes together all combinations of subpopulation and health
# state names. Indexing is done using modular arithmetic, so x %% y is remainder
# in x/y and x %/% y is x/y without the remainder.
full_names <- sapply(1:(n_subpops*n_states),
                     FUN = function(i){
                       paste(subpop_names[(i+n_states-1) %/% n_states],
                             markov_states[((i+n_states-1) %% n_states)+1],
                             sep="_")
                     })

#### Now build the decision tree model ####

parameters_STEMI <- read_xlsx("data-inputs/masterfile_100625.xlsx",
                              sheet = "Parameters.STEMI") %>% # Skip row 1 since this doesn't match the format of other rows
  rename_all(make.names) %>% # Converts parameter names in valid R names
  rename_all(.funs = function(name){
    name %>%
      make.names() %>%
      str_replace("parmeter..STEMI.",
                  "parameter.list") %>%
      str_replace("description..STEMI.",
                  "parameter.list") %>%
      str_replace_all("\\.\\.",
                      ".")}) %>% # Converts parameter names in valid R names and corrects typos
  type.convert(as.is = TRUE) %>% # Make sure numbers are numbers, not characters
  filter(!is.na(parameter.list)) %>% # Remove empty lines
  filter(!grepl("_sa", variable.name)) %>% # Drop sensitivity analysis values
  filter(!(parameter.list == "CE threshold")) %>% # Remove this since a single value isn't provided
  mutate(variable.name = variable.name %>%
           str_to_lower() %>% # Standardise parameter names for use with other objects
           str_replace_all(" ", "_") %>%
           str_replace_all("-", "_") %>%
           str_replace_all("__", "_") %>%
           str_replace_all("with_", "") %>%
           str_replace_all("in_", "") %>%
           str_replace_all("nolof", "no_lof") %>%
           str_replace_all("hazard_ratio", "hr") %>%
           str_replace_all("standard_care", "sc") %>%
           str_replace_all("reinfarction", "mi") %>%
           str_replace_all("tica_vs_clop", "at") %>%
           str_replace_all("pras_vs_clop", "ap") %>%
           str_replace_all("ac$", "ac_lof") %>%
           str_replace_all("_vs_ac_standard", "") %>%
           str_replace_all("_vs_at_standard", "") %>%
           str_replace_all("mbleed", "minor_bleed") %>%
           str_replace_all("maj_bleed", "major_bleed") %>%
           str_replace_all("bleeding", "bleed")) %>%
  mutate(value = as.numeric(value),
         sample.size = as.numeric(sample.size)) # Convert values from characters to numbers

# Function for fixing variable names for tree utilities:
rename_utility_variables <- function(name){
  if (str_detect(name, "^utility_")){
    return(paste0(name, "_tree"))
  }else{
    return(name)
  }
} %>% Vectorize()

uncertainty_df <- read_xlsx("data-inputs/masterfile_100625.xlsx",
                            sheet = "random") %>%
  rename_all(.funs = function(name){
    name %>%
      make.names() %>%
      str_replace("alpha..natural.mean",
                  "par1") %>%
      str_replace("beta..natural.SE",
                  "par2")}) %>%
  select(c(variable.name,
           Value,
           par1,
           par2,
           distribution)) %>%
  drop_na() %>%
  mutate(variable.name = variable.name %>%
           str_to_lower() %>% # Standardise parameter names for use with other objects
           str_replace_all(" ", "_") %>%
           str_replace_all("-", "_") %>%
           str_replace_all("__", "_") %>%
           str_replace_all("with_", "") %>%
           str_replace_all("in_", "") %>%
           str_replace_all("nolof", "no_lof") %>%
           str_replace_all("hazard_ratio", "hr") %>%
           str_replace_all("standard_care", "sc") %>%
           str_replace_all("reinfarction", "mi") %>%
           str_replace_all("tica_vs_clop", "at") %>%
           str_replace_all("pras_vs_clop", "ap") %>%
           str_replace_all("ac$", "ac_lof") %>%
           str_replace_all("_vs_ac_standard", "") %>%
           str_replace_all("_vs_at_standard", "") %>%
           str_replace_all("mbleed", "minor_bleed") %>%
           str_replace_all("maj_bleed", "major_bleed") %>%
           str_replace_all("bleeding", "bleed") %>%
           str_replace_all("nofurther", "no")) %>%
  # mutate(variable.name = map(variable.name, rename_utility_variables)) %>%
  mutate(par1 = as.numeric(par1),
         par2 = as.numeric(par2)) # Convert values from characters to numbers

# Reorder so that when we take common lines with parameter dataframe everything
# is in the correct row:
uncertainty_df <- uncertainty_df[order(match(uncertainty_df$variable.name,
                                             parameters_STEMI$variable.name)), ]

draw_df <- do_PSA_draw(uncertainty_df)
rewrite <- rewrite_pars_from_draw(parameters_STEMI, draw_df)
par_df <- rewrite[[1]]
markov_df <- rewrite[[2]]
