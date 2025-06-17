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

# Function for building a new probability dataframe based on updates to baseline
# risks or odds/hazard ratios. A set of baseline risks needs to be provided
rescale_probs <- function(baseline_prob_df,
                          at_ratio_df,
                          ap_ratio_df,
                          ac_no_lof_ratio_df
){
  # Extract probabilities and odds/hazard ratios for each drug/genotype:
  ac_lof_probs <- baseline_prob_df
  
  # Work out Ticagrelor probabilities
  at_probs <- at_ratio_df %>%
    mutate(value = prob_from_or(value, baseline_prob_df$odds))
  
  # Now do Prasagruel:
  ap_probs <- ap_ratio_df %>%
    mutate(value = prob_from_or(value, baseline_prob_df$odds))
  
  # And clopidogrel with no loss of function:
  
  if (UNIQUE_NO_LOF_PROBS){
    ac_no_lof_probs <- ac_no_lof_ratio_df %>%
      mutate(value = value * baseline_prob_df$value)
  }else{
    ac_no_lof_probs <- ac_no_lof_ratio_df %>%
      mutate(value = value * baseline_prob_df$value) %>%
      add_column(at_value = at_probs$value) %>%
      mutate(value = ifelse(grepl("death", variable.name)|grepl("_mi_", variable.name),
                            yes = at_value,
                            no = value)) %>%
      select(c(variable.name, value))
  }
  
  # Assemble into single probability dataframe:
  prob_df <- rbind(ac_lof_probs %>% select(-odds),
                   ac_no_lof_probs,
                   at_probs,
                   ap_probs) %>%
    mutate(variable.name = variable.name %>%
             str_replace("_vs_ac_lof", "") %>%
             str_replace("^or_", "prob_") %>%
             str_replace("^hr_", "prob_") %>%
             str_replace("^rr_", "prob_"))
  return(prob_df)
}


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
  mutate(value = as.numeric(value)) # Convert values from characters to numbers

# Extract probabilities and odds/hazard ratios for each drug/genotype:
baseline_prob_df <- parameters_STEMI %>%
  filter(grepl("prob", variable.name)) %>%
  filter(grepl("_ac_", variable.name)) %>%
  filter(!grepl("no_lof", variable.name)) %>%
  select(c(variable.name, value)) %>%
  mutate(odds = value/(1-value))

# Work out Ticagrelor probabilities
at_ratio_df <- parameters_STEMI %>%
  filter(grepl("^or_", variable.name)) %>%
  filter(grepl("_at", variable.name)) %>%
  select(c(variable.name, value))

# Now do Prasagruel:
ap_ratio_df <- parameters_STEMI %>%
  filter(grepl("^or_", variable.name)) %>%
  filter(grepl("_ap", variable.name)) %>%
  select(c(variable.name, value)) %>%
  add_row(variable.name = "or_dyspnoea_ap", value = 1)

# And clopidogrel with no loss of function:

if (UNIQUE_NO_LOF_PROBS){
  ac_no_lof_ratio_df <- parameters_STEMI %>%
    filter(grepl("^hr_", variable.name)|grepl("^rr_", variable.name)) %>%
    filter(grepl("_ac_no_lof", variable.name)) %>%
    select(c(variable.name, value)) %>%
    add_row(variable.name = "prob_dyspnoea_ac_no_lof", value = 1)
}else{
  ac_no_lof_ratio_df <- parameters_STEMI %>%
    filter(grepl("^hr_", variable.name)|grepl("^rr_", variable.name)) %>%
    filter(grepl("_ac_no_lof", variable.name)) %>%
    add_row(variable.name = "prob_dyspnoea_ac_no_lof", value = 1)
}

# Assemble into single probability dataframe:
prob_df <- rescale_probs(baseline_prob_df,
                         at_ratio_df,
                         ap_ratio_df,
                         ac_no_lof_ratio_df)

# Set up scaling for discounting by time step
discount_by_cycle <- (1 / (1 + parameters_STEMI$value[(
  parameters_STEMI$variable.name=="qaly_discount_rate")|
    (parameters_STEMI$variable.name=="disc_effect_b")])^seq(
    1.0, time_hor-1, by = time_step))

# Note: following appears to be unnecessary
# # Convert hazard rates to probabilities:
# parameters_STEMI$value[which(parameters_STEMI$type=="hazard rate")] <-
#   1 - exp( - parameters_STEMI$value[which(parameters_STEMI$type=="hazard rate")])


# Start ages specified on row 1 of parameter matrices:
start_age_female <- parameters_STEMI$value[parameters_STEMI$variable.name=="start_age_female"]
start_age_male <- parameters_STEMI$value[parameters_STEMI$variable.name=="start_age_male"]

# List parameters for decision tree
pc_uptake <- parameters_STEMI$value[parameters_STEMI$variable.name=="prob_test_order"]
pc_test_cost <- parameters_STEMI$value[parameters_STEMI$variable.name=="poct_cost"]
pc_sens <- 1.
pc_spec <- 1.

# In practice we don't model lab testing. I've left outlines of this in the code
# in case we later decide we do want to make this comparison.
# l_uptake <- parameters_STEMI$value[parameters_STEMI$variable.name=="probability_of_ordering_test"]
# l_test_cost <- parameters_STEMI$value[parameters_STEMI$variable.name=="gendrive_test_price"]
# l_resource_cost <- 100.
# l_sens <- .95
# l_spec <- .95

p_test_followed = parameters_STEMI$value[parameters_STEMI$variable.name=="prob_test_followed"]

p_ac_sc <- parameters_STEMI$value[parameters_STEMI$variable.name=="proportion_ac_standard"] # Proportion prescriped clopidogrel under standard care
p_at_sc <- parameters_STEMI$value[parameters_STEMI$variable.name=="proportion_at_standard"]  # Proportion prescriped ticagrelor under standard care
p_ap_sc <- parameters_STEMI$value[parameters_STEMI$variable.name=="proportion_ap_standard"]  # Proportion prescriped prasugrel under standard care

# Assign probabilities of AT and AP prescription following testing, assuming
# same relative proportions as under standard care.
p_at_A <- p_at_sc / (p_at_sc + p_ap_sc) # Proportion prescriped ticagrelor during subroutine A
p_ap_A <- p_ap_sc / (p_at_sc + p_ap_sc) # Proportion prescriped prasugrel during subroutine A

# Assign prevalences
lof_prev_eu <- parameters_STEMI$value[parameters_STEMI$variable.name=="prevalence_lof_base_case"]
lof_prev_as <- parameters_STEMI$value[parameters_STEMI$variable.name=="prevalence_lof_sensitivity"]
lof_prev <- lof_prev_eu

# Get cost_pci, baseline cost applied in all courses
cost_pci <- parameters_STEMI$value[parameters_STEMI$variable.name=="cost_pci"]

# Get loading doses for drugs. Note that ac_lof and ac_no_lof are identical, but
# for formatting purposes it's more convenient to use the same subpopulations as
# the model.

# Assign waiting period between initial admission and acting on PGX test result:
test_waiting_period <- 0.

ld_costs_sc <- data.frame(drug = c("ac_lof",
                                   "ac_no_lof",
                                   "at",
                                   "ap"),
                          cost = c(parameters_STEMI$value[parameters_STEMI$variable.name=="ld_clop_600"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_clop_600"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_tica_180"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_pras_60"]))

ld_costs_pc <- data.frame(drug = c("ac_lof",
                                   "ac_no_lof",
                                   "at",
                                   "ap"),
                          cost = c(parameters_STEMI$value[parameters_STEMI$variable.name=="ld_clop_600"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_clop_600"] +
                                     test_waiting_period * parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_tica"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_tica_180"],
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_pras_60"]))

# Get daily costs for drugs
daily_costs <- data.frame(drug = c("ac_lof",
                                   "ac_no_lof",
                                   "at",
                                   "ap"),
                          cost = parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_asa"] +
                            c(parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_clop"],
                              parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_clop"],
                              parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_tica"],
                              parameters_STEMI$value[parameters_STEMI$variable.name=="day_cost_pras"]))

# Get expected durations of courses by event
course_dur_by_event_sc <- data.frame(event = c("no_event",
                                               "stroke",
                                               "mi",
                                               "death"),
                                     duration = c(parameters_STEMI$value[(parameters_STEMI$variable.name=="md_no_event")|(parameters_STEMI$variable.name=="md_328")],
                                                  parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_sc")|(parameters_STEMI$variable.name=="md_365_day")],
                                                  parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_sc")|(parameters_STEMI$variable.name=="md_365_day")],
                                                  AVE_TIME_TO_EVENT * parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_sc")|(parameters_STEMI$variable.name=="md_365_day")])) # Assume death occurs half way through course
course_dur_by_event_pc <- data.frame(event = c("no_event",
                                               "stroke",
                                               "mi",
                                               "death"),
                                     duration = c(parameters_STEMI$value[(parameters_STEMI$variable.name=="md_no_event")|(parameters_STEMI$variable.name=="md_328")],
                                                  parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_pc")|(parameters_STEMI$variable.name=="md_365_day")],
                                                  parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_pc")|(parameters_STEMI$variable.name=="md_365_day")],
                                                  AVE_TIME_TO_EVENT * parameters_STEMI$value[(parameters_STEMI$variable.name=="md_mi_stroke_pc")|(parameters_STEMI$variable.name=="md_365_day")])) # Assume death occurs half way through course

# Now get course costs by combination of drug and event (bleeds do not affect
# costs here)
drug_course_costs_sc <-  merge(daily_costs, course_dur_by_event_sc) %>% # For now just do standard care - I don't see why times should be different under PC
  mutate(drug_event = paste(drug, event, sep="_"),
         cost = cost * duration) %>%
  select(-c(drug, event, duration))

drug_course_costs_pc <-  merge(daily_costs, course_dur_by_event_sc) %>% # For now just do standard care - I don't see why times should be different under PC
  mutate(drug_event = paste(drug, event, sep="_"),
         cost = ifelse(drug=="ac_no_lof",
                       cost * (duration - test_waiting_period),
                       cost * duration)) %>%
  select(-c(drug, event, duration))

# Utilities associated with events:
event_utilities <- parameters_STEMI %>%
  filter(grepl("utility", parameter.list)) %>%
  select(c(variable.name, value)) %>%
  mutate(variable.name = gsub("utility_", "", variable.name) %>%
           str_replace_all("_tree", "")) %>%
  filter(!grepl("_ac", variable.name)) %>% # Drop any derived drug-specific values
  filter(!grepl("_at", variable.name)) %>%
  filter(!grepl("_ap", variable.name)) %>%
  filter(!grepl("annual", variable.name)) %>% # Drop any annual-scale values
  add_row(variable.name = "death",
          value = 0) # Apply no_event utility before death, 0 afterwards

# Add dyspnoea utilities, converting from whole-year to 30 day values, and
# convert bleed decrements to utilities
# NOTE: Utility decrements for bleed events are loaded in in units of years so
# do not need to be converted.
duration_dyspnoea <- parameters_STEMI$value[
  which(parameters_STEMI$variable.name=="duration_dyspnoea")] / 365
event_utilities <- event_utilities %>%
  add_row(variable.name = "dyspnoea",
          value = 1 - duration_dyspnoea * event_utilities$value[
            which(event_utilities$variable.name=="dec_dyspnoea")]) %>%
  add_row(variable.name = "major_bleed",
          value = 1 - event_utilities$value[
            which(event_utilities$variable.name=="dec_major_bleed")]) %>%
  add_row(variable.name = "minor_bleed",
          value = 1 - event_utilities$value[
            which(event_utilities$variable.name=="dec_minor_bleed")]) %>%
  filter(!grepl("dec_", variable.name))

# Similar formula to get costs for DT model:
event_costs <- parameters_STEMI %>%
  filter(grepl("cost", parameter.list)) %>%
  filter(grepl("tree|bleed|dysp", parameter.list)) %>%
  filter(!grepl("_ac", variable.name)) %>% # Drop any derived drug-specific values
  filter(!grepl("_at", variable.name)) %>%
  filter(!grepl("_ap", variable.name)) %>%
  select(c(variable.name, value)) %>%
  mutate(variable.name = gsub("cost_", "", variable.name)) %>%
  mutate(variable.name = gsub("_tree", "", variable.name)) %>%
  mutate(variable.name = gsub("_pp", "", variable.name))

### Parameters for Markov cohort model ####

prop_male <- parameters_STEMI$value[
  parameters_STEMI$variable.name=="proportion_male"]

# Read in standardised mortality ratios
smr_df <- read_xlsx("data-inputs/masterfile_100625.xlsx",
                    sheet = "age_sex_dependant_mortality",
                    range = "A3:E8") %>%
  mutate(SMRs = SMRs %>%
           str_to_lower() %>%
           str_replace("smr_", "") %>%
           str_replace("no further event", "no_event") %>%
           str_replace("reinfarction", "mi") %>%
           str_replace("post-", "post_"))

# Make copy containing only values - this distinction should be useful for PSA
# purposes
smr_vals <- smr_df %>%
  select(c(SMRs,
           value)) %>%
  spread(SMRs, value)

# Function for calculating mortality by age for a cohort given an SMR
apply_smr <- function(mort_female,
                      mort_male,
                      smr,
                      prop_male){
  adj_female <- pmin(1, smr * mort_female)
  adj_male <- pmin(1, smr * mort_male)
  return((1 - prop_male) * adj_female +
           prop_male * adj_male)
} 

# Read in life table for healthy individuals
mortality_prob_by_age <- read_xlsx("data-inputs/masterfile_100625.xlsx",
                        sheet = "age_sex_dependant_mortality",
                        range = "A12:D52") %>%
  mutate(no_event = apply_smr(mortality_female,
                              mortality_male,
                              smr_vals$no_event,
                              prop_male),
         mi = apply_smr(mortality_female,
                        mortality_male,
                        smr_vals$mi,
                        prop_male),
         post_mi = apply_smr(mortality_female,
                             mortality_male,
                             smr_vals$post_mi,
                             prop_male),
         stroke = apply_smr(mortality_female,
                            mortality_male,
                            smr_vals$stroke,
                            prop_male),
         post_stroke = apply_smr(mortality_female,
                                 mortality_male,
                                 smr_vals$post_stroke,
                                 prop_male)) %>%
  select(-c(mortality_female,
            mortality_male,
            age_female,
            age_male))

# Write in Markov parameters directly (these come from the literature)
base_markov_pars <- data.frame(parameter.list = c("mi",
                                             "stroke"),
                          value = c(.043,
                                    .0112))

# Load utilities and add zeros for death
markov_utils <- read_xlsx("data-inputs/masterfile_240325.xlsx",
                          sheet = "time_event_utility") %>%
  mutate(death = 0) %>%
  select(all_of(markov_states)) # Last step just makes sure ordering matches model
markov_utils <- markov_utils[1:39, ] # Don't end up using last entry

# Extract costs from main parameter table
markov_costs <- parameters_STEMI %>%
  filter(grepl("cost", parameter.list)) %>%
  filter(grepl("markov", parameter.list)) %>%
  select(c(variable.name, value)) %>%
  mutate(variable.name = gsub("cost_", "", variable.name)) %>%
  mutate(variable.name = gsub("_markov", "", variable.name)) %>%
  add_row(variable.name = "death",
          value = 0) %>% # Assign 0 cost to deaths in this section of model
  arrange(factor(variable.name, levels = markov_states))