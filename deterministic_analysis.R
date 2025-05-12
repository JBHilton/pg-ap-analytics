# In this script we carry out a deterministic simulation of the decision process
# and long-term cohort dynamics
#
# Subpopulations to model are:
# 1: Clopidogrel, LOF
# 2: Clopidogrel, no LOF
# 3: Ticagrelor
# 4: Prasugrel
#
# Here "LOF" and "no LOF" refer to true genotype, NOT to test results
#
# Basic approach: Work through decision tree to assign subpopulation sizes, then
# use these to initialise a Markov model with a structured population cohort.
#
# Note: The decision "tree" has a very non-treelike structure, with all routes
# passing through some combination of two subroutines, A and B. This limits the
# effectiveness of existing R packages for decision tree modelling, and so we
# perform the decision analysis using direct calculation.
#
# All checks for individual genotype/test results only test whether or not the
# recorded genotype is "no_lof". Any other string will be interpreted as
# indicating the wild type loss of function genotype.

# Flag to determine whether to convert rates in the data to probabilities
CONVERSIONS_REQUIRED <- FALSE

# Next flag is here to suppress some code which isn't currently useful but might
# be at some point
SUBPOP_PLOTTING <- FALSE

# Step correction for death in decision tree and Markov model, i.e. average time
# to death given death occurs within a given year.
AVE_TIME_TO_EVENT <- 0.5

# Set to true to save the results of the arm comparison as a .csv file. The file
# path can be set on the following line:
SAVE_ARM_COMPARISON <- TRUE
SAVE_FILEPATH <- ""

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

parameters_STEMI <- read_xlsx("data-inputs/masterfile_270425.xlsx",
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
           str_replace_all("with_", "") %>%
           str_replace_all("in_", "") %>%
           str_replace_all("hazard_ratio", "hr") %>%
           str_replace_all("standard_care", "sc") %>%
           str_replace_all("reinfarction", "mi") %>%
           str_replace_all("tica_vs_clop", "at") %>%
           str_replace_all("pras_vs_clop", "ap") %>%
           str_replace_all("ac$", "ac_lof") %>%
           str_replace_all("mbleed", "minor_bleed") %>%
           str_replace_all("maj_bleed", "major_bleed") %>%
           str_replace_all("bleeding", "bleed")) %>%
  mutate(value = as.numeric(value)) # Convert values from characters to numbers

# Extract probabilities and odds/hazard ratios for each drug/genotype:
ac_lof_probs <- parameters_STEMI %>%
  filter(grepl("prob", variable.name)) %>%
  filter(grepl("_ac_", variable.name)) %>%
  filter(!grepl("no_lof", variable.name)) %>%
  select(c(variable.name, value)) %>%
  mutate(odds = value/(1-value))

prob_from_or <- function(or, base_odds){
  return (or * base_odds / (1 + or * base_odds))
} %>% Vectorize()

# Work out Ticagrelor probabilities
at_probs <- parameters_STEMI %>%
  filter(grepl("^or_", variable.name)) %>%
  filter(grepl("_at", variable.name)) %>%
  select(c(variable.name, value)) %>%
  mutate(value = prob_from_or(value, ac_lof_probs$odds))

# We can inspect the following to make sure the odds ratio conversion is
# consistent with the Excel workbook:
at_probs_load <- parameters_STEMI %>%
  filter(grepl("prob", variable.name)) %>%
  filter(grepl("_at", variable.name)) %>%
  select(c(variable.name, value))

# Now do Prasagruel:
ap_probs <- parameters_STEMI %>%
  filter(grepl("^or_", variable.name)) %>%
  filter(grepl("_ap", variable.name)) %>%
  select(c(variable.name, value)) %>%
  add_row(variable.name = "or_dyspnoea_ap", value = 1) %>%
  mutate(value = prob_from_or(value, ac_lof_probs$odds))

ap_probs_load <- parameters_STEMI %>%
  filter(grepl("prob", variable.name)) %>%
  filter(grepl("_ap", variable.name)) %>%
  select(c(variable.name, value))

# And clopidogrel with no loss of function:
ac_no_lof_probs <- parameters_STEMI %>%
  filter(grepl("^hr_", variable.name)) %>%
  filter(grepl("_ac_no_lof", variable.name)) %>%
  select(c(variable.name, value)) %>%
  add_row(variable.name = "prob_dyspnoea_ac_no_lof", value = 1) %>%
  mutate(value = value * ac_lof_probs$value)

ac_no_lof_probs_load <- parameters_STEMI %>%
  filter(grepl("prob", variable.name)) %>%
  filter(grepl("_ac_no_lof", variable.name)) %>%
  select(c(variable.name, value))

# Assemble into single probability dataframe:
prob_df <- rbind(ac_lof_probs %>% select(-odds),
                 ac_no_lof_probs,
                 at_probs,
                 ap_probs) %>%
  mutate(variable.name = variable.name %>%
           str_replace("_vs_ac_lof", "") %>%
           str_replace("^or_", "prob_") %>%
           str_replace("^hr_", "prob_"))

# Set up scaling for discounting by time step
discount_by_cycle <- as.matrix((1 / (1 + parameters_STEMI$value[
  parameters_STEMI$variable.name=="qaly_discount_rate"])^seq(
    1.0, (time_hor - 1), by = time_step)))

# Note: following appears to be unnecessary
# # Convert hazard rates to probabilities:
# parameters_STEMI$value[which(parameters_STEMI$type=="hazard rate")] <-
#   1 - exp( - parameters_STEMI$value[which(parameters_STEMI$type=="hazard rate")])


# Start ages specified on row 1 of parameter matrices:
start_age_female <- parameters_STEMI$value[parameters_STEMI$variable.name=="start_age_female"]
start_age_male <- parameters_STEMI$value[parameters_STEMI$variable.name=="start_age_male"]

# List parameters for decision tree
pc_uptake <- parameters_STEMI$value[parameters_STEMI$variable.name=="prob_test_order"]
pc_test_cost <- parameters_STEMI$value[parameters_STEMI$variable.name=="gendrive_test_price"]
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
                                   parameters_STEMI$value[parameters_STEMI$variable.name=="ld_tica_180"] +
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
                                     duration = c(parameters_STEMI$value[parameters_STEMI$variable.name=="md_no_event"],
                                                  parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_sc"],
                                                  parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_sc"],
                                                  AVE_TIME_TO_EVENT * parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_sc"])) # Assume death occurs half way through course
course_dur_by_event_pc <- data.frame(event = c("no_event",
                                               "stroke",
                                               "mi",
                                               "death"),
                                     duration = c(parameters_STEMI$value[parameters_STEMI$variable.name=="md_no_event"],
                                                  parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_pc"],
                                                  parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_pc"],
                                                  AVE_TIME_TO_EVENT * parameters_STEMI$value[parameters_STEMI$variable.name=="md_mi_stroke_pc"])) # Assume death occurs half way through course

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
duration_dyspnoea <- parameters_STEMI$value[
  which(parameters_STEMI$variable.name=="duration_dyspnoea")]
event_utilities <- event_utilities %>%
  add_row(variable.name = "dyspnoea",
          value = 1 - event_utilities$value[
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

# Put some filler values in for utilities of bleeds
minor_bleed_utility <- event_utilities$value[event_utilities$variable.name=="minor_bleed"]
major_bleed_utility <- event_utilities$value[event_utilities$variable.name=="major_bleed"]

# Basic patient status dataframe

patient_status <- data.frame(subpop = subpop_names,
                             prob = rep(0, length(subpop_names)),
                             exp_cost = rep(0, length(subpop_names)),
                             exp_utility = rep(1, length(subpop_names)))

# Define functions to calculate end state probablities after passing through
# some combination of standard care, subroutine A, and subroutine B.

# Standard care
implement_sc <- function(true_genotype){
  
  # If prescribed clopidogrel then subpopulation depends on genotype, otherwise
  # just used presciption proportions.
  p_ac_lof <- p_ac_sc
  p_ac_no_lof <- 0
  p_at <- p_at_sc
  p_ap <- p_ap_sc
  
  sc_results <- data.frame(subpops = subpop_names,
                                     prob = c(p_ac_lof,
                                              p_ac_no_lof,
                                              p_at,
                                              p_ap),
                           cost = c(ld_costs_sc$cost[ld_costs_sc$drug=="ac"],
                                    ld_costs_sc$cost[ld_costs_sc$drug=="ac"],
                                    ld_costs_sc$cost[ld_costs_sc$drug=="at"],
                                    ld_costs_sc$cost[ld_costs_sc$drug=="ap"]),
                           utility = c(1,
                                       1,
                                       1,
                                       1))
  return(sc_results)
}

implement_A <- function(true_genotype,
                        test_result
                        ){
  
  # If test is followed, test says no lof, and true genotype is lof, then
  # patient goes to clopidogrel_lof
  p_ac_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 0, 1)
  
  # If test is followed, test says no lof, and true genotype is no_lof, then
  # patient goes to clopidogrel_no_lof
  p_ac_no_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 1, 0)
  
  # If test is followed and test says lof, and true genotype is lof, then
  # patient is prescribed one of the other drugs
  p_at <- p_test_followed *
    ifelse(test_result=="no_lof", 0, 1) * 
    p_at_A
  p_ap <- p_test_followed *
    ifelse(test_result=="no_lof", 0, 1) * 
    p_ap_A
  
  subroutine_A_results <- data.frame(subpop = c("ac_lof",
                                                   "ac_no_lof",
                                                   "at",
                                                   "ap",
                                                   "sc"),
             prob = c(p_ac_lof,
                      p_ac_no_lof,
                      p_at,
                      p_ap,
                      1 - p_test_followed), # Go to standard care if test not followed
             cost = c(ld_costs_pc$cost[ld_costs_pc$drug=="ac_lof"],
                      ld_costs_pc$cost[ld_costs_pc$drug=="ac_no_lof"],
                      ld_costs_pc$cost[ld_costs_pc$drug=="at"],
                      ld_costs_pc$cost[ld_costs_pc$drug=="ap"],
                      0),
             utility = c(1,
                         1,
                         1,
                         1,
                         1)) # costs assigned in standard care
  return(subroutine_A_results)
}

implement_B <- function(subpop_id,
                        test){
  # subpop_pars <- parameters_STEMI[grepl(paste(subpop_id, "$", sep=""), # Gather subpopulation-specific parameters
  #                                       parameters_STEMI$variable.name),] %>%
  #   filter(!grepl("hr_", variable.name)) %>%
  #   filter(!grepl("hazard_ratio_", variable.name)) %>%
  #   filter(!(variable.name == "or_dyspnoea_at")) %>% # Drop stray odds ratio that appears for Ticagrelor
  #   filter(!grepl("utility_decrement_", variable.name)) %>% # Drop derived drug-specific utility decrements
  #   filter(!grepl("u_dec_", variable.name)) %>% # Drop derived drug-specific utility decrements
  #   filter(!grepl("cost_", variable.name)) %>%
  #   mutate(variable.name = variable.name %>% str_replace_all(paste("_", subpop_id, sep=""), ""))
  
  subpop_pars <- prob_df[grepl(paste(subpop_id, "$", sep=""), # Gather subpopulation-specific parameters
                                        prob_df$variable.name),] %>%
      mutate(variable.name = variable.name %>% str_replace_all(paste("_", subpop_id, sep=""), ""))
  
  # Get probabilities of each event from dataframe:
  minor_bleed_prob <- subpop_pars$value[
    grep("minor_bleed", subpop_pars$variable.name)]
  major_bleed_prob <- subpop_pars$value[
    grep("major_bleed", subpop_pars$variable.name)]
  # For ac_no_lof, just use ac_lof value of dyspnoea probability since LOF gene does not affect this
  if (subpop_id == "ac_no_lof"){
    dysp_prob <- prob_df$value[which(
      prob_df$variable.name == "prob_dyspnoea_ac_lof"
    )]
  }else{
    dysp_prob <- subpop_pars$value[
      grep("dyspnoea", subpop_pars$variable.name)]
  }
  mi_prob <- subpop_pars$value[
    grep("mi$", subpop_pars$variable.name)]
  stroke_prob <- subpop_pars$value[
    grep("stroke", subpop_pars$variable.name)]
  death_prob <- subpop_pars$value[
    grep("death", subpop_pars$variable.name)]
  
  if (test=="pc"){
    drug_course_costs = drug_course_costs_pc
  }else{
    drug_course_costs = drug_course_costs_sc
  }
  
  mi_cost <- event_costs$value[event_costs$variable.name=="mi"] +
    drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                               "mi",
                                                               sep="_")]
  stroke_cost <- event_costs$value[event_costs$variable.name=="stroke"] +
    drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                               "stroke",
                                                               sep="_")]
  death_cost <- event_costs$value[event_costs$variable.name=="death"] +
    drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                               "death",
                                                               sep="_")]
  no_event_cost <- event_costs$value[event_costs$variable.name=="no_event"] +
    drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                               "no_event",
                                                               sep="_")]
  
  bleed_results <- data.frame(event = c("minor_bleed",
                                     "major_bleed",
                                     "no_event"),
                           prob = c(minor_bleed_prob,
                                    major_bleed_prob,
                                    1 - minor_bleed_prob - major_bleed_prob),
                           cost = c(event_costs$value[
                             event_costs$variable.name=="minor_bleed"],
                             event_costs$value[
                               event_costs$variable.name=="major_bleed"],
                             0),
                           utility = c(event_utilities$value[
                             event_utilities$variable.name=="minor_bleed"],
                             event_utilities$value[
                               event_utilities$variable.name=="major_bleed"],
                                       1))
  exp_cost_bleed <- sum(bleed_results$prob * bleed_results$cost)
  exp_util_bleed <- sum(bleed_results$prob * bleed_results$utility)
  
  dysp_results <- data.frame(event = c("dysp",
                                       "no_event"),
                              prob = c(dysp_prob,
                                       1 - dysp_prob),
                             cost = c(event_costs$value[
                               event_costs$variable.name=="dyspnoea"],
                               0),
                              utility = c(event_utilities$value[
                                event_utilities$variable.name=="dyspnoea"],
                                1))
  exp_cost_dysp <- sum(dysp_results$prob * dysp_results$cost)
  exp_util_dysp <- sum(dysp_results$prob * dysp_results$utility)
  
  # Note on results: we currently assign zero cost to each of these outcomes
  subroutine_B_results <- data.frame(event = c("reinfarction",
                                     "stroke",
                                     "death",
                                     "no_event"),
                              prob = c(mi_prob,
                                       stroke_prob,
                                       death_prob,
                                       1 - mi_prob - stroke_prob - death_prob),
                              cost = cost_pci + exp_cost_bleed + exp_cost_dysp +
                                c(mi_cost,
                                       stroke_cost,
                                       death_cost,
                                       no_event_cost),
                           utility = -(1 - exp_util_bleed) + # Subtract utility decrement due to bleed
                             -(1 - exp_util_dysp) + # Subtract utility decrement due to dyspnoea
                             AVE_TIME_TO_EVENT * # Apply step correction; events happen part way through year
                             event_utilities$value[event_utilities$variable.name == "no_event"] +
                             AVE_TIME_TO_EVENT *
                             c(event_utilities$value[event_utilities$variable.name == "mi"],
                                       event_utilities$value[event_utilities$variable.name == "stroke"],
                                       event_utilities$value[event_utilities$variable.name == "death"],
                                       event_utilities$value[event_utilities$variable.name == "no_event"]))
  return(subroutine_B_results)
}

# Patient steps through subroutine B regardless of testing method, so we
# calculate it here
B_result_list_sc <- lapply(subpop_names,
                        implement_B,
                        "sc")

B_result_list_pc <- lapply(subpop_names,
                           implement_B,
                           "pc")

# Make a dataframe to match up names of events in B with the Markov model states
# patients end up in afterwards. Note leading underscore "_" used in string
# concatentation later.
post_event_match <- data.frame(event = c("reinfarction",
                                         "stroke",
                                         "death",
                                         "no_event"),
                               state = c("_post_mi",
                                         "_post_stroke",
                                         "_death",
                                         "_no_event"))


# Function to do entire decision tree
# Test should be one of "pc", "l", or "sc" (standard care).
run_forward <- function(test = "sc"){
  patient_status <- data.frame(subpop = subpop_names,
                               prob = rep(1, length(subpop_names)),
                               exp_cost = rep(0, length(subpop_names)),
                               exp_utility = rep(1, length(subpop_names)))
  
  # Always need to calculate standard care results as all routes have a chance
  # of going to it
  sc_results_lof <- implement_sc("lof")
  sc_results_no_lof <- implement_sc("no_lof")
  
  if (test == "sc"){
    patient_status$prob <- patient_status$prob *
      (lof_prev * sc_results_lof$prob +
      (1 - lof_prev) * sc_results_no_lof$prob)
    patient_status$exp_cost <- patient_status$exp_cost +
      lof_prev * sc_results_lof$prob * sc_results_lof$cost +
      (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost
  }
  
  # Only need results of subroutine A if we actually do testing. Also decide
  # here which version of B_result_list to use
  if((test == "pc") | (test == "l")){
    
    # Assign performance metrics for relevant test
    if (test == "pc"){
      sens <- pc_sens
      spec <- pc_spec
      test_uptake <- pc_uptake
    }
    if (test == "l"){
      sens <- l_sens
      spec <- l_spec
      test_uptake <- l_uptake
    }
    A_results_lof_lof <- implement_A(true_genotype = "lof",
                                 test_result = "lof")
    A_results_lof_no_lof <- implement_A(true_genotype = "lof",
                                     test_result = "no_lof")
    A_results_no_lof_lof <- implement_A(true_genotype = "no_lof",
                                     test_result = "lof")
    A_results_no_lof_no_lof <- implement_A(true_genotype = "no_lof",
                                     test_result = "no_lof")
    
    # Take weighted average of these based on expected frequency of genotype-
    # test result combinations. This ignores any potential variance which would
    # need to be incorporated into a fully stochastic model.
    A_results_ave <- data.frame(subpop = A_results_lof_lof$subpop,
                                 prob = rep(1, length(subpop_names) + 1),
                                 cost = rep(0, length(subpop_names) + 1),
                                 utility = rep(1, length(subpop_names) + 1))
    A_results_ave[, -1] <- lof_prev * spec * A_results_lof_lof[, -1] + # True positive
      lof_prev * (1 - spec) * A_results_lof_no_lof[, -1] + # False negative
      (1 - lof_prev) * (1 - sens) * A_results_no_lof_lof[, -1] + # False positive
      (1 - lof_prev) * sens * A_results_no_lof_no_lof[, -1] # True negative
    
    # Probability of going into standard care is probability of not testing plus
    # probability of testing and then ignoring. Note that this is a bit messy
    # because test uptake is currently implemented outside of the A subroutine
    # but whether the clinician follows the test is implemented within it.
    prob_sc <- (1 - test_uptake) + test_uptake * A_results_ave$prob[5]
    
    # Now add results to patient status
    patient_status$prob <- test_uptake *
      patient_status$prob *
      A_results_ave$prob[1:4]  +
      prob_sc * (lof_prev * sc_results_lof$prob + # This adds possibility of reverting to standard care
         (1 - lof_prev) * sc_results_no_lof$prob)
    patient_status$exp_cost <- patient_status$exp_cost +
      test_uptake * parameters_STEMI$value[parameters_STEMI$variable.name=="poct_cost"] + # Account for cost of testing
      test_uptake * A_results_ave$cost[1:4] +
      prob_sc *
    (lof_prev * sc_results_lof$prob * sc_results_lof$cost +
      (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost)
    print(test)
    print(patient_status$exp_cost)
    B_result_list <- B_result_list_pc
  }else{
    # If we aren't testing we use the standard care version of B results
    B_result_list <- B_result_list_sc
  }
  
  # Now extend to combined subpopulation-health state space by implementing
  # the outcomes of subroutine B
  extended_patient_status <- data.frame(subpop = full_names,
                                        prob = rep(0, length(full_names)),
                                        exp_cost = rep(0, length(full_names)),
                                        exp_utility = rep(1, length(full_names)))
  for (i in (1:nrow(patient_status))){
    name <- patient_status$subpop[i]
    
    # Get total probability of being in subpopulation, to divide between subpop-
    # event combinations, and utility associated with subpopulation.
    p_subgroup <- patient_status$prob[which(patient_status$subpop==name)]
    c_subgroup <- patient_status$exp_cost[which(patient_status$subpop==name)]
    u_subgroup <- patient_status$exp_utility[which(patient_status$subpop==name)]
    
    # Safer to bring in costs for all Markov states, not just the ones we have
    # nonzero probability of starting the Markov chain in. For probability and
    # utility, I think we want to leave them at zero before implementing B.
    extended_patient_status$exp_cost[
      grep(paste(name, "_",sep=""),
           extended_patient_status$subpop)] <- c_subgroup
    
    # Assign probabilities and update utilities based on subroutine B.
    # Note: this assigns utility 1 to anyone in 
    for (j in (1:nrow(post_event_match))){
      # print(paste(name, post_event_match$state[j], sep=""))
      # print(which(grepl(paste(name, post_event_match$state[j], sep=""),
      #                  extended_patient_status$subpop)))
      extended_patient_status$prob[
        grepl(paste(name, post_event_match$state[j], sep=""),
                                        extended_patient_status$subpop)] <-
        B_result_list[[i]]$prob[
          which(B_result_list[[i]]$event==post_event_match$event[j])] *
        p_subgroup
      # print(paste("ping", extended_patient_status$exp_cost[6]))
      extended_patient_status$exp_cost[
        grepl(paste(name, post_event_match$state[j], sep=""),
             extended_patient_status$subpop)] <-
        extended_patient_status$exp_cost[
          grepl(paste(name, post_event_match$state[j], sep=""),
               extended_patient_status$subpop)] +
        B_result_list[[i]]$cost[
          which(B_result_list[[i]]$event==post_event_match$event[j])]
      # print(extended_patient_status$exp_cost[6])
      extended_patient_status$exp_utility[
        grep(paste(name, post_event_match$state[j], sep=""),
             extended_patient_status$subpop)] <-
        B_result_list[[i]]$utility[
          which(B_result_list[[i]]$event==post_event_match$event[j])] *
        u_subgroup
    }
  }
  return(extended_patient_status)
}

### Set up Markov cohort model ####

# Assemble generic transition matrix without accounting for sex, genotype, or drug
base_transition_matrix <- matrix(0,
                                 nrow = n_states,
                            ncol = n_states,
                            dimnames = list(markov_states,
                                            markov_states))

build_markov_submodel <- function(age,
                               mort_prob_by_age = mort_prob_by_age,
                               stroke_prob = .043,
                               mi_prob = .01){
  
  transition_matrix <- matrix(0,
                              nrow = n_states,
                              ncol = n_states,
                              dimnames = list(markov_states,
                                              markov_states))
  
  # Get probabilities of each event from dataframe:
  
  # Outflow from no_event
  transition_matrix["no_event", "mi"] <- subpop_pars$value[
    grep("reinfarction", subpop_pars$parameter.list)]
  transition_matrix["no_event", "stroke"] <- subpop_pars$value[
    grep("^stroke", subpop_pars$parameter.list)]
  transition_matrix["no_event", "death"] <- mortality_prob_by_age[
    which(age==age)]
  
  # Outflow from mi
  transition_matrix["mi", "post_mi"] <- 1. - mortality_prob_by_age[
    which(age==age)]
  transition_matrix["mi", "death"] <- mortality_prob_by_age[
    which(age==age)]
  
  # Outflow from post_mi
  transition_matrix["post_mi", "death"] <- mortality_prob_by_age[
    which(age==age)]
  
  # Outflow from stroke
  transition_matrix["stroke", "post_stroke"] <- 1. - mortality_prob_by_age[
    which(age==age)]
  transition_matrix["stroke", "death"] <- mortality_prob_by_age[
    which(age==age)]
  
  # Outflow from post_stroke
  transition_matrix["post_stroke", "death"] <- mortality_prob_by_age[
    which(age==age)]
  
  diag(transition_matrix) <- 1 - rowSums(transition_matrix)
  
  return(transition_matrix)
}

# Read in age/health state-dependent mortality rates:
mortality_prob_by_age <- read_xlsx("data-inputs/masterfile_240325.xlsx",
                                   sheet = "time_event_mortality")

# Write in Markov parameters directly (these come from the literature)
markov_pars <- data.frame(parameter.list = c("mi",
                                             "stroke"),
                          value = c(.043,
                                    .0112))

# Load utilities and add zeros for death
markov_utils <- read_xlsx("data-inputs/masterfile_240325.xlsx",
                          sheet = "time_event_utility") %>%
  mutate(death = 0) %>%
  slice(-1) %>% # Remove first row since this is covered by decision tree
  select(all_of(markov_states)) # Last step just makes sure ordering matches model

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
  

# Function for building single-population Markov model
# This is almost identical to build_markov_submodel, but without string
# matching etc for choosing the subpopulation.
build_markov_model <- function(tstep){
  
  transition_matrix <- matrix(0,
                              nrow = n_states,
                              ncol = n_states,
                              dimnames = list(markov_states,
                                              markov_states))
  
  # Get probabilities of each event from dataframe:
  
  # Outflow from no_event
  transition_matrix["no_event", "mi"] <- markov_pars$value[
    grep("mi", markov_pars$parameter.list)]
  transition_matrix["no_event", "stroke"] <- markov_pars$value[
    grep("^stroke", markov_pars$parameter.list)]
  transition_matrix["no_event", "death"] <- mortality_prob_by_age$no_event[tstep]
  
  # Outflow from mi
  transition_matrix["mi", "post_mi"] <- 1. - mortality_prob_by_age$mi[tstep]
  transition_matrix["mi", "death"] <- mortality_prob_by_age$mi[tstep]
  
  # Outflow from post_mi
  transition_matrix["post_mi", "death"] <- mortality_prob_by_age$post_mi[tstep]
  
  # Outflow from stroke
  transition_matrix["stroke", "post_stroke"] <- 1. - mortality_prob_by_age$stroke[tstep]
  transition_matrix["stroke", "death"] <- mortality_prob_by_age$stroke[tstep]
  
  # Outflow from post_stroke
  transition_matrix["post_stroke", "death"] <- mortality_prob_by_age$post_stroke[tstep]
  
  diag(transition_matrix) <- 1 - rowSums(transition_matrix)
  
  
  return(transition_matrix)
}


#### Now run model ####
# The prob column from the dataframe outputted by the run_forward function is
# the initial condition for the Markov model.

# Choose time horizon
n_tsteps <- 39

# Choose test to analyse
test <- "pc"

# Run decision tree analysis
pc_dt_results <- run_forward(test)

# Marginalise over drug subpopulations to get an initial condition for Markov
# model
grouped_results <- pc_dt_results %>%
  mutate(event = subpop %>%
           str_remove_all(paste(paste(subpop_names, collapse = "_|"),
                                "_",
                                sep = "")))
P0 <- sapply(markov_states,
                   FUN=function(event){
                     sum(grouped_results$prob[grouped_results$event == event])}
                   )

markov_trace <- sapply(1:n_tsteps,
                        FUN = function(t){
                          P0 %*% (build_markov_model(t+1) %^% t)
                          }
                        ) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 1:n_tsteps)

# We can get, for example, a plot of expected strokes per capita per year, split
# into subgroups. Note that each curve is scaled by proportion of population in
# subgroup - to get probability by subgroup you'd need to divide through by this
# proportion.
stroke_plot_df <- markov_trace %>%
  pivot_longer(grep("stroke",colnames(markov_trace)),
               names_to = "drug",
               values_to = "risk_of_stroke") %>%
  filter(!grepl("post", drug)) %>%
  mutate(drug = gsub("_stroke", "", drug)) %>%
  select(time_step, drug, risk_of_stroke)

stroke_plot_df %>%
  ggplot(aes(x = time_step, y = risk_of_stroke, colour = drug)) +
  geom_line()

# Calculate expected utility over time:
utility_df <- data.frame(time_step = 1:n_tsteps,
                         utility = rowSums(as.matrix(markov_trace %>%
                                                select(-time_step)) *
                           as.matrix(markov_utils)))

#### ICER calculation ####
# Now that we've seen how the modelling workflow works, let's calculate an ICER
# for point of care vs standard care:


# Run decision tree analysis
dt_results_pc <- run_forward("pc") %>%
  mutate(event = subpop %>%
           str_remove_all(paste(paste(subpop_names, collapse = "_|"),
                                "_",
                                sep = "")))

# Expected costs and utilities at DT stage:
dt_pc_cost <- sum(dt_results_pc$prob * dt_results_pc$exp_cost)
dt_pc_util <- sum(dt_results_pc$prob * dt_results_pc$exp_utility)

P0_pc <- sapply(markov_states,
             FUN=function(event){
               sum(dt_results_pc$prob[dt_results_pc$event == event])}
)

MT_pc <- sapply(1:n_tsteps,
                       FUN = function(t){
                         P0_pc %*% (build_markov_model(t) %^% t)
                       }
) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 1:n_tsteps)


# Calculate expected utility over time:
utility_pc <- data.frame(time_step = 1:n_tsteps,
                         utility = rowSums(as.matrix(MT_pc %>%
                                                           select(-time_step)) *
                                                 as.matrix(markov_utils))) %>%
  mutate(discounted_utility =
           utility * discount_by_cycle)
# Expected costs:
MC_costs_pc <- data.frame(time_step = 1:n_tsteps,
                         cost = rowSums(as.matrix(MT_pc %>%
                                                       select(-time_step)) %*%
                                             as.matrix(markov_costs$value))) %>%
  mutate(discounted_cost =
           cost * discount_by_cycle)


# Mean life years:
life_years_pc <- (1-P0_pc["death"]) + sum(1 - MT_pc$death)

# Now do sc

# Run decision tree analysis
dt_results_sc <- run_forward("sc") %>%
  mutate(event = subpop %>%
           str_remove_all(paste(paste(subpop_names, collapse = "_|"),
                                "_",
                                sep = "")))

# Expected costs at DT stage:
dt_sc_cost <- sum(dt_results_sc$prob * dt_results_sc$exp_cost)
dt_sc_util <- sum(dt_results_sc$prob * dt_results_sc$exp_utility)

P0_sc <- sapply(markov_states,
                FUN=function(event){
                  sum(dt_results_sc$prob[dt_results_sc$event == event])}
)

MT_sc <- sapply(1:n_tsteps,
                FUN = function(t){
                  P0_sc %*% (build_markov_model(t) %^% t)
                }
) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 1:n_tsteps)


# Calculate expected utility over time:
utility_sc <- data.frame(time_step = 1:n_tsteps,
                         utility = rowSums(as.matrix(MT_sc %>%
                                                           select(-time_step)) *
                                                 as.matrix(markov_utils))) %>%
  mutate(discounted_utility =
           utility * discount_by_cycle)
# Expected costs:
MC_costs_sc <- data.frame(time_step = 1:n_tsteps,
                          cost = rowSums(as.matrix(MT_sc %>%
                                                     select(-time_step)) %*%
                                           as.matrix(markov_costs$value))) %>%
  mutate(discounted_cost =
           cost * discount_by_cycle)

# Mean life years:
life_years_sc <- (1-P0_sc["death"]) + sum(1 - MT_sc$death)

# Now calculate ICER
ICER <- ((dt_pc_cost + sum(MC_costs_pc$cost)) - (dt_sc_cost + sum(MC_costs_sc$cost))) /
  ((dt_pc_util + sum(utility_pc$utility)) - (dt_sc_util + sum(utility_sc$utility)))


# Now calculate ICER
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
{
  print(paste("Mean life years under SC = ",
              (life_years_sc)))
  print(paste("Mean life years under PC = ",
              (life_years_pc)))
  print("")
  print("With discounting:")
  print(paste("Mean cost under SC = ",
              (dt_sc_cost + sum(MC_costs_sc$discounted_cost))))
  print(paste("Mean cost under PC = ",
        (dt_pc_cost + sum(MC_costs_pc$discounted_cost))))
  print(paste("Mean utility under PC = ",
              (dt_pc_util + sum(utility_pc$discounted_utility))))
  print(paste("Mean utility under SC = ",
              (dt_sc_util + sum(utility_sc$discounted_utility))))
  print(paste("Estimated incremental cost is",
        (dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))))
  print(paste("Estimated incremental utility is",
              (dt_pc_util + sum(utility_pc$discounted_utility)) -
                (dt_sc_util + sum(utility_sc$discounted_utility))))
  print(paste("Estimated ICER is",
              ICER_disc))
  print("")
  print("Without discounting:")
  print(paste("Mean cost under PC = ",
              (dt_pc_cost + sum(MC_costs_pc$cost))))
  print(paste("Mean cost under SC = ",
              (dt_sc_cost + sum(MC_costs_sc$cost))))
  print(paste("Mean utility under PC = ",
              (dt_pc_util + sum(utility_pc$utility))))
  print(paste("Mean utility under SC = ",
              (dt_sc_util + sum(utility_sc$utility))))
  print(paste("Estimated incremental cost is",
              (dt_pc_cost + sum(MC_costs_pc$cost)) - (dt_sc_cost + sum(MC_costs_sc$cost))))
  print(paste("Estimated incremental utility is",
              (dt_pc_util + sum(utility_pc$utility)) -
                (dt_sc_util + sum(utility_sc$utility))))
  print(paste("Estimated ICER is",
              ICER))
  }

# Create a dataframe storing outputs by case:
arm_comparison <- data.frame(arm = c("SC", "PC", "Increment"),
                             utility_udc = c(dt_sc_util + sum(utility_sc$utility),
                                             dt_pc_util + sum(utility_pc$utility),
                                             (dt_pc_util + sum(utility_pc$utility)) -
                                               (dt_sc_util + sum(utility_sc$utility))),
                             cost_udc = c(dt_sc_cost + sum(MC_costs_sc$cost),
                                             dt_pc_cost + sum(MC_costs_pc$cost),
                                             (dt_pc_cost + sum(MC_costs_pc$cost)) -
                                               (dt_sc_cost + sum(MC_costs_sc$cost))),
                             ratio_udc = c(NA,
                                           NA,
                                           ICER),
                             utility = c(dt_sc_util + sum(utility_sc$discounted_utility),
                                             dt_pc_util + sum(utility_pc$discounted_utility),
                                             (dt_pc_util + sum(utility_pc$discounted_utility)) -
                                               (dt_sc_util + sum(utility_sc$discounted_utility))),
                             cost = c(dt_sc_cost + sum(MC_costs_sc$discounted_cost),
                                          dt_pc_cost + sum(MC_costs_pc$discounted_cost),
                                          (dt_pc_cost + sum(MC_costs_pc$discounted_cost)) -
                                            (dt_sc_cost + sum(MC_costs_sc$discounted_cost))),
                             ratio = c(NA,
                                       NA,
                                       ICER_disc))

if (SAVE_ARM_COMPARISON){
  fwrite(arm_comparison,
         file = paste(SAVE_FILEPATH,
                      "arm_comparison.csv",
                      sep = ""))
}

#### Further exploration of decision tree results ####
{
print("Decision tree results:")
print(paste("Expected cost under SC =",
      dt_sc_cost))
print(paste("Expected cost under PC =",
      dt_pc_cost))
print(paste("Expected utility under SC =",
            dt_sc_util))
print(paste("Expected utility under PC =",
            dt_pc_util))
print(paste("Expected cost difference =",
            dt_pc_cost-dt_sc_cost))
print(paste("Expected utility difference =",
            dt_pc_util-dt_sc_util))
print(paste("Post-DT ICER =",
            (dt_pc_cost-dt_sc_cost)/(dt_pc_util-dt_sc_util)))
}

# For clarity, make decision tree states comparable with those from the Excel
# workbook:

dt_results_sc <- dt_results_sc %>%
  filter(grepl("event|death|post", subpop)) %>%
  filter(grepl("ac_no_lof", subpop))
dt_states <- c("no_event",
               "post_stroke",
               "post_mi",
               "death")

# Merge ac_lof and ac_no_lof rows; there is definitely a better way to do this
# using dplyr but I can't quite figure it out!
# for (state in dt_states){
#   ac_lof_row = dt_results_sc[dt_results_sc$subpop==paste("ac_lof_", state, sep=""), ]
#   ac_no_lof_row = dt_results_sc[dt_results_sc$subpop==paste("ac_no_lof_", state, sep=""), ]
#   ac_row <- ac_lof_row %>% mutate(subpop = paste("ac_", state, sep=""))
#   ac_row$prob <- ac_lof_row$prob + ac_no_lof_row$prob
#   
#   #Indexing here extracts expected utility and cost:
#   ac_row[c(3,4)] <- lof_prev * ac_lof_row[c(3,4)] + 
#     (1-lof_prev) * ac_no_lof_row[c(3,4)]
#   dt_results_sc <- dt_results_sc %>%
#     add_row(ac_row)
# }


dt_results_pc <- dt_results_pc %>%
  filter(grepl("event|death|post", subpop)) %>%
  filter(grepl("ac_lof", subpop))
dt_states <- c("no_event",
               "post_stroke",
               "post_mi",
               "death")

# # Merge ac_lof and ac_no_lof rows; there is definitely a better way to do this
# # using dplyr but I can't quite figure it out!
# for (state in dt_states){
#   ac_lof_row = dt_results_pc[dt_results_pc$subpop==paste("ac_lof_", state, sep=""), ]
#   ac_no_lof_row = dt_results_pc[dt_results_pc$subpop==paste("ac_no_lof_", state, sep=""), ]
#   ac_row <- ac_lof_row %>% mutate(subpop = paste("ac_", state, sep=""))
#   ac_row$prob <- ac_lof_row$prob + ac_no_lof_row$prob
#   
#   #Indexing here extracts expected utility and cost:
#   ac_row[c(3,4)] <- lof_prev * ac_lof_row[c(3,4)] + 
#     (1-lof_prev) * ac_no_lof_row[c(3,4)]
#   dt_results_pc <- dt_results_pc %>%
#     add_row(ac_row)
# }
