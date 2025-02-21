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

library("data.tree")
library("dplyr")
library("readxl")
library("stringr")
library("tidyverse")

# Set time horizon
time_hor <- 10

# Define subpopulations
subpop_names <- c("clo_lof",
                  "clo_no_lof",
                  "tic",
                  "pra")

parameters_STEMI <- read_xlsx("data-inputs/MasterFile-PGXACS.xlsx",
                              sheet = "Parameters.STEMI")[-1,] %>% # Skip row 1 since this doesn't match the format of other rows
  rename_all(make.names) %>% # Converts parameter names in valid R names
  type.convert(as.is = TRUE) %>% # Make sure numbers are numbers, not characters
  filter(!is.na(parameter.list)) %>% # Remove empty lines
  mutate(parameter.list = parameter.list %>%
           str_to_lower() %>% # Standardise parameter names for use with other objects
           str_replace_all(" ", "_") %>%
           str_replace_all("-", "_") %>%
           str_replace_all("with_", "") %>%
           str_replace_all("in_", "") %>%
           str_replace_all("standard_care", "sc") %>%
           str_replace_all("clopidogrel_no_lof", "clo_no_lof") %>%
           str_replace_all("clopidogrel", "clo_lof") %>%
           str_replace_all("ticagrelor", "tic") %>%
           str_replace_all("prasugrel", "pra"))

# Function to convert rate parameters to probability of happening within 1 time unit
rates_to_probs <- function(rate, unit){
  return(ifelse(grepl("rate", unit),
         yes = 1 - exp(-rate),
         no = NA))
}

# Add these probabilities to parameter table
parameters_STEMI <- parameters_STEMI %>%
  mutate(prob = rates_to_probs(mean, unit))


# Start ages specified on row 1 of parameter matrices:
start_age_male <- 62
start_age_female <- 69

# List parameters for decision tree
pc_test_cost <- 100.
pc_resource_cost <- 100.
pc_sens <- .95
pc_spec <- .95

l_test_cost <- 100.
l_resource_cost <- 100.
l_sens <- .95
l_spec <- .95

p_test_followed = .9

p_clo_sc <- parameters_STEMI["proportion_of_clopidogrel_in_standard_care", "mean"] # Proportion prescriped clopidogrel under standard care
p_tic_sc <- parameters_STEMI["proportion_of_ticagrelor_in_standard_care", "mean"]  # Proportion prescriped ticagrelor under standard care
p_pra_sc <- parameters_STEMI["proportion_of_prasugrel_in_standard_care", "mean"]  # Proportion prescriped prasugrel under standard care

p_tic_A <- .5 # Proportion prescriped ticagrelor during subroutine A
p_pra_A <- .5 # Proportion prescriped prasugrel during subroutine A

lof_prev_eu <- .30
lof_prev_as <- .58

# Define functions to calculate end state probablities after passing through
# some combination of standard care, subroutine A, and subroutine B.

# Standard care
implement_sc <- function(true_genotype){
  
  # If prescribed clopidogrel then subpopulation depends on genotype, otherwise
  # just used presciption proportions.
  p_clo_lof <- ifelse(true_genotype=="no_lof", 0, p_clo_sc)
  p_clo_no_lof <- ifelse(true_genotype=="no_lof", p_clo_sc, 0)
  p_tri <- p_test_followed * p_tic_sc
  p_pra <- p_test_followed * p_pra_sc
  
  sc_results <- data.frame(row.names = c("clo_lof",
                                         "clo_no_lof",
                                         "ticr",
                                         "pra"),
                                     prob = c(p_clo_lof,
                                              p_clo_no_lof,
                                              p_tri,
                                              p_pra),
                                     cost = c(clopidogrel_cost,
                                              clopidogrel_cost,
                                              ticagrelor_cost,
                                              prasugrel_cost))
  return(sc_results)
}

implement_A <- function(true_genotype,
                        test_result
                        ){
  
  # If test is followed, test says no lof, and true genotype is lof, then
  # patient goes to clopidogrel_lof
  p_clo_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 1, 0)
  
  # If test is followed, test says no lof, and true genotype is no_lof, then
  # patient goes to clopidogrel_no_lof
  p_clo_no_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 0, 1)
  
  # If test is followed and test says lof, and true genotype is lof, then
  # patient is prescribed one of the other drugs
  p_tri <- p_test_followed * p_tic_A
  p_pra <- p_test_followed * p_pra_A
  
  subroutine_A_results <- data.frame(subpop = c("clo_lof",
                                                   "clo_no_lof",
                                                   "ticr",
                                                   "pra",
                                                   "standard_care"),
             prob = c(p_clo_lof,
                      p_clo_no_lof,
                      p_tri,
                      p_pra,
                      1 - p_test_followed), # Go to standard care if test not followed
             cost = c(clopidogrel_cost,
                      clopidogrel_cost,
                      ticagrelor_cost,
                      prasugrel_cost,
                      0)) # costs assigned in standard care
  return(subroutine_A_results)
}

implement_B <- function(subpop_id){
  subpop_pars <- parameters_STEMI[grepl(paste(subpop_id, "|utility", sep=""), # Gather subpopulation-specific parameters plus utility weights
                                        parameters_STEMI$parameter.list),] %>%
    mutate(parameter.list = gsub(paste("_", subpop_id, sep=""), "", parameter.list))
  
  # Get probabilities of each event from dataframe:
  minor_bleed_prob <- subpop_pars$prob[
    grep("minor_bleeding", subpop_pars$parameter.list)]
  major_bleed_prob <- subpop_pars$prob[
    grep("minor_bleeding", subpop_pars$parameter.list)]
  reinfarction_prob <- subpop_pars$prob[
    grep("reinfarction", subpop_pars$parameter.list)]
  stroke_prob <- subpop_pars$prob[
    grep("stroke", subpop_pars$parameter.list)]
  death_prob <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  # Get utilities attached to each event
  reinfarction_util <- subpop_pars$prob[
    grep("reinfarction", subpop_pars$parameter.list)]
  stroke_util <- subpop_pars$prob[
    grep("stroke", subpop_pars$parameter.list)]
  death_util <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  subroutine_B_results <- data.frame(event = c("no_bleed",
                                                   "minor_bleed",
                                                   "major_bleed",
                                                   "no_event",
                                                   "reinfarction",
                                                   "stroke",
                                                   "death"),
                                     prob = c(1 - minor_bleed_prob - major_bleed_prob,
                                              minor_bleed_prob,
                                              major_bleed_prob,
                                              1 - reinfarction_prob - stroke_prob - death_prob,
                                              reinfarction_prob,
                                              stroke_prob,
                                              death_prob
                                              ))
  return(subroutine_B_results)
}


### Set up Markov cohort model ####

# Get list of states for Markov cohort model
markov_states <- c("no_event",
                   "stroke",
                   "post_stroke",
                   "mi",
                   "post_mi",
                   "death")
n_states <- length(markov_states)

# Assemble generic transition matrix without accounting for sex, genotype, or drug
base_transition_matrix <- matrix(0,
                                 nrow = n_states,
                            ncol = n_states,
                            dimnames = list(markov_states,
                                            markov_states))

build_markov_submodel <- function(subpop_id){
  subpop_pars <- parameters_STEMI[grepl(paste(subpop_id, "|utility", sep=""), # Gather subpopulation-specific parameters plus utility weights
                                        parameters_STEMI$parameter.list),] %>%
    mutate(parameter.list = gsub(paste("_", subpop_id, sep=""), "", parameter.list))
  
  subpop_transition_matrix <- base_transition_matrix
  
  # Get probabilities of each event from dataframe:
  
  # Outflow from no_event
  subpop_transition_matrix["no_event", "mi"] <- subpop_pars$prob[
    grep("reinfarction", subpop_pars$parameter.list)]
  subpop_transition_matrix["no_event", "stroke"] <- subpop_pars$prob[
    grep("^stroke", subpop_pars$parameter.list)]
  subpop_transition_matrix["no_event", "death"] <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  # Outflow from mi
  subpop_transition_matrix["mi", "post_mi"] <- 1. - subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  subpop_transition_matrix["mi", "death"] <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  # Outflow from post_mi
  subpop_transition_matrix["post_mi", "death"] <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  # Outflow from stroke
  subpop_transition_matrix["stroke", "post_stroke"] <- 1. - subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  subpop_transition_matrix["stroke", "death"] <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  # Outflow from post_stroke
  subpop_transition_matrix["post_stroke", "death"] <- subpop_pars$prob[
    grep("all_cause_mortality", subpop_pars$parameter.list)]
  
  diag(subpop_transition_matrix) <- 1 - rowSums(subpop_transition_matrix)
  
  # Get utilities attached to each event
  no_event_util <- subpop_pars$mean[
    grep("utility_of_no_further_event", subpop_pars$parameter.list)]
  mi_util <- subpop_pars$mean[
    grep("utility_of_mi", subpop_pars$parameter.list)]
  post_mi_util <- subpop_pars$mean[
    grep("utility_of_post_mi", subpop_pars$parameter.list)]
  stroke_util <- subpop_pars$mean[
    grep("utility_of_mi", subpop_pars$parameter.list)]
  post_stroke_util <- subpop_pars$mean[
    grep("utility_of_post_stroke", subpop_pars$parameter.list)]
  death_util <- 0
  
  subpop_util <- data.frame(event = markov_states,
                            utility = c(no_event_util,
                                        stroke_util,
                                        post_stroke_util,
                                        mi_util,
                                        post_mi_util,
                                        death_util))
  return(list(subpop_transition_matrix,
              subpop_util))
}

n_subpops <- length(subpop_names)

# The following pastes together all combinations of subpopulation and health
# state names. Indexing is done using modular arithmetic, so x %% y is remainder
# in x/y and x %/% y is x/y without the remainder.
full_names <- sapply(1:(n_subpops*n_states),
                     FUN = function(i){
                       paste(subpop_names[(i+n_states-1) %/% n_states],
                             markov_states[((i+n_states-1) %% n_states)+1],
                             sep="_")
                     })

# Build blocks containing transition matrices for each subpopulation
matrix_blocks <- lapply(subpop_names,
                        FUN = function(subpop_id){
                          A = build_markov_submodel(subpop_id)[[1]]
                          return(A)})

# Assemble into block matrix containing all events for all subpopulations
full_transition_matrix <- bdiag(matrix_blocks)
rownames(full_transition_matrix) <- full_names
colnames(full_transition_matrix) <- full_names

full_utility_df <- data.frame(event = 0,
                            utility = )

#### Now run model ####
