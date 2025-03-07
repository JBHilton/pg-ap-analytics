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
library("ggplot2")
library("readxl")
library("stringr")
library("tidyverse")

# Set time horizon
time_hor <- 10

#### Start by setting up names ####
# There are four subpopulations based on genotype and drug, and six health
# states in the Markov cohort model.

# Define subpopulations
subpop_names <- c("clo_lof",
                  "clo_no_lof",
                  "tic",
                  "pra")

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
pc_uptake <- .75
pc_test_cost <- 100.
pc_resource_cost <- 100.
pc_sens <- .95
pc_spec <- .95

l_uptake <- .75
l_test_cost <- 100.
l_resource_cost <- 100.
l_sens <- .95
l_spec <- .95

p_test_followed = .9

p_clo_sc <- parameters_STEMI$mean[parameters_STEMI$parameter.list=="proportion_of_clo_lof_sc"] # Proportion prescriped clopidogrel under standard care
p_tic_sc <- parameters_STEMI$mean[parameters_STEMI$parameter.list=="proportion_of_tic_sc"]  # Proportion prescriped ticagrelor under standard care
p_pra_sc <- parameters_STEMI$mean[parameters_STEMI$parameter.list=="proportion_of_pra_sc"]  # Proportion prescriped prasugrel under standard care

p_tic_A <- .5 # Proportion prescriped ticagrelor during subroutine A
p_pra_A <- .5 # Proportion prescriped prasugrel during subroutine A

lof_prev_eu <- .30
lof_prev_as <- .58

lof_prev <- lof_prev_eu

drug_costs <- data.frame(drug = c("clo", "tic", "pra"),
                         cost = c(100, 100, 100))

event_utilities <- parameters_STEMI %>%
  filter(grepl("utility", parameter.list)) %>%
  select(c(parameter.list, mean)) %>%
  mutate(parameter.list = gsub("utility_of_", "", parameter.list)) %>%
  add_row(parameter.list = "death",
          mean = 0)

# Put some filler values in for utilities of bleeds
minor_bleed_utility <- 0.9
major_bleed_utility <- 0.9

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
  p_clo_lof <- ifelse(true_genotype=="no_lof", 0, p_clo_sc)
  p_clo_no_lof <- ifelse(true_genotype=="no_lof", p_clo_sc, 0)
  p_tri <- p_tic_sc
  p_pra <- p_pra_sc
  
  sc_results <- data.frame(subpops = subpop_names,
                                     prob = c(p_clo_lof,
                                              p_clo_no_lof,
                                              p_tri,
                                              p_pra),
                           cost = c(drug_costs$cost[drug_costs$drug=="clo"],
                                    drug_costs$cost[drug_costs$drug=="clo"],
                                    drug_costs$cost[drug_costs$drug=="tic"],
                                    drug_costs$cost[drug_costs$drug=="pra"]),
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
  p_clo_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 0, 1)
  
  # If test is followed, test says no lof, and true genotype is no_lof, then
  # patient goes to clopidogrel_no_lof
  p_clo_no_lof <- p_test_followed *
    ifelse(test_result=="no_lof", 1, 0) *  
    ifelse(true_genotype=="no_lof", 1, 0)
  
  # If test is followed and test says lof, and true genotype is lof, then
  # patient is prescribed one of the other drugs
  p_tri <- p_test_followed *
    ifelse(test_result=="no_lof", 0, 1) * 
    p_tic_A
  p_pra <- p_test_followed *
    ifelse(test_result=="no_lof", 0, 1) * 
    p_pra_A
  
  subroutine_A_results <- data.frame(subpop = c("clo_lof",
                                                   "clo_no_lof",
                                                   "tic",
                                                   "pra",
                                                   "sc"),
             prob = c(p_clo_lof,
                      p_clo_no_lof,
                      p_tri,
                      p_pra,
                      1 - p_test_followed), # Go to standard care if test not followed
             cost = c(drug_costs$cost[drug_costs$drug=="clo"],
                      drug_costs$cost[drug_costs$drug=="clo"],
                      drug_costs$cost[drug_costs$drug=="tic"],
                      drug_costs$cost[drug_costs$drug=="pra"],
                      0),
             utility = c(1,
                         1,
                         1,
                         1,
                         1)) # costs assigned in standard care
  return(subroutine_A_results)
}

implement_B <- function(subpop_id){
  subpop_pars <- parameters_STEMI[grepl(paste(subpop_id, sep=""), # Gather subpopulation-specific parameters plus utility weights
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
  
  bleed_results <- data.frame(event = c("minor_bleed",
                                     "major_bleed",
                                     "no_event"),
                           prob = c(minor_bleed_prob,
                                    major_bleed_prob,
                                    1 - minor_bleed_prob - major_bleed_prob),
                           utility = c(minor_bleed_utility,
                                       major_bleed_utility,
                                       1))
  exp_util_bleed <- sum(bleed_results$prob * bleed_results$utility)
  
  # Note on results: we currently assign zero cost to each of these outcomes
  subroutine_B_results <- data.frame(event = c("reinfarction",
                                     "stroke",
                                     "death",
                                     "no_event"),
                              prob = c(reinfarction_prob,
                                       stroke_prob,
                                       death_prob,
                                       1 - reinfarction_prob - stroke_prob - death_prob),
                              cost = c(0,
                                       0,
                                       0,
                                       0),
                           utility = exp_util_bleed * c(event_utilities$mean[event_utilities$parameter.list == "utility_of_mi"],
                                       event_utilities$mean[event_utilities$parameter.list == "utility_of_stroke"],
                                       0,
                                       event_utilities$mean[event_utilities$parameter.list == "utility_of_no_further_event"]))
  return(subroutine_B_results)
}

# Patient steps through subroutine B regardless of testing method, so we
# calculate it here
B_result_list <- lapply(subpop_names,
                        implement_B)

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
  
  # Only need results of subroutine A if we actually do testing
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
      test_uptake * A_results_ave$cost[1:4] +
      prob_sc *
    (lof_prev * sc_results_lof$prob * sc_results_lof$cost +
      (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost)
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
      grep(name,
           extended_patient_status$subpop)] <- c_subgroup
    
    # Assign probabilities and update utilities based on subroutine B.
    # Note: this assigns utility 1 to anyone in 
    for (j in (1:nrow(post_event_match))){
      extended_patient_status$prob[
        grep(paste(name, post_event_match$state[j], sep=""),
                                        extended_patient_status$subpop)] <-
        B_result_list[[i]]$prob[
          which(B_result_list[[i]]$event==post_event_match$event[j])] *
        p_subgroup
      
      extended_patient_status$exp_cost[
        grep(paste(name, post_event_match$state[j], sep=""),
             extended_patient_status$subpop)] <-
        extended_patient_status$exp_cost[
          grep(paste(name, post_event_match$state[j], sep=""),
               extended_patient_status$subpop)] +
        B_result_list[[i]]$cost[
          which(B_result_list[[i]]$event==post_event_match$event[j])]
      
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



# Build blocks containing transition matrices for each subpopulation
matrix_blocks <- lapply(subpop_names,
                        FUN = function(subpop_id){
                          A = build_markov_submodel(subpop_id)[[1]]
                          return(A)})

# Assemble into block matrix containing all events for all subpopulations
full_transition_matrix <- bdiag(matrix_blocks) %>%
  as.matrix()
rownames(full_transition_matrix) <- full_names
colnames(full_transition_matrix) <- full_names

# Add utility by state:
util_by_state <- replicate(n_subpops,
                           event_utilities,
                           simplify = FALSE) %>%
  bind_rows()

#### Now run model ####
# The prob column from the dataframe outputted by the run_forward function is
# the initial condition for the Markov model.

# Choose time horizon
n_tsteps <- 10

# Choose test to analyse
test <- "pc"

pc_dt_results <- run_forward(test)

P0 <- as.vector(pc_dt_results$prob)
markov_trace <- sapply(1:n_tsteps,
                        FUN = function(t){
                          P0 %*% (full_transition_matrix %^% t)
                          }
                        ) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(full_names) %>%
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
                         allpop_util = as.matrix(markov_trace %>%
                                                select(-time_step)) %*%
                           as.matrix(util_by_state$mean))

# Slightly hacky piece of code to add subpopulation-specific utilities. The
# function adds the expected utility over time for a given subpopulation, then
# we add them by cycling over the indices:
add_subpop <- function(df, i) {
  cols_to_select <- grepl(subpop_names[i],
                          colnames(markov_trace))
  varname <- paste(subpop_names[i],
                   "_util",
                   sep = "")[1]
  mutate(df, !!varname :=
           c((1 / rowSums(as.matrix(markov_trace[1, cols_to_select]))) *
           as.matrix(markov_trace[, cols_to_select]) %*%
           as.vector(util_by_state$mean[cols_to_select])))
}
for (i in (1:length(subpop_names))){
  utility_df <- add_subpop(utility_df, i)
}


utility_df %>%
  pivot_longer(grep("util",colnames(utility_df)),
               names_to = "subpop",
               values_to = "exp_util") %>%
  mutate(subpop = gsub("_util", "", subpop)) %>%
  ggplot(aes(x = time_step,
             y = exp_util,
             colour = subpop)) +
  geom_line() +
  labs(title = paste("Test strategy", test))
