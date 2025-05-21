# This file contains the functions for running the decision tree analyses and
# the Markov model.

library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

# Define functions to calculate end state probablities after passing through
# some combination of standard care, subroutine A, and subroutine B.

# Standard care
implement_sc <- function(true_genotype,
                         sc_props = c(p_ac_sc,
                                      p_at_sc,
                                      p_ap_sc)){
  
  # If prescribed clopidogrel then subpopulation depends on genotype, otherwise
  # just used presciption proportions.
  p_ac_lof <- sc_props[1]
  p_ac_no_lof <- 0
  p_at <- sc_props[2]
  p_ap <- sc_props[3]
  
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
                                       (1-AVE_TIME_TO_EVENT) *
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
run_forward <- function(test = "sc",
                        sc_props = c(p_ac_sc,
                                     p_at_sc,
                                     p_ap_sc),
                        poct_cost = pc_test_cost){
  patient_status <- data.frame(subpop = subpop_names,
                               prob = rep(1, length(subpop_names)),
                               exp_cost = rep(0, length(subpop_names)),
                               exp_utility = rep(1, length(subpop_names)))
  
  # Always need to calculate standard care results as all routes have a chance
  # of going to it
  sc_results_lof <- implement_sc("lof", sc_props)
  sc_results_no_lof <- implement_sc("no_lof", sc_props)
  
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
      test_uptake * poct_cost + # Account for cost of testing
      test_uptake * A_results_ave$cost[1:4] +
      prob_sc *
      (lof_prev * sc_results_lof$prob * sc_results_lof$cost +
         (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost)
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


# Function for building Markov model
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

run_arm_comparison <- function(
    sc_props = c(p_ac_sc,
                 p_at_sc,
                 p_ap_sc),
    poct_cost = pc_test_cost){
  
  # Run decision tree analysis
  dt_results_pc <- run_forward("pc",
                               sc_props = sc_props,
                               poct_cost = poct_cost) %>%
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
  
  MT_pc <- sapply(1:(time_hor-1),
                  FUN = function(t){
                    P0_pc %*% (build_markov_model(t) %^% t)
                  }
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(markov_states) %>%
    mutate(time_step = 1:(time_hor-1))
  
  
  # Calculate expected utility over time:
  utility_pc <- data.frame(time_step = 1:(time_hor-1),
                           utility = rowSums(as.matrix(MT_pc %>%
                                                         select(-time_step)) *
                                               as.matrix(markov_utils))) %>%
    mutate(discounted_utility =
             utility * discount_by_cycle)
  # Expected costs:
  MC_costs_pc <- data.frame(time_step = 1:(time_hor-1),
                            cost = rowSums(as.matrix(MT_pc %>%
                                                       select(-time_step)) %*%
                                             as.matrix(markov_costs$value))) %>%
    mutate(discounted_cost =
             cost * discount_by_cycle)
  
  
  # Mean life years:
  life_years_pc <- (1-P0_pc["death"]) + sum(1 - MT_pc$death)
  
  # Now do sc
  
  # Run decision tree analysis
  dt_results_sc <- run_forward("sc", sc_props) %>%
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
  
  MT_sc <- sapply(1:(time_hor-1),
                  FUN = function(t){
                    P0_sc %*% (build_markov_model(t) %^% t)
                  }
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(markov_states) %>%
    mutate(time_step = 1:(time_hor-1))
  
  
  # Calculate expected utility over time:
  utility_sc <- data.frame(time_step = 1:(time_hor-1),
                           utility = rowSums(as.matrix(MT_sc %>%
                                                         select(-time_step)) *
                                               as.matrix(markov_utils))) %>%
    mutate(discounted_utility =
             utility * discount_by_cycle)
  # Expected costs:
  MC_costs_sc <- data.frame(time_step = 1:(time_hor-1),
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
  
  return(arm_comparison)
}