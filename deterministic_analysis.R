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

# Set to true to directly calculate probabilities for the no-LOF subpopulation,
# or false to use values from other subpopulations
UNIQUE_NO_LOF_PROBS <- FALSE

# Step correction for death in decision tree and Markov model, i.e. average time
# to death given death occurs within a given year.
AVE_TIME_TO_EVENT <- 0.5

# Set to true to save the results of the arm comparison as a .csv file. The file
# path can be set on the following line:
SAVE_ARM_COMPARISON <- TRUE
SAVE_FILEPATH <- ""

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

source("make_parameters.R")
source("simulation_functions.R")

# Basic patient status dataframe

patient_status <- data.frame(subpop = subpop_names,
                             prob = rep(0, length(subpop_names)),
                             exp_cost = rep(0, length(subpop_names)),
                             exp_utility = rep(1, length(subpop_names)))


  


#### ICER calculation ####
# Now that we've defined the modelling workflow, let's calculate an ICER
# for point of care vs standard care:

n_tsteps <- time_hor - 1


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

markov_update <- function(p, t){
  print(t)
  p <- p %*% (build_markov_model(t))
  print(p)
  return(p)
}
MT_pc <- sapply(Reduce("%*%", lapply(1:n_tsteps,
                              build_markov_model), accumulate = TRUE),
                FUN = function(A){
                  P0_pc %*% A
                }
                ) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 1:n_tsteps)

MT_pc_alt <- sapply(1:n_tsteps,
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

MT_sc <- sapply(Reduce("%*%", lapply(1:n_tsteps,
                                     build_markov_model), accumulate = TRUE),
                FUN = function(A){
                  P0_sc %*% A
                }
) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 1:n_tsteps)

MT_sc_alt <- sapply(1:n_tsteps,
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
# 
dt_results_sc <- dt_results_sc %>%
  filter(grepl("event|death|post", subpop)) %>%
  filter(!grepl("ac_no_lof", subpop))


dt_results_pc <- dt_results_pc %>%
  filter(grepl("event|death|post", subpop)) %>%
  filter(!grepl("ac_lof", subpop))
