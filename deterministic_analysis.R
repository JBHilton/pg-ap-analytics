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
SAVE_FILEPATH <- "stemi_"

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("janitor")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

source("make_parameters.R")
source("simulation_functions.R")

#### ICER calculation ####

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
MT_pc <- sapply(Reduce("%*%", lapply(-1:(n_tsteps-1),
                                     FUN = function(t){
                                       if (t==-1){
                                         return(diag(nrow = n_states))
                                       }else{
                                         if (t==0){
                                           return(build_markov_model(1))
                                         }else{
                                           return(build_markov_model(t))
                                         }
                                       }
                                     }), accumulate = TRUE),
                FUN = function(A){
                  P0_pc %*% A
                }
) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 0:n_tsteps)


# Calculate expected utility over time:
utility_pc <- utils_from_markov_trace(MT_pc,
                                      markov_utils)
# Expected costs:
MC_costs_pc <- costs_from_markov_trace(MT_pc,
                                       markov_costs)


# Mean life years:
life_years_pc <- sum(1 - c(MT_pc$death[1],
                           0.5 * (MT_pc$death[2:39] + MT_pc$death[3:40])))

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

MT_sc <- sapply(Reduce("%*%", lapply(-1:(n_tsteps-1),
                                     FUN = function(t){
                                       if (t==-1){
                                         return(diag(nrow = n_states))
                                       }else{
                                         if (t==0){
                                           return(build_markov_model(1))
                                         }else{
                                           return(build_markov_model(t))
                                         }
                                       }
                                       }), accumulate = TRUE),
                FUN = function(A){
                  P0_sc %*% A
                }
) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(markov_states) %>%
  mutate(time_step = 0:n_tsteps)


# Calculate expected utility over time:

utility_sc <- utils_from_markov_trace(MT_sc,
                                      markov_utils)

# Expected costs:
MC_costs_sc <- costs_from_markov_trace(MT_sc,
                                       markov_costs)

# Mean life years:
life_years_sc <- sum(1 - c(MT_sc$death[1],
                           0.5 * (MT_sc$death[2:39] + MT_sc$death[3:40])))

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$halfstep)) - (dt_sc_cost + sum(MC_costs_sc$halfstep))) /
  ((dt_pc_util + sum(utility_pc$halfstep)) - (dt_sc_util + sum(utility_sc$halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
{
  print(paste("Mean life years under SC = ",
              (life_years_sc)))
  print(paste("Mean life years under PC = ",
              (life_years_pc)))
  print("")
  print("With discounting:")
  print(paste("Mean cost under PC = ",
        (dt_pc_cost + sum(MC_costs_pc$discounted_cost))))
  print(paste("Mean cost under SC = ",
              (dt_sc_cost + sum(MC_costs_sc$discounted_cost))))
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
              (dt_pc_cost + sum(MC_costs_pc$undiscounted_cost))))
  print(paste("Mean cost under SC = ",
              (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))))
  print(paste("Mean utility under PC = ",
              (dt_pc_util + sum(utility_pc$undiscounted_utility))))
  print(paste("Mean utility under SC = ",
              (dt_sc_util + sum(utility_sc$undiscounted_utility))))
  print(paste("Estimated incremental cost is",
              (dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))))
  print(paste("Estimated incremental utility is",
              (dt_pc_util + sum(utility_pc$undiscounted_utility)) -
                (dt_sc_util + sum(utility_sc$undiscounted_utility))))
  print(paste("Estimated ICER is",
              ICER_undisc))
  }

# Create a dataframe storing outputs by case:
arm_comparison <- data.frame(arm = c("sc", "pc", "inc"),
                             lifeyears = c(life_years_sc,
                                            life_years_pc,
                                            life_years_pc - life_years_sc),
                             util_udc = c(dt_sc_util + sum(utility_sc$halfstep),
                                                dt_pc_util + sum(utility_pc$halfstep),
                                                (dt_pc_util + sum(utility_pc$halfstep)) -
                                                  (dt_sc_util + sum(utility_sc$halfstep))),
                             cost_udc = c(dt_sc_cost + sum(MC_costs_sc$halfstep),
                                             dt_pc_cost + sum(MC_costs_pc$halfstep),
                                             (dt_pc_cost + sum(MC_costs_pc$halfstep)) -
                                               (dt_sc_cost + sum(MC_costs_sc$halfstep))),
                             ratio_udc = c(NA,
                                              NA,
                                              ICER_undisc),
                             util = c(dt_sc_util + sum(utility_sc$discounted_halfstep),
                                               dt_pc_util + sum(utility_pc$discounted_halfstep),
                                               (dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                                 (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                             cost = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep),
                                            dt_pc_cost + sum(MC_costs_pc$discounted_halfstep),
                                            (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                              (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                             icer = c(NA,
                                             NA,
                                             ICER_disc),
                             nmb = 20000 * c(dt_sc_util + sum(utility_sc$discounted_halfstep),
                                                dt_pc_util + sum(utility_pc$discounted_halfstep),
                                                (dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                                  (dt_sc_util + sum(utility_sc$discounted_halfstep))) -
                               c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep),
                                 dt_pc_cost + sum(MC_costs_pc$discounted_halfstep),
                                 (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)))) %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column(var = "output_name")

if (SAVE_ARM_COMPARISON){
  fwrite(arm_comparison,
         file = paste(SAVE_FILEPATH,
                      "arm_comparison.csv",
                      sep = ""))
}

bleed_results_sc <- aggregate_bleed_results("sc")
bleed_results_pc <- aggregate_bleed_results("pc")

# Count events per 1000:
events_df <- data.frame(arm = c("sc",
                                "pc"),
                        mi_dt = 1000 * c(sum(dt_results_sc$prob[grepl("mi", dt_results_sc$event)]),
                                         sum(dt_results_pc$prob[grepl("mi", dt_results_pc$event)])),
                        mi_mc = 1000 * c(sum(MT_sc$mi),
                                         sum(MT_pc$mi)),
                        stroke_dt = 1000 * c(sum(dt_results_sc$prob[grepl("stroke", dt_results_sc$event)]),
                                         sum(dt_results_pc$prob[grepl("stroke", dt_results_pc$event)])),
                        stroke_mc = 1000 * c(sum(MT_sc$stroke),
                                         sum(MT_pc$stroke)),
                        major_bleed_dt = 1000 * c(bleed_results_sc$major,
                                                  bleed_results_pc$major),
                        minor_bleed_mc = 1000 * c(bleed_results_sc$minor,
                                                  bleed_results_pc$minor),
                        death_dt = 1000 * c(sum(dt_results_sc$prob[grepl("death", dt_results_sc$event)]),
                                         sum(dt_results_pc$prob[grepl("death", dt_results_pc$event)])),
                        death_mc = 1000 * c(sum(MT_sc$death),
                                         sum(MT_pc$death)))
if (SAVE_ARM_COMPARISON){
  fwrite(events_df,
         file = paste(SAVE_FILEPATH,
                      "event_counts.csv",
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

# Add extra column to PC arm results which subtracts probability mass
# corresponding to individuals who enter standard care, for ease of comparison
# with outputs in workbook
dt_results_pc <- dt_results_pc %>%
  filter(grepl("event|death|post", subpop)) %>%
  filter(!grepl("ac_lof", subpop)) %>%
  mutate(prob_sc = dt_results_sc$prob) %>%
  mutate(prob_minus_sc = ifelse(grepl("^ac", subpop),
                                prob,
                                prob - ((1 - pc_uptake) + pc_uptake * (1 - p_test_followed)) * prob_sc)) %>%
  select(-prob_sc)
