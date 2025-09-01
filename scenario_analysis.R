# In this script we perform scenario analyses on of the deterministic model
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
SAVE_SCENARIO_ANALYSIS <- FALSE
SAVE_FILEPATH <- "stemi_"

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

ce_thresh <- 20000

#### Baseline scenario ####

source("make_parameters.R")

source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- data.frame(scenario = c("base"),
                          lygs_standard = c(life_years_sc),
                          lygs_poc = c(life_years_pc),
                          cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                          cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                          cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                          (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                          utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                          utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                          utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                                  (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                          icer = c(ICER_disc_hs)
                          )

#### Scenario 1: full test uptake ####

source("make_parameters.R")
pc_uptake <- 1.
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("prob test ordered = 100%"),
                          lygs_standard = c(life_years_sc),
                          lygs_poc = c(life_years_pc),
                          cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                          cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                          cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                          (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                          utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                          utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                          utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                             (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                          icer = c(ICER_disc_hs)
  )
)


#### Scenario 2: full test uptake ####

source("make_parameters.R")
p_test_followed <- .699
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("prob test followed = 69.9%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 3A: drug proportions ####

source("make_parameters.R")
p_ac_sc <- .1
p_at_sc <- .6
p_ap_sc <- .3
p_at_A <- p_at_sc / (1 - p_ac_sc)
p_ap_A <- p_ap_sc / (1 - p_ac_sc)
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("AC=10%, AT=60%, AP=30%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 3B: drug proportions ####

source("make_parameters.R")
p_ac_sc <- .01
p_at_sc <- .01
p_ap_sc <- .98
p_at_A <- p_at_sc / (1 - p_ac_sc)
p_ap_A <- p_ap_sc / (1 - p_ac_sc)
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("AC=1%, AT=1%, AP=98%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 3C: drug proportions ####

source("make_parameters.R")
p_ac_sc <- .01
p_at_sc <- .98
p_ap_sc <- .01
p_at_A <- p_at_sc / (1 - p_ac_sc)
p_ap_A <- p_ap_sc / (1 - p_ac_sc)
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("AC=1%, AT=98%, AP=1%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 4: LOF prevalence ####

source("make_parameters.R")
lof_prev <- lof_prev_as
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("LOF prev = 56.8%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 5: Baseline risk of stroke ####

source("make_parameters.R")
baseline_prob_df$value[
  grepl("stroke", baseline_prob_df$variable.name)] <- 0.0072
baseline_prob_df <- baseline_prob_df %>%
  mutate(odds = odds_from_prob(value))
prob_df <- rescale_probs(baseline_prob_df,
                         at_ratio_df,
                         ap_ratio_df,
                         ac_no_lof_ratio_df)
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Baseline risk of stroke = 0.0072"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )
#### Scenario 6: Baseline risk of MI ####

source("make_parameters.R")
baseline_prob_df$value[
  grepl("_mi_", baseline_prob_df$variable.name)] <- 0.054
baseline_prob_df <- baseline_prob_df %>%
  mutate(odds = odds_from_prob(value))
prob_df <- rescale_probs(baseline_prob_df,
                         at_ratio_df,
                         ap_ratio_df,
                         ac_no_lof_ratio_df)
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Baseline risk of MI = 0.054"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )
#### Scenario 7: Reducing SMRs ####

source("make_parameters.R")

smr_vals$no_event <- 1.6
smr_vals$mi <- 3.6
smr_vals$post_mi <- 2.4

mortality_prob_by_age <- read_xlsx("data-inputs/masterfile_070725.xlsx",
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
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Lower SMRs"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )
#### Scenario 8: Cost of minor bleed ####

source("make_parameters.R")
event_costs$value[event_costs$variable.name=="minor_bleed"] <- 893
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of minor bleed = £893"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 9A: discount_rate ####

source("make_parameters.R")
discount_by_cycle <- (1 / (1 + 0.015)^seq(
    1.0, time_hor-1, by = time_step))
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Discount rate = 1.5%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )

#### Scenario 9B: discount_rate ####

source("make_parameters.R")
discount_by_cycle <- (1 / (1 + 0.0)^seq(
  1.0, time_hor-1, by = time_step))
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Discount rate = 0%"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 10A: Cost of major bleed ####

source("make_parameters.R")
event_costs$value[event_costs$variable.name=="major_bleed"] <- 3528.
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of major bleed = £3528"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 10B: Cost of major bleed ####

source("make_parameters.R")
event_costs$value[event_costs$variable.name=="major_bleed"] <- 3849
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of major bleed = £3849"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 10C: Cost of major bleed ####

source("make_parameters.R")
event_costs$value[event_costs$variable.name=="major_bleed"] <- 4170
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of major bleed = £4170"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 10C: Cost of major bleed ####

source("make_parameters.R")
event_costs$value[event_costs$variable.name=="major_bleed"] <- 4490
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of major bleed = £4490"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Scenario 11: cost of POC test ####

source("make_parameters.R")
pc_test_cost <- 250.
source("simulation_functions.R")

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
life_years_pc <- sum(1 - MT_pc$death)

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
life_years_sc <- sum(1 - MT_sc$death)

# Now calculate ICER
ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
  ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
  ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
  ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))

# Create a dataframe storing outputs by case:
scenario_df <- scenario_df %>%
  rbind(data.frame(scenario = c("Cost of POC test = £250"),
                   lygs_standard = c(life_years_sc),
                   lygs_poc = c(life_years_pc),
                   cost_standard = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep)),
                   cost_poc = c(dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)),
                   cost_diff = c((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                   (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                   utility_standard = c(dt_sc_util + sum(utility_sc$discounted_halfstep)),
                   utility_poc = c(dt_pc_util + sum(utility_pc$discounted_halfstep)),
                   utility_diff = c((dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                      (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                   icer = c(ICER_disc_hs)
  )
  )


#### Save results ####

if (SAVE_SCENARIO_ANALYSIS){
  fwrite(scenario_df,
         file = paste(SAVE_FILEPATH,
                      "scenario_analysis.csv",
                      sep = ""))
}

