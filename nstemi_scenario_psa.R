# In this script we load parameters and assemble all the parameter objects we
# will need for the analysis.

# Flag to determine whether to convert rates in the data to probabilities
CONVERSIONS_REQUIRED <- FALSE

# Set to true to directly calculate probabilities for the no-LOF subpopulation,
# or false to use values from other subpopulations
UNIQUE_NO_LOF_PROBS <- TRUE

# Step correction for death in decision tree and Markov model, i.e. average time
# to death given death occurs within a given year.
AVE_TIME_TO_EVENT <- 0.5

CASE <- "NSTEMI"

SAVE_OUTPUTS <- TRUE

n_sample <- 1e4
SAVE_FILEPATH <- paste("nstemi_n_", n_sample, sep="")

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("matrixStats")
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

parameters_STEMI <- read_xlsx("data-inputs/NSTEMI_masterfile_160725.xlsm",
                               sheet = "Parameters.NSTEMI") %>%
  rename_all(make.names) %>% # Converts parameter names in valid R names
  rename_all(.funs = function(name){
    name %>%
      make.names() %>%
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

source("psa_functions.R")

# Load utilities and add zeros for death
base_markov_utils <- read_xlsx("data-inputs/NSTEMI_masterfile_160725.xlsm",
                               sheet = "markov",
                               range = "AH7:AM48") %>% # Skip row 1 since this doesn't match the format of other rows
  rename_all(make.names) %>% # Converts parameter names in valid R names
  rename_all(.funs = function(name){
    name %>%
      make.names() %>%
      str_replace("\\.", "_") %>%
      str_replace("reinfarction", "mi") %>%
      str_replace("dead", "death")})  %>%
  filter(!row_number() %in% c(1,2)) %>%
  # mutate(death = 0) %>%
  select(all_of(markov_states)) %>% # Last step just makes sure ordering matches model
  mutate_at(markov_states, as.numeric)
base_markov_utils <- base_markov_utils[1:39, ] # Don't end up using last entry

# Function for fixing variable names for tree utilities:
rename_utility_variables <- function(name){
  if (str_detect(name, "^utility_")){
    return(paste0(name, "_tree"))
  }else{
    return(name)
  }
} %>% Vectorize()

uncertainty_df <- read_xlsx("data-inputs/NSTEMI_masterfile_160725.xlsm",
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
           SE,
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
  mutate(distribution = str_to_lower(distribution)) %>%
  # mutate(variable.name = map(variable.name, rename_utility_variables)) %>%
  mutate(SE = as.numeric(SE),
         par1 = as.numeric(par1),
         par2 = as.numeric(par2)) # Convert values from characters to numbers

# Reorder so that when we take common lines with parameter dataframe everything
# is in the correct row:
uncertainty_df <- uncertainty_df[order(match(uncertainty_df$variable.name,
                                             parameters_STEMI$variable.name)), ]

#### Do PSA for each scenario ####

n_tsteps <- time_hor - 1
scenarios <- c("SA1",
               "SA2",
               "SA3",
               "SA4",
               "SA5",
               "SA6",
               "SA7",
               "SA8",
               "SA9",
               "SA10",
               "SA11",
               "SA12")

# Get central estimates
baseline_results <- lapply(scenarios,
                           FUN = function(scenario){
                             draw_df <- uncertainty_df %>% mutate(draw = Value)
                             psa_results <- run_PSA_arm_comparison(parameters_STEMI,
                                                                   draw_df,
                                                                   scenario = scenario)
                             res_df <- psa_results[[1]] %>%
                               cbind(psa_results[[3]] %>% # Attach event counts
                                       pivot_wider(names_from = arm, 
                                                   values_from = c(-arm)))%>%
                               mutate(icer = inc_cost_dc_hs / inc_util_dc_hs) %>%
                               mutate(scenario = scenario)
                           }) %>%
  bind_rows() %>%
  group_by(scenario) %>%
  select(-matches("dc_(.)c$")) %>% # Drop all cost/util columns without half step correction
  select(-matches("^inc_(\\w+)dc$")) %>% # Drop all cost/util columns without half step correction
  rename_all(.funs = function(name){ # Remove half step correction indicator
    name %>%
      str_replace_all("_hs", "") %>%
      str_replace_all("_dc", "")
  })

if (SAVE_OUTPUTS){
  fwrite(baseline_results,
         file = "nstemi_scenario_central_estimate.csv")
}

start_time <- Sys.time()
multi_results <- lapply(scenarios,
                        FUN = function(scenario){
                          lapply(1:n_sample,
                          FUN = function(i){
                            draw_i <- do_PSA_draw(uncertainty_df)
                            psa_results <- run_PSA_arm_comparison(parameters_STEMI,
                                                                  draw_i,
                                                                  scenario = scenario)
                            res_i <- psa_results[[1]] %>%
                              mutate(sample_id = as.character(i)) %>%
                              cbind(psa_results[[3]] %>% # Attach event counts
                                      pivot_wider(names_from = arm, 
                                                  values_from = c(-arm)))
                            return(res_i)
                          }) %>%
    bind_rows() %>%
    mutate(icer = inc_cost_dc_hs / inc_util_dc_hs) %>%
    mutate(scenario = scenario)
                          }
    ) %>%
  bind_rows() %>%
  group_by(scenario)
end_time <- Sys.time()
print(paste("PSA for",
            n_sample,
            "samples conducted in",
            difftime(end_time,
                     start_time,
                     units = "secs"),
            "seconds."))

# Add indicator for whether ICER passes threshold for each sample
ce_threshold = 20000


# Function to calculate the acceptance probability for a given cost
# effectiveness threshold given a set of ICERs from PSA
calculate_acc_prob <- function(utils,
                               icers){
  length(which((icers < ce_threshold) & (utils > 0))) /
    length(icers)
}

ce_thresh_df <- multi_results %>%
  summarise(ce_prob = calculate_acc_prob(inc_util_dc_hs, icer))

mr_split <- split(multi_results,
                  multi_results$scenario)

# Add bootstrap samples:
n_bootstrap <- 1000
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:n_sample,
                       n_sample,
                       replace = TRUE)
  ce_by_scenario <- sapply(mr_split,
                           FUN = function(df){
                             sample_utils <- df$inc_util_dc_hs[sample_ids]
                             sample_icers <- df$icer[sample_ids]
                             calculate_acc_prob(sample_utils, sample_icers)
                           })
  ce_thresh_df[[varname]] <- ce_by_scenario
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                 probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))]

if (SAVE_OUTPUTS){
  fwrite(ce_thresh_df,
         file = paste(SAVE_FILEPATH,
                      "_psa_acceptance_probability.csv",
                      sep=""))
}


ci_df <- data.frame(scenario = scenarios,
                    life_years_sc_mean = summarise(multi_results, mean = mean(life_years_sc))$mean,
                    life_years_sc_L = summarise(multi_results, q = quantile((life_years_sc), 0.025))$q,
                    life_years_sc_U = summarise(multi_results, q = quantile((life_years_sc), 0.975))$q,
                    cost_udc_sc_mean = summarise(multi_results, mean = mean(cost_udc_sc_hs))$mean,
                    cost_udc_sc_L = summarise(multi_results, q = quantile((cost_udc_sc_hs), 0.025))$q,
                    cost_udc_sc_U = summarise(multi_results, q = quantile((cost_udc_sc_hs), 0.975))$q,
                    util_udc_sc_mean = summarise(multi_results, mean = mean(util_udc_sc_hs))$mean,
                    util_udc_sc_L = summarise(multi_results, q = quantile((util_udc_sc_hs), 0.025))$q,
                    util_udc_sc_U = summarise(multi_results, q = quantile((util_udc_sc_hs), 0.975))$q,
                    cost_sc_mean = summarise(multi_results, mean = mean(cost_dc_sc_hs))$mean,
                    cost_sc_L = summarise(multi_results, q = quantile((cost_dc_sc_hs), 0.025))$q,
                    cost_sc_U = summarise(multi_results, q = quantile((cost_dc_sc_hs), 0.975))$q,
                    util_sc_mean = summarise(multi_results, mean = mean(util_dc_sc_hs))$mean,
                    util_sc_L = summarise(multi_results, q = quantile((util_dc_sc_hs), 0.025))$q,
                    util_sc_U = summarise(multi_results, q = quantile((util_dc_sc_hs), 0.975))$q,
                    nmb_sc_mean = summarise(multi_results, mean = mean(nmb_sc))$mean,
                    nmb_sc_L = summarise(multi_results, q = quantile((nmb_sc), 0.025))$q,
                    nmb_sc_U = summarise(multi_results, q = quantile((nmb_sc), 0.975))$q,
                    life_years_pc_mean = summarise(multi_results, mean = mean(life_years_pc))$mean,
                    life_years_pc_L = summarise(multi_results, q = quantile((life_years_pc), 0.025))$q,
                    life_years_pc_U = summarise(multi_results, q = quantile((life_years_pc), 0.975))$q,
                    cost_udc_pc_mean = summarise(multi_results, mean = mean(cost_udc_pc_hs))$mean,
                    cost_udc_pc_L = summarise(multi_results, q = quantile((cost_udc_pc_hs), 0.025))$q,
                    cost_udc_pc_U = summarise(multi_results, q = quantile((cost_udc_pc_hs), 0.975))$q,
                    util_udc_pc_mean = summarise(multi_results, mean = mean(util_udc_pc_hs))$mean,
                    util_udc_pc_L = summarise(multi_results, q = quantile((util_udc_pc_hs), 0.025))$q,
                    util_udc_pc_U = summarise(multi_results, q = quantile((util_udc_pc_hs), 0.975))$q,
                    cost_pc_mean = summarise(multi_results, mean = mean(cost_dc_pc_hs))$mean,
                    cost_pc_L = summarise(multi_results, q = quantile((cost_dc_pc_hs), 0.025))$q,
                    cost_pc_U = summarise(multi_results, q = quantile((cost_dc_pc_hs), 0.975))$q,
                    util_pc_mean = summarise(multi_results, mean = mean(util_dc_pc_hs))$mean,
                    util_pc_L = summarise(multi_results, q = quantile((util_dc_pc_hs), 0.025))$q,
                    util_pc_U = summarise(multi_results, q = quantile((util_dc_pc_hs), 0.975))$q,
                    nmb_pc_mean = summarise(multi_results, mean = mean(nmb_pc))$mean,
                    nmb_pc_L = summarise(multi_results, q = quantile((nmb_pc), 0.025))$q,
                    nmb_pc_U = summarise(multi_results, q = quantile((nmb_pc), 0.975))$q,
                    life_years_inc_mean = summarise(multi_results, mean = mean(life_years_inc))$mean,
                    life_years_inc_L = summarise(multi_results, q = quantile((life_years_inc), 0.025))$q,
                    life_years_inc_U = summarise(multi_results, q = quantile((life_years_inc), 0.975))$q,
                    cost_udc_inc_mean = summarise(multi_results, mean = mean(inc_cost_udc_hs))$mean,
                    cost_udc_inc_L = summarise(multi_results, q = quantile((inc_cost_udc_hs), 0.025))$q,
                    cost_udc_inc_U = summarise(multi_results, q = quantile((inc_cost_udc_hs), 0.975))$q,
                    util_udc_inc_mean = summarise(multi_results, mean = mean(inc_util_udc_hs))$mean,
                    util_udc_inc_L = summarise(multi_results, q = quantile((inc_util_udc_hs), 0.025))$q,
                    util_udc_inc_U = summarise(multi_results, q = quantile((inc_util_udc_hs), 0.975))$q,
                    cost_inc_mean = summarise(multi_results, mean = mean(inc_cost_dc_hs))$mean,
                    cost_inc_L = summarise(multi_results, q = quantile((inc_cost_dc_hs), 0.025))$q,
                    cost_inc_U = summarise(multi_results, q = quantile((inc_cost_dc_hs), 0.975))$q,
                    util_inc_mean = summarise(multi_results, mean = mean(inc_util_dc_hs))$mean,
                    util_inc_L = summarise(multi_results, q = quantile((inc_util_dc_hs), 0.025))$q,
                    util_inc_U = summarise(multi_results, q = quantile((inc_util_dc_hs), 0.975))$q,
                    nmb_inc_mean = summarise(multi_results, mean = mean(inc_nmb))$mean,
                    nmb_inc_L = summarise(multi_results, q = quantile((inc_nmb), 0.025))$q,
                    nmb_inc_U = summarise(multi_results, q = quantile((inc_nmb), 0.975))$q,
                    icer_mean = summarise(multi_results, mean = mean(icer))$mean,
                    icer_L = summarise(multi_results, q = quantile((icer), 0.025))$q,
                    icer_U = summarise(multi_results, q = quantile((icer), 0.975))$q
)
if (SAVE_OUTPUTS){
  fwrite(ci_df,
         file = paste(SAVE_FILEPATH,
                      "_scenario_psa_stats.csv",
                      sep = ""))
}