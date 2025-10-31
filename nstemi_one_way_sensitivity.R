# In this script we perform a 1-way deterministic sensitivity analysis


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

dir.create("outputs",
           showWarnings = FALSE)
SAVE_FILEPATH <- "outputs/nstemi_"

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
                              sheet = "Parameters.NSTEMI") %>% # Skip row 1 since this doesn't match the format of other rows
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
  mutate(SE = as.numeric(SE),
         par1 = as.numeric(par1),
         par2 = as.numeric(par2)) # Convert values from characters to numbers

# Reorder so that when we take common lines with parameter dataframe everything
# is in the correct row:
uncertainty_df <- uncertainty_df[order(match(uncertainty_df$variable.name,
                                             parameters_STEMI$variable.name)), ]

# Set number of time steps
n_tsteps <- time_hor - 1

# Check what happens when we do arm comparison with baseline parameters
det_df <- uncertainty_df %>%
  mutate(draw = Value)
det_results <- run_PSA_arm_comparison(parameters_STEMI,
                                      det_df)
sc_comparison <- det_results[[2]] %>%
  filter(arm == "SC") %>%
  select(-c(arm, ratio_udc, ratio_dc, ratio_udc_hs, ratio_dc_hs))
pc_comparison <- det_results[[2]] %>%
  filter(arm == "PC") %>%
  select(-c(arm, ratio_udc, ratio_dc, ratio_udc_hs, ratio_dc_hs))
inc_comparison <- det_results[[2]] %>%
  filter(arm == "Increment") %>%
  select(-c(arm, ratio_udc, ratio_dc, ratio_udc_hs, ratio_dc_hs))

quantile_from_keyword <- function(dist_name,
                                  par1,
                                  par2,
                                  direction,
                                  base_cost = 0., # Fill this in if working with costs
                                  print_errors = TRUE){
  if (is.na(dist_name)){
    if (direction=="high"){
      p <- 1.2
      return(p * base_cost)
    }
    if (direction=="low"){
      p <- .8
      return(p * base_cost)
    }
    
  }else{
    if (direction=="high"){
      p <- .975
    }
    if (direction=="low"){
      p <- .025
    }
    if (!(dist_name %in% c("beta",
                           "gamma",
                           "lognormal"))){
      if (print_errors){
        print(paste(dist_name, 
                    "is not a valid distribution name.",
                    sep = " "))
      }
      return(NA)
    }
    if (dist_name == "beta"){
      return(qbeta(p,
                   par1,
                   par2))
    }
    if (dist_name == "gamma"){
      return(qgamma(p,
                    shape = par1,
                    scale = par2))
    }
    if (dist_name == "lognormal"){
      return(qlnorm(p,
                    par1,
                    par2))
    }
  }
} %>% Vectorize()

make_one_way_pars <- function(u_df,
                              direction,
                              varname){
  var_idx <- which(u_df$variable.name==varname)
  u_df$draw <- u_df$Value
  u_df$draw[var_idx] <- quantile_from_keyword(u_df$distribution[var_idx],
                                              u_df$par1[var_idx],
                                              u_df$par2[var_idx],
                                              direction,
                                              base_cost = ifelse(is.na(u_df$distribution[var_idx]),
                                                                 yes = u_df$Value[var_idx],
                                                                 no = 0))
  return(u_df)
}

pars_to_vary <- c(uncertainty_df$variable.name,
                  "poct_cost",
                  "day_cost_tica",
                  "day_cost_clop")

dsa_options <- data.frame(varname = pars_to_vary) %>%
  crossing(direction = c("high",
                         "low"))

missing_sa_pars <- base::setdiff(pars_to_vary,
                                 uncertainty_df$variable.name)

uncertainty_df <- uncertainty_df %>%
  rbind.fill(data.frame(variable.name = missing_sa_pars,
                        Value = parameters_STEMI %>%
                          filter(variable.name %in% missing_sa_pars) %>%
                          arrange(match(variable.name, missing_sa_pars)) %>%
                          select(value) %>%
                          unlist() )) %>%
  rowwise()

# Reorder so that when we take common lines with parameter dataframe everything
# is in the correct row:
uncertainty_df <- uncertainty_df[order(match(uncertainty_df$variable.name,
                                             parameters_STEMI$variable.name)), ]

start_time <- Sys.time()
dsa_results <- lapply(1:nrow(dsa_options),
                      FUN = function(i){
                        draw_i <- make_one_way_pars(uncertainty_df,
                                                    dsa_options$direction[i],
                                                    dsa_options$varname[i]
                        )
                        arm_results <- run_PSA_arm_comparison(parameters_STEMI,
                                                              draw_i)
                        res_i <- arm_results[[1]] %>%
                          mutate(scenario = paste(dsa_options$direction[i],
                                                  "_",
                                                  dsa_options$varname[i],
                                                  sep = "")
                          ) %>%
                          relocate(scenario, .before = life_years_sc) %>%
                          cbind(arm_results[[3]] %>% # Attach event counts
                                  pivot_wider(names_from = arm, 
                                              values_from = c(-arm)))
                        return(res_i)
                      }) %>%
  bind_rows() %>%
  mutate(icer = inc_cost_dc_hs / inc_util_dc_hs)
end_time <- Sys.time()
print(paste("One-way sensitivity analysis conducted in",
            difftime(end_time,
                     start_time,
                     units = "secs"),
            "seconds."))
if (SAVE_OUTPUTS){
  fwrite(dsa_results,
         file = paste(SAVE_FILEPATH,
                      "one_way_sensitivity_analysis.csv",
                      sep = ""))
}