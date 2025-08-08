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

CASE <- "STEMI"

SAVE_OUTPUTS <- FALSE

n_sample <- 1e2
SAVE_FILEPATH <- paste("stemi_n_", n_sample, sep="")

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

parameters_STEMI <- read_xlsx("data-inputs/masterfile_070725.xlsx",
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
  mutate(value = as.numeric(value),
         sample.size = as.numeric(sample.size)) # Convert values from characters to numbers

source("psa_functions.R")

# Load utilities and add zeros for death
base_markov_utils <- read_xlsx("data-inputs/masterfile_240325.xlsx",
                               sheet = "time_event_utility") %>%
  mutate(death = 0) %>%
  select(all_of(markov_states)) # Last step just makes sure ordering matches model
base_markov_utils <- base_markov_utils[1:39, ] # Don't end up using last entry

# Function for fixing variable names for tree utilities:
rename_utility_variables <- function(name){
  if (str_detect(name, "^utility_")){
    return(paste0(name, "_tree"))
  }else{
    return(name)
  }
} %>% Vectorize()

uncertainty_df <- read_xlsx("data-inputs/masterfile_070725.xlsx",
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

draw_df <- do_PSA_draw(uncertainty_df)
rewrite <- rewrite_pars_from_draw(parameters_STEMI, draw_df)
par_df <- rewrite[[1]]
markov_df <- rewrite[[2]]

# Try multiple draws at once:
tall_df <- do_tall_PSA_draw(uncertainty_df,
                            10)
n_tsteps <- time_hor - 1
example_results <- run_PSA_arm_comparison(parameters_STEMI,
                                          draw_df)

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

start_time <- Sys.time()
multi_results <- lapply(1:n_sample,
                        FUN = function(i){
                          draw_i <- do_PSA_draw(uncertainty_df)
                          res_i <- run_PSA_arm_comparison(parameters_STEMI,
                                                                    draw_i)[[1]] %>%
                            mutate(sample_id = as.character(i))
                          return(res_i)
                        }) %>%
  bind_rows() %>%
  mutate(icer = inc_cost_dc_hs / inc_util_dc_hs) %>% 
  summarise(sample_id = c(sample_id, 'mean'),
            across(where(is.numeric), ~ c(., mean(.))))
end_time <- Sys.time()
print(paste("PSA for",
            n_sample,
            "samples conducted in",
            difftime(end_time,
                     start_time,
                     units = "secs"),
            "seconds."))

ce_line_df <- data.frame(x = c(min(multi_results$inc_util_dc_hs),
                               max(multi_results$inc_util_dc_hs),
                               min(multi_results$inc_util_dc_hs),
                               max(multi_results$inc_util_dc_hs)),
                         ce_threshold = c(2*10e4,
                                          2*10e4,
                                          3*10e4,
                                          3*10e4)) %>%
  mutate(y = x * ce_threshold)

ce_plane_plot_with_mean <- ggplot(multi_results[1:n_sample, ],
                                  aes(x = inc_util_dc_hs,
                                      y = inc_cost_dc_hs)) +
  geom_point() +
  geom_point(data = multi_results[n_sample+1, ],
             aes(x = inc_util_dc_hs,
                 y = inc_cost_dc_hs,
                 color = "red"),
             size = 5.,
             shape = 18) +
  theme(legend.position = "none") +
  xlab("Incremental utility") +
  ylab("Incremental cost")
ce_plane_plot_with_lines <- ggplot(multi_results[1:n_sample, ],
                                  aes(x = inc_util_dc_hs,
                                      y = inc_cost_dc_hs)) +
  geom_point() +
  geom_line(data = ce_line_df,
            aes(x=x,
                y=y,
                color = factor(ce_threshold))) +
  scale_color_manual(labels = c("£20,000",
                                "£30,000"),
                     values = c("orange",
                                "green")) +
  labs(x = "Incremental utility (QALY's)",
       y = "Incremental cost (£)",
       color = "Cost-effectiveness threshold")

# Add indicator for whether ICER passes threshold for each sample
ce_threshold = seq(from = 0.,
                   to = 30000.,
                   by = 500.)
for (i in 1:length(ce_threshold)){
  varname <- paste("ce_threshold_", ce_threshold[[i]])
  multi_results[[varname]] <- ifelse(
    multi_results$icer < ce_threshold[[i]],
    1,
    0)
}

# Function to calculate the acceptance probability for a given cost
# effectiveness threshold given a set of ICERs from PSA
calculate_acc_prob <- function(ce_threshold,
                               utils,
                               icers){
  length(which((icers < ce_threshold) & (which(utils > 0)))) /
    n_sample
} %>% Vectorize()

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:n_sample],
                              multi_results$icer[1:n_sample]))

# Add bootstrap samples:
n_bootstrap <- 1000
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:n_sample,
                       n_sample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                 probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))]
if (SAVE_OUTPUTS){
  fwrite(ce_thresh_df,
         file = paste(SAVE_FILEPATH,
                      "_acceptance_probability.csv",
                      sep=""))
}

ce_thresh_plot <- ggplot(ce_thresh_df,
                         aes(x = ce_threshold,
                             y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
              ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")


ci_df <- data.frame(lifeyears_sc_mean = mean(multi_results$life_years_sc[1:n_sample]),
                    lifeyears_sc_L = quantile(multi_results$life_years_sc[1:n_sample], .025),
                    lifeyears_sc_U = quantile(multi_results$life_years_sc[1:n_sample], .975),
                    cost_udc_sc_mean = mean(multi_results$sc_cost_udc_hs[1:n_sample]),
                    cost_udc_sc_L = quantile(multi_results$sc_cost_udc_hs[1:n_sample], .025),
                    cost_udc_sc_U = quantile(multi_results$sc_cost_udc_hs[1:n_sample], .975),
                    util_udc_sc_mean = mean(multi_results$sc_util_udc_hs[1:n_sample]),
                    util_udc_sc_L = quantile(multi_results$sc_util_udc_hs[1:n_sample], .025),
                    util_udc_sc_U = quantile(multi_results$sc_util_udc_hs[1:n_sample], .975),
                    cost_sc_mean = mean(multi_results$sc_cost_dc_hs[1:n_sample]),
                    cost_sc_L = quantile(multi_results$sc_cost_dc_hs[1:n_sample], .025),
                    cost_sc_U = quantile(multi_results$sc_cost_dc_hs[1:n_sample], .975),
                    util_sc_mean = mean(multi_results$sc_util_dc_hs[1:n_sample]),
                    util_sc_L = quantile(multi_results$sc_util_dc_hs[1:n_sample], .025),
                    util_sc_U = quantile(multi_results$sc_util_dc_hs[1:n_sample], .975),
                    nmb_sc_mean = mean(multi_results$sc_nmb[1:n_sample]),
                    nmb_sc_L = quantile(multi_results$sc_nmb[1:n_sample], .025),
                    nmb_sc_U = quantile(multi_results$sc_nmb[1:n_sample], .975),
                    lifeyears_pc_mean = mean(multi_results$life_years_pc[1:n_sample]),
                    lifeyears_pc_L = quantile(multi_results$life_years_pc[1:n_sample], .025),
                    lifeyears_pc_U = quantile(multi_results$life_years_pc[1:n_sample], .975),
                    cost_udc_pc_mean = mean(multi_results$pc_cost_udc_hs[1:n_sample]),
                    cost_udc_pc_L = quantile(multi_results$pc_cost_udc_hs[1:n_sample], .025),
                    cost_udc_pc_U = quantile(multi_results$pc_cost_udc_hs[1:n_sample], .975),
                    util_udc_pc_mean = mean(multi_results$pc_util_udc_hs[1:n_sample]),
                    util_udc_pc_L = quantile(multi_results$pc_util_udc_hs[1:n_sample], .025),
                    util_udc_pc_U = quantile(multi_results$pc_util_udc_hs[1:n_sample], .975),
                    cost_pc_mean = mean(multi_results$pc_cost_dc_hs[1:n_sample]),
                    cost_pc_L = quantile(multi_results$pc_cost_dc_hs[1:n_sample], .025),
                    cost_pc_U = quantile(multi_results$pc_cost_dc_hs[1:n_sample], .975),
                    util_pc_mean = mean(multi_results$pc_util_dc_hs[1:n_sample]),
                    util_pc_L = quantile(multi_results$pc_util_dc_hs[1:n_sample], .025),
                    util_pc_U = quantile(multi_results$pc_util_dc_hs[1:n_sample], .975),
                    nmb_pc_mean = mean(multi_results$pc_nmb[1:n_sample]),
                    nmb_pc_L = quantile(multi_results$pc_nmb[1:n_sample], .025),
                    nmb_pc_U = quantile(multi_results$pc_nmb[1:n_sample], .975),
                    lifeyears_inc_mean = mean(multi_results$life_years_inc[1:n_sample]),
                    lifeyears_inc_L = quantile(multi_results$life_years_inc[1:n_sample], .025),
                    lifeyears_inc_U = quantile(multi_results$life_years_inc[1:n_sample], .975),
                    cost_udc_inc_mean = mean(multi_results$inc_cost_udc_hs[1:n_sample]),
                    cost_udc_inc_L = quantile(multi_results$inc_cost_udc_hs[1:n_sample], .025),
                    cost_udc_inc_U = quantile(multi_results$inc_cost_udc_hs[1:n_sample], .975),
                    util_udc_inc_mean = mean(multi_results$inc_util_udc_hs[1:n_sample]),
                    util_udc_inc_L = quantile(multi_results$inc_util_udc_hs[1:n_sample], .025),
                    util_udc_inc_U = quantile(multi_results$inc_util_udc_hs[1:n_sample], .975),
                    cost_inc_mean = mean(multi_results$inc_cost_dc_hs[1:n_sample]),
                    cost_inc_L = quantile(multi_results$inc_cost_dc_hs[1:n_sample], .025),
                    cost_inc_U = quantile(multi_results$inc_cost_dc_hs[1:n_sample], .975),
                    util_inc_mean = mean(multi_results$inc_util_dc_hs[1:n_sample]),
                    util_inc_L = quantile(multi_results$inc_util_dc_hs[1:n_sample], .025),
                    util_inc_U = quantile(multi_results$inc_util_dc_hs[1:n_sample], .975),
                    nmb_inc_mean = mean(multi_results$inc_nmb[1:n_sample]),
                    nmb_inc_L = quantile(multi_results$inc_nmb[1:n_sample], .025),
                    nmb_inc_U = quantile(multi_results$inc_nmb[1:n_sample], .975),
                    icer_mean = mean(multi_results$icer[1:n_sample]),
                    icer_L = quantile(multi_results$icer[1:n_sample], .025),
                    icer_U = quantile(multi_results$icer[1:n_sample], .975),
                    row.names = ""
)

sc_cols <- colnames(ci_df) %>% grep(pattern = "sc", value = TRUE)
pc_cols <- colnames(ci_df) %>% grep(pattern = "pc", value = TRUE)
inc_cols <- colnames(ci_df) %>% grep(pattern = "inc", value = TRUE)
icer_cols <- colnames(ci_df) %>% grep(pattern = "icer", value = TRUE)

mean_cols <- colnames(ci_df) %>%
  grep(pattern = "mean", value = TRUE) %>%
  sub(pattern = "sc_|pc_|inc_",
      replace = "") %>%
  unique()
L_cols <- colnames(ci_df) %>%
  grep(pattern = "_L", value = TRUE) %>%
  sub(pattern = "sc_|pc_|inc_",
      replace = "") %>%
  unique()
U_cols <- colnames(ci_df) %>%
  grep(pattern = "_U", value = TRUE) %>%
  sub(pattern = "sc_|pc_|inc_",
      replace = "") %>%
  unique()

base_cols <- sub("_mean",
                 "",
                 mean_cols)

output_df <- ci_df %>%
  select(-all_of(icer_cols)) %>%
  pivot_longer(cols = all_of(sc_cols),
               values_to = "sc",
               names_to = "sc_names") %>%
  pivot_longer(cols = all_of(pc_cols),
               values_to = "pc",
               names_to = "pc_names") %>%
  pivot_longer(cols = all_of(inc_cols),
               values_to = "inc",
               names_to = "inc_names") %>%
  mutate(sc_names = sub("sc_", "", sc_names),
         pc_names = sub("pc_", "", pc_names),
         inc_names = sub("inc_", "", inc_names)) %>%
  filter((sc_names==pc_names) & (sc_names==inc_names)) %>%
  relocate(output_name = sc_names) %>%
  select(-c(pc_names, inc_names)) %>%
  mutate(sc_mean = ifelse(output_name %in% mean_cols,
                          sc,
                          NA),
         sc_L = ifelse(output_name %in% L_cols,
                          sc,
                          NA),
         sc_U = ifelse(output_name %in% U_cols,
                          sc,
                          NA),
         pc_mean = ifelse(output_name %in% mean_cols,
                          pc,
                          NA),
         pc_L = ifelse(output_name %in% L_cols,
                       pc,
                       NA),
         pc_U = ifelse(output_name %in% U_cols,
                       pc,
                       NA),
         inc_mean = ifelse(output_name %in% mean_cols,
                          inc,
                          NA),
         inc_L = ifelse(output_name %in% L_cols,
                       inc,
                       NA),
         inc_U = ifelse(output_name %in% U_cols,
                       inc,
                       NA)) %>%
  select(-c(sc, pc, inc)) %>%
  mutate(output_name = sub("_mean|_L|_U",
                          "",
                          output_name)) %>%
  as.list() %>%
  lapply(FUN = function(x){unique(x[!is.na(x)])}) %>%
  as.data.frame() %>%
  add_row(output_name = "icer",
          sc_mean = NA,
          sc_L = NA,
          sc_U = NA,
          pc_mean = NA,
          pc_L = NA,
          pc_U = NA,
          inc_mean = ci_df$icer_mean,
          inc_L = ci_df$icer_L,
          inc_U = ci_df$icer_U)

if (SAVE_OUTPUTS){
  fwrite(output_df,
         file = paste(SAVE_FILEPATH,
                      "_psa_stats.csv",
                      sep = ""))
}