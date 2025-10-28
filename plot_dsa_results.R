# In this script we generate a tornado plot visualising the deterministic
# sensitivity analysis results

dir.create("plots", showWarnings = FALSE)

library("data.table")
library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("matrixStats")
library("readxl")
library("stringr")
library("tibble")
library("tidyverse")

# varnames <- data.frame(short = c("annual_utility_dec_major_bleed",
#                                  "annual_utility_dec_minor_bleed",
#                                  "day_cost_clop",
#                                  "day_cost_tica",
#                                  "poct_cost",
#                                  "prevalence_lof_base_case",
#                                  "prob_major_bleed_ac_lof",
#                                  "prob_mi_ac_lof",
#                                  "prob_minor_bleed_ac_lof",
#                                  "prob_stroke_ac_lof",
#                                  "prob_test_followed",
#                                  "prob_test_order",
#                                  "rr_death_ac_no_lof",
#                                  "rr_mi_ac_no_lof",
#                                  "u_dec_dyspnoea",
#                                  "utility_no_event"),
#                        long = c("Annual utility dec., major bleed",
#                                 "Annual utility dec., minor bleed",
#                                 "Daily cost, clopidogrel",
#                                 "Daily cost, ticagrelor",
#                                 "Cost of POC testing",
#                                 "LOF allele prevalence",
#                                 "Prob. major bleed in AC general population",
#                                 "Prob. reinfarction in AC general population",
#                                 "Prob. minor bleed in AC general population",
#                                 "Prob. stroke in AC general population",
#                                 "Prob. clinician follows test results",
#                                 "Prob. test ordered",
#                                 "Risk ratio for death, AC no-LOF to LOF",
#                                 "Risk ratio for reinfarction, AC no-LOF to LOF",
#                                 "Utility dec., dyspnoea",
#                                 "Utility of no event state"))


### First do STEMI population ###

varnames <- read_xlsx("data-inputs/masterfile_111025.xlsx",
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
  select(variable.name, parameter.list) %>%
  dplyr::rename(scenario = variable.name) %>%
  add_row(scenario = c("prob_nevent_to_stk",
                       "prob_nevent_to_mi"),
          parameter.list = c("Stroke probability in Markov model",
                             "Reinfarction probability in Markov model"))

stemi_arm_comparison <- read.csv("outputs/stemi_arm_comparison.csv")
baseline_inc_nmb <- stemi_arm_comparison$inc[
  which(stemi_arm_comparison$output_name=="nmb")]

stemi_dsa_results <- read.csv(
  file = "outputs/stemi_one_way_sensitivity_analysis.csv") %>%
  select(c(scenario,
           inc_nmb)) %>%
  mutate(direction = ifelse(grepl("high", scenario),
                            "High",
                            "Low")) %>%
  mutate(scenario = str_replace(scenario,
                                "high_|low_",
                                "")) %>%
  pivot_wider(names_from = direction,
              values_from = inc_nmb) %>%
  mutate(width = abs(High - Low)) %>%
  pivot_longer(cols = c(High, Low),
               names_to = "direction",
               values_to = "inc_nmb") %>%
  filter(width > 20) %>%
  mutate(diff = inc_nmb - baseline_inc_nmb) %>%
  left_join(varnames, by = c("scenario" = "scenario")) %>%
  mutate(scenario = factor(parameter.list) %>%
           fct_reorder(width))

p_stemi <-stemi_dsa_results %>%
  ggplot(aes(x = inc_nmb,
             y = scenario,
             fill = direction)) +
  geom_col(orientation = "y") +
  geom_vline(aes(xintercept = baseline_inc_nmb),
             alpha = .5) +
  scale_x_continuous(trans = scales::trans_new("shift",
                                               transform = function(x) {x - baseline_inc_nmb},
                                               inverse = function(x) {x + baseline_inc_nmb})) +
  xlab("Incremental net monetary benefit") +
  ylab("")

ggsave(paste("plots/stemi_dsa_tornado.png", sep=""),
       plot = p_stemi,
       width = 10,
       height = 10,
       dpi = 1200)



### Now do NSTEMI population ###


varnames <- read_xlsx("data-inputs/NSTEMI_masterfile_160725.xlsm",
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
  select(variable.name, parameter.list) %>%
  dplyr::rename(scenario = variable.name) %>%
  add_row(scenario = c("prob_nevent_to_stk",
                       "prob_nevent_to_mi",
                       "utility_post_mi"),
          parameter.list = c("Stroke probability in Markov model",
                             "Reinfarction probability in Markov model",
                             "Utiliy post-reinfarction"))

nstemi_arm_comparison <- read.csv("outputs/nstemi_arm_comparison.csv")
baseline_inc_nmb <- nstemi_arm_comparison$inc[
  which(nstemi_arm_comparison$output_name=="nmb")]

nstemi_dsa_results <- read.csv(
  file = "outputs/nstemi_one_way_sensitivity_analysis.csv") %>%
  select(c(scenario,
           inc_nmb)) %>%
  mutate(direction = ifelse(grepl("high", scenario),
                            "High",
                            "Low")) %>%
  mutate(scenario = str_replace(scenario,
                                "high_|low_",
                                "")) %>%
  pivot_wider(names_from = direction,
              values_from = inc_nmb) %>%
  mutate(width = abs(High - Low)) %>%
  pivot_longer(cols = c(High, Low),
               names_to = "direction",
               values_to = "inc_nmb") %>%
  filter(width > 20) %>%
  mutate(diff = inc_nmb - baseline_inc_nmb) %>%
  left_join(varnames, by = c("scenario" = "scenario")) %>%
  mutate(scenario = factor(parameter.list) %>%
           fct_reorder(width))

p_nstemi <- nstemi_dsa_results %>%
  ggplot(aes(x = inc_nmb,
             y = scenario,
             fill = direction)) +
  geom_col(orientation = "y") +
  geom_vline(aes(xintercept = baseline_inc_nmb),
             alpha = .5) +
  scale_x_continuous(trans = scales::trans_new("shift",
                                               transform = function(x) {x - baseline_inc_nmb},
                                               inverse = function(x) {x + baseline_inc_nmb})) +
  xlab("Incremental net monetary benefit") +
  ylab("")

ggsave(paste("plots/nstemi_dsa_tornado.png", sep=""),
       plot = p_nstemi,
       width = 10,
       height = 10,
       dpi = 1200)

