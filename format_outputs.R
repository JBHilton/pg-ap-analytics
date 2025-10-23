# In this script we load and reformat the different bits of output data into
# tables of base cases and confidence intervals.

SAVE_OUTPUTS <- FALSE
dir.create("formatted_outputs",
           showWarnings = FALSE)

library(stringr)
library(tidyverse)

stemi_arm_comparison <- read.csv("outputs/stemi_arm_comparison.csv")
stemi_psa <- read.csv("outputs/stemi_n_10000_psa_stats.csv")
stemi_prob_ce <- read.csv("outputs/stemi_n_10000_acceptance_probability.csv")
nstemi_arm_comparison <- read.csv("outputs/nstemi_arm_comparison.csv")
nstemi_psa <- read.csv("outputs/nstemi_n_10000_psa_stats.csv")
nstemi_prob_ce <- read.csv("outputs/nstemi_n_10000_acceptance_probability.csv")

# Extra stuff for checking
stemi_ac_short <- stemi_arm_comparison %>% filter(output_name %in% stemi_psa$output_name) %>% arrange(output_name)
stemi_psa <- stemi_psa %>% arrange(output_name)
nstemi_ac_short <- nstemi_arm_comparison %>% filter(output_name %in% nstemi_psa$output_name) %>% arrange(output_name)
nstemi_psa <- nstemi_psa %>% arrange(output_name)

# Population sizes for total NMB:
stemi_size = 41690
nstemi_size = 61248

stemi_nmb_df <- stemi_arm_comparison %>% filter(output_name == "nmb")
stemi_arm_comparison <- stemi_arm_comparison %>%
  add_row(output_name = "total_nmb",
          sc = stemi_size * stemi_nmb_df$sc,
          pc = stemi_size * stemi_nmb_df$pc,
          inc = stemi_size * stemi_nmb_df$inc)
stemi_psa_nmb <- stemi_psa %>%
  filter(output_name == "nmb") %>%
  select(-output_name)
stemi_psa <- stemi_psa %>%
  rbind(data.frame(output_name="total_nmb",
                   stemi_size * stemi_psa_nmb))

nstemi_nmb_df <- nstemi_arm_comparison %>% filter(output_name == "nmb")
nstemi_arm_comparison <- nstemi_arm_comparison %>%
  add_row(output_name = "total_nmb",
          sc = nstemi_size * nstemi_nmb_df$sc,
          pc = nstemi_size * nstemi_nmb_df$pc,
          inc = nstemi_size * nstemi_nmb_df$inc)
nstemi_psa_nmb <- nstemi_psa %>%
  filter(output_name == "nmb") %>%
  select(-output_name)
nstemi_psa <- nstemi_psa %>%
  rbind(data.frame(output_name="total_nmb",
                   nstemi_size * nstemi_psa_nmb))

print_stemi_output <- function(arm,
                         output,
                         digits=3){
  L_diff <- stemi_arm_comparison[which(stemi_arm_comparison$output_name==output),
                                 arm] - stemi_psa[which(stemi_psa$output_name==output),
                                                  paste(arm,
                                                        "_L",
                                                        sep="")]
  U_diff <- stemi_arm_comparison[which(stemi_arm_comparison$output_name==output),
                                 arm] - stemi_psa[which(stemi_psa$output_name==output),
                                                  paste(arm,
                                                        "_U",
                                                        sep="")]
  if ((L_diff<0)|(U_diff>0)){
    print(paste("Central estimate out of bounds for",
                output,
                "in STEMI",
                arm,
                "arm."))
  }
  paste(format(stemi_arm_comparison[which(stemi_arm_comparison$output_name==output),
                                    arm],
               digits = digits),
        " (",
        format(stemi_psa[which(stemi_psa$output_name==output),
                         paste(arm,
                               "_L",
                               sep="")],
               digits = digits),
        ", ",
        format(stemi_psa[which(stemi_psa$output_name==output),
                         paste(arm,
                               "_U",
                               sep="")],
               digits = digits),
        ")",
        sep="")
}

print_nstemi_output <- function(arm,
                               output,
                               digits=3){
  L_diff <- nstemi_arm_comparison[which(nstemi_arm_comparison$output_name==output),
                                 arm] - nstemi_psa[which(nstemi_psa$output_name==output),
                                                  paste(arm,
                                                        "_L",
                                                        sep="")]
  U_diff <- nstemi_arm_comparison[which(nstemi_arm_comparison$output_name==output),
                                 arm] - nstemi_psa[which(nstemi_psa$output_name==output),
                                                  paste(arm,
                                                        "_U",
                                                        sep="")]
  if ((L_diff<0)|(U_diff>0)){
    print(paste("Central estimate out of bounds for",
                output,
                "in NSTEMI",
                arm,
                "arm."))
  }
  paste(format(nstemi_arm_comparison[which(nstemi_arm_comparison$output_name==output),
                                    arm],
               digits = digits),
        " (",
        format(nstemi_psa[which(nstemi_psa$output_name==output),
                         paste(arm,
                               "_L",
                               sep="")],
               digits = digits),
        ", ",
        format(nstemi_psa[which(nstemi_psa$output_name==output),
                         paste(arm,
                               "_U",
                               sep="")],
               digits = digits),
        ")",
        sep="")
}

print_ce_prob <- function(ce_df,
                          digits){
  idx <- which.min(abs(ce_df$ce_threshold - 2*1e4))
  paste(format(ce_df$acceptance_prob[idx],
               digits = digits),
        ", (",
        format(ce_df$lower_95[idx],
               digits = digits),
        ", ",
        format(ce_df$upper_95[idx],
               digits = digits),
        ")",
        sep="")
}

main_df <- data.frame(group = rep("STEMI",
                                  2),
                      intervention = c("standard DAPT",
                                       "genotype-guided DAPT"),
                      life_years = sapply(c("sc", "pc"),
                                          print_stemi_output,
                                          output = "life_years",
                                          digits=5),
                      costs_udc = sapply(c("sc", "pc"),
                                         print_stemi_output,
                                         output = "cost_udc",
                                         digits=5),
                      costs_dc = sapply(c("sc", "pc"),
                                        print_stemi_output,
                                        output = "cost",
                                        digits=5),
                      qalys_udc = sapply(c("sc", "pc"),
                                         print_stemi_output,
                                         output = "util_udc",
                                         digits=5),
                      qalys_dc = sapply(c("sc", "pc"),
                                        print_stemi_output,
                                        output = "util",
                                        digits=5),
                      inc_costs = c(0,
                                    print_stemi_output("inc",
                                                       "cost",
                                                       digits=5)),
                      inc_qalys = c(0,
                                    print_stemi_output("inc",
                                                       "util",
                                                       digits=5)),
                      icer = c(0,
                               print_stemi_output("inc",
                                                  "icer",
                                                  digits=5)),
                      nmb_20k = sapply(c("sc", "pc"),
                                       print_stemi_output,
                                       output = "nmb",
                                       digits=5),
                      total_nmb = sapply(c("sc", "pc"),
                                       print_stemi_output,
                                       output = "total_nmb",
                                       digits=5),
                      inc_nmb = c(0,
                                  print_stemi_output("inc",
                                                     "nmb",
                                                     digits=5)),
                      percent_ce = c(0,
                                     print_ce_prob(stemi_prob_ce,
                                                   digits=5))
                      ) %>%
  rbind(
    data.frame(group = rep("NSTEMI",
                           2),
               intervention = c("standard DAPT",
                                "genotype-guided DAPT"),
               life_years = sapply(c("sc", "pc"),
                                   print_nstemi_output,
                                   output = "life_years",
                                   digits=5),
               costs_udc = sapply(c("sc", "pc"),
                                  print_nstemi_output,
                                  output = "cost_udc",
                                  digits=5),
               costs_dc = sapply(c("sc", "pc"),
                                 print_nstemi_output,
                                 output = "cost",
                                 digits=5),
               qalys_udc = sapply(c("sc", "pc"),
                                  print_nstemi_output,
                                  output = "util_udc",
                                  digits=5),
               qalys_dc = sapply(c("sc", "pc"),
                                 print_nstemi_output,
                                 output = "util",
                                 digits=5),
               inc_costs = c(0,
                             print_nstemi_output("inc",
                                                "cost",
                                                digits=5)),
               inc_qalys = c(0,
                             print_nstemi_output("inc",
                                                "util",
                                                digits=5)),
               icer = c(0,
                        print_nstemi_output("inc",
                                           "icer",
                                           digits=5)),
               nmb_20k = sapply(c("sc", "pc"),
                                print_nstemi_output,
                                output = "nmb",
                                digits=5),
               total_nmb = sapply(c("sc", "pc"),
                                  print_stemi_output,
                                  output = "total_nmb",
                                  digits=5),
               inc_nmb = c(0,
                           print_nstemi_output("inc",
                                              "nmb",
                                              digits=5)),
               percent_ce = c(0,
                              print_ce_prob(nstemi_prob_ce,
                                            digits=5))
    )
  ) %>% 
  mutate(grp_split = factor(group, levels = c("STEMI", "NSTEMI"))) %>%
  group_split(grp_split) %>% 
  map_dfr(~ add_row(.x, .after = Inf)) %>%
  select(-grp_split)

if (SAVE_OUTPUTS){
  write.csv(main_df,
            "formatted_outputs/base_case_and_PSA.csv")
}
# Now load in event count results

stemi_event_counts <- read.csv("outputs/stemi_event_counts.csv")
stemi_event_psa <- read.csv("outputs/stemi_n_1e+05_event_stats.csv")
nstemi_event_counts <- read.csv("outputs/nstemi_event_counts.csv")
nstemi_event_psa <- read.csv("outputs/nstemi_n_1e+05_event_stats.csv")

arm_pos <- function(arm) ifelse(arm=="sc", 1, 2)

print_stemi_events <- function(arm,
                               output,
                               digits=3){
  paste(format(stemi_event_counts[arm_pos(arm), output],
               digits = digits),
        " (",
        format(stemi_event_psa[1, paste(output, "_",  arm, "_L", sep="")],
               digits = digits),
        ", ",
        format(stemi_event_psa[1, paste(output, "_",  arm, "_U", sep="")],
               digits = digits),
        ")",
        sep="")
}

print_nstemi_events <- function(arm,
                                output,
                                digits=3){
  paste(format(nstemi_event_counts[arm_pos(arm), output],
               digits = digits),
        " (",
        format(nstemi_event_psa[1, paste(output, "_",  arm, "_L", sep="")],
               digits = digits),
        ", ",
        format(nstemi_event_psa[1, paste(output, "_",  arm, "_U", sep="")],
               digits = digits),
        ")",
        sep="")
}

events_df <- data.frame(group = rep("STEMI",
                                  2),
                      intervention = c("standard DAPT",
                                       "genotype-guided DAPT"),
                      mi_y1 = sapply(c("sc", "pc"),
                                     print_stemi_events,
                                     output = "mi_dt",
                                     digits=5),
                      mi_post_y1 = sapply(c("sc", "pc"),
                                          print_stemi_events,
                                          output = "mi_mc",
                                          digits=5),
                      stroke_y1 = sapply(c("sc", "pc"),
                                     print_stemi_events,
                                     output = "stroke_dt",
                                     digits=5),
                      stroke_post_y1 = sapply(c("sc", "pc"),
                                          print_stemi_events,
                                          output = "stroke_mc",
                                          digits=5),
                      major_bleed_y1 = sapply(c("sc", "pc"),
                                     print_stemi_events,
                                     output = "major_bleed_dt",
                                     digits=5),
                      minor_bleed_y1 = sapply(c("sc", "pc"),
                                              print_stemi_events,
                                              output = "minor_bleed_dt",
                                              digits=5),
                      death_y1 = sapply(c("sc", "pc"),
                                              print_stemi_events,
                                              output = "death_dt",
                                              digits=5),
                      death_post_y1 = sapply(c("sc", "pc"),
                                                   print_stemi_events,
                                                   output = "death_mc",
                                                   digits=5)
) %>% rbind(
  data.frame(group = rep("NSTEMI",
                         2),
             intervention = c("standard DAPT",
                              "genotype-guided DAPT"),
             mi_y1 = sapply(c("sc", "pc"),
                            print_nstemi_events,
                            output = "mi_dt",
                            digits=5),
             mi_post_y1 = sapply(c("sc", "pc"),
                                 print_nstemi_events,
                                 output = "mi_mc",
                                 digits=5),
             stroke_y1 = sapply(c("sc", "pc"),
                                print_nstemi_events,
                                output = "stroke_dt",
                                digits=5),
             stroke_post_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "stroke_mc",
                                     digits=5),
             major_bleed_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "major_bleed_dt",
                                     digits=5),
             minor_bleed_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "minor_bleed_dt",
                                     digits=5),
             death_y1 = sapply(c("sc", "pc"),
                               print_nstemi_events,
                               output = "death_dt",
                               digits=5),
             death_post_y1 = sapply(c("sc", "pc"),
                                    print_nstemi_events,
                                    output = "death_mc",
                                    digits=5)
  )
) %>% 
  mutate(grp_split = factor(group, levels = c("STEMI", "NSTEMI"))) %>%
  group_split(grp_split) %>% 
  map_dfr(~ add_row(.x, .after = Inf)) %>%
  select(-grp_split)

if (SAVE_OUTPUTS){
  write.csv(events_df,
            "formatted_outputs/base_case_event_counts.csv")
}

# Now load in scenario analysis results

stemi_sa_baseline <- read.csv("outputs/stemi_scenario_central_estimate.csv") %>%
  mutate(total_nmb_sc = stemi_size * nmb_sc,
         total_nmb_pc = stemi_size * nmb_pc)
stemi_sa_psa <- read.csv("outputs/stemi_n_10000_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = stemi_size * nmb_sc_mean,
         total_nmb_sc_L = stemi_size * nmb_sc_L,
         total_nmb_sc_U = stemi_size * nmb_sc_U,
         total_nmb_pc_mean = stemi_size * nmb_pc_mean,
         total_nmb_pc_L = stemi_size * nmb_pc_L,
         total_nmb_pc_U = stemi_size * nmb_pc_U)
stemi_scenario_prob_ce <- read.csv("outputs/stemi_n_10000_psa_acceptance_probability.csv")
nstemi_sa_baseline <- read.csv("outputs/nstemi_scenario_central_estimate.csv") %>%
  mutate(total_nmb_sc = stemi_size * nmb_sc,
         total_nmb_pc = stemi_size * nmb_pc)
nstemi_sa_psa <- read.csv("outputs/nstemi_n_10000_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = stemi_size * nmb_sc_mean,
         total_nmb_sc_L = stemi_size * nmb_sc_L,
         total_nmb_sc_U = stemi_size * nmb_sc_U,
         total_nmb_pc_mean = stemi_size * nmb_pc_mean,
         total_nmb_pc_L = stemi_size * nmb_pc_L,
         total_nmb_pc_U = stemi_size * nmb_pc_U)
nstemi_scenario_prob_ce <- read.csv("outputs/nstemi_n_10000_psa_acceptance_probability.csv")

n_stemi_scenario <- nrow(stemi_sa_baseline)
print_stemi_scenarios <- function(arm,
                               output,
                               digits=3){
  for (i in 1:n_stemi_scenario){
    L_diff <- stemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      stemi_sa_psa[1, paste(output, "_",  arm, "_L", sep="")]
    U_diff <- stemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      stemi_sa_psa[1, paste(output, "_",  arm, "_U", sep="")]
    if ((L_diff<0)|(U_diff>0)){
      print(paste("Central estimate out of bounds for",
                  output,
                  "in STEMI",
                  arm,
                  "arm."))
    }
  }
  sapply(1:n_stemi_scenario,
         FUN = function(i)
           paste(format(stemi_sa_baseline[i, paste(output, "_",  arm, sep="")],
                       digits = digits),
                " (",
                format(stemi_sa_psa[i, paste(output, "_",  arm, "_L", sep="")],
                       digits = digits),
                ", ",
                format(stemi_sa_psa[i, paste(output, "_",  arm, "_U", sep="")],
                       digits = digits),
                ")",
                sep="")
  )
    
}

n_nstemi_scenario <- nrow(nstemi_sa_baseline)
print_nstemi_scenarios <- function(arm,
                                  output,
                                  digits=3){
  for (i in 1:n_nstemi_scenario){
    L_diff <- nstemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      nstemi_sa_psa[1, paste(output, "_",  arm, "_L", sep="")]
    U_diff <- nstemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      nstemi_sa_psa[1, paste(output, "_",  arm, "_U", sep="")]
    if ((L_diff<0)|(U_diff>0)){
      print(paste("Central estimate out of bounds for",
                  output,
                  "in NSTEMI",
                  arm,
                  "arm."))
    }
  }
  sapply(1:n_nstemi_scenario,
         FUN = function(i)
           paste(format(nstemi_sa_baseline[i, paste(output, "_",  arm, sep="")],
                        digits = digits),
                 " (",
                 format(nstemi_sa_psa[i, paste(output, "_",  arm, "_L", sep="")],
                        digits = digits),
                 ", ",
                 format(nstemi_sa_psa[i, paste(output, "_",  arm, "_U", sep="")],
                        digits = digits),
                 ")",
                 sep="")
  )
  
}

all_outputs <- stemi_sa_psa %>%
  select(contains("_sc_mean")) %>%
  colnames() %>%
  str_replace("_sc_mean", "") %>%
  unique()

stemi_scenario_df <- data.frame(group = rep("STEMI",
                                    2 * n_stemi_scenario),
                                scenario = rep(1:n_stemi_scenario, 2),
                        intervention = c(rep("standard DAPT", n_stemi_scenario),
                                         rep("genotype-guided DAPT", n_stemi_scenario)))
for (output in all_outputs){
  stemi_scenario_df <- stemi_scenario_df %>%
    mutate("{output}" := sapply(c("sc", "pc"),
                                print_stemi_scenarios,
                                output = output,
                                digits=5) %>% c())
}
  
stemi_scenario_df <- stemi_scenario_df %>%
  relocate(cost, .after = cost_udc) %>%
  mutate(inc_cost = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(stemi_sa_baseline$inc_cost[scenario],
                                          digits = 5),
                                   " (",
                                   format(stemi_sa_psa$cost_inc_L[scenario],
                                          digits = 5),
                                   ", ",
                                   format(stemi_sa_psa$cost_inc_U[scenario],
                                          digits = 5),
                                   ")",
                                   sep="")
  )) %>%
  mutate(inc_util = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(stemi_sa_baseline$inc_util[scenario],
                                        digits = 5),
                                 " (",
                                 format(stemi_sa_psa$util_inc_L[scenario],
                                        digits = 5),
                                 ", ",
                                 format(stemi_sa_psa$util_inc_U[scenario],
                                        digits = 5),
                                 ")",
                                 sep="")
  )) %>%
  mutate(icer = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(stemi_sa_baseline$icer[scenario],
                                       digits = 5),
                                " (",
                                format(stemi_sa_psa$icer_L[scenario],
                                       digits = 5),
                                ", ",
                                format(stemi_sa_psa$icer_U[scenario],
                                       digits = 5),
                                ")",
                                sep="")
  )) %>%
  mutate(inc_nmb = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(stemi_sa_baseline$inc_nmb[scenario],
                                        digits = 5),
                                 " (",
                                 format(stemi_sa_psa$nmb_inc_L[scenario],
                                        digits = 5),
                                 ", ",
                                 format(stemi_sa_psa$nmb_inc_U[scenario],
                                        digits = 5),
                                 ")",
                                 sep="")
  )) %>%
  mutate(prob_ce = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(stemi_scenario_prob_ce$ce_prob[scenario],
                                       digits = 5),
                                " (",
                                format(stemi_scenario_prob_ce$lower_95[scenario],
                                       digits = 5),
                                ", ",
                                format(stemi_scenario_prob_ce$upper_95[scenario],
                                       digits = 5),
                                ")",
                                sep="")
  )) %>%
  relocate(nmb, .before = inc_nmb) %>%
  relocate(total_nmb, .before = inc_nmb) %>%
  arrange(scenario) %>% 
  group_split(scenario) %>% 
  map_dfr(~ add_row(.x, .after = Inf))

if (SAVE_OUTPUTS){
  write.csv(stemi_scenario_df,
            "formatted_outputs/stemi_scenario_psa.csv")
}

nstemi_scenario_df <- data.frame(group = rep("NSTEMI",
                                            2 * n_nstemi_scenario),
                                scenario = rep(1:n_nstemi_scenario, 2),
                                intervention = c(rep("standard DAPT", n_nstemi_scenario),
                                                 rep("genotype-guided DAPT", n_nstemi_scenario)))
for (output in all_outputs){
  nstemi_scenario_df <- nstemi_scenario_df %>%
    mutate("{output}" := sapply(c("sc", "pc"),
                                print_nstemi_scenarios,
                                output = output,
                                digits=5) %>% c())
}

nstemi_scenario_df <- nstemi_scenario_df %>%
  relocate(cost, .after = cost_udc) %>%
  mutate(inc_cost = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(nstemi_sa_baseline$inc_cost[scenario],
                                        digits = 5),
                                 " (",
                                 format(nstemi_sa_psa$cost_inc_L[scenario],
                                        digits = 5),
                                 ", ",
                                 format(nstemi_sa_psa$cost_inc_U[scenario],
                                        digits = 5),
                                 ")",
                                 sep="")
  )) %>%
  mutate(inc_util = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(nstemi_sa_baseline$inc_util[scenario],
                                        digits = 5),
                                 " (",
                                 format(nstemi_sa_psa$util_inc_L[scenario],
                                        digits = 5),
                                 ", ",
                                 format(nstemi_sa_psa$util_inc_U[scenario],
                                        digits = 5),
                                 ")",
                                 sep="")
  )) %>%
  mutate(icer = ifelse(intervention == "standard DAPT",
                       NA,
                       paste(format(nstemi_sa_baseline$icer[scenario],
                                    digits = 5),
                             " (",
                             format(nstemi_sa_psa$icer_L[scenario],
                                    digits = 5),
                             ", ",
                             format(nstemi_sa_psa$icer_U[scenario],
                                    digits = 5),
                             ")",
                             sep="")
  )) %>%
  mutate(inc_nmb = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(nstemi_sa_baseline$inc_nmb[scenario],
                                       digits = 5),
                                " (",
                                format(nstemi_sa_psa$nmb_inc_L[scenario],
                                       digits = 5),
                                ", ",
                                format(nstemi_sa_psa$nmb_inc_U[scenario],
                                       digits = 5),
                                ")",
                                sep="")
  )) %>%
  mutate(prob_ce = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(nstemi_scenario_prob_ce$ce_prob[scenario],
                                       digits = 5),
                                " (",
                                format(nstemi_scenario_prob_ce$lower_95[scenario],
                                       digits = 5),
                                ", ",
                                format(nstemi_scenario_prob_ce$upper_95[scenario],
                                       digits = 5),
                                ")",
                                sep="")
  )) %>%
  relocate(nmb, .before = inc_nmb) %>%
  relocate(total_nmb, .before = inc_nmb) %>%
  arrange(scenario) %>% 
  group_split(scenario) %>% 
  map_dfr(~ add_row(.x, .after = Inf))

if (SAVE_OUTPUTS){
  write.csv(nstemi_scenario_df,
            "formatted_outputs/nstemi_scenario_psa.csv")
}