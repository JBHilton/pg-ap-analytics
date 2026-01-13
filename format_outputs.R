# In this script we load and reformat the different bits of output data into
# tables of base cases and confidence intervals.

SAVE_OUTPUTS <- TRUE
dir.create("formatted_outputs",
           showWarnings = FALSE)

library(stringr)
library(tidyverse)

stemi_arm_comparison <- read.csv("outputs/stemi_arm_comparison.csv") %>%
  mutate(output_name = str_replace_all(output_name, "ratio_udc", "icer_udc"))
stemi_psa <- read.csv("outputs/stemi_n_1e+05_psa_stats.csv")
stemi_psa_samples <- read.csv("outputs/stemi_n_1e+05_psa_samples.csv")
stemi_prob_ce <- read.csv("outputs/stemi_n_1e+05_acceptance_probability.csv")
nstemi_arm_comparison <- read.csv("outputs/nstemi_arm_comparison.csv") %>%
  mutate(output_name = str_replace_all(output_name, "ratio_udc", "icer_udc"))
nstemi_psa <- read.csv("outputs/nstemi_n_1e+05_psa_stats.csv")
nstemi_psa_samples <- read.csv("outputs/nstemi_n_1e+05_psa_samples.csv")
nstemi_prob_ce <- read.csv("outputs/nstemi_n_1e+05_acceptance_probability.csv")

# Extra stuff for checking
stemi_ac_short <- stemi_arm_comparison %>% filter(output_name %in% stemi_psa$output_name) %>% arrange(output_name)
stemi_psa <- stemi_psa %>% arrange(output_name)
nstemi_ac_short <- nstemi_arm_comparison %>% filter(output_name %in% nstemi_psa$output_name) %>% arrange(output_name)
nstemi_psa <- nstemi_psa %>% arrange(output_name)

# Get deterministic undiscounted NMBs
threshold <- 20000
stemi_nmb_udc_sc <- threshold * (
  stemi_arm_comparison[which(stemi_arm_comparison$output_name=="util_udc"), 2:4]) -
  stemi_arm_comparison[which(stemi_arm_comparison$output_name=="cost_udc"), 2:4]
stemi_arm_comparison <- stemi_arm_comparison %>%
  add_row(output_name = "nmb_udc",
          sc = stemi_nmb_udc_sc$sc,
          pc = stemi_nmb_udc_sc$pc,
          inc = stemi_nmb_udc_sc$inc)

nstemi_nmb_udc_sc <- threshold * (
  nstemi_arm_comparison[which(nstemi_arm_comparison$output_name=="util_udc"), 2:4]) -
  nstemi_arm_comparison[which(nstemi_arm_comparison$output_name=="cost_udc"), 2:4]
nstemi_arm_comparison <- nstemi_arm_comparison %>%
  add_row(output_name = "nmb_udc",
          sc = nstemi_nmb_udc_sc$sc,
          pc = nstemi_nmb_udc_sc$pc,
          inc = nstemi_nmb_udc_sc$inc)

# Add undiscounted ICER and NMB to PSA samples
stemi_psa_samples <- stemi_psa_samples %>%
  mutate(icer_udc = inc_cost_udc_hs / inc_util_udc_hs,
         nmb_sc_udc = threshold * util_udc_sc_hs - cost_udc_sc_hs,
         nmb_pc_udc = threshold * util_udc_pc_hs - cost_udc_pc_hs) %>%
  mutate(inc_nmb_udc = nmb_pc_udc - nmb_sc_udc)

nstemi_psa_samples <- nstemi_psa_samples %>%
  mutate(icer_udc = inc_cost_udc_hs / inc_util_udc_hs,
         nmb_sc_udc = threshold * util_udc_sc_hs - cost_udc_sc_hs,
         nmb_pc_udc = threshold * util_udc_pc_hs - cost_udc_pc_hs) %>%
  mutate(inc_nmb_udc = nmb_pc_udc - nmb_sc_udc)

# Add these to stats
stemi_psa <- stemi_psa %>%
  add_row(output_name = "icer_udc",
         sc_mean = NA,
         sc_L = NA,
         sc_U = NA,
         pc_mean = NA,
         pc_L = NA,
         pc_U = NA,
         inc_mean = mean(stemi_psa_samples$icer_udc),
         inc_L = quantile((stemi_psa_samples$icer_udc), 0.025),
         inc_U = quantile((stemi_psa_samples$icer_udc), 0.975)) %>%
  add_row(output_name = "nmb_udc",
          sc_mean = mean(stemi_psa_samples$nmb_sc_udc),
          sc_L = quantile((stemi_psa_samples$nmb_sc_udc), 0.025),
          sc_U = quantile((stemi_psa_samples$nmb_sc_udc), 0.975),
          pc_mean = mean(stemi_psa_samples$nmb_pc_udc),
          pc_L = quantile((stemi_psa_samples$nmb_pc_udc), 0.025),
          pc_U = quantile((stemi_psa_samples$nmb_pc_udc), 0.975),
          inc_mean = mean(stemi_psa_samples$inc_nmb_udc),
          inc_L = quantile((stemi_psa_samples$inc_nmb_udc), 0.025),
          inc_U = quantile((stemi_psa_samples$inc_nmb_udc), 0.975))

nstemi_psa <- nstemi_psa %>%
  add_row(output_name = "icer_udc",
          sc_mean = NA,
          sc_L = NA,
          sc_U = NA,
          pc_mean = NA,
          pc_L = NA,
          pc_U = NA,
          inc_mean = mean(nstemi_psa_samples$icer_udc),
          inc_L = quantile((nstemi_psa_samples$icer_udc), 0.025),
          inc_U = quantile((nstemi_psa_samples$icer_udc), 0.975)) %>%
  add_row(output_name = "nmb_udc",
          sc_mean = mean(nstemi_psa_samples$nmb_sc_udc),
          sc_L = quantile((nstemi_psa_samples$nmb_sc_udc), 0.025),
          sc_U = quantile((nstemi_psa_samples$nmb_sc_udc), 0.975),
          pc_mean = mean(nstemi_psa_samples$nmb_pc_udc),
          pc_L = quantile((nstemi_psa_samples$nmb_pc_udc), 0.025),
          pc_U = quantile((nstemi_psa_samples$nmb_pc_udc), 0.975),
          inc_mean = mean(nstemi_psa_samples$inc_nmb_udc),
          inc_L = quantile((nstemi_psa_samples$inc_nmb_udc), 0.025),
          inc_U = quantile((nstemi_psa_samples$inc_nmb_udc), 0.975))

# Population sizes for total NMB:
stemi_size = 41690
nstemi_size = 61248

stemi_nmb_df <- stemi_arm_comparison %>% filter(grepl("nmb", output_name)) %>%
  mutate(output_name = paste("total_", output_name, sep = ""),
         sc = stemi_size * sc,
         pc = stemi_size * pc,
         inc= stemi_size * inc)
stemi_arm_comparison <- stemi_arm_comparison %>%
  rbind(stemi_nmb_df)
stemi_psa_nmb <- stemi_psa %>%
  filter(grepl("nmb", output_name)) %>%
  select(-output_name)
stemi_psa <- stemi_psa %>%
  rbind(data.frame(output_name=c("total_nmb",
                                 "total_nmb_udc"),
                   stemi_size * stemi_psa_nmb))

nstemi_nmb_df <- nstemi_arm_comparison %>% filter(grepl("nmb", output_name)) %>%
  mutate(output_name = paste("total_", output_name, sep = ""),
         sc = nstemi_size * sc,
         pc = nstemi_size * pc,
         inc= nstemi_size * inc)
nstemi_arm_comparison <- nstemi_arm_comparison %>%
  rbind(nstemi_nmb_df)
nstemi_psa_nmb <- nstemi_psa %>%
  filter(grepl("nmb", output_name)) %>%
  select(-output_name)
nstemi_psa <- nstemi_psa %>%
  rbind(data.frame(output_name=c("total_nmb",
                                 "total_nmb_udc"),
                   nstemi_size * nstemi_psa_nmb))

print_stemi_det_output <- function(arm,
                                   output,
                                   decimals = 2){

  paste(format(round(stemi_arm_comparison[which(stemi_arm_comparison$output_name==output),
                                          arm], decimals),
               nsmall = decimals))
}

print_nstemi_det_output <- function(arm,
                                    output,
                                    decimals = 2){
  
  paste(format(round(nstemi_arm_comparison[which(nstemi_arm_comparison$output_name==output),
                                          arm], decimals),
               nsmall = decimals))
}

print_stemi_psa_output <- function(arm,
                         output,
                         decimals = 2){
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
  paste(format(round(stemi_psa[which(stemi_psa$output_name==output),
                               paste(arm,
                                     "_mean",
                                     sep="")], decimals),
               nsmall = decimals),
        " (",
        format(round(stemi_psa[which(stemi_psa$output_name==output),
                         paste(arm,
                               "_L",
                               sep="")], decimals),
               nsmall = decimals),
        ", ",
        format(round(stemi_psa[which(stemi_psa$output_name==output),
                         paste(arm,
                               "_U",
                               sep="")], decimals),
               nsmall = decimals),
        ")",
        sep="")
}

print_nstemi_psa_output <- function(arm,
                               output,
                               decimals = 2){
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
  paste(format(round(nstemi_psa[which(nstemi_psa$output_name==output),
                                paste(arm,
                                      "_mean",
                                      sep="")], decimals),
               nsmall = decimals),
        " (",
        format(round(nstemi_psa[which(nstemi_psa$output_name==output),
                         paste(arm,
                               "_L",
                               sep="")], decimals),
               nsmall = decimals),
        ", ",
        format(round(nstemi_psa[which(nstemi_psa$output_name==output),
                         paste(arm,
                               "_U",
                               sep="")], decimals),
               nsmall = decimals),
        ")",
        sep="")
}

print_ce_prob <- function(ce_df,
                          decimals){
  idx <- which.min(abs(ce_df$ce_threshold - 2*1e4))
  paste(format(round(100 * ce_df$acceptance_prob[idx],
                     decimals),
               nsmall = decimals),
        ", (",
        format(round(100 * ce_df$lower_95[idx],
                     decimals),
               nsmall = decimals),
        ", ",
        format(round(100 * ce_df$upper_95[idx],
                     decimals),
               nsmall = decimals),
        ")",
        sep="")
}

det_df <- data.frame(group = rep("STEMI",
                                  2),
                      intervention = c("standard DAPT",
                                       "genotype-guided DAPT"),
                      life_years = sapply(c("sc", "pc"),
                                          print_stemi_det_output,
                                          output = "life_years",
                                          decimals = 2),
                      costs_udc = sapply(c("sc", "pc"),
                                         print_stemi_det_output,
                                         output = "cost_udc",
                                         decimals = 0),
                      costs_dc = sapply(c("sc", "pc"),
                                        print_stemi_det_output,
                                        output = "cost",
                                        decimals = 0),
                      qalys_udc = sapply(c("sc", "pc"),
                                         print_stemi_det_output,
                                         output = "util_udc",
                                         decimals = 4),
                      qalys_dc = sapply(c("sc", "pc"),
                                        print_stemi_det_output,
                                        output = "util",
                                        decimals = 4),
                      inc_costs_dc = c(0,
                                    print_stemi_det_output("inc",
                                                           "cost",
                                                           decimals = 0)),
                      inc_costs_udc = c(0,
                                    print_stemi_det_output("inc",
                                                           "cost_udc",
                                                           decimals = 0)),
                      inc_qalys_dc = c(0,
                                    print_stemi_det_output("inc",
                                                           "util",
                                                           decimals = 4)),
                      inc_qalys_udc = c(0,
                                    print_stemi_det_output("inc",
                                                           "util_udc",
                                                           decimals = 4)),
                      icer = c(0,
                               print_stemi_det_output("inc",
                                                      "icer",
                                                      decimals = 0)),
                     icer_udc = c(0,
                              print_stemi_det_output("inc",
                                                     "icer_udc",
                                                     decimals = 0)),
                      total_nmb = sapply(c("sc", "pc"),
                                         print_stemi_det_output,
                                         output = "total_nmb",
                                         decimals = 0),
                     total_nmb_udc = sapply(c("sc", "pc"),
                                        print_stemi_det_output,
                                        output = "total_nmb_udc",
                                        decimals = 0),
                      inc_nmb = c(0,
                                  print_stemi_det_output("inc",
                                                         "total_nmb",
                                                         decimals = 0)),
                     inc_nmb_udc = c(0,
                                 print_stemi_det_output("inc",
                                                        "total_nmb_udc",
                                                        decimals = 0))
) %>%
  rbind(
    data.frame(group = rep("NSTEMI",
                           2),
               intervention = c("standard DAPT",
                                "genotype-guided DAPT"),
               life_years = sapply(c("sc", "pc"),
                                   print_nstemi_det_output,
                                   output = "life_years",
                                   decimals = 2),
               costs_udc = sapply(c("sc", "pc"),
                                  print_nstemi_det_output,
                                  output = "cost_udc",
                                  decimals = 0),
               costs_dc = sapply(c("sc", "pc"),
                                 print_nstemi_det_output,
                                 output = "cost",
                                 decimals = 0),
               qalys_udc = sapply(c("sc", "pc"),
                                  print_nstemi_det_output,
                                  output = "util_udc",
                                  decimals = 4),
               qalys_dc = sapply(c("sc", "pc"),
                                 print_nstemi_det_output,
                                 output = "util",
                                 decimals = 4),
               inc_costs_dc = c(0,
                                print_nstemi_det_output("inc",
                                                       "cost",
                                                       decimals = 0)),
               inc_costs_udc = c(0,
                                 print_nstemi_det_output("inc",
                                                        "cost_udc",
                                                        decimals = 0)),
               inc_qalys_dc = c(0,
                                print_nstemi_det_output("inc",
                                                       "util",
                                                       decimals = 4)),
               inc_qalys_udc = c(0,
                                 print_nstemi_det_output("inc",
                                                        "util_udc",
                                                        decimals = 4)),
               icer = c(0,
                        print_nstemi_det_output("inc",
                                                "icer",
                                                decimals = 0)),
               icer_udc = c(0,
                            print_nstemi_det_output("inc",
                                                   "icer_udc",
                                                   decimals = 0)),
               total_nmb = sapply(c("sc", "pc"),
                                  print_nstemi_det_output,
                                  output = "total_nmb",
                                  decimals = 0),
               total_nmb_udc = sapply(c("sc", "pc"),
                                      print_nstemi_det_output,
                                      output = "total_nmb_udc",
                                      decimals = 0),
               inc_nmb = c(0,
                           print_nstemi_det_output("inc",
                                                  "total_nmb",
                                                  decimals = 0)),
               inc_nmb_udc = c(0,
                               print_nstemi_det_output("inc",
                                                      "total_nmb_udc",
                                                      decimals = 0))
    )
  ) %>% 
  mutate(grp_split = factor(group, levels = c("STEMI", "NSTEMI"))) %>%
  group_split(grp_split) %>% 
  map_dfr(~ add_row(.x, .after = Inf)) %>%
  select(-grp_split)

main_df <- data.frame(group = rep("STEMI",
                                  2),
                      intervention = c("standard DAPT",
                                       "genotype-guided DAPT"),
                      life_years = sapply(c("sc", "pc"),
                                          print_stemi_psa_output,
                                          output = "life_years",
                                          decimals = 2),
                      costs_udc = sapply(c("sc", "pc"),
                                         print_stemi_psa_output,
                                         output = "cost_udc",
                                         decimals = 0),
                      costs_dc = sapply(c("sc", "pc"),
                                        print_stemi_psa_output,
                                        output = "cost",
                                        decimals = 0),
                      qalys_udc = sapply(c("sc", "pc"),
                                         print_stemi_psa_output,
                                         output = "util_udc",
                                         decimals = 4),
                      qalys_dc = sapply(c("sc", "pc"),
                                        print_stemi_psa_output,
                                        output = "util",
                                        decimals = 4),
                      inc_costs_dc = c(0,
                                       print_stemi_psa_output("inc",
                                                              "cost",
                                                              decimals = 0)),
                      inc_costs_udc = c(0,
                                        print_stemi_psa_output("inc",
                                                               "cost_udc",
                                                               decimals = 0)),
                      inc_qalys_dc = c(0,
                                       print_stemi_psa_output("inc",
                                                              "util",
                                                              decimals = 4)),
                      inc_qalys_udc = c(0,
                                        print_stemi_psa_output("inc",
                                                               "util_udc",
                                                               decimals = 4)),
                      icer = c(0,
                               print_stemi_psa_output("inc",
                                                      "icer",
                                                      decimals = 0)),
                      icer_udc = c(0,
                                   print_stemi_psa_output("inc",
                                                          "icer_udc",
                                                          decimals = 0)),
                      total_nmb = sapply(c("sc", "pc"),
                                         print_stemi_psa_output,
                                         output = "total_nmb",
                                         decimals = 0),
                      total_nmb_udc = sapply(c("sc", "pc"),
                                             print_stemi_psa_output,
                                             output = "total_nmb_udc",
                                             decimals = 0),
                      inc_nmb = c(0,
                                  print_stemi_psa_output("inc",
                                                         "total_nmb",
                                                         decimals = 0)),
                      inc_nmb_udc = c(0,
                                      print_stemi_psa_output("inc",
                                                             "total_nmb_udc",
                                                             decimals = 0))
) %>%
  rbind(
    data.frame(group = rep("NSTEMI",
                           2),
               intervention = c("standard DAPT",
                                "genotype-guided DAPT"),
               life_years = sapply(c("sc", "pc"),
                                   print_nstemi_psa_output,
                                   output = "life_years",
                                   decimals = 2),
               costs_udc = sapply(c("sc", "pc"),
                                  print_nstemi_psa_output,
                                  output = "cost_udc",
                                  decimals = 0),
               costs_dc = sapply(c("sc", "pc"),
                                 print_nstemi_psa_output,
                                 output = "cost",
                                 decimals = 0),
               qalys_udc = sapply(c("sc", "pc"),
                                  print_nstemi_psa_output,
                                  output = "util_udc",
                                  decimals = 4),
               qalys_dc = sapply(c("sc", "pc"),
                                 print_nstemi_psa_output,
                                 output = "util",
                                 decimals = 4),
               inc_costs_dc = c(0,
                                print_nstemi_psa_output("inc",
                                                        "cost",
                                                        decimals = 0)),
               inc_costs_udc = c(0,
                                 print_nstemi_psa_output("inc",
                                                         "cost_udc",
                                                         decimals = 0)),
               inc_qalys_dc = c(0,
                                print_nstemi_psa_output("inc",
                                                        "util",
                                                        decimals = 4)),
               inc_qalys_udc = c(0,
                                 print_nstemi_psa_output("inc",
                                                         "util_udc",
                                                         decimals = 4)),
               icer = c(0,
                        print_nstemi_psa_output("inc",
                                                "icer",
                                                decimals = 0)),
               icer_udc = c(0,
                            print_nstemi_psa_output("inc",
                                                    "icer_udc",
                                                    decimals = 0)),
               total_nmb = sapply(c("sc", "pc"),
                                  print_nstemi_psa_output,
                                  output = "total_nmb",
                                  decimals = 0),
               total_nmb_udc = sapply(c("sc", "pc"),
                                      print_nstemi_psa_output,
                                      output = "total_nmb_udc",
                                      decimals = 0),
               inc_nmb = c(0,
                           print_nstemi_psa_output("inc",
                                                   "total_nmb",
                                                   decimals = 0)),
               inc_nmb_udc = c(0,
                               print_nstemi_psa_output("inc",
                                                       "total_nmb_udc",
                                                       decimals = 0))
    )
  ) %>% 
  mutate(grp_split = factor(group, levels = c("STEMI", "NSTEMI"))) %>%
  group_split(grp_split) %>% 
  map_dfr(~ add_row(.x, .after = Inf)) %>%
  select(-grp_split)

det_df_dc <- det_df %>%
  select(-c(grep("_udc", colnames(main_df))))

det_df_udc <- det_df %>%
  select(c("group",
           "intervention",
           grep("_udc", colnames(main_df))))

psa_df_dc <- main_df %>%
  select(-c(grep("_udc", colnames(main_df))))

psa_df_udc <- main_df %>%
  select(c("group",
           "intervention",
           grep("_udc", colnames(main_df))))

if (SAVE_OUTPUTS){
  write.csv(det_df,
            "formatted_outputs/deterministic_outcomes.csv")
  write.csv(det_df_dc,
            "formatted_outputs/deterministic_with_discounting.csv")
  write.csv(det_df_udc,
            "formatted_outputs/deterministic_undiscounted.csv")
  write.csv(main_df,
            "formatted_outputs/full_psa_outcomes.csv")
  write.csv(psa_df_dc,
            "formatted_outputs/PSA_with_discounting.csv")
  write.csv(psa_df_udc,
            "formatted_outputs/PSA_undiscounted.csv")
}
# Now load in event count results

stemi_event_counts <- read.csv("outputs/stemi_event_counts.csv")
stemi_event_psa <- read.csv("outputs/stemi_n_1e+05_event_stats.csv")
nstemi_event_counts <- read.csv("outputs/nstemi_event_counts.csv")
nstemi_event_psa <- read.csv("outputs/nstemi_n_1e+05_event_stats.csv")

arm_pos <- function(arm) ifelse(arm=="sc", 1, 2)

print_stemi_events <- function(arm,
                               output,
                               decimals = 2){
  paste(format(round(stemi_event_counts[arm_pos(arm), output],
                     decimals),
               nsmall = decimals),
        " (",
        format(round(stemi_event_psa[1, paste(output, "_",  arm, "_L", sep="")],
               decimals),
               nsmall = decimals),
        ", ",
        format(round(stemi_event_psa[1, paste(output, "_",  arm, "_U", sep="")],
               decimals),
               nsmall = decimals),
        ")",
        sep="")
}

print_nstemi_events <- function(arm,
                               output,
                               decimals = 2){
  paste(format(round(nstemi_event_counts[arm_pos(arm), output],
                     decimals),
               nsmall = decimals),
        " (",
        format(round(nstemi_event_psa[1, paste(output, "_",  arm, "_L", sep="")],
                     decimals),
               nsmall = decimals),
        ", ",
        format(round(nstemi_event_psa[1, paste(output, "_",  arm, "_U", sep="")],
                     decimals),
               nsmall = decimals),
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
                                     decimals = 1),
                      mi_post_y1 = sapply(c("sc", "pc"),
                                          print_stemi_events,
                                          output = "mi_mc",
                                          decimals = 1),
                      stroke_y1 = sapply(c("sc", "pc"),
                                     print_stemi_events,
                                     output = "stroke_dt",
                                     decimals = 1),
                      stroke_post_y1 = sapply(c("sc", "pc"),
                                          print_stemi_events,
                                          output = "stroke_mc",
                                          decimals = 1),
                      major_bleed_y1 = sapply(c("sc", "pc"),
                                     print_stemi_events,
                                     output = "major_bleed_dt",
                                     decimals = 1),
                      minor_bleed_y1 = sapply(c("sc", "pc"),
                                              print_stemi_events,
                                              output = "minor_bleed_dt",
                                              decimals = 1),
                      death_y1 = sapply(c("sc", "pc"),
                                              print_stemi_events,
                                              output = "death_dt",
                                              decimals = 1),
                      death_post_y1 = sapply(c("sc", "pc"),
                                                   print_stemi_events,
                                                   output = "death_mc",
                                                   decimals = 1)
) %>% rbind(
  data.frame(group = rep("NSTEMI",
                         2),
             intervention = c("standard DAPT",
                              "genotype-guided DAPT"),
             mi_y1 = sapply(c("sc", "pc"),
                            print_nstemi_events,
                            output = "mi_dt",
                            decimals = 1),
             mi_post_y1 = sapply(c("sc", "pc"),
                                 print_nstemi_events,
                                 output = "mi_mc",
                                 decimals = 1),
             stroke_y1 = sapply(c("sc", "pc"),
                                print_nstemi_events,
                                output = "stroke_dt",
                                decimals = 1),
             stroke_post_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "stroke_mc",
                                     decimals = 1),
             major_bleed_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "major_bleed_dt",
                                     decimals = ),
             minor_bleed_y1 = sapply(c("sc", "pc"),
                                     print_nstemi_events,
                                     output = "minor_bleed_dt",
                                     decimals = ),
             death_y1 = sapply(c("sc", "pc"),
                               print_nstemi_events,
                               output = "death_dt",
                               decimals = 1),
             death_post_y1 = sapply(c("sc", "pc"),
                                    print_nstemi_events,
                                    output = "death_mc",
                                    decimals = 1)
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
stemi_sa_psa <- read.csv("outputs/stemi_n_1e+05_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = stemi_size * nmb_sc_mean,
         total_nmb_sc_L = stemi_size * nmb_sc_L,
         total_nmb_sc_U = stemi_size * nmb_sc_U,
         total_nmb_pc_mean = stemi_size * nmb_pc_mean,
         total_nmb_pc_L = stemi_size * nmb_pc_L,
         total_nmb_pc_U = stemi_size * nmb_pc_U,
         total_nmb_inc_mean = stemi_size * nmb_inc_mean,
         total_nmb_inc_L = stemi_size * nmb_inc_L,
         total_nmb_inc_U = stemi_size * nmb_inc_U)
stemi_scenario_prob_ce <- read.csv("outputs/stemi_n_1e+05_psa_acceptance_probability.csv")
nstemi_sa_baseline <- read.csv("outputs/nstemi_scenario_central_estimate.csv") %>%
  mutate(total_nmb_sc = stemi_size * nmb_sc,
         total_nmb_pc = stemi_size * nmb_pc)
nstemi_sa_psa <- read.csv("outputs/nstemi_n_1e+05_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = stemi_size * nmb_sc_mean,
         total_nmb_sc_L = stemi_size * nmb_sc_L,
         total_nmb_sc_U = stemi_size * nmb_sc_U,
         total_nmb_pc_mean = stemi_size * nmb_pc_mean,
         total_nmb_pc_L = stemi_size * nmb_pc_L,
         total_nmb_pc_U = stemi_size * nmb_pc_U,
         total_nmb_inc_mean = stemi_size * nmb_inc_mean,
         total_nmb_inc_L = stemi_size * nmb_inc_L,
         total_nmb_inc_U = stemi_size * nmb_inc_U)
nstemi_scenario_prob_ce <- read.csv("outputs/nstemi_n_1e+05_psa_acceptance_probability.csv")

n_stemi_scenario <- nrow(stemi_sa_baseline)
print_stemi_scenarios <- function(arm,
                               output,
                               decimals = 2){
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
           paste(format(round(stemi_sa_psa[i, paste(output, "_",  arm, "_mean", sep="")],
                              decimals),
                       nsmall = decimals),
                " (",
                format(round(stemi_sa_psa[i, paste(output, "_",  arm, "_L", sep="")],
                       decimals),
                       nsmall = decimals),
                ", ",
                format(round(stemi_sa_psa[i, paste(output, "_",  arm, "_U", sep="")],
                       decimals),
                       nsmall = decimals),
                ")",
                sep="")
  )
    
}

n_nstemi_scenario <- nrow(nstemi_sa_baseline)
print_nstemi_scenarios <- function(arm,
                                  output,
                                  decimals = 2){
  for (i in 1:n_nstemi_scenario){
    L_diff <- nstemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      nstemi_sa_psa[1, paste(output, "_",  arm, "_L", sep="")]
    U_diff <- nstemi_sa_baseline[1, paste(output, "_",  arm, sep="")] -
      nstemi_sa_psa[1, paste(output, "_",  arm, "_U", sep="")]
    if ((L_diff<0)|(U_diff>0)){
      print(paste("Central estimate out of bounds for",
                  output,
                  "in STEMI",
                  arm,
                  "arm."))
    }
  }
  sapply(1:n_nstemi_scenario,
         FUN = function(i)
           paste(format(round(nstemi_sa_psa[i, paste(output, "_",  arm, "_mean", sep="")],
                              decimals),
                        nsmall = decimals),
                 " (",
                 format(round(nstemi_sa_psa[i, paste(output, "_",  arm, "_L", sep="")],
                              decimals),
                        nsmall = decimals),
                 ", ",
                 format(round(nstemi_sa_psa[i, paste(output, "_",  arm, "_U", sep="")],
                              decimals),
                        nsmall = decimals),
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
  decimals <- ifelse(grepl("cost|icer|nmb", output),
                     0,
                     no = ifelse(grepl("qalys", output),
                                 yes = 4,
                                 no = 2))
  stemi_scenario_df[[output]] <- sapply(c("sc", "pc"),
                                print_stemi_scenarios,
                                output = output,
                                decimals = decimals) %>% c()
}
  
stemi_scenario_df <- stemi_scenario_df %>%
  relocate(cost, .after = cost_udc) %>%
  mutate(inc_cost = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(stemi_sa_psa$cost_inc_mean[scenario],
                                              0),
                                          nsmall = 0),
                                   " (",
                                   format(round(stemi_sa_psa$cost_inc_L[scenario],
                                                0),
                                          nsmall = 0),
                                   ", ",
                                   format(round(stemi_sa_psa$cost_inc_U[scenario],
                                                0),
                                          nsmall = 0),
                                   ")",
                                   sep="")
  )) %>%
  mutate(inc_cost_udc = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(stemi_sa_psa$cost_udc_inc_mean[scenario],
                                              0),
                                        nsmall = 0),
                                 " (",
                                 format(round(stemi_sa_psa$cost_udc_inc_L[scenario],
                                              0),
                                        nsmall = 0),
                                 ", ",
                                 format(round(stemi_sa_psa$cost_udc_inc_U[scenario],
                                              0),
                                        nsmall = 0),
                                 ")",
                                 sep="")
  )) %>%
  mutate(inc_util = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(stemi_sa_psa$util_inc_mean[scenario],
                                              3),
                                        nsmall = 3),
                                 " (",
                                 format(round(stemi_sa_psa$util_inc_L[scenario],
                                              3),
                                        nsmall = 3),
                                 ", ",
                                 format(round(stemi_sa_psa$util_inc_U[scenario],
                                              3),
                                        nsmall = 3),
                                 ")",
                                 sep="")
  )) %>%
  mutate(inc_util_udc = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(stemi_sa_psa$util_udc_inc_mean[scenario],
                                              3),
                                        nsmall = 3),
                                 " (",
                                 format(round(stemi_sa_psa$util_udc_inc_L[scenario],
                                              3),
                                        nsmall = 3),
                                 ", ",
                                 format(round(stemi_sa_psa$util_udc_inc_U[scenario],
                                              3),
                                        nsmall = 3),
                                 ")",
                                 sep="")
  )) %>%
  mutate(icer = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(round(stemi_sa_psa$icer_mean[scenario],
                                             0),
                                       nsmall = 0),
                                " (",
                                format(round(stemi_sa_psa$icer_L[scenario],
                                             0),
                                       nsmall = 0),
                                ", ",
                                format(round(stemi_sa_psa$icer_U[scenario],
                                             0),
                                       nsmall = 0),
                                ")",
                                sep="")
  )) %>%
  mutate(inc_nmb = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(stemi_sa_psa$total_nmb_inc_mean[scenario],
                                              0),
                                        nsmall = 0),
                                 " (",
                                 format(round(stemi_sa_psa$total_nmb_inc_L[scenario],
                                              0),
                                        nsmall = 0),
                                 ", ",
                                 format(round(stemi_sa_psa$total_nmb_inc_U[scenario],
                                              0),
                                        nsmall = 0),
                                 ")",
                                 sep="")
  )) %>%
  mutate(prob_ce = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(round(100 * stemi_scenario_prob_ce$ce_prob[scenario],
                                             1),
                                       nsmall = 0),
                                " (",
                                format(round(100 * stemi_scenario_prob_ce$lower_95[scenario],
                                             1),
                                       nsmall = 0),
                                ", ",
                                format(round(100 * stemi_scenario_prob_ce$upper_95[scenario],
                                             1),
                                       nsmall = 0),
                                ")",
                                sep="")
  )) %>%
  relocate(nmb, .before = inc_nmb) %>%
  relocate(total_nmb, .before = inc_nmb) %>%
  arrange(scenario) %>% 
  group_split(scenario) %>% 
  map_dfr(~ add_row(.x, .after = Inf))

scenario_df_dc <- stemi_scenario_df %>%
  select(-c(grep("_udc", colnames(stemi_scenario_df)),
            "nmb"))

scenario_df_udc <- stemi_scenario_df %>%
  select(c("group",
           "scenario",
           "intervention",
           grep("_udc", colnames(stemi_scenario_df))))

if (SAVE_OUTPUTS){
  write.csv(scenario_df_dc,
            "formatted_outputs/stemi_dc_scenario_psa.csv")
  write.csv(scenario_df_udc,
            "formatted_outputs/stemi_udc_scenario_psa.csv")
}

# Now do nstemi scenarios
nstemi_scenario_df <- data.frame(group = rep("NSTEMI",
                                            2 * n_nstemi_scenario),
                                scenario = rep(1:n_nstemi_scenario, 2),
                                intervention = c(rep("standard DAPT", n_nstemi_scenario),
                                                 rep("genotype-guided DAPT", n_nstemi_scenario)))

for (output in all_outputs){
  decimals <- ifelse(grepl("cost|icer|nmb", output),
                     0,
                     no = ifelse(grepl("qalys", output),
                                 yes = 4,
                                 no = 2))
  nstemi_scenario_df[[output]] <- sapply(c("sc", "pc"),
                                        print_nstemi_scenarios,
                                        output = output,
                                        decimals = decimals) %>% c()
}

nstemi_scenario_df <- nstemi_scenario_df %>%
  relocate(cost, .after = cost_udc) %>%
  mutate(inc_cost = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(nstemi_sa_psa$cost_inc_mean[scenario],
                                              0),
                                        nsmall = 0),
                                 " (",
                                 format(round(nstemi_sa_psa$cost_inc_L[scenario],
                                              0),
                                        nsmall = 0),
                                 ", ",
                                 format(round(nstemi_sa_psa$cost_inc_U[scenario],
                                              0),
                                        nsmall = 0),
                                 ")",
                                 sep="")
  )) %>%
  mutate(inc_cost_udc = ifelse(intervention == "standard DAPT",
                               NA,
                               paste(format(round(nstemi_sa_psa$cost_udc_inc_mean[scenario],
                                                  0),
                                            nsmall = 0),
                                     " (",
                                     format(round(nstemi_sa_psa$cost_udc_inc_L[scenario],
                                                  0),
                                            nsmall = 0),
                                     ", ",
                                     format(round(nstemi_sa_psa$cost_udc_inc_U[scenario],
                                                  0),
                                            nsmall = 0),
                                     ")",
                                     sep="")
  )) %>%
  mutate(inc_util = ifelse(intervention == "standard DAPT",
                           NA,
                           paste(format(round(nstemi_sa_psa$util_inc_mean[scenario],
                                              3),
                                        nsmall = 3),
                                 " (",
                                 format(round(nstemi_sa_psa$util_inc_L[scenario],
                                              3),
                                        nsmall = 3),
                                 ", ",
                                 format(round(nstemi_sa_psa$util_inc_U[scenario],
                                              3),
                                        nsmall = 3),
                                 ")",
                                 sep="")
  )) %>%
  mutate(inc_util_udc = ifelse(intervention == "standard DAPT",
                               NA,
                               paste(format(round(nstemi_sa_psa$util_udc_inc_mean[scenario],
                                                  3),
                                            nsmall = 3),
                                     " (",
                                     format(round(nstemi_sa_psa$util_udc_inc_L[scenario],
                                                  3),
                                            nsmall = 3),
                                     ", ",
                                     format(round(nstemi_sa_psa$util_udc_inc_U[scenario],
                                                  3),
                                            nsmall = 3),
                                     ")",
                                     sep="")
  )) %>%
  mutate(icer = ifelse(intervention == "standard DAPT",
                       NA,
                       paste(format(round(nstemi_sa_psa$icer_mean[scenario],
                                          0),
                                    nsmall = 0),
                             " (",
                             format(round(nstemi_sa_psa$icer_L[scenario],
                                          0),
                                    nsmall = 0),
                             ", ",
                             format(round(nstemi_sa_psa$icer_U[scenario],
                                          0),
                                    nsmall = 0),
                             ")",
                             sep="")
  )) %>%
  mutate(inc_nmb = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(round(nstemi_sa_psa$total_nmb_inc_mean[scenario],
                                             0),
                                       nsmall = 0),
                                " (",
                                format(round(nstemi_sa_psa$total_nmb_inc_L[scenario],
                                             0),
                                       nsmall = 0),
                                ", ",
                                format(round(nstemi_sa_psa$total_nmb_inc_U[scenario],
                                             0),
                                       nsmall = 0),
                                ")",
                                sep="")
  )) %>%
  mutate(prob_ce = ifelse(intervention == "standard DAPT",
                          NA,
                          paste(format(round(100 * nstemi_scenario_prob_ce$ce_prob[scenario],
                                             1),
                                       nsmall = 0),
                                " (",
                                format(round(100 * nstemi_scenario_prob_ce$lower_95[scenario],
                                             1),
                                       nsmall = 0),
                                ", ",
                                format(round(100 * nstemi_scenario_prob_ce$upper_95[scenario],
                                             1),
                                       nsmall = 0),
                                ")",
                                sep="")
  )) %>%
  relocate(nmb, .before = inc_nmb) %>%
  relocate(total_nmb, .before = inc_nmb) %>%
  arrange(scenario) %>% 
  group_split(scenario) %>% 
  map_dfr(~ add_row(.x, .after = Inf))

scenario_df_dc <- nstemi_scenario_df %>%
  select(-c(grep("_udc", colnames(nstemi_scenario_df)),
            "nmb"))

scenario_df_udc <- nstemi_scenario_df %>%
  select(c("group",
           "scenario",
           "intervention",
           grep("_udc", colnames(nstemi_scenario_df))))

if (SAVE_OUTPUTS){
  write.csv(scenario_df_dc,
            "formatted_outputs/nstemi_dc_scenario_psa.csv")
  write.csv(scenario_df_udc,
            "formatted_outputs/nstemi_udc_scenario_psa.csv")
}