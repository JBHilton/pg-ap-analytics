# In this script we load and reformat the different bits of output data into
# tables of base cases and confidence intervals.

library(tidyverse)

stemi_arm_comparison <- read.csv("stemi_arm_comparison.csv")
stemi_psa <- read.csv("stemi_n_1e+05_psa_stats.csv")
stemi_prob_ce <- read.csv("stemi_n_1e+05_acceptance_probability.csv")
nstemi_arm_comparison <- read.csv("nstemi_arm_comparison.csv")
nstemi_psa <- read.csv("nstemi_n_1e+05_psa_stats.csv")
nstemi_prob_ce <- read.csv("nstemi_n_1e+05_acceptance_probability.csv")

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
                                          output = "lifeyears",
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
                                   output = "lifeyears",
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
               inc_nmb = c(0,
                           print_nstemi_output("inc",
                                              "nmb",
                                              digits=5)),
               percent_ce = c(0,
                              print_ce_prob(nstemi_prob_ce,
                                            digits=5))
    )
  )
