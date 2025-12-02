# Script for gathering extra outputs from PSA

SAVE_OUTPUTS <- TRUE
dir.create("additional_outputs",
           showWarnings = FALSE)

library(stringr)
library(tidyverse)

# Useful formatting function
print_ci <- function(m,
                     l,
                     u,
                     decimals=2){
  paste(format(round(m, decimals),
               nsmall = decimals),
        " (",
        format(round(l, decimals),
               nsmall = decimals),
        ", ",
        format(round(u, decimals),
               nsmall = decimals),
        ")",
        sep="")
}

n_sample <- 1e5

stemi_samples <- read.csv("outputs/stemi_n_1e+05_psa_samples.csv")
nstemi_samples <- read.csv("outputs/nstemi_n_1e+05_psa_samples.csv")

stemi_size = 41690
nstemi_size = 61248

extra_stemi_outputs <- data.frame(pop_inc_life_years_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(stemi_size * (life_years_pc - life_years_sc)))$mean,
                                  pop_inc_life_years_L = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * (life_years_pc - life_years_sc)), 0.025))$q,
                                  pop_inc_life_years_U = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * (life_years_pc - life_years_sc)), 0.975))$q,
                                  pop_inc_qalys_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(stemi_size * (util_dc_pc_hs - util_dc_sc_hs)))$mean,
                                  pop_inc_qalys_L = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * (util_dc_pc_hs - util_dc_sc_hs)), 0.025))$q,
                                  pop_inc_qalys_U = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * (util_dc_pc_hs - util_dc_sc_hs)), 0.975))$q,
                                  pop_nmb_inc_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(stemi_size * inc_nmb))$mean,
                                  pop_nmb_inc_L = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * inc_nmb), 0.025))$q,
                                  pop_nmb_inc_U = summarise(stemi_samples[1:n_sample, ], q = quantile((stemi_size * inc_nmb), 0.975))$q,
                                  inc_mi_dt_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(mi_dt_sc - mi_dt_pc))$mean,
                                  inc_mi_dt_L = summarise(stemi_samples[1:n_sample, ], q = quantile((mi_dt_sc - mi_dt_pc), 0.025))$q,
                                  inc_mi_dt_U = summarise(stemi_samples[1:n_sample, ], q = quantile((mi_dt_sc - mi_dt_pc), 0.975))$q,
                                  inc_mi_mc_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(mi_mc_sc - mi_mc_pc))$mean,
                                  inc_mi_mc_L = summarise(stemi_samples[1:n_sample, ], q = quantile((mi_mc_sc - mi_mc_pc), 0.025))$q,
                                  inc_mi_mc_U = summarise(stemi_samples[1:n_sample, ], q = quantile((mi_mc_sc - mi_mc_pc), 0.975))$q,
                                  inc_stroke_dt_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(stroke_dt_sc - stroke_dt_pc))$mean,
                                  inc_stroke_dt_L = summarise(stemi_samples[1:n_sample, ], q = quantile((stroke_dt_sc - stroke_dt_pc), 0.025))$q,
                                  inc_stroke_dt_U = summarise(stemi_samples[1:n_sample, ], q = quantile((stroke_dt_sc - stroke_dt_pc), 0.975))$q,
                                  inc_stroke_mc_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(stroke_mc_sc - stroke_mc_pc))$mean,
                                  inc_stroke_mc_L = summarise(stemi_samples[1:n_sample, ], q = quantile((stroke_mc_sc - stroke_mc_pc), 0.025))$q,
                                  inc_stroke_mc_U = summarise(stemi_samples[1:n_sample, ], q = quantile((stroke_mc_sc - stroke_mc_pc), 0.975))$q,
                                  inc_major_bleed_dt_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(major_bleed_dt_sc - major_bleed_dt_pc))$mean,
                                  inc_major_bleed_dt_L = summarise(stemi_samples[1:n_sample, ], q = quantile((major_bleed_dt_sc - major_bleed_dt_pc), 0.025))$q,
                                  inc_major_bleed_dt_U = summarise(stemi_samples[1:n_sample, ], q = quantile((major_bleed_dt_sc - major_bleed_dt_pc), 0.975))$q,
                                  inc_minor_bleed_dt_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(minor_bleed_dt_sc - minor_bleed_dt_pc))$mean,
                                  inc_minor_bleed_dt_L = summarise(stemi_samples[1:n_sample, ], q = quantile((minor_bleed_dt_sc - minor_bleed_dt_pc), 0.025))$q,
                                  inc_minor_bleed_dt_U = summarise(stemi_samples[1:n_sample, ], q = quantile((minor_bleed_dt_sc - minor_bleed_dt_pc), 0.975))$q,
                                  inc_death_dt_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(death_dt_sc - death_dt_pc))$mean,
                                  inc_death_dt_L = summarise(stemi_samples[1:n_sample, ], q = quantile((death_dt_sc - death_dt_pc), 0.025))$q,
                                  inc_death_dt_U = summarise(stemi_samples[1:n_sample, ], q = quantile((death_dt_sc - death_dt_pc), 0.975))$q,
                                  inc_death_mc_mean = summarise(stemi_samples[1:n_sample, ], mean = mean(death_mc_sc - death_mc_pc))$mean,
                                  inc_death_mc_L = summarise(stemi_samples[1:n_sample, ], q = quantile((death_mc_sc - death_mc_pc), 0.025))$q,
                                  inc_death_mc_U = summarise(stemi_samples[1:n_sample, ], q = quantile((death_mc_sc - death_mc_pc), 0.975))$q)

stemi_df <- data.frame(var = colnames(extra_stemi_outputs),
                       value = as.numeric(as.vector(extra_stemi_outputs[1, ]))) %>%
  mutate(direction = ifelse(grepl("_mean", var),
                            yes = "mean",
                            no = ifelse(grepl("_L", var),
                                        yes = "L",
                                        no = "U"))) %>%
  mutate(var = var %>%
           str_replace_all("_mean|_L|_U",
                           "")) %>%
  pivot_wider(names_from = direction,
              values_from = value) %>%
  mutate(decimals = ifelse(grepl("nmb", var),
                           yes = 0,
                           no = 2)) %>%
  mutate(estimate = sapply(1:length(var), FUN = function(i) print_ci(mean[i], L[i], U[i], decimals[i]))) %>%
  select(var, estimate)
write.csv(stemi_df,
          file = "formatted_outputs/additional_outputs_stemi.csv")


extra_nstemi_outputs <- data.frame(pop_inc_life_years_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(nstemi_size * (life_years_pc - life_years_sc)))$mean,
                                  pop_inc_life_years_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * (life_years_pc - life_years_sc)), 0.025))$q,
                                  pop_inc_life_years_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * (life_years_pc - life_years_sc)), 0.975))$q,
                                  pop_inc_qalys_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(nstemi_size * (util_dc_pc_hs - util_dc_sc_hs)))$mean,
                                  pop_inc_qalys_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * (util_dc_pc_hs - util_dc_sc_hs)), 0.025))$q,
                                  pop_inc_qalys_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * (util_dc_pc_hs - util_dc_sc_hs)), 0.975))$q,
                                  pop_nmb_inc_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(nstemi_size * inc_nmb))$mean,
                                  pop_nmb_inc_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * inc_nmb), 0.025))$q,
                                  pop_nmb_inc_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((nstemi_size * inc_nmb), 0.975))$q,
                                  inc_mi_dt_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(mi_dt_sc - mi_dt_pc))$mean,
                                  inc_mi_dt_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((mi_dt_sc - mi_dt_pc), 0.025))$q,
                                  inc_mi_dt_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((mi_dt_sc - mi_dt_pc), 0.975))$q,
                                  inc_mi_mc_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(mi_mc_sc - mi_mc_pc))$mean,
                                  inc_mi_mc_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((mi_mc_sc - mi_mc_pc), 0.025))$q,
                                  inc_mi_mc_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((mi_mc_sc - mi_mc_pc), 0.975))$q,
                                  inc_stroke_dt_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(stroke_dt_sc - stroke_dt_pc))$mean,
                                  inc_stroke_dt_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((stroke_dt_sc - stroke_dt_pc), 0.025))$q,
                                  inc_stroke_dt_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((stroke_dt_sc - stroke_dt_pc), 0.975))$q,
                                  inc_stroke_mc_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(stroke_mc_sc - stroke_mc_pc))$mean,
                                  inc_stroke_mc_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((stroke_mc_sc - stroke_mc_pc), 0.025))$q,
                                  inc_stroke_mc_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((stroke_mc_sc - stroke_mc_pc), 0.975))$q,
                                  inc_major_bleed_dt_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(major_bleed_dt_sc - major_bleed_dt_pc))$mean,
                                  inc_major_bleed_dt_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((major_bleed_dt_sc - major_bleed_dt_pc), 0.025))$q,
                                  inc_major_bleed_dt_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((major_bleed_dt_sc - major_bleed_dt_pc), 0.975))$q,
                                  inc_minor_bleed_dt_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(minor_bleed_dt_sc - minor_bleed_dt_pc))$mean,
                                  inc_minor_bleed_dt_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((minor_bleed_dt_sc - minor_bleed_dt_pc), 0.025))$q,
                                  inc_minor_bleed_dt_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((minor_bleed_dt_sc - minor_bleed_dt_pc), 0.975))$q,
                                  inc_death_dt_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(death_dt_sc - death_dt_pc))$mean,
                                  inc_death_dt_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((death_dt_sc - death_dt_pc), 0.025))$q,
                                  inc_death_dt_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((death_dt_sc - death_dt_pc), 0.975))$q,
                                  inc_death_mc_mean = summarise(nstemi_samples[1:n_sample, ], mean = mean(death_mc_sc - death_mc_pc))$mean,
                                  inc_death_mc_L = summarise(nstemi_samples[1:n_sample, ], q = quantile((death_mc_sc - death_mc_pc), 0.025))$q,
                                  inc_death_mc_U = summarise(nstemi_samples[1:n_sample, ], q = quantile((death_mc_sc - death_mc_pc), 0.975))$q)

nstemi_df <- data.frame(var = colnames(extra_nstemi_outputs),
                       value = as.numeric(as.vector(extra_nstemi_outputs[1, ]))) %>%
  mutate(direction = ifelse(grepl("_mean", var),
                            yes = "mean",
                            no = ifelse(grepl("_L", var),
                                        yes = "L",
                                        no = "U"))) %>%
  mutate(var = var %>%
           str_replace_all("_mean|_L|_U",
                           "")) %>%
  pivot_wider(names_from = direction,
              values_from = value) %>%
  mutate(decimals = ifelse(grepl("nmb", var),
                           yes = 0,
                           no = 2)) %>%
  mutate(estimate = sapply(1:length(var), FUN = function(i) print_ci(mean[i], L[i], U[i], decimals[i]))) %>%
  select(var, estimate)
write.csv(stemi_df,
          file = "formatted_outputs/additional_outputs_nstemi.csv")