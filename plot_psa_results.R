# In this script we generate plots of the PSA results for both the STEMI and
# NSTEMI subgroups.

dir.create("plots", showWarnings = FALSE)

nsample <- 1e5
figure_dpi <- 1200

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


#### Do STEMI first #####
SAVE_FILEPATH <- paste("stemi_n_", nsample, sep="")

det_results <- read_csv(file = paste("outputs/stemi_arm_comparison.csv",
                                     sep = ""))
multi_results <- read_csv(file = paste("outputs/",
                                       SAVE_FILEPATH,
                                       "_psa_samples.csv",
                                       sep = ""))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))
# ylim <- c(-10000,
#           10000)

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

stemi_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                              aes(x = inc_util_dc_hs,
                                  y = inc_cost_dc_hs)) +
  geom_point(size = .5, alpha = .75) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0),
             alpha = 1.) +
  geom_vline(aes(xintercept = 0),
             alpha = 1.) +
  geom_line(data = ce_line_df,
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "CE threshold £30,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#D56C25",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane.png", sep=""),
       plot = stemi_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

stemi_ce_plane_plot_20k <- ggplot(multi_results[1:nsample, ],
                              aes(x = inc_util_dc_hs,
                                  y = inc_cost_dc_hs)) +
  geom_point(size = .5, alpha = .75) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0),
             alpha = 1.) +
  geom_vline(aes(xintercept = 0),
             alpha = 1.) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = "#156082"),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("Willingness to pay threshold",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane_20k_only.png", sep=""),
       plot = stemi_ce_plane_plot_20k,
       width = 10,
       height = 10,
       dpi = figure_dpi)

ce_thresh_df <- read_csv(file = paste("outputs/",
                                      SAVE_FILEPATH,
                                      "_acceptance_probability.csv",
                                      sep="")) %>%
  filter(ce_threshold <= 18000)

stemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_probability.png", sep=""),
       plot = stemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

#### Now do NSTEMI ####

SAVE_FILEPATH <- paste("nstemi_n_", nsample, sep="")

det_results <- read_csv(file = paste("outputs/nstemi_arm_comparison.csv",
                                     sep = ""))
multi_results <- read_csv(file = paste("outputs/",
                                       SAVE_FILEPATH,
                                       "_psa_samples.csv",
                                       sep = ""))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

nstemi_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                               aes(x = inc_util_dc_hs,
                                   y = inc_cost_dc_hs)) +
  geom_point(size = .5, alpha = .75) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  geom_line(data = ce_line_df,
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  geom_hline(aes(yintercept = 0),
             alpha = 1.) +
  geom_vline(aes(xintercept = 0),
             alpha = 1.) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "CE threshold £30,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#D56C25",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane.png", sep=""),
       plot = nstemi_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

nstemi_ce_plane_plot_20k <- ggplot(multi_results[1:nsample, ],
                                  aes(x = inc_util_dc_hs,
                                      y = inc_cost_dc_hs)) +
  geom_point(size = .5, alpha = .75) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0),
             alpha = 1.) +
  geom_vline(aes(xintercept = 0),
             alpha = 1.) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = "#156082"),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("Willingness to pay threshold",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane_20k_only.png", sep=""),
       plot = nstemi_ce_plane_plot_20k,
       width = 10,
       height = 10,
       dpi = figure_dpi)

ce_thresh_df <- read_csv(file = paste("outputs/",
                                      SAVE_FILEPATH,
                                      "_acceptance_probability.csv",
                                      sep="")) %>%
  filter(ce_threshold < 18000)

nstemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                                aes(x = ce_threshold,
                                    y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_probability.png", sep=""),
       plot = nstemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

stemi_ce_thresh_df <- read_csv(file = "outputs/stemi_n_1e+05_acceptance_probability.csv") %>%
  filter(ce_threshold < 18000) %>%
  mutate(Population = "STEMI")
nstemi_ce_thresh_df <- read_csv(file = "outputs/nstemi_n_1e+05_acceptance_probability.csv") %>%
  filter(ce_threshold < 18000) %>%
  mutate(Population = "NSTEMI")
ce_thresh_df <- rbind(stemi_ce_thresh_df,
                      nstemi_ce_thresh_df)

both_ce_thresh_plot <- ggplot(ce_thresh_df,
                                aes(x = ce_threshold,
                                    y = acceptance_prob,
                                    colour = Population)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  scale_color_manual(values = c("#156082",
                                "#D56C25")) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/both_groups_ce_probability.png", sep=""),
       plot = both_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

#### Now look at scenario analyses ####

# Want CE planes and acceptance probabilities for:
# - Scenario 1
# - Scenario 11 (cost poct doubled)
# - Scenario 13 (Ticagrelor off-patent)
# Need these for both STEMI and NSTEMI

# Function to calculate the acceptance probability for a given cost
# effectiveness threshold given a set of ICERs from PSA
calculate_acc_prob <- function(ce_threshold,
                               utils,
                               icers,
                               nsample){
  length(which((icers < ce_threshold) & (utils > 0))) /
    nsample
} %>% Vectorize()

nsample <- 1e4

# Population sizes for total NMB:
stemi_size = 41690
nstemi_size = 61248

stemi_sa_psa <- read.csv("outputs/stemi_n_10000_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = stemi_size * nmb_sc_mean,
         total_nmb_sc_L = stemi_size * nmb_sc_L,
         total_nmb_sc_U = stemi_size * nmb_sc_U,
         total_nmb_pc_mean = stemi_size * nmb_pc_mean,
         total_nmb_pc_L = stemi_size * nmb_pc_L,
         total_nmb_pc_U = stemi_size * nmb_pc_U)
stemi_scenario_prob_ce <- read.csv("outputs/stemi_n_1000_psa_acceptance_probability.csv")
stemi_scenario_samples <- read.csv('outputs/stemi_n_10000_scenario_samples.csv')

# Plots for SA1:

scenario_to_plot <- 'SA1'
multi_results <- stemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario)) %>%
  filter(!grepl("SA10|SA11|SA12|SA13", scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

stemi_SA1_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                              aes(x = inc_util_dc_hs,
                                  y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/stemi_SA1_ce_plane.png", sep=""),
       plot = stemi_SA1_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

stemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/stemi_SA1_ce_probability.png", sep=""),
       plot = stemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Now do SA11

scenario_to_plot <- 'SA11'
multi_results <- stemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

stemi_SA11_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                                  aes(x = inc_util_dc_hs,
                                      y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/stemi_SA11_ce_plane.png", sep=""),
       plot = stemi_SA11_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

stemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/stemi_SA11_ce_probability.png", sep=""),
       plot = stemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Now do SA13

scenario_to_plot <- 'SA13'
multi_results <- stemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

stemi_SA13_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                                   aes(x = inc_util_dc_hs,
                                       y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/stemi_SA13_ce_plane.png", sep=""),
       plot = stemi_SA13_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

stemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/stemi_SA13_ce_probability.png", sep=""),
       plot = stemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

### Now do same for NSTEMI ###
nstemi_sa_psa <- read.csv("outputs/nstemi_n_10000_scenario_psa_stats.csv") %>%
  mutate(total_nmb_sc_mean = nstemi_size * nmb_sc_mean,
         total_nmb_sc_L = nstemi_size * nmb_sc_L,
         total_nmb_sc_U = nstemi_size * nmb_sc_U,
         total_nmb_pc_mean = nstemi_size * nmb_pc_mean,
         total_nmb_pc_L = nstemi_size * nmb_pc_L,
         total_nmb_pc_U = nstemi_size * nmb_pc_U)
nstemi_scenario_prob_ce <- read.csv("outputs/nstemi_n_1000_psa_acceptance_probability.csv")
nstemi_scenario_samples <- read.csv('outputs/nstemi_n_10000_scenario_samples.csv')

# Plots for SA1:

scenario_to_plot <- 'SA1'
multi_results <- nstemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario)) %>%
  filter(!grepl("SA10|SA11|SA12|SA13", scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

nstemi_SA1_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                                  aes(x = inc_util_dc_hs,
                                      y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/nstemi_SA1_ce_plane.png", sep=""),
       plot = nstemi_SA1_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

nstemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/nstemi_SA1_ce_probability.png", sep=""),
       plot = nstemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Now do SA11

scenario_to_plot <- 'SA11'
multi_results <- nstemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

nstemi_SA11_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                                   aes(x = inc_util_dc_hs,
                                       y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/nstemi_SA11_ce_plane.png", sep=""),
       plot = nstemi_SA11_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

nstemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/nstemi_SA11_ce_probability.png", sep=""),
       plot = nstemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Now do SA13

scenario_to_plot <- 'SA13'
multi_results <- nstemi_scenario_samples %>%
  filter(grepl(scenario_to_plot, scenario))

# Calculate axis limits; utility cuts off at nearest .1, cost to nearest £200
xlim <- c(0.1 * floor(min(multi_results$inc_util_dc_hs) / 0.1),
          0.1 * ceiling(max(multi_results$inc_util_dc_hs) / 0.1))
ylim <- c(200 * floor(min(multi_results$inc_cost_dc_hs) / 200),
          200 * ceiling(max(multi_results$inc_cost_dc) / 200))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 1e4
ce_thresh_2 = 3 * 1e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

nstemi_SA13_ce_plane_plot <- ggplot(multi_results[1:nsample, ],
                                   aes(x = inc_util_dc_hs,
                                       y = inc_cost_dc_hs)) +
  geom_point(size = 1.) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "#C00000"),
             size = 5.,
             shape = 18) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_line(data = ce_line_df %>% filter(ce_threshold < 25000),
            aes(x=x,
                y=y,
                color = factor(ce_threshold)),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "Point estimate"),
                     values = c("#156082",
                                "#C00000")) +
  theme(legend.title=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .1),
        text = element_text(size = 20),
        axis.text = element_text(colour = "black",
                                 size = 24)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/nstemi_SA13_ce_plane.png", sep=""),
       plot = nstemi_SA13_ce_plane_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)

# Acceptance probability:

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

ce_thresh_df <- data.frame(ce_threshold = ce_threshold) %>%
  rowwise() %>%
  mutate(acceptance_prob =
           calculate_acc_prob(ce_threshold,
                              multi_results$inc_util_dc_hs[1:nsample],
                              multi_results$icer[1:nsample],
                              nsample))

# Add bootstrap samples:
n_bootstrap <- 100
for (i in 1:n_bootstrap){
  varname <- paste("bootstrap_sample", i)
  sample_ids <- sample(1:nsample,
                       nsample,
                       replace = TRUE)
  sample_utils <- multi_results$inc_util_dc_hs[sample_ids]
  sample_icers <- multi_results$icer[sample_ids]
  temp_df <- data.frame(ce_threshold = ce_threshold) %>%
    rowwise() %>%
    mutate(acceptance_prob =
             calculate_acc_prob(ce_threshold, sample_utils, sample_icers, nsample))
  ce_thresh_df[[varname]] <- temp_df$acceptance_prob
}
ce_thresh_df$lower_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.025)
ce_thresh_df$upper_95 = rowQuantiles(as.matrix(ce_thresh_df[, -c(1, 2)]),
                                     probs = 0.975)
ce_thresh_df <- ce_thresh_df[, -which(grepl("bootstrap", colnames(ce_thresh_df)))] %>%
  filter(ce_threshold <= 18000)

nstemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 18000, 2000)) +
  xlab("Cost-effectiveness threshold (£)") +
  ylab("Cost-effectiveness probability")

ggsave(paste("plots/nstemi_SA13_ce_probability.png", sep=""),
       plot = nstemi_ce_thresh_plot,
       width = 10,
       height = 10,
       dpi = figure_dpi)