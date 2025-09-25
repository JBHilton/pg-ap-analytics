# In this script we generate plots of the PSA results for both the STEMI and
# NSTEMI subgroups.

dir.create("plots", showWarnings = FALSE)

n_sample <- 1e5

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
SAVE_FILEPATH <- paste("stemi_n_", n_sample, sep="")

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

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 10e4
ce_thresh_2 = 3 * 10e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

stemi_ce_plane_plot <- ggplot(multi_results[1:n_sample, ],
                              aes(x = inc_util_dc_hs,
                                  y = inc_cost_dc_hs)) +
  geom_point(size = .1) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "red"),
             size = 5.,
             shape = 18) +
  geom_line(data = ce_line_df,
            aes(x=x,
                y=y,
                color = factor(ce_threshold))) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "CE threshold £30,000",
                                "Point estimate"),
                     values = c("green",
                                "orange",
                                "red")) +
  theme(legend.title=element_blank(),
        text = element_text(size = 20)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane.png", sep=""),
       plot = stemi_ce_plane_plot,
       width = 12,
       height = 10)

ce_thresh_df <- read_csv(file = paste("outputs/",
                                      SAVE_FILEPATH,
                                      "_acceptance_probability.csv",
                                      sep="")) %>%
  filter(ce_threshold < 18000)

stemi_ce_thresh_plot <- ggplot(ce_thresh_df,
                               aes(x = ce_threshold,
                                   y = acceptance_prob)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_95,
                  ymax = upper_95),
              alpha = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
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
       height = 10)

#### Now do NSTEMI ####

SAVE_FILEPATH <- paste("nstemi_n_", n_sample, sep="")

det_results <- read_csv(file = paste("outputs/nstemi_arm_comparison.csv",
                                     sep = ""))
multi_results <- read_csv(file = paste("outputs/",
                                       SAVE_FILEPATH,
                                       "_psa_samples.csv",
                                       sep = ""))

det_cost_util <- det_results %>%
  filter(output_name %in% c("cost", "util")) %>%
  select(c(output_name, inc)) %>%
  column_to_rownames("output_name") %>%
  t() %>%
  as.data.frame()

# Choose limits of CE threshold lines so that there are minimal changes to axes
ce_thresh_1 = 2 * 10e4
ce_thresh_2 = 3 * 10e4
ce_line_df <- data.frame(x = c(.99 * ylim / ce_thresh_1,
                               .99 * ylim / ce_thresh_2),
                         ce_threshold = c(ce_thresh_1,
                                          ce_thresh_1,
                                          ce_thresh_2,
                                          ce_thresh_2)) %>%
  mutate(y = x * ce_threshold)

nstemi_ce_plane_plot <- ggplot(multi_results[1:n_sample, ],
                               aes(x = inc_util_dc_hs,
                                   y = inc_cost_dc_hs)) +
  geom_point(size = .1) +
  geom_point(data = det_cost_util,
             aes(x = util,
                 y = cost,
                 color = "red"),
             size = 5.,
             shape = 18) +
  geom_line(data = ce_line_df,
            aes(x=x,
                y=y,
                color = factor(ce_threshold))) +
  scale_x_continuous(breaks = seq(xlim[1], xlim[2], .05),
                     limits = xlim) +
  scale_y_continuous(breaks = seq(ylim[1], ylim[2], 200),
                     limits = ylim) +
  scale_color_manual(labels = c("CE threshold £20,000",
                                "CE threshold £30,000",
                                "Point estimate"),
                     values = c("green",
                                "orange",
                                "red")) +
  theme(legend.title=element_blank(),
        text = element_text(size = 20)) +
  xlab("Incremental utility (QALYs)") +
  ylab("Incremental cost (£s)")

ggsave(paste("plots/",
             SAVE_FILEPATH,
             "ce_plane.png", sep=""),
       plot = nstemi_ce_plane_plot,
       width = 12,
       height = 10)

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
  theme(text = element_text(size = 20)) +
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
       height = 10)