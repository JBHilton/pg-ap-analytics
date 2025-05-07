# In this script we run the deterministic analysis over a range of values for
# the cost of PGX testing to identify the cost-effectiveness threshold.

# Flag to determine whether to convert rates in the data to probabilities
CONVERSIONS_REQUIRED <- FALSE

# Next flag is here to suppress some code which isn't currently useful but might
# be at some point
SUBPOP_PLOTTING <- FALSE

# Step correction for death in decision tree and Markov model, i.e. average time
# to death given death occurs within a given year.
AVE_TIME_TO_EVENT <- 0.5

library("dplyr")
library("expm")
library("ggplot2")
library("Matrix")
library("readxl")
library("stringr")
library("tidyverse")

source(file = "simulation_functions.R")

icer_df <- data.frame(poct_cost = c(0),
                      ICER_udc = c(0),
                      ICER_dc = c(0))

for (poct_cost in seq(125., 1000., 25.)){
  arm_df <- run_arm_comparison(poct_cost = poct_cost)
  ICER_udc <- arm_df$ratio_udc[arm_df$arm=="Increment"]
  ICER_dc <- arm_df$ratio[arm_df$arm=="Increment"]
  icer_df <- icer_df %>% add_row(poct_cost = poct_cost,
                                 ICER_udc = ICER_udc,
                                 ICER_dc = ICER_dc)
}
icer_df <- icer_df[2:nrow(icer_df), ]

p <- ggplot(icer_df, aes(poct_cost, ICER_dc)) +
  geom_line()
p

p <- ggplot(icer_df, aes(poct_cost, ICER_udc)) +
  geom_line()
p
