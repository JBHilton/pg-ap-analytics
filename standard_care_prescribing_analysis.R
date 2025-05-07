# In this script we run the deterministic analysis over a range of values for
# the proportions of different drugs prescribed under standard care.

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

icer_df <- data.frame(p_ac = c(0),
                      p_at = c(0),
                      p_ap = c(0),
                      ICER_udc = c(0),
                      ICER_dc = c(0))

for (p_ac in seq(0., 1., 0.01)){
  for (p_at in seq(0., 1. - p_ac, 0.01)){
    p_ap = 1. - p_ac - p_at
    arm_df <- run_arm_comparison(c(p_ac,
                                   p_at,
                                   p_ap))
    ICER_udc <- arm_df$ratio_udc[arm_df$arm=="Increment"]
    ICER_dc <- arm_df$ratio[arm_df$arm=="Increment"]
    icer_df <- icer_df %>% add_row(p_ac = p_ac,
                                   p_at = p_at,
                                   p_ap = p_ap,
                                   ICER_udc = ICER_udc,
                                   ICER_dc = ICER_dc)
  }
}
icer_df <- icer_df[2:nrow(icer_df), ]

p <- ggplot(icer_df, aes(p_at, p_ap, fill = ICER_dc)) +
  geom_tile()
p

p <- ggplot(icer_df, aes(p_at, p_ap, fill = ICER_udc)) +
  geom_tile()
p
