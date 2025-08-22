# This script contains functions needed for carrying out probabilistic
# sensitivity analysis.

# Start ages specified on row 1 of parameter matrices:
start_age_female <- parameters_STEMI$value[(parameters_STEMI$variable.name=="start_age_female")|(parameters_STEMI$variable.name=="start_age_f")]
start_age_male <- parameters_STEMI$value[(parameters_STEMI$variable.name=="start_age_male")|(parameters_STEMI$variable.name=="start_age_m")]

prop_male <- parameters_STEMI$value[
  parameters_STEMI$variable.name=="proportion_male"]

# Need to load in mortality table
if (CASE=="STEMI"){
  mort_table <- read_xlsx("data-inputs/masterfile_070725.xlsx",
                         sheet = "age_sex_dependant_mortality",
                         range = "A12:D52")
}else{
  mort_table <- read_xlsx("data-inputs/NSTEMI_masterfile_160725.xlsm",
                          sheet = "age_sex_dependant_mortality",
                          range = "A12:D51")
}
# Function for calculating mortality by age for a cohort given an SMR
apply_smr <- function(mort_female,
                      mort_male,
                      smr,
                      prop_male){
  adj_female <- pmin(1, smr * mort_female)
  adj_male <- pmin(1, smr * mort_male)
  return((1 - prop_male) * adj_female +
           prop_male * adj_male)
}



# Function for calculating utility of a given event in a given year
utility_formula <- function(util_val,
                            tstep,
                            prop_male){
  u_male <- util_val * (0.9508566 +
                     (0.0212126 * 1) -
                     (0.0002587 * (tstep + start_age_male)) -
                     (0.0000332 * (tstep + start_age_male)^2))
  u_female <- util_val * (0.9508566 +
                       (0.0212126 * 0) -
                       (0.0002587 * (tstep + start_age_female)) -
                       (0.0000332 * (tstep + start_age_female)^2))
  return(prop_male * u_male + (1 - prop_male) * u_female)
}

# Function for calculating a utility over time matrix from a PSA draw
rewrite_markov_utils <- function(draw_df,
                                 time_hor,
                                 prop_male){
  markov_utils <- data.frame(no_event = sapply(1:time_hor,
                                               FUN = function(t){
                                                 utility_formula(draw_df$draw[
                                                   which(draw_df$variable.name=="utility_no_event")],
                                                   t,
                                                   prop_male)
                                               }),
                             stroke = sapply(1:time_hor,
                                               FUN = function(t){
                                                 utility_formula(draw_df$draw[
                                                   which(draw_df$variable.name=="utility_stroke")],
                                                   t,
                                                   prop_male)
                                               }),
                             post_stroke = sapply(1:time_hor,
                                               FUN = function(t){
                                                 utility_formula(draw_df$draw[
                                                   which(draw_df$variable.name=="utility_post_stroke")],
                                                   t,
                                                   prop_male)
                                               }),
                             mi = sapply(1:time_hor,
                                               FUN = function(t){
                                                 utility_formula(draw_df$draw[
                                                   which(draw_df$variable.name=="utility_mi")],
                                                   t,
                                                   prop_male)
                                               }),
                             post_mi = sapply(1:time_hor,
                                               FUN = function(t){
                                                 utility_formula(draw_df$draw[
                                                   which(draw_df$variable.name=="utility_post_mi")],
                                                   t,
                                                   prop_male)
                                               }))
  return(markov_utils)
}

# Function for building a new probability dataframe based on updates to baseline
# risks or odds/hazard ratios. A set of baseline risks needs to be provided
rescale_probs <- function(baseline_prob_df,
                          at_ratio_df,
                          ap_ratio_df,
                          ac_no_lof_ratio_df
){
  # Extract probabilities and odds/hazard ratios for each drug/genotype:
  ac_lof_probs <- baseline_prob_df
  
  # Work out Ticagrelor probabilities
  at_probs <- at_ratio_df %>%
    mutate(value = prob_from_or(value, baseline_prob_df$odds))
  
  # Now do Prasagruel:
  ap_probs <- ap_ratio_df %>%
    mutate(value = prob_from_or(value, baseline_prob_df$odds))
  
  # And clopidogrel with no loss of function:
  
  if (UNIQUE_NO_LOF_PROBS){
    ac_no_lof_probs <- ac_no_lof_ratio_df %>%
      mutate(value = ifelse((grepl("death", variable.name)|
                              grepl("_mi_", variable.name)),
                            value * at_probs$value,
                            value * baseline_prob_df$value))
  }else{
    ac_no_lof_probs <- ac_no_lof_ratio_df %>%
      mutate(value = value * baseline_prob_df$value) %>%
      add_column(at_value = at_probs$value) %>%
      mutate(value = ifelse(grepl("death", variable.name)|grepl("_mi_", variable.name),
                            yes = at_value,
                            no = value)) %>%
      select(c(variable.name, value))
  }
  
  # Assemble into single probability dataframe:
  prob_df <- rbind(ac_lof_probs %>% select(-odds),
                   ac_no_lof_probs,
                   at_probs,
                   ap_probs) %>%
    mutate(variable.name = variable.name %>%
             str_replace("_vs_ac_lof", "") %>%
             str_replace("^or_", "prob_") %>%
             str_replace("^hr_", "prob_") %>%
             str_replace("^rr_", "prob_"))
  return(prob_df)
}

rewrite_dt_utilities <- function(baseline_util_df){
  # Utilities associated with events:
  event_utilities <- par_df %>%
    filter(grepl("utility", parameter.list)) %>%
    select(c(variable.name, value)) %>%
    mutate(variable.name = gsub("utility_", "", variable.name) %>%
             str_replace_all("_tree", "")) %>%
    filter(!grepl("_ac", variable.name)) %>% # Drop any derived drug-specific values
    filter(!grepl("_at", variable.name)) %>%
    filter(!grepl("_ap", variable.name)) %>%
    filter(!grepl("annual", variable.name)) %>% # Drop any annual-scale values
    add_row(variable.name = "death",
            value = 0) # Apply no_event utility before death, 0 afterwards
  
  # Add dyspnoea utilities, converting from whole-year to 30 day values, and
  # convert bleed decrements to utilities
  # NOTE: Utility decrements for bleed events are loaded in in units of years so
  # do not need to be converted.
  duration_dyspnoea <- par_df$value[
    which(par_df$variable.name=="duration_dyspnoea")] / 365
  event_utilities <- event_utilities %>%
    add_row(variable.name = "dyspnoea",
            value = 1 - duration_dyspnoea * event_utilities$value[
              which(event_utilities$variable.name=="u_dec_dyspnoea")]) %>%
    add_row(variable.name = "major_bleed",
            value = 1 - event_utilities$value[
              which(event_utilities$variable.name=="dec_major_bleed_*_(duration_major_bleed/365)")]) %>%
    add_row(variable.name = "minor_bleed",
            value = 1 - event_utilities$value[
              which(event_utilities$variable.name=="dec_minor_bleed_*_(duration_minor_bleed/365)")]) %>%
    filter(!grepl("dec_", variable.name))
  return(event_utilities)
}


do_PSA_draw <- function(uc_df){
  uc_df$draw <- 0
  uc_df$draw[uc_df$distribution=="beta"] <-
    rbeta(length(which(uc_df$distribution=="beta")),
          uc_df$par1[uc_df$distribution=="beta"],
          uc_df$par2[uc_df$distribution=="beta"])
  uc_df$draw[uc_df$distribution=="gamma"] <-
    rgamma(length(which(uc_df$distribution=="gamma")),
           shape = uc_df$par1[uc_df$distribution=="gamma"],
           scale = uc_df$par2[uc_df$distribution=="gamma"])
  uc_df$draw[uc_df$distribution=="lognormal"] <-
    rlnorm(length(which(uc_df$distribution=="lognormal")),
           uc_df$par1[uc_df$distribution=="lognormal"],
           uc_df$par2[uc_df$distribution=="lognormal"])
  uc_df$draw[uc_df$variable.name=="prevalence_lof_base_case"] <- max(0.26,
    uc_df$draw[uc_df$variable.name=="prevalence_lof_base_case"])
  return(uc_df)
}

# Optional function to generate multiple parameter samples in a single
# dataframe. There may be a way to speed up the parameter draws using this basic
# idea, and it could be convenient for situations where we want to save the
# parameters used in a specific simulation.
do_tall_PSA_draw <- function(uc_df,
                            n_samples){
  tall_df <- lapply(1:n_samples,
                    FUN = function(i){
                    do_PSA_draw(uc_df) %>%
                        mutate(sample_id = i)}) %>%
    bind_rows()
  return(tall_df)
}

# Function to update large parameter dataframe with new parameters from random
# draw.
rewrite_pars_from_draw <- function(par_df, draw_df){
  par_df$value[which(par_df$variable.name %in% draw_df$variable.name)] <-
    draw_df$draw[which(draw_df$variable.name %in% par_df$variable.name)]
  par_df$value[par_df$variable.name == "utility_no_event_tree"] <-
    utility_formula(draw_df$draw[
      which(draw_df$variable.name=="utility_no_event")],
      0,
      prop_male)
  par_df$value[par_df$variable.name == "utility_mi_tree"] <-
    utility_formula(draw_df$draw[
      which(draw_df$variable.name=="utility_mi")],
      0,
      prop_male)
  par_df$value[par_df$variable.name == "utility_stroke_tree"] <-
    utility_formula(draw_df$draw[
      which(draw_df$variable.name=="utility_stroke")],
      0,
      prop_male)
  
  markov_pars <- data.frame(parameter.list = c("mi",
                                               "stroke"),
                            value = c(draw_df$draw[draw_df$variable.name=="prob_nevent_to_rinfarc"],
                                      draw_df$draw[draw_df$variable.name=="prob_nevent_to_stk"]))
  return(list(par_df, markov_pars))
}

# Simulation function carrying out a two-arm comparison given a parameter
# dataframe
run_PSA_arm_comparison <- function(par_df,
                                   draw_df,
                                   pop = "STEMI",
                                   scenario = ""){
  rewrite_list <- rewrite_pars_from_draw(par_df,
                                         draw_df)
  par_df <- rewrite_list[[1]]
  markov_pars <- rewrite_list[[2]]
  
  # Extract probabilities and odds/hazard ratios for each drug/genotype:
  baseline_prob_df <- par_df %>%
    filter(grepl("prob", variable.name)) %>%
    filter(grepl("_ac_", variable.name)) %>%
    filter(!grepl("no_lof", variable.name)) %>%
    select(c(variable.name, value)) %>%
    mutate(odds = value/(1-value))
  
  # SA5: adjusted baseline risk of stroke
  if (scenario == "SA5"){
    if (pop == "STEMI"){
      baseline_prob_df$value[
        grepl("stroke", baseline_prob_df$variable.name)] <- 0.0072
    }else{
      baseline_prob_df$value[
        grepl("stroke", baseline_prob_df$variable.name)] <- 0.003
    }
  }
  
  # SA6: adjusted baseline risk of reinfarction
  if (scenario == "SA6"){
    if (pop == "STEMI"){
      baseline_prob_df$value[
        grepl("mi", baseline_prob_df$variable.name)] <- 0.054
    }else{
      baseline_prob_df$value[
        grepl("mi", baseline_prob_df$variable.name)] <- 0.0338
    }
  }
  
  # Work out Ticagrelor probabilities
  at_ratio_df <- par_df %>%
    filter(grepl("^or_", variable.name)) %>%
    filter(grepl("_at", variable.name)) %>%
    select(c(variable.name, value))
  
  # Now do Prasagruel:
  ap_ratio_df <- par_df %>%
    filter(grepl("^or_", variable.name)) %>%
    filter(grepl("_ap", variable.name)) %>%
    select(c(variable.name, value)) %>%
    add_row(variable.name = "or_dyspnoea_ap", value = 1)
  
  # And clopidogrel with no loss of function:
  
  if (UNIQUE_NO_LOF_PROBS){
    ac_no_lof_ratio_df <- par_df %>%
      filter(grepl("^hr_", variable.name)|grepl("^rr_", variable.name)) %>%
      filter(grepl("_ac_no_lof", variable.name)) %>%
      select(c(variable.name, value)) %>%
      add_row(variable.name = "prob_dyspnoea_ac_no_lof", value = 1)
  }else{
    ac_no_lof_ratio_df <- par_df %>%
      filter(grepl("^hr_", variable.name)|grepl("^rr_", variable.name)) %>%
      filter(grepl("_ac_no_lof", variable.name)) %>%
      add_row(variable.name = "prob_dyspnoea_ac_no_lof", value = 1)
  }
  
  # Assemble into single probability dataframe:
  prob_df <- rescale_probs(baseline_prob_df,
                           at_ratio_df,
                           ap_ratio_df,
                           ac_no_lof_ratio_df)
  
  # Set up scaling for discounting by time step
  
  # SA12: 1.5% discount rate
  if (scenario == "SA12"){
    discount_by_cycle <- (1 / (1 + 0.015)^seq(
      1.0, time_hor-1, by = time_step))
  }else{
    discount_by_cycle <- (1 / (1 + par_df$value[(
      par_df$variable.name=="qaly_discount_rate")|
        (par_df$variable.name=="disc_effect_b")])^seq(
          1.0, time_hor-1, by = time_step))
  }
  
  # List parameters for decision tree
  
  # SA2: full uptake
  if (scenario == "SA2"){
    pc_uptake = 1.
  }else{
    pc_uptake <- par_df$value[par_df$variable.name=="prob_test_order"]
  }
  
  # SA10: Cost of POCT test doubled:
  if (scenario == "SA11"){
    pc_test_cost <- 250
  }else{
    pc_test_cost <- par_df$value[par_df$variable.name=="poct_cost"]
  }
  pc_sens <- 1.
  pc_spec <- 1.
  
  # SA3: always follow test
  if (scenario == "SA3"){
    p_test_followed = 1.
  }else{
    p_test_followed = par_df$value[par_df$variable.name=="prob_test_followed"]
  }
  
  
  
  # SA1: proportion of DAPT
  if (scenario == "SA1"){
    p_ac_sc = .1
    p_at_sc = .45
    p_ap_sc = .45
  }else{
    p_ac_sc <- par_df$value[par_df$variable.name=="proportion_ac_standard"] # Proportion prescriped clopidogrel under standard care
    p_at_sc <- par_df$value[par_df$variable.name=="proportion_at_standard"]  # Proportion prescriped ticagrelor under standard care
    p_ap_sc <- par_df$value[par_df$variable.name=="proportion_ap_standard"]  # Proportion prescriped prasugrel under standard care
  }
  
  # Assign probabilities of AT and AP prescription following testing, assuming
  # same relative proportions as under standard care.
  p_at_A <- p_at_sc / (p_at_sc + p_ap_sc) # Proportion prescriped ticagrelor during subroutine A
  p_ap_A <- p_ap_sc / (p_at_sc + p_ap_sc) # Proportion prescriped prasugrel during subroutine A
  
  # Assign prevalences
  lof_prev_eu <- par_df$value[par_df$variable.name=="prevalence_lof_base_case"]
  lof_prev_as <- par_df$value[par_df$variable.name=="prevalence_lof_sensitivity"]
  
  # SA4: prevalence=56.8%
  if (scenario == "SA4"){
    lof_prev = lof_prev_as
  }else{
    lof_prev <- lof_prev_eu
  }
  
  # Get cost_pci, baseline cost applied in all courses
  cost_pci <- par_df$value[par_df$variable.name=="cost_pci"]
  
  # Get loading doses for drugs. Note that ac_lof and ac_no_lof are identical, but
  # for formatting purposes it's more convenient to use the same subpopulations as
  # the model.
  
  # Assign waiting period between initial admission and acting on PGX test result:
  test_waiting_period <- 0.
  
  ld_costs_sc <- data.frame(drug = c("ac_lof",
                                     "ac_no_lof",
                                     "at",
                                     "ap"),
                            cost = c(par_df$value[par_df$variable.name=="ld_clop_600"],
                                     par_df$value[par_df$variable.name=="ld_clop_600"],
                                     par_df$value[par_df$variable.name=="ld_tica_180"],
                                     par_df$value[par_df$variable.name=="ld_pras_60"]))
  
  ld_costs_pc <- data.frame(drug = c("ac_lof",
                                     "ac_no_lof",
                                     "at",
                                     "ap"),
                            cost = c(par_df$value[par_df$variable.name=="ld_clop_600"],
                                     par_df$value[par_df$variable.name=="ld_clop_600"] +
                                       test_waiting_period * par_df$value[par_df$variable.name=="day_cost_tica"],
                                     par_df$value[par_df$variable.name=="ld_tica_180"],
                                     par_df$value[par_df$variable.name=="ld_pras_60"]))
  
  # Get daily costs for drugs
  daily_costs <- data.frame(drug = c("ac_lof",
                                     "ac_no_lof",
                                     "at",
                                     "ap"),
                            cost = par_df$value[par_df$variable.name=="day_cost_asa"] +
                              c(par_df$value[par_df$variable.name=="day_cost_clop"],
                                par_df$value[par_df$variable.name=="day_cost_clop"],
                                par_df$value[par_df$variable.name=="day_cost_tica"],
                                par_df$value[par_df$variable.name=="day_cost_pras"]))
  
  # Get expected durations of courses by event
  course_dur_by_event_sc <- data.frame(event = c("no_event",
                                                 "stroke",
                                                 "mi",
                                                 "death"),
                                       duration = c(par_df$value[(par_df$variable.name=="md_no_event")|(par_df$variable.name=="md_328")],
                                                    par_df$value[(par_df$variable.name=="md_mi_stroke_sc")|(par_df$variable.name=="md_365_day")],
                                                    par_df$value[(par_df$variable.name=="md_mi_stroke_sc")|(par_df$variable.name=="md_365_day")],
                                                    AVE_TIME_TO_EVENT * par_df$value[(par_df$variable.name=="md_mi_stroke_sc")|(par_df$variable.name=="md_365_day")])) # Assume death occurs half way through course
  course_dur_by_event_pc <- data.frame(event = c("no_event",
                                                 "stroke",
                                                 "mi",
                                                 "death"),
                                       duration = c(par_df$value[(par_df$variable.name=="md_no_event")|(par_df$variable.name=="md_328")],
                                                    par_df$value[(par_df$variable.name=="md_mi_stroke_pc")|(par_df$variable.name=="md_365_day")],
                                                    par_df$value[(par_df$variable.name=="md_mi_stroke_pc")|(par_df$variable.name=="md_365_day")],
                                                    AVE_TIME_TO_EVENT * par_df$value[(par_df$variable.name=="md_mi_stroke_pc")|(par_df$variable.name=="md_365_day")])) # Assume death occurs half way through course
  
  # Now get course costs by combination of drug and event (bleeds do not affect
  # costs here)
  drug_course_costs_sc <-  merge(daily_costs, course_dur_by_event_sc) %>% # For now just do standard care - I don't see why times should be different under PC
    mutate(drug_event = paste(drug, event, sep="_"),
           cost = cost * duration) %>%
    select(-c(drug, event, duration))
  
  drug_course_costs_pc <-  merge(daily_costs, course_dur_by_event_sc) %>% # For now just do standard care - I don't see why times should be different under PC
    mutate(drug_event = paste(drug, event, sep="_"),
           cost = ifelse(drug=="ac_no_lof",
                         cost * (duration - test_waiting_period),
                         cost * duration)) %>%
    select(-c(drug, event, duration))
  
  # Utilities associated with events:
  event_utilities <- par_df %>%
    filter(grepl("utility", parameter.list)) %>%
    select(c(variable.name, value)) %>%
    mutate(variable.name = gsub("utility_", "", variable.name) %>%
             str_replace_all("_tree", "")) %>%
    filter(!grepl("_ac", variable.name)) %>% # Drop any derived drug-specific values
    filter(!grepl("_at", variable.name)) %>%
    filter(!grepl("_ap", variable.name)) %>%
    filter(!grepl("annual", variable.name)) %>% # Drop any annual-scale values
    add_row(variable.name = "death",
            value = 0) # Apply no_event utility before death, 0 afterwards
  
  # Add dyspnoea utilities, converting from whole-year to 30 day values, and
  # convert bleed decrements to utilities
  # NOTE: Utility decrements for bleed events are loaded in in units of years so
  # do not need to be converted.
  duration_dyspnoea <- par_df$value[
    which(par_df$variable.name=="duration_dyspnoea")] / 365
  event_utilities <- event_utilities %>%
    add_row(variable.name = "dyspnoea",
            value = 1 - duration_dyspnoea * event_utilities$value[
              which(event_utilities$variable.name=="u_dec_dyspnoea")]) %>%
    add_row(variable.name = "major_bleed",
            value = 1 - event_utilities$value[
              which(event_utilities$variable.name=="dec_major_bleed_*_(duration_major_bleed/365)")]) %>%
    add_row(variable.name = "minor_bleed",
            value = 1 - event_utilities$value[
              which(event_utilities$variable.name=="dec_minor_bleed_*_(duration_minor_bleed/365)")]) %>%
    filter(!grepl("dec_", variable.name))
  
  # Similar formula to get costs for DT model:
  event_costs <- par_df %>%
    filter(grepl("cost", parameter.list)) %>%
    filter(grepl("tree|bleed|dysp", parameter.list)) %>%
    filter(!grepl("_ac", variable.name)) %>% # Drop any derived drug-specific values
    filter(!grepl("_at", variable.name)) %>%
    filter(!grepl("_ap", variable.name)) %>%
    select(c(variable.name, value)) %>%
    mutate(variable.name = gsub("cost_", "", variable.name)) %>%
    mutate(variable.name = gsub("_tree", "", variable.name)) %>%
    mutate(variable.name = gsub("_pp", "", variable.name))
  
  # SA8: Minor bleeding costs set to GI bleed: £893
  if (scenario == "SA8"){
    event_costs$value[event_costs$variable.name=="minor_bleed"] <- 893
  }
  
  # SA9: Major bleeding costs (10%): £2502
  if (scenario == "SA9"){
    event_costs$value[event_costs$variable.name=="major_bleed"] <- 2502
  }
  
  # SA10: Major bleeding costs (40%): £2843
  if (scenario == "SA10"){
    event_costs$value[event_costs$variable.name=="major_bleed"] <- 2843
  }
  
  ### Parameters for Markov cohort model ####
  
  # Get standardised mortality ratios from PSA draw
  smr_df <- draw_df %>%
    filter(grepl("smr", variable.name))
  
  # Make copy containing only values - this distinction should be useful for PSA
  # purposes
  smr_vals <- smr_df %>%
    select(c(variable.name,
             Value)) %>%
    mutate(variable.name = variable.name %>%
             str_replace_all("smr_",
                             "") %>%
             str_replace_all("_further",
                             "")) %>%
    spread(variable.name, Value)
  
  # SA7: baseline SMR for ACS/reinfarction reduced by 2-%
  if (scenario == "SA7"){
    smr_vals$no_event = 0.8 * smr_vals$no_event
    smr_vals$mi = 0.8 * smr_vals$mi
    smr_vals$post_mi = 0.8 * smr_vals$post_mi
  }
  
  # Read in life table for healthy individuals
  mortality_prob_by_age <- mort_table %>%
    mutate(no_event = apply_smr(mortality_female,
                                mortality_male,
                                smr_vals$no_event,
                                prop_male),
           mi = apply_smr(mortality_female,
                          mortality_male,
                          smr_vals$mi,
                          prop_male),
           post_mi = apply_smr(mortality_female,
                               mortality_male,
                               smr_vals$post_mi,
                               prop_male),
           stroke = apply_smr(mortality_female,
                              mortality_male,
                              smr_vals$stroke,
                              prop_male),
           post_stroke = apply_smr(mortality_female,
                                   mortality_male,
                                   smr_vals$post_stroke,
                                   prop_male)) %>%
    select(-c(mortality_female,
              mortality_male,
              age_female,
              age_male))
  
  # Load utilities and add zeros for death  
  markov_utils <- rewrite_markov_utils(draw_df,
                                       time_hor,
                                       prop_male) %>%
    mutate(death = 0) %>%
    select(all_of(markov_states)) # Last step just makes sure ordering matches model
  markov_utils <- markov_utils[1:39, ] # Don't end up using last entry

  
  # Extract costs from main parameter table
  markov_costs <- par_df %>%
    filter(grepl("cost", parameter.list)) %>%
    filter(grepl("markov", parameter.list)) %>%
    select(c(variable.name, value)) %>%
    mutate(variable.name = gsub("cost_", "", variable.name)) %>%
    mutate(variable.name = gsub("_markov", "", variable.name)) %>%
    add_row(variable.name = "death",
            value = 0) %>% # Assign 0 cost to deaths in this section of model
    arrange(factor(variable.name, levels = markov_states))
  
  # Redefine simulation functions with new parameters
  {
    # Standard care
    implement_sc <- function(true_genotype,
                             sc_props = c(p_ac_sc,
                                          p_at_sc,
                                          p_ap_sc)){
      
      # If prescribed clopidogrel then subpopulation depends on genotype, otherwise
      # just used presciption proportions.
      p_ac_lof <- sc_props[1]
      p_ac_no_lof <- 0
      p_at <- sc_props[2]
      p_ap <- sc_props[3]
      
      sc_results <- data.frame(subpops = subpop_names,
                               prob = c(p_ac_lof,
                                        p_ac_no_lof,
                                        p_at,
                                        p_ap),
                               cost = c(ld_costs_sc$cost[ld_costs_sc$drug=="ac"],
                                        ld_costs_sc$cost[ld_costs_sc$drug=="ac"],
                                        ld_costs_sc$cost[ld_costs_sc$drug=="at"],
                                        ld_costs_sc$cost[ld_costs_sc$drug=="ap"]),
                               utility = c(1,
                                           1,
                                           1,
                                           1))
      return(sc_results)
    }
    
    implement_A <- function(true_genotype,
                            test_result
    ){
      
      # If test is followed, test says no lof, and true genotype is lof, then
      # patient goes to clopidogrel_lof
      p_ac_lof <- p_test_followed *
        ifelse(test_result=="no_lof", 1, 0) *  
        ifelse(true_genotype=="no_lof", 0, 1)
      
      # If test is followed, test says no lof, and true genotype is no_lof, then
      # patient goes to clopidogrel_no_lof
      p_ac_no_lof <- p_test_followed *
        ifelse(test_result=="no_lof", 1, 0) *  
        ifelse(true_genotype=="no_lof", 1, 0)
      
      # If test is followed and test says lof, and true genotype is lof, then
      # patient is prescribed one of the other drugs
      p_at <- p_test_followed *
        ifelse(test_result=="no_lof", 0, 1) * 
        p_at_A
      p_ap <- p_test_followed *
        ifelse(test_result=="no_lof", 0, 1) * 
        p_ap_A
      
      subroutine_A_results <- data.frame(subpop = c("ac_lof",
                                                    "ac_no_lof",
                                                    "at",
                                                    "ap",
                                                    "sc"),
                                         prob = c(p_ac_lof,
                                                  p_ac_no_lof,
                                                  p_at,
                                                  p_ap,
                                                  1 - p_test_followed), # Go to standard care if test not followed
                                         cost = c(ld_costs_pc$cost[ld_costs_pc$drug=="ac_lof"],
                                                  ld_costs_pc$cost[ld_costs_pc$drug=="ac_no_lof"],
                                                  ld_costs_pc$cost[ld_costs_pc$drug=="at"],
                                                  ld_costs_pc$cost[ld_costs_pc$drug=="ap"],
                                                  0),
                                         utility = c(1,
                                                     1,
                                                     1,
                                                     1,
                                                     1)) # costs assigned in standard care
      return(subroutine_A_results)
    }
    
    implement_B <- function(subpop_id,
                            test){
      
      subpop_pars <- prob_df[grepl(paste(subpop_id, "$", sep=""), # Gather subpopulation-specific parameters
                                   prob_df$variable.name),] %>%
        mutate(variable.name = variable.name %>% str_replace_all(paste("_", subpop_id, sep=""), ""))
      
      # Get probabilities of each event from dataframe:
      minor_bleed_prob <- subpop_pars$value[
        grep("minor_bleed", subpop_pars$variable.name)]
      major_bleed_prob <- subpop_pars$value[
        grep("major_bleed", subpop_pars$variable.name)]
      # For ac_no_lof, just use ac_lof value of dyspnoea probability since LOF gene does not affect this
      if (subpop_id == "ac_no_lof"){
        dysp_prob <- prob_df$value[which(
          prob_df$variable.name == "prob_dyspnoea_ac_lof"
        )]
      }else{
        dysp_prob <- subpop_pars$value[
          grep("dyspnoea", subpop_pars$variable.name)]
      }
      mi_prob <- subpop_pars$value[
        grep("mi$", subpop_pars$variable.name)]
      stroke_prob <- subpop_pars$value[
        grep("stroke", subpop_pars$variable.name)]
      death_prob <- subpop_pars$value[
        grep("death", subpop_pars$variable.name)]
      
      if (test=="pc"){
        drug_course_costs = drug_course_costs_pc
      }else{
        drug_course_costs = drug_course_costs_sc
      }
      
      mi_cost <- event_costs$value[event_costs$variable.name=="mi"] +
        drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                                   "mi",
                                                                   sep="_")]
      stroke_cost <- event_costs$value[event_costs$variable.name=="stroke"] +
        drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                                   "stroke",
                                                                   sep="_")]
      death_cost <- event_costs$value[event_costs$variable.name=="death"] +
        drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                                   "death",
                                                                   sep="_")]
      no_event_cost <- event_costs$value[event_costs$variable.name=="no_event"] +
        drug_course_costs$cost[drug_course_costs$drug_event==paste(subpop_id,
                                                                   "no_event",
                                                                   sep="_")]
      
      bleed_results <- data.frame(event = c("minor_bleed",
                                            "major_bleed",
                                            "no_event"),
                                  prob = c(minor_bleed_prob,
                                           major_bleed_prob,
                                           1 - minor_bleed_prob - major_bleed_prob),
                                  cost = c(event_costs$value[
                                    event_costs$variable.name=="minor_bleed"],
                                    event_costs$value[
                                      event_costs$variable.name=="major_bleed"],
                                    0),
                                  utility = c(event_utilities$value[
                                    event_utilities$variable.name=="minor_bleed"],
                                    event_utilities$value[
                                      event_utilities$variable.name=="major_bleed"],
                                    1))
      exp_cost_bleed <- sum(bleed_results$prob * bleed_results$cost)
      exp_util_bleed <- sum(bleed_results$prob * bleed_results$utility)
      
      dysp_results <- data.frame(event = c("dysp",
                                           "no_event"),
                                 prob = c(dysp_prob,
                                          1 - dysp_prob),
                                 cost = c(event_costs$value[
                                   event_costs$variable.name=="dyspnoea"],
                                   0),
                                 utility = c(event_utilities$value[
                                   event_utilities$variable.name=="dyspnoea"],
                                   1))
      exp_cost_dysp <- sum(dysp_results$prob * dysp_results$cost)
      exp_util_dysp <- sum(dysp_results$prob * dysp_results$utility)
      
      # Note on results: we currently assign zero cost to each of these outcomes
      subroutine_B_results <- data.frame(event = c("reinfarction",
                                                   "stroke",
                                                   "death",
                                                   "no_event"),
                                         prob = c(mi_prob,
                                                  stroke_prob,
                                                  death_prob,
                                                  1 - mi_prob - stroke_prob - death_prob),
                                         cost = cost_pci + exp_cost_bleed + exp_cost_dysp +
                                           c(mi_cost,
                                             stroke_cost,
                                             death_cost,
                                             no_event_cost),
                                         utility = -(1 - exp_util_bleed) + # Subtract utility decrement due to bleed
                                           -(1 - exp_util_dysp) + # Subtract utility decrement due to dyspnoea
                                           AVE_TIME_TO_EVENT * # Apply step correction; events happen part way through year
                                           event_utilities$value[event_utilities$variable.name == "no_event"] +
                                           (1-AVE_TIME_TO_EVENT) *
                                           c(event_utilities$value[event_utilities$variable.name == "mi"],
                                             event_utilities$value[event_utilities$variable.name == "stroke"],
                                             event_utilities$value[event_utilities$variable.name == "death"],
                                             event_utilities$value[event_utilities$variable.name == "no_event"]))
      return(subroutine_B_results)
    }
    
    # Patient steps through subroutine B regardless of testing method, so we
    # calculate it here
    B_result_list_sc <- lapply(subpop_names,
                               implement_B,
                               "sc")
    
    B_result_list_pc <- lapply(subpop_names,
                               implement_B,
                               "pc")
    
    # Make a dataframe to match up names of events in B with the Markov model states
    # patients end up in afterwards. Note leading underscore "_" used in string
    # concatentation later.
    post_event_match <- data.frame(event = c("reinfarction",
                                             "stroke",
                                             "death",
                                             "no_event"),
                                   state = c("_post_mi",
                                             "_post_stroke",
                                             "_death",
                                             "_no_event"))
    
    
    # Function to do entire decision tree
    # Test should be one of "pc", "l", or "sc" (standard care).
    run_forward <- function(test = "sc",
                            sc_props = c(p_ac_sc,
                                         p_at_sc,
                                         p_ap_sc),
                            poct_cost = pc_test_cost){
      patient_status <- data.frame(subpop = subpop_names,
                                   prob = rep(1, length(subpop_names)),
                                   exp_cost = rep(0, length(subpop_names)),
                                   exp_utility = rep(1, length(subpop_names)))
      
      # Always need to calculate standard care results as all routes have a chance
      # of going to it
      sc_results_lof <- implement_sc("lof", sc_props)
      sc_results_no_lof <- implement_sc("no_lof", sc_props)
      
      if (test == "sc"){
        patient_status$prob <- patient_status$prob *
          (lof_prev * sc_results_lof$prob +
             (1 - lof_prev) * sc_results_no_lof$prob)
        patient_status$exp_cost <- patient_status$exp_cost +
          lof_prev * sc_results_lof$prob * sc_results_lof$cost +
          (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost
      }
      
      # Only need results of subroutine A if we actually do testing. Also decide
      # here which version of B_result_list to use
      if((test == "pc") | (test == "l")){
        
        # Assign performance metrics for relevant test
        if (test == "pc"){
          sens <- pc_sens
          spec <- pc_spec
          test_uptake <- pc_uptake
        }
        if (test == "l"){
          sens <- l_sens
          spec <- l_spec
          test_uptake <- l_uptake
        }
        A_results_lof_lof <- implement_A(true_genotype = "lof",
                                         test_result = "lof")
        A_results_lof_no_lof <- implement_A(true_genotype = "lof",
                                            test_result = "no_lof")
        A_results_no_lof_lof <- implement_A(true_genotype = "no_lof",
                                            test_result = "lof")
        A_results_no_lof_no_lof <- implement_A(true_genotype = "no_lof",
                                               test_result = "no_lof")
        
        # Take weighted average of these based on expected frequency of genotype-
        # test result combinations. This ignores any potential variance which would
        # need to be incorporated into a fully stochastic model.
        A_results_ave <- data.frame(subpop = A_results_lof_lof$subpop,
                                    prob = rep(1, length(subpop_names) + 1),
                                    cost = rep(0, length(subpop_names) + 1),
                                    utility = rep(1, length(subpop_names) + 1))
        A_results_ave[, -1] <- lof_prev * spec * A_results_lof_lof[, -1] + # True positive
          lof_prev * (1 - spec) * A_results_lof_no_lof[, -1] + # False negative
          (1 - lof_prev) * (1 - sens) * A_results_no_lof_lof[, -1] + # False positive
          (1 - lof_prev) * sens * A_results_no_lof_no_lof[, -1] # True negative
        
        # Probability of going into standard care is probability of not testing plus
        # probability of testing and then ignoring. Note that this is a bit messy
        # because test uptake is currently implemented outside of the A subroutine
        # but whether the clinician follows the test is implemented within it.
        prob_sc <- (1 - test_uptake) + test_uptake * A_results_ave$prob[5]
        
        # Now add results to patient status
        patient_status$prob <- test_uptake *
          patient_status$prob *
          A_results_ave$prob[1:4]  +
          prob_sc * (lof_prev * sc_results_lof$prob + # This adds possibility of reverting to standard care
                       (1 - lof_prev) * sc_results_no_lof$prob)
        patient_status$exp_cost <- patient_status$exp_cost +
          test_uptake * poct_cost + # Account for cost of testing
          test_uptake * A_results_ave$cost[1:4] +
          prob_sc *
          (lof_prev * sc_results_lof$prob * sc_results_lof$cost +
             (1 - lof_prev) * sc_results_no_lof$prob * sc_results_no_lof$cost)
        B_result_list <- B_result_list_pc
      }else{
        # If we aren't testing we use the standard care version of B results
        B_result_list <- B_result_list_sc
      }
      
      # Now extend to combined subpopulation-health state space by implementing
      # the outcomes of subroutine B
      extended_patient_status <- data.frame(subpop = full_names,
                                            prob = rep(0, length(full_names)),
                                            exp_cost = rep(0, length(full_names)),
                                            exp_utility = rep(1, length(full_names)))
      for (i in (1:nrow(patient_status))){
        name <- patient_status$subpop[i]
        
        # Get total probability of being in subpopulation, to divide between subpop-
        # event combinations, and utility associated with subpopulation.
        p_subgroup <- patient_status$prob[which(patient_status$subpop==name)]
        c_subgroup <- patient_status$exp_cost[which(patient_status$subpop==name)]
        u_subgroup <- patient_status$exp_utility[which(patient_status$subpop==name)]
        
        # Safer to bring in costs for all Markov states, not just the ones we have
        # nonzero probability of starting the Markov chain in. For probability and
        # utility, I think we want to leave them at zero before implementing B.
        extended_patient_status$exp_cost[
          grep(paste(name, "_",sep=""),
               extended_patient_status$subpop)] <- c_subgroup
        
        # Assign probabilities and update utilities based on subroutine B.
        # Note: this assigns utility 1 to anyone in 
        for (j in (1:nrow(post_event_match))){
          # print(paste(name, post_event_match$state[j], sep=""))
          # print(which(grepl(paste(name, post_event_match$state[j], sep=""),
          #                  extended_patient_status$subpop)))
          extended_patient_status$prob[
            grepl(paste(name, post_event_match$state[j], sep=""),
                  extended_patient_status$subpop)] <-
            B_result_list[[i]]$prob[
              which(B_result_list[[i]]$event==post_event_match$event[j])] *
            p_subgroup
          # print(paste("ping", extended_patient_status$exp_cost[6]))
          extended_patient_status$exp_cost[
            grepl(paste(name, post_event_match$state[j], sep=""),
                  extended_patient_status$subpop)] <-
            extended_patient_status$exp_cost[
              grepl(paste(name, post_event_match$state[j], sep=""),
                    extended_patient_status$subpop)] +
            B_result_list[[i]]$cost[
              which(B_result_list[[i]]$event==post_event_match$event[j])]
          # print(extended_patient_status$exp_cost[6])
          extended_patient_status$exp_utility[
            grep(paste(name, post_event_match$state[j], sep=""),
                 extended_patient_status$subpop)] <-
            B_result_list[[i]]$utility[
              which(B_result_list[[i]]$event==post_event_match$event[j])] *
            u_subgroup
        }
      }
      return(extended_patient_status)
    }
    
    
    # Function for building Markov model
    build_markov_model <- function(tstep){
      
      transition_matrix <- matrix(0,
                                  nrow = n_states,
                                  ncol = n_states,
                                  dimnames = list(markov_states,
                                                  markov_states))
      
      # Get probabilities of each event from dataframe:
      
      # Outflow from no_event
      transition_matrix["no_event", "mi"] <- markov_pars$value[
        grep("mi", markov_pars$parameter.list)]
      transition_matrix["no_event", "stroke"] <- markov_pars$value[
        grep("^stroke", markov_pars$parameter.list)]
      transition_matrix["no_event", "death"] <- mortality_prob_by_age$no_event[tstep]
      
      # Outflow from mi
      transition_matrix["mi", "post_mi"] <- 1. - mortality_prob_by_age$mi[tstep]
      transition_matrix["mi", "death"] <- mortality_prob_by_age$mi[tstep]
      
      # Outflow from post_mi
      transition_matrix["post_mi", "death"] <- mortality_prob_by_age$post_mi[tstep]
      
      # Outflow from stroke
      transition_matrix["stroke", "post_stroke"] <- 1. - mortality_prob_by_age$stroke[tstep]
      transition_matrix["stroke", "death"] <- mortality_prob_by_age$stroke[tstep]
      
      # Outflow from post_stroke
      transition_matrix["post_stroke", "death"] <- mortality_prob_by_age$post_stroke[tstep]
      
      diag(transition_matrix) <- 1 - rowSums(transition_matrix)
      
      
      return(transition_matrix)
    }
    
    utils_from_markov_trace <- function(MT,
                                        value_by_state){
      MT_utils <- MT[-1,] %>%
        mutate(no_event =
                 no_event * value_by_state$no_event ) %>%
        mutate(stroke =
                 stroke * value_by_state$stroke ) %>%
        mutate(post_stroke =
                 post_stroke * value_by_state$post_stroke ) %>%
        mutate(mi =
                 mi * value_by_state$mi ) %>%
        mutate(post_mi =
                 post_mi * value_by_state$post_mi ) %>%
        mutate(death =
                 death * value_by_state$death ) %>%
        mutate(undiscounted_utility = no_event +
                 stroke +
                 post_stroke +
                 mi +
                 post_mi +
                 death) %>%
        mutate(discounted_utility = undiscounted_utility * discount_by_cycle) %>%
        mutate(no_event =
                 no_event * discount_by_cycle ) %>%
        mutate(stroke =
                 stroke * discount_by_cycle ) %>%
        mutate(post_stroke =
                 post_stroke * discount_by_cycle ) %>%
        mutate(mi =
                 mi * discount_by_cycle ) %>%
        mutate(post_mi =
                 post_mi * discount_by_cycle ) %>%
        mutate(death =
                 death * discount_by_cycle )
      halfstep_utils <- (0.5 * (MT[2:(time_hor), ] + rbind(MT[3:(time_hor), ], 0))) %>%
        mutate(no_event =
                 no_event * value_by_state$no_event ) %>%
        mutate(stroke =
                 stroke * value_by_state$stroke ) %>%
        mutate(post_stroke =
                 post_stroke * value_by_state$post_stroke ) %>%
        mutate(mi =
                 mi * value_by_state$mi ) %>%
        mutate(post_mi =
                 post_mi * value_by_state$post_mi ) %>%
        mutate(death =
                 death * value_by_state$death ) %>%
        mutate(undiscounted_utility = no_event +
                 stroke +
                 post_stroke +
                 mi +
                 post_mi +
                 death) %>%
        mutate(discounted_utility = undiscounted_utility * discount_by_cycle) %>%
        mutate(no_event =
                 no_event * discount_by_cycle ) %>%
        mutate(stroke =
                 stroke * discount_by_cycle ) %>%
        mutate(post_stroke =
                 post_stroke * discount_by_cycle ) %>%
        mutate(mi =
                 mi * discount_by_cycle ) %>%
        mutate(post_mi =
                 post_mi * discount_by_cycle ) %>%
        mutate(death =
                 death * discount_by_cycle )
      MT_utils$halfstep <- halfstep_utils$undiscounted_utility
      MT_utils$discounted_halfstep <- halfstep_utils$discounted_utility
      return(MT_utils)
    }
    
    costs_from_markov_trace <- function(MT,
                                        value_by_state){
      MT_costs <- MT[-1,] %>%
        mutate(no_event =
                 no_event * value_by_state$value[
                   value_by_state$variable.name=="no_event"] ) %>%
        mutate(stroke =
                 stroke * value_by_state$value[
                   value_by_state$variable.name=="stroke"] ) %>%
        mutate(post_stroke =
                 post_stroke * value_by_state$value[
                   value_by_state$variable.name=="post_stroke"] ) %>%
        mutate(mi =
                 mi * value_by_state$value[
                   value_by_state$variable.name=="mi"] ) %>%
        mutate(post_mi =
                 post_mi * value_by_state$value[
                   value_by_state$variable.name=="post_mi"] ) %>%
        mutate(death =
                 death * value_by_state$value[
                   value_by_state$variable.name=="death"] ) %>%
        mutate(undiscounted_cost = no_event +
                 stroke +
                 post_stroke +
                 mi +
                 post_mi +
                 death) %>%
        mutate(discounted_cost = undiscounted_cost * discount_by_cycle) %>%
        mutate(no_event =
                 no_event * discount_by_cycle ) %>%
        mutate(stroke =
                 stroke * discount_by_cycle ) %>%
        mutate(post_stroke =
                 post_stroke * discount_by_cycle ) %>%
        mutate(mi =
                 mi * discount_by_cycle ) %>%
        mutate(post_mi =
                 post_mi * discount_by_cycle ) %>%
        mutate(death =
                 death * discount_by_cycle )
      MT_costs$halfstep <- 0.5 * (MT_costs$undiscounted_cost[1:39] +
                                    c(MT_costs$undiscounted_cost[2:39], 0))
      MT_costs$discounted_halfstep <- MT_costs$halfstep * discount_by_cycle
      return(MT_costs)
    }
    }
  
  # Run decision tree analysis
  dt_results_pc <- run_forward("pc") %>%
    mutate(event = subpop %>%
             str_remove_all(paste(paste(subpop_names, collapse = "_|"),
                                  "_",
                                  sep = "")))
  
  # Expected costs and utilities at DT stage:
  dt_pc_cost <- sum(dt_results_pc$prob * dt_results_pc$exp_cost)
  dt_pc_util <- sum(dt_results_pc$prob * dt_results_pc$exp_utility)
  
  P0_pc <- sapply(markov_states,
                  FUN=function(event){
                    sum(dt_results_pc$prob[dt_results_pc$event == event])}
  )
  
  
  MT_pc <- sapply(Reduce("%*%", lapply(-1:(n_tsteps-1),
                                       FUN = function(t){
                                         if (t==-1){
                                           return(diag(nrow = n_states))
                                         }else{
                                           if (t==0){
                                             return(build_markov_model(1))
                                           }else{
                                             return(build_markov_model(t+1))
                                           }
                                         }
                                       }), accumulate = TRUE),
                  FUN = function(A){
                    P0_pc %*% A
                  }
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(markov_states) %>%
    mutate(time_step = 0:n_tsteps)
  
  
  # Calculate expected utility over time:
  utility_pc <- utils_from_markov_trace(MT_pc,
                                        markov_utils)
  # Expected costs:
  MC_costs_pc <- costs_from_markov_trace(MT_pc,
                                         markov_costs)
  
  
  # Mean life years:
  life_years_pc <- sum(1 - c(MT_pc$death[1],
                             0.5 * (MT_pc$death[2:39] + MT_pc$death[3:40])))
  
  # Now do sc
  
  # Run decision tree analysis
  dt_results_sc <- run_forward("sc") %>%
    mutate(event = subpop %>%
             str_remove_all(paste(paste(subpop_names, collapse = "_|"),
                                  "_",
                                  sep = "")))
  
  # Expected costs at DT stage:
  dt_sc_cost <- sum(dt_results_sc$prob * dt_results_sc$exp_cost)
  dt_sc_util <- sum(dt_results_sc$prob * dt_results_sc$exp_utility)
  
  P0_sc <- sapply(markov_states,
                  FUN=function(event){
                    sum(dt_results_sc$prob[dt_results_sc$event == event])}
  )
  
  MT_sc <- sapply(Reduce("%*%", lapply(-1:(n_tsteps-1),
                                       FUN = function(t){
                                         if (t==-1){
                                           return(diag(nrow = n_states))
                                         }else{
                                           if (t==0){
                                             return(build_markov_model(1))
                                           }else{
                                             return(build_markov_model(t+1))
                                           }
                                         }
                                       }), accumulate = TRUE),
                  FUN = function(A){
                    P0_sc %*% A
                  }
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(markov_states) %>%
    mutate(time_step = 0:n_tsteps)
  
  
  # Calculate expected utility over time:
  
  utility_sc <- utils_from_markov_trace(MT_sc,
                                        markov_utils)
  
  # Expected costs:
  MC_costs_sc <- costs_from_markov_trace(MT_sc,
                                         markov_costs)
  
  # Mean life years:
  life_years_sc <- sum(1 - c(MT_sc$death[1],
                             0.5 * (MT_sc$death[2:39] + MT_sc$death[3:40])))
  
  # Now calculate ICER
  ICER_undisc <- ((dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))) /
    ((dt_pc_util + sum(utility_pc$undiscounted_utility)) - (dt_sc_util + sum(utility_sc$undiscounted_utility)))
  ICER_disc <- ((dt_pc_cost + sum(MC_costs_pc$discounted_cost)) - (dt_sc_cost + sum(MC_costs_sc$discounted_cost))) /
    ((dt_pc_util + sum(utility_pc$discounted_utility)) - (dt_sc_util + sum(utility_sc$discounted_utility)))
  ICER_undisc_hs <- ((dt_pc_cost + sum(MC_costs_pc$halfstep)) - (dt_sc_cost + sum(MC_costs_sc$halfstep))) /
    ((dt_pc_util + sum(utility_pc$halfstep)) - (dt_sc_util + sum(utility_sc$halfstep)))
  ICER_disc_hs <- ((dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) - (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))) /
    ((dt_pc_util + sum(utility_pc$discounted_halfstep)) - (dt_sc_util + sum(utility_sc$discounted_halfstep)))
  
  # Create a dataframe storing outputs by case:
  arm_comparison <- data.frame(arm = c("SC", "PC", "Increment"),
                               life_years = c(life_years_sc,
                                              life_years_pc,
                                              life_years_pc - life_years_sc),
                               utility_udc = c(dt_sc_util + sum(utility_sc$undiscounted_utility),
                                               dt_pc_util + sum(utility_pc$undiscounted_utility),
                                               (dt_pc_util + sum(utility_pc$undiscounted_utility)) -
                                                 (dt_sc_util + sum(utility_sc$undiscounted_utility))),
                               cost_udc = c(dt_sc_cost + sum(MC_costs_sc$undiscounted_cost),
                                            dt_pc_cost + sum(MC_costs_pc$undiscounted_cost),
                                            (dt_pc_cost + sum(MC_costs_pc$undiscounted_cost)) -
                                              (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))),
                               ratio_udc = c(NA,
                                             NA,
                                             ICER_undisc),
                               utility_dc = c(dt_sc_util + sum(utility_sc$discounted_utility),
                                              dt_pc_util + sum(utility_pc$discounted_utility),
                                              (dt_pc_util + sum(utility_pc$discounted_utility)) -
                                                (dt_sc_util + sum(utility_sc$discounted_utility))),
                               cost_dc = c(dt_sc_cost + sum(MC_costs_sc$discounted_cost),
                                           dt_pc_cost + sum(MC_costs_pc$discounted_cost),
                                           (dt_pc_cost + sum(MC_costs_pc$discounted_cost)) -
                                             (dt_sc_cost + sum(MC_costs_sc$discounted_cost))),
                               ratio_dc = c(NA,
                                            NA,
                                            ICER_disc),
                               utility_udc_hs = c(dt_sc_util + sum(utility_sc$halfstep),
                                                  dt_pc_util + sum(utility_pc$halfstep),
                                                  (dt_pc_util + sum(utility_pc$halfstep)) -
                                                    (dt_sc_util + sum(utility_sc$halfstep))),
                               cost_udc_hs = c(dt_sc_cost + sum(MC_costs_sc$halfstep),
                                               dt_pc_cost + sum(MC_costs_pc$halfstep),
                                               (dt_pc_cost + sum(MC_costs_pc$halfstep)) -
                                                 (dt_sc_cost + sum(MC_costs_sc$halfstep))),
                               ratio_udc_hs = c(NA,
                                                NA,
                                                ICER_undisc_hs),
                               utility_dc_hs = c(dt_sc_util + sum(utility_sc$discounted_halfstep),
                                                 dt_pc_util + sum(utility_pc$discounted_halfstep),
                                                 (dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                                   (dt_sc_util + sum(utility_sc$discounted_halfstep))),
                               cost_dc_hs = c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep),
                                              dt_pc_cost + sum(MC_costs_pc$discounted_halfstep),
                                              (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                                (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))),
                               ratio_dc_hs = c(NA,
                                               NA,
                                               ICER_disc_hs),
                               nmb_per_capita = 20000 * c(dt_sc_util + sum(utility_sc$discounted_halfstep),
                                                          dt_pc_util + sum(utility_pc$discounted_halfstep),
                                                          (dt_pc_util + sum(utility_pc$discounted_halfstep)) -
                                                            (dt_sc_util + sum(utility_sc$discounted_halfstep))) -
                                 c(dt_sc_cost + sum(MC_costs_sc$discounted_halfstep),
                                   dt_pc_cost + sum(MC_costs_pc$discounted_halfstep),
                                   (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep)) -
                                     (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))))
  
  # Mean life years:
  life_years_sc <- sum(1 - 0.5 * (MT_sc$death[1:39] + (MT_sc$death[2:40])))
  
  # Build output lists and calculate incremental values
  sc_cost_udc <- (dt_sc_cost + sum(MC_costs_sc$undiscounted_cost))
  sc_util_udc <- (dt_sc_util + sum(utility_sc$undiscounted_utility))
  sc_cost_dc <- (dt_sc_cost + sum(MC_costs_sc$discounted_cost))
  sc_util_dc <- (dt_sc_util + sum(utility_sc$discounted_utility))
  sc_cost_udc_hs <- (dt_sc_cost + sum(MC_costs_sc$halfstep))
  sc_util_udc_hs <- (dt_sc_util + sum(utility_sc$halfstep))
  sc_cost_dc_hs <- (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))
  sc_util_dc_hs <- (dt_sc_util + sum(utility_sc$discounted_halfstep))
  sc_nmb <- 20000 * (dt_sc_util + sum(utility_sc$discounted_halfstep)) -
    (dt_sc_cost + sum(MC_costs_sc$discounted_halfstep))
  
  
  pc_cost_udc <- (dt_pc_cost + sum(MC_costs_pc$undiscounted_cost))
  pc_util_udc <- (dt_pc_util + sum(utility_pc$undiscounted_utility))
  pc_cost_dc <- (dt_pc_cost + sum(MC_costs_pc$discounted_cost))
  pc_util_dc <- (dt_pc_util + sum(utility_pc$discounted_utility))
  pc_cost_udc_hs <- (dt_pc_cost + sum(MC_costs_pc$halfstep))
  pc_util_udc_hs <- (dt_pc_util + sum(utility_pc$halfstep))
  pc_cost_dc_hs <- (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep))
  pc_util_dc_hs <- (dt_pc_util + sum(utility_pc$discounted_halfstep))
  pc_nmb <- 20000 * (dt_pc_util + sum(utility_pc$discounted_halfstep)) -
    (dt_pc_cost + sum(MC_costs_pc$discounted_halfstep))
  
  life_years_inc <- life_years_pc - life_years_sc
  inc_util_udc <- pc_util_udc - sc_util_udc
  inc_cost_udc <- pc_cost_udc - sc_cost_udc
  inc_util_dc <- pc_util_dc - sc_util_dc
  inc_cost_dc <- pc_cost_dc - sc_cost_dc
  inc_util_udc_hs <- pc_util_udc_hs - sc_util_udc_hs
  inc_cost_udc_hs <- pc_cost_udc_hs - sc_cost_udc_hs
  inc_util_dc_hs <- pc_util_dc_hs - sc_util_dc_hs
  inc_cost_dc_hs <- pc_cost_dc_hs - sc_cost_dc_hs
  inc_nmb <- pc_nmb - sc_nmb
  
  # Create a one-line dataframe containing outputs under different discounting assumptions:
  outcome_df <- data.frame(life_years_sc,
                           sc_util_udc,
                           sc_cost_udc,
                           sc_util_dc,
                           sc_cost_dc,
                           sc_util_udc_hs,
                           sc_cost_udc_hs,
                           sc_util_dc_hs,
                           sc_cost_dc_hs,
                           sc_nmb,
                           life_years_pc,
                           pc_util_udc,
                           pc_cost_udc,
                           pc_util_dc,
                           pc_cost_dc,
                           pc_util_udc_hs,
                           pc_cost_udc_hs,
                           pc_util_dc_hs,
                           pc_cost_dc_hs,
                           pc_nmb,
                           life_years_inc,
                           inc_util_udc,
                           inc_cost_udc,
                           inc_util_dc,
                           inc_cost_dc,
                           inc_util_udc_hs,
                           inc_cost_udc_hs,
                           inc_util_dc_hs,
                           inc_cost_dc_hs,
                           inc_nmb)
  
  # Extra function for returning number of bleed events
  get_bleed_counts <- function(subpop_id,
                               test){
    
    subpop_pars <- prob_df[grepl(paste(subpop_id, "$", sep=""), # Gather subpopulation-specific parameters
                                 prob_df$variable.name),] %>%
      mutate(variable.name = variable.name %>% str_replace_all(paste("_", subpop_id, sep=""), ""))
    
    # Get probabilities of each event from dataframe:
    minor_bleed_prob <- subpop_pars$value[
      grep("minor_bleed", subpop_pars$variable.name)]
    major_bleed_prob <- subpop_pars$value[
      grep("major_bleed", subpop_pars$variable.name)]
    
    bleed_df <- data.frame(subpop = subpop_id,
                           minor_bleed = minor_bleed_prob,
                           major_bleed = major_bleed_prob)
    return(bleed_df)
  }
  
  bleeds_by_subpop_sc <- lapply(subpop_names,
                                get_bleed_counts,
                                "sc") %>%
    bind_rows()
  bleeds_by_subpop_pc <- lapply(subpop_names,
                                get_bleed_counts,
                                "pc") %>%
    bind_rows
  
  aggregate_bleed_results <- function(test = "sc",
                                      sc_props = c(p_ac_sc,
                                                   p_at_sc,
                                                   p_ap_sc)){
    patient_status <- data.frame(subpop = subpop_names,
                                 prob = rep(1, length(subpop_names)))
    
    # Always need to calculate standard care results as all routes have a chance
    # of going to it
    sc_results_lof <- implement_sc("lof", sc_props)
    sc_results_no_lof <- implement_sc("no_lof", sc_props)
    
    if (test == "sc"){
      patient_status$prob <- patient_status$prob *
        (lof_prev * sc_results_lof$prob +
           (1 - lof_prev) * sc_results_no_lof$prob)
    }
    
    # Only need results of subroutine A if we actually do testing. Also decide
    # here which version of B_result_list to use
    if((test == "pc") | (test == "l")){
      
      # Assign performance metrics for relevant test
      if (test == "pc"){
        sens <- pc_sens
        spec <- pc_spec
        test_uptake <- pc_uptake
      }
      if (test == "l"){
        sens <- l_sens
        spec <- l_spec
        test_uptake <- l_uptake
      }
      A_results_lof_lof <- implement_A(true_genotype = "lof",
                                       test_result = "lof")
      A_results_lof_no_lof <- implement_A(true_genotype = "lof",
                                          test_result = "no_lof")
      A_results_no_lof_lof <- implement_A(true_genotype = "no_lof",
                                          test_result = "lof")
      A_results_no_lof_no_lof <- implement_A(true_genotype = "no_lof",
                                             test_result = "no_lof")
      
      # Take weighted average of these based on expected frequency of genotype-
      # test result combinations. This ignores any potential variance which would
      # need to be incorporated into a fully stochastic model.
      A_results_ave <- data.frame(subpop = A_results_lof_lof$subpop,
                                  prob = rep(1, length(subpop_names) + 1))
      A_results_ave$prob <- lof_prev * spec * A_results_lof_lof$prob + # True positive
        lof_prev * (1 - spec) * A_results_lof_no_lof$prob + # False negative
        (1 - lof_prev) * (1 - sens) * A_results_no_lof_lof$prob + # False positive
        (1 - lof_prev) * sens * A_results_no_lof_no_lof$prob # True negative
      
      # Probability of going into standard care is probability of not testing plus
      # probability of testing and then ignoring. Note that this is a bit messy
      # because test uptake is currently implemented outside of the A subroutine
      # but whether the clinician follows the test is implemented within it.
      prob_sc <- (1 - test_uptake) + test_uptake * A_results_ave$prob[5]
      
      # Now add results to patient status
      patient_status$prob <- test_uptake *
        patient_status$prob *
        A_results_ave$prob[1:4]  +
        prob_sc * (lof_prev * sc_results_lof$prob + # This adds possibility of reverting to standard care
                     (1 - lof_prev) * sc_results_no_lof$prob)
      bleeds_by_subpop <- bleeds_by_subpop_pc
    }else{
      # If we aren't testing we use the standard care version of B results
      bleeds_by_subpop <- bleeds_by_subpop_sc
    }
    
    bleed_df <- data.frame(minor = sum(bleeds_by_subpop$minor * patient_status$prob),
                           major = sum(bleeds_by_subpop$major * patient_status$prob))
    
    return(bleed_df)
  }
  
  bleed_results_sc <- aggregate_bleed_results("sc")
  bleed_results_pc <- aggregate_bleed_results("pc")
  
  # Count events per 1000:
  events_df <- data.frame(arm = c("sc",
                                  "pc"),
                          mi_dt = 1000 * c(sum(dt_results_sc$prob[grepl("mi", dt_results_sc$event)]),
                                           sum(dt_results_pc$prob[grepl("mi", dt_results_pc$event)])),
                          mi_mc = 1000 * c(sum(MT_sc$mi),
                                           sum(MT_pc$mi)),
                          stroke_dt = 1000 * c(sum(dt_results_sc$prob[grepl("stroke", dt_results_sc$event)]),
                                               sum(dt_results_pc$prob[grepl("stroke", dt_results_pc$event)])),
                          stroke_mc = 1000 * c(sum(MT_sc$stroke),
                                               sum(MT_pc$stroke)),
                          major_bleed_dt = 1000 * c(bleed_results_sc$major,
                                                    bleed_results_pc$major),
                          minor_bleed_dt = 1000 * c(bleed_results_sc$minor,
                                                    bleed_results_pc$minor),
                          death_dt = 1000 * c(sum(dt_results_sc$prob[grepl("death", dt_results_sc$event)]),
                                              sum(dt_results_pc$prob[grepl("death", dt_results_pc$event)])),
                          death_mc = 1000 * c(MT_sc$death[nrow(MT_sc)],
                                              MT_pc$death[nrow(MT_pc)]))
  
  return(list(outcome_df,
              arm_comparison,
              events_df,
              mortality_prob_by_age,
              markov_utils,
              markov_costs,
              dt_results_sc,
              dt_results_pc,
              event_costs,
              event_utilities,
              MT_sc,
              MT_pc))
}