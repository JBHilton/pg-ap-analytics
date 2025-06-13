# This script contains functions needed for carrying out probabilistic
# sensitivity analysis.

do_PSA_draw <- function(uc_df){
  uc_df$draw <- 0
  uc_df$draw[uc_df$distribution=="Beta"] <-
    rbeta(length(which(uc_df$distribution=="Beta")),
          uc_df$par1[uc_df$distribution=="Beta"],
          uc_df$par2[uc_df$distribution=="Beta"])
  uc_df$draw[uc_df$distribution=="Gamma"] <-
    rgamma(length(which(uc_df$distribution=="Gamma")),
           shape = uc_df$par1[uc_df$distribution=="Gamma"],
           scale = uc_df$par2[uc_df$distribution=="Gamma"])
  uc_df$draw[uc_df$distribution=="lognormal"] <-
    rlnorm(length(which(uc_df$distribution=="lognormal")),
           uc_df$par1[uc_df$distribution=="lognormal"],
           uc_df$par2[uc_df$distribution=="lognormal"])
  return(uc_df)
}

rewrite_pars_from_draw <- function(par_df, draw_df){
  par_df$value[which(par_df$variable.name %in% draw_df$variable.name)] <-
    draw_df$draw[which(draw_df$variable.name %in% par_df$variable.name)]
  markov_pars <- data.frame(parameter.list = c("mi",
                                               "stroke"),
                            value = c(draw_df$draw[draw_df$variable.name=="prob_nevent_to_rinfarc"],
                                      draw_df$draw[draw_df$variable.name=="prob_nevent_to_rinfarc"]))
  return(list(par_df, markov_pars))
}
