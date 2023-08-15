################################
##
## Generating uncertainty parameters
##
#################################
library(malariasimulation)

setwd("C:/Users/esherrar/Documents/Rprojects/uganda_covid_impact")

# mosquito nets

## Now run Mosq_Net_efficacy (Net Parameters.R to produce matrix)
## This is for 2014 nets (ALL PYRETHROIDS)
matrix_dn0 = as.data.frame(matrix_dn0)
matrix_rn0 = as.data.frame(matrix_rn0)
matrix_halflife = as.data.frame(matrix_halflife)

matrix_dn0$res = seq(0,1,by=0.01)

## Similarly for PBOs
matrix_dn0_pbo = as.data.frame(matrix_dn0_pbo)
matrix_rn0_pbo = as.data.frame(matrix_rn0_pbo)
matrix_halflife_pbo = as.data.frame(matrix_halflife_pbo)
matrix_dn0_pbo$res = seq(0,1,by=0.01)

dt2 = read.csv("raw data/dt2.csv",header=TRUE)
## determine the resistance level to explore
which(matrix_dn0$res == dt2$res_1[1])

## sample 20 random draws
cols_to_sample = sample(1:1000, size = 20, replace = FALSE)

## confirm these are representative
matrix_dn0[which(matrix_dn0$res == dt2$res_1[1]),cols_to_sample]
hist(as.numeric(c(matrix_dn0[which(matrix_dn0$res == dt2$res_1[1]),])))
hist(as.numeric(c(matrix_rn0[which(matrix_dn0$res == dt2$res_1[1]),cols_to_sample])))

which(dt2$net_for_2017_trial == "Standard")

uncertainty_1 = list()

for(i in 1:nrow(dt2)){
  
  parms1 = expand.grid(
    dn0_2014 =         as.numeric(c(matrix_dn0[which(round(matrix_dn0$res,2) == round(dt2$res_1[i],2)),cols_to_sample])))
  parms1$rn_2014 =     as.numeric(c(matrix_rn0[which(round(matrix_dn0$res,2) == round(dt2$res_1[i],2)),cols_to_sample]))
  parms1$gamman_2014 = as.numeric(c(matrix_halflife[which(round(matrix_dn0$res,2) == round(dt2$res_1[i],2)),cols_to_sample]))
  parms1$itn_use_baseline = rbinom(n = 20,prob = dt2$ITN_Use_prior[i],size = 1000)/1000 
  parms1$ITN_Use_A = rbinom(n = 20,prob = dt2$ITN_Use_A[i],size = 1000)/1000
  parms1$itn_leave_dur =  rnorm(n = 20, mean = dt2$itn_leave_durMean[i], sd = 1)
  parms1$itn_leave_dur = ifelse(parms1$itn_leave_dur < 0, 1,parms1$itn_leave_dur)

  if (dt2$net_for_2017_trial[i] == "Standard") { 
    parms1$dn0_2017 = as.numeric(c(matrix_dn0[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample])) 
  } else if (dt2$net_for_2017_trial[i] == "PBO") {
    parms1$dn0_2017 = as.numeric(c(matrix_dn0_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample]))
  }
  
  if (dt2$net_for_2017_trial[i] == "Standard") { 
    parms1$rn_2017 = as.numeric(c(matrix_rn0[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample])) 
  } else if (dt2$net_for_2017_trial[i] == "PBO"){
    parms1$rn_2017 = as.numeric(c(matrix_rn0_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample]))
  }
  
  if (dt2$net_for_2017_trial[i] == "Standard") { 
    parms1$gamman_2017 = as.numeric(c(matrix_halflife[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample])) 
  } else if (dt2$net_for_2017_trial[i] == "PBO"){
    parms1$gamman_2017 = as.numeric(c(matrix_halflife_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_2[i],2)),cols_to_sample]))
  }
  
  ## And 2020 nets
  if (dt2$net_type[i] == "Standard") { 
    parms1$dn0_2020 = as.numeric(c(matrix_dn0[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample])) 
  } else if (dt2$net_type[i] == "PBO") {
    parms1$dn0_2020 = as.numeric(c(matrix_dn0_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample]))
  }
  
  if (dt2$net_type[i] == "Standard") { 
    parms1$rn_2020 = as.numeric(c(matrix_rn0[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample])) 
  } else if (dt2$net_type[i] == "PBO"){
    parms1$rn_2020 = as.numeric(c(matrix_rn0_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample]))
  }
  
  if (dt2$net_type[i] == "Standard") { 
    parms1$gamman_2020 = as.numeric(c(matrix_halflife[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample])) 
  } else if (dt2$net_type[i] == "PBO"){
    parms1$gamman_2020 = as.numeric(c(matrix_halflife_pbo[which(round(matrix_dn0$res,2) == round(dt2$res_3[i],2)),cols_to_sample]))
  }
  
  uncertainty_1[[i]] = parms1  
}




## now we want the function to draw from the uncertainty parameters
malsim_actual_f = function(test_data, ## global dataset (dt2)
                           dat_row,   ## row from dt2
                           uncertainty, ## uncertainty around each row in dt2
                           uncert_row){  ## 1 - 20 uncertainty draws
  
  year <- 365
  month <- 30
  sim_length <- 11 * year ## initially just til 2020 as need extra resistance data!
  ## This is spanning Jan 2014 - Dec 2025
  ## Assume that all places have received nets at start of 2014
  ## then relative to jan 2017 (see uganda_what_if...)
  ## then as noted for 2020
  
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 10000
  starting_EIR <- test_data$eir_preds[dat_row]
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      
      prevalence_rendering_min_ages = 2 * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = 10 * 365,
      
      clinical_incidence_rendering_min_ages = 0 * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = 100 * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      g0 = test_data$seasonal_a0[dat_row],
      g = c(test_data$seasonal_a1[dat_row], test_data$seasonal_a2[dat_row], test_data$seasonal_a3[dat_row]),
      h = c(test_data$seasonal_b1[dat_row], test_data$seasonal_b2[dat_row], test_data$seasonal_b3[dat_row]),
      
      individual_mosquitoes = FALSE ## Update next
    )
  )
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set species
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           
                           proportions=c(test_data$prop_gam[dat_row],
                                         test_data$prop_ara[dat_row],
                                         1-test_data$prop_gam[dat_row]-test_data$prop_ara[dat_row]))
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(100),
                                      coverage=c(0.1424563))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(100),
                                      coverage=c(0.409862))  
  
  ## Set up bed nets
  
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = c(0, 3, 6) * year + c(0,
                                     test_data$days_after_jan_2017[dat_row],
                                     test_data$net_delivered_2020[dat_row]),
    name=c("background", 
           "trial_nets",
           "2020_nets")
  )
  
  
  # Will run over a loop so this is for either net
  # this is using what was done i.e. actually used
  bednetparams_2 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(uncertainty[[dat_row]]$itn_use_baseline[uncert_row], ## historic prior to RCT
                  uncertainty[[dat_row]]$ITN_Use_A[uncert_row], ## historic during RCT (on deployment) 
                  uncertainty[[dat_row]]$ITN_Use_A[uncert_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = uncertainty[[dat_row]]$itn_leave_dur[uncert_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = t(matrix(as.numeric(c(uncertainty[[dat_row]]$dn0_2014[uncert_row], uncertainty[[dat_row]]$dn0_2014[uncert_row], uncertainty[[dat_row]]$dn0_2014[uncert_row],
                                uncertainty[[dat_row]]$dn0_2017[uncert_row], uncertainty[[dat_row]]$dn0_2017[uncert_row], uncertainty[[dat_row]]$dn0_2017[uncert_row],
                                uncertainty[[dat_row]]$dn0_2020[uncert_row], uncertainty[[dat_row]]$dn0_2020[uncert_row], uncertainty[[dat_row]]$dn0_2020[uncert_row])), nrow=3, ncol=3)),
    rn = t(matrix(as.numeric(c(uncertainty[[dat_row]]$rn_2014[uncert_row], uncertainty[[dat_row]]$rn_2014[uncert_row], uncertainty[[dat_row]]$rn_2014[uncert_row],
                               uncertainty[[dat_row]]$rn_2017[uncert_row], uncertainty[[dat_row]]$rn_2017[uncert_row], uncertainty[[dat_row]]$rn_2017[uncert_row],
                               uncertainty[[dat_row]]$rn_2020[uncert_row], uncertainty[[dat_row]]$rn_2020[uncert_row], uncertainty[[dat_row]]$rn_2020[uncert_row])), nrow=3, ncol=3)),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = as.numeric(c(uncertainty[[dat_row]]$gamman_2014[uncert_row], uncertainty[[dat_row]]$gamman_2017[uncert_row], uncertainty[[dat_row]]$gamman_2020[uncert_row])) * 365
  )
  
  
  ##scheduled
  ##including the 2020 just for saving time later
  bednet_events_scheduled = data.frame(
    timestep = c(0, 3, 6) * year + c(0,
                                     test_data$days_after_jan_2017[dat_row],
                                     test_data$net_schedule_2020[dat_row]),
    name=c("background", 
           "trial_nets",
           "2020_nets")
  )
  
  bednetparams_2b <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events_scheduled$timestep,
    
    coverages = c(uncertainty[[dat_row]]$itn_use_baseline[uncert_row], ## historic prior to RCT
                  uncertainty[[dat_row]]$ITN_Use_A[uncert_row], ## historic during RCT (on deployment) 
                  uncertainty[[dat_row]]$ITN_Use_A[uncert_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = uncertainty[[dat_row]]$itn_leave_dur[uncert_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = t(matrix(as.numeric(c(uncertainty[[dat_row]]$dn0_2014[uncert_row], uncertainty[[dat_row]]$dn0_2014[uncert_row], uncertainty[[dat_row]]$dn0_2014[uncert_row],
                                uncertainty[[dat_row]]$dn0_2017[uncert_row], uncertainty[[dat_row]]$dn0_2017[uncert_row], uncertainty[[dat_row]]$dn0_2017[uncert_row],
                                uncertainty[[dat_row]]$dn0_2020[uncert_row], uncertainty[[dat_row]]$dn0_2020[uncert_row], uncertainty[[dat_row]]$dn0_2020[uncert_row])), nrow=3, ncol=3)),
    rn = t(matrix(as.numeric(c(uncertainty[[dat_row]]$rn_2014[uncert_row], uncertainty[[dat_row]]$rn_2014[uncert_row], uncertainty[[dat_row]]$rn_2014[uncert_row],
                               uncertainty[[dat_row]]$rn_2017[uncert_row], uncertainty[[dat_row]]$rn_2017[uncert_row], uncertainty[[dat_row]]$rn_2017[uncert_row],
                               uncertainty[[dat_row]]$rn_2020[uncert_row], uncertainty[[dat_row]]$rn_2020[uncert_row], uncertainty[[dat_row]]$rn_2020[uncert_row])), nrow=3, ncol=3)),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = as.numeric(c(uncertainty[[dat_row]]$gamman_2014[uncert_row], uncertainty[[dat_row]]$gamman_2017[uncert_row], uncertainty[[dat_row]]$gamman_2020[uncert_row])) * 365
  )
  
  ## assume the same people are getting nets each round
  correlationsb1 <- get_correlation_parameters(bednetparams_2)
  correlationsb1$inter_round_rho('bednets', 1)
  
  correlationsb2 <- get_correlation_parameters(bednetparams_2b)
  correlationsb2$inter_round_rho('bednets', 1)
  
  
  ## Run the simulations
  output1_matching_rct_and_2020 <- run_simulation(sim_length, bednetparams_2,correlationsb1)
  output2_matching_rct_and_sched <- run_simulation(sim_length, bednetparams_2b,correlationsb2)
  
  output1_matching_rct_and_2020$pv_730_3650 = output1_matching_rct_and_2020$n_detect_730_3650/output1_matching_rct_and_2020$n_730_3650
  output2_matching_rct_and_sched$pv_730_3650 = output2_matching_rct_and_sched$n_detect_730_3650/output2_matching_rct_and_sched$n_730_3650
  
  output1_matching_rct_and_2020$clin_inc_0_36500 = output1_matching_rct_and_2020$n_detect_730_3650
  output2_matching_rct_and_sched$clin_inc_0_36500 = output2_matching_rct_and_sched$n_detect_730_3650
  
  return(data.frame(timestep = output1_matching_rct_and_2020$timestep,
                    
                    prev_match = output1_matching_rct_and_2020$pv_730_3650,
                    clin_match = output1_matching_rct_and_2020$clin_inc_0_36500,
                    
                    prev_scheduled = output2_matching_rct_and_sched$pv_730_3650,
                    clin_scheduled = output2_matching_rct_and_sched$clin_inc_0_36500
                    
  ))
}

uncertainty_data = list()

for(j in 1:nrow(dt2)){
  
  for(i in 1:20){
    test = malsim_actual_f(test_data = dt2,dat_row = j, 
                           uncertainty = uncertainty_1, uncert_row = i)
    
    uncertainty_data[[i]] = test
  }
  saveRDS(uncertainty_data, file=paste0("simulations/uncertainty_outputs_cluster_",print(dt2$cluster[j]),".RData"))
  
}

