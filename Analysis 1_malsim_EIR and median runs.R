########################################################################
##
## Task: Two consecutive randomised control trials were conducted
## in Uganda looking at mosquito nets. Subsequently, the Covid-19 pandemic
## hit and the cases of malaria have since increased in 2021 and 2022 
## in many places across the country. Using the trial data and modelling
## we aim to explore whether we can answer if the timing of the distribution
## of campaign nets by National Programs, that was disrupted by the pandemic
## could be one reason for the impacts.

## 104 Clusters were completed in the first trial with less conducted later
## but in over lapping areas. 

## In the second trial different areas were considered 
## but these regions are grouped to learn from the original trial, calibrate
## the model to a reasonable degree and then project forward.

## In the original trial, 3 clusters were dropped due to checks showing that
## a non-trial planned net was dominant in these clusters. And 2 clusters 
## were re-designated based on high use of a different trial net. 

## subsequent deployments were principly pyrethroid nets or pyrethroid-PBO
## nets with one exception (cluster 104 received a pyrethroid-pyriproxyfen 
## net). We cannot simulate this yet so this was also dropped.

## We can create timelines from these data for the remaining 100 clusters

## Input data
dt = read.csv("raw data/collated_trial_data.csv",header=TRUE)

## drop the clusters we are not using
dt1 = subset(dt,dt$include == 1)
dim(dt1)

## add the level of resistance and net information
pyr_net_efficacy = read.csv("raw data/pyrethroid_only_nets.csv",header=TRUE)
pbo_net_efficacy = read.csv("raw data/pyrethroid_pbo_nets.csv",header=TRUE)

## align but with uncertainty....
## initial median runs first though and for calibrations
stor = expand.grid(cluster = dt1$cluster)
for(i in 1:nrow(dt1)){
  stor[i,2] = pyr_net_efficacy$dn0_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_1[i])]
  stor[i,3] = pyr_net_efficacy$rn0_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_1[i])]
  stor[i,4] = pyr_net_efficacy$gamman_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_1[i])]
  }
colnames(stor) = c("cluster","pyr_net_dn0_2014","pyr_net_rn_2014","pyr_net_gamman_2014")
  
stor2017 = expand.grid(cluster = dt1$cluster)
for(i in 1:nrow(dt1)){
  stor2017[i,2] = pyr_net_efficacy$dn0_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_2[i])]
  stor2017[i,3] = pyr_net_efficacy$rn0_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_2[i])]
  stor2017[i,4] = pyr_net_efficacy$gamman_med[which(round(pyr_net_efficacy$resistance,2) == dt1$res_2[i])]

  stor2017[i,5] = pbo_net_efficacy$dn0_med[which(round(pbo_net_efficacy$resistance,2) == dt1$res_2[i])]
  stor2017[i,6] = pbo_net_efficacy$rn0_med[which(round(pbo_net_efficacy$resistance,2) == dt1$res_2[i])]
  stor2017[i,7] = pbo_net_efficacy$gamman_med[which(round(pbo_net_efficacy$resistance,2) == dt1$res_2[i])]
}
colnames(stor2017) = c("cluster","pyr_net_dn0_2017","pyr_net_rn_2017","pyr_net_gamman_2017",
                   "pbo_net_dn0_2017","pbo_net_rn_2017","pbo_net_gamman_2017")

dt1$dn0_med_2014 = stor$pyr_net_dn0_2014
dt1$rn0_med_2014 = stor$pyr_net_rn_2014
dt1$gamman_med_2014 = stor$pyr_net_gamman_2014

dt1$net_for_2017_trial = ifelse(dt1$Net_Type == "PermaNet 2.0", "Standard",
                                ifelse(dt1$Net_Type == "Olyset Net", "Standard",
                                       "PBO"))

dt1$dn0_med_2017 = ifelse(dt1$net_for_2017_trial == "Standard", stor2017$pyr_net_dn0_2017,stor2017$pbo_net_dn0_2017)
dt1$rn0_med_2017 = ifelse(dt1$net_for_2017_trial == "Standard", stor2017$pyr_net_rn_2017,stor2017$pbo_net_rn_2017)
dt1$gamman_med_2017 = ifelse(dt1$net_for_2017_trial == "Standard", stor2017$pyr_net_gamman_2017,stor2017$pbo_net_gamman_2017)


dt1$dn0_med_2020 = ifelse(dt1$net_type == "Standard", stor2017$pyr_net_dn0_2017,stor2017$pbo_net_dn0_2017)
dt1$rn0_med_2020 = ifelse(dt1$net_type == "Standard", stor2017$pyr_net_rn_2017,stor2017$pbo_net_rn_2017)
dt1$gamman_med_2020 = ifelse(dt1$net_type == "Standard", stor2017$pyr_net_gamman_2017,stor2017$pbo_net_gamman_2017)

# write.csv(dt1,"raw data/collated_trial_data.csv")
######################
##
## Data visualisation for paper and supplement figures
##
######################

## quick check that PBO kills better and increased resistance decreases mortality
plot(dt1$dn0_med_2014 ~ dt1$res_1,ylim=c(0,1),xlim=c(0,1),
     ylab = "Net mortality efficacy",
     xlab = "Resistance")
points(dt1$dn0_med_2017[dt1$net_for_2017_trial == "Standard"] ~ dt1$res_2[dt1$net_for_2017_trial == "Standard"],pch=19,col="red")
points(dt1$dn0_med_2017[dt1$net_for_2017_trial == "PBO"] ~ dt1$res_2[dt1$net_for_2017_trial == "PBO"],pch=19,col="blue")
points(dt1$dn0_med_2020[dt1$net_type == "PBO"] ~ dt1$res_3[dt1$net_type == "PBO"],pch=19,col="blue")
points(dt1$dn0_med_2020[dt1$net_type == "Standard"] ~ dt1$res_3[dt1$net_type == "Standard"],pch=19,col="red")

######################
##
## Model code using malariasimulation
## first to calibrate clusters to prevalence at 
## trial baseline
##
######################

## Calibration to baseline prev
library(malariasimulation)
malsim_calib_f = function(test_data,dat_row){
  
  year <- 365
  month <- 30
  sim_length <- 7 * year
  ## This is spanning Jan 2014 - Dec 2025
  ## Assume that all places have received nets at start of 2014
  ## then relative to jan 2017 (see input file (collated_trial_data)...)
  ## then as noted for 2020
  
  
  ## and we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 10000
  starting_EIR = test_data$eir_est[dat_row]
    # starting_EIR <- ifelse(test_data$Prevalence_baseline_2_10_yrs[dat_row] < 0.2,
  #                        0.005 + test_data$total_M.x[dat_row]*0.7,
  #                        ifelse(test_data$Prevalence_baseline_2_10_yrs[dat_row] > 0.2 & test_data$Prevalence_baseline_2_10_yrs[dat_row] < 0.35,
  #                               0.005 + test_data$total_M.x[dat_row]*0.8,0.005 + test_data$total_M.x[dat_row]))
  ## upper estimate
  
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      
      prevalence_rendering_min_ages = 2 * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = 10 * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      g0 = test_data$seasonal_a0[dat_row],
      g = c(test_data$seasonal_a1[dat_row], test_data$seasonal_a2[dat_row], test_data$seasonal_a3[dat_row]),
      h = c(test_data$seasonal_b1[dat_row], test_data$seasonal_b2[dat_row], test_data$seasonal_b3[dat_row]),
      
      individual_mosquitoes = FALSE ## Update next
    )
  )
  
  # simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set species
  # if these were unknown, we are assuming the mean situation across known clusters
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(test_data$prop_gam[dat_row],
                                         test_data$prop_ara[dat_row],
                                         1 - test_data$prop_gam[dat_row] - test_data$prop_ara[dat_row]))
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
  # and this is what was done
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(test_data$ITN_Use_prior[dat_row], ## historic prior to RCT
                  test_data$ITN_Use_A[dat_row], ## historic during RCT (on deployment) 
                  test_data$ITN_Use_A[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = test_data$itn_leave_durMean[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = t(matrix(as.numeric(c(test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row],
                                test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row],
                                test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rn = t(matrix(as.numeric(c(test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row],
                               test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row],
                               test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = as.numeric(c(test_data$gamman_med_2014[dat_row], test_data$gamman_med_2017[dat_row], test_data$gamman_med_2020[dat_row])) * 365
  )
  
  init_EIR <- c(0.01, 0.1, 1, 5, 10, 15, 20, 25, 50, 80, 100, 150, 200)
  
  # For each initial EIR, calculate equilibrium parameter set and run the simulation
  
  malSim_outs <- lapply(
    init_EIR,
    function(init) {
      p_i <- set_equilibrium(bednetparams_1, init)
      run_simulation(5 * year, p_i)
    }
  )
  
  # Convert the default EIR output (per vector species, per timestep, across 
  # the entire human population) to a cross-vector species average EIR per
  # person per year across the final year of the simulation:
  malSim_EIR <- lapply(
    malSim_outs,
    function(output) {
      mean(
        rowSums(
          output[
            output$timestep %in% seq(bednet_events$timestep[2]-365, bednet_events$timestep[2]),
            grepl('EIR_', names(output))
          ] / human_population * year
        )
      )
    }
  )
  
  # Calculate the average PfPR2-10 value across the final year for each initial
  # EIR value:
  malSim_prev <- lapply(
    malSim_outs,
    function(output) {
      mean(
        output[
          output$timestep %in% seq(bednet_events$timestep[2]-365, bednet_events$timestep[2]),
          'n_detect_730_3650'
        ] / output[
          output$timestep %in% seq(bednet_events$timestep[2]-365, bednet_events$timestep[2]),
          'n_730_3650'
        ]
      )
    }
  )
  
  # Create dataframe of initial EIR, output EIR, and PfPR2-10 results:
  malSim_P2E <- cbind.data.frame(init_EIR, EIR = unlist(malSim_EIR), prev = unlist(malSim_prev))
  return(malSim_P2E)

}

row1test = malsim_calib_f(dt1,1)
EIR_estimates = list()

for(i in 1:5){
  
  EIR_estimates[[i]] = malsim_calib_f(dt1,i)
    
}

plot(EIR_estimates[[1]]$EIR ~ EIR_estimates[[1]]$prev, pch=19,ylim = c(0,200),xlim = c(0,1),
     ylab = "Estimated EIR",xlab = "Prevalence 2 to 10 years")
lines(EIR_estimates[[1]]$EIR ~ EIR_estimates[[1]]$prev,lty=1,col="darkred",cex=1.5)
abline(v=dt1$Prevalence_baseline_2_10_yrs[1],col="darkred")

for(i in 1:5){
  lines(EIR_estimates[[i]]$EIR ~ EIR_estimates[[i]]$prev,lty=1,col="blue",cex=1.5)
  abline(v=dt1$Prevalence_baseline_2_10_yrs[i],col="blue")
  
}
DT1_PREV_ORDER = dt1[order(dt1$Prevalence_baseline_2_10_yrs), ]

dt_eir = data.frame(y_eir = EIR_estimates[[1]]$EIR,
                    prev = EIR_estimates[[1]]$prev)

y_EIR <- loess(y_eir ~ prev, dt_eir)
eir_preds = predict(y_EIR, data.frame(prev = DT1_PREV_ORDER$Prevalence_baseline_2_10_yrs), se = TRUE)
DT1_PREV_ORDER$eir_preds = c(rep(0.1,6),eir_preds$fit[7:97],160,185,190)

dt2 = DT1_PREV_ORDER[order(DT1_PREV_ORDER$cluster),]
head(dt2)

# write.csv(dt1,"raw data/collated_trial_data_updateEIR.csv")

plot(dt2$eir_preds ~ dt2$Prevalence_baseline_2_10_yrs,pch=19,cex=1.4)
lines(EIR_estimates[[1]]$EIR ~ EIR_estimates[[1]]$prev,lty=1,col="blue",cex=1)


## visual checks need to increase EIR for some clusters:
# dt2$eir_preds[c(1,5,6,7,10,11,12,13,15,24,34,38,46,48,54,66,69,71,75,78,84)] = 
#   1.4 * dt2$eir_preds[c(1,5,6,7,10,11,12,13,15,24,34,38,46,48,54,66,69,71,75,78,84)]

################################
##
## For Uncertainty introduced through net parameters
## net use, adherence, efficacy
## See Analysis 2
##
################################

## Create input data files
## using uncertainty from Sherrard-Smith et al, Nat Communs, 2022
## using 10% range around median estimates for net use and adherrence

##########################################
##
## Function to run all sims

malsim_actual_f = function(test_data,dat_row){
  
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
    
    coverages = c(test_data$ITN_Use_prior[dat_row], ## historic prior to RCT
                  test_data$ITN_Use_A[dat_row], ## historic during RCT (on deployment) 
                  test_data$ITN_Use_A[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = test_data$itn_leave_durMean[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = t(matrix(as.numeric(c(test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row],
                                test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row],
                                test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rn = t(matrix(as.numeric(c(test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row],
                               test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row],
                               test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = as.numeric(c(test_data$gamman_med_2014[dat_row], test_data$gamman_med_2017[dat_row], test_data$gamman_med_2020[dat_row])) * 365
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
    
    coverages = c(test_data$ITN_Use_prior[dat_row], ## historic prior to RCT
                  test_data$ITN_Use_A[dat_row], ## historic during RCT (on deployment) 
                  test_data$ITN_Use_A[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = test_data$itn_leave_dur[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = t(matrix(as.numeric(c(test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row], test_data$dn0_med_2014[dat_row],
                                test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row], test_data$dn0_med_2017[dat_row],
                                test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row], test_data$dn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rn = t(matrix(as.numeric(c(test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row], test_data$rn0_med_2014[dat_row],
                               test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row], test_data$rn0_med_2017[dat_row],
                               test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row], test_data$rn0_med_2020[dat_row])), nrow=3, ncol=3)),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = as.numeric(c(test_data$gamman_med_2014[dat_row], test_data$gamman_med_2017[dat_row], test_data$gamman_med_2020[dat_row])) * 365
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

  output1_matching_rct_and_2020$n_clin_inc_0_36500 = output1_matching_rct_and_2020$n_inc_clinical_0_36500
  output2_matching_rct_and_sched$n_clin_inc_0_36500 = output2_matching_rct_and_sched$n_inc_clinical_0_36500
  
  return(data.frame(timestep = output1_matching_rct_and_2020$timestep,
                    
                    prev_match = output1_matching_rct_and_2020$pv_730_3650,
                    clin_match = output1_matching_rct_and_2020$clin_inc_0_36500,
                    n_clin_match = output1_matching_rct_and_2020$n_clin_inc_0_36500,
                    
                    prev_scheduled = output2_matching_rct_and_sched$pv_730_3650,
                    clin_scheduled = output2_matching_rct_and_sched$clin_inc_0_36500,
                    n_clin_scheduled = output2_matching_rct_and_sched$n_clin_inc_0_36500
                    
  ))
}


#####################
##
## Loop through the clusters to produce median simulations for each one
test1 = list() ## create an object store model outputs in
check1 = matrix(nrow=nrow(dt2),ncol=3,data = NA) ## create a data.frame to store delivery of nets in for checks

for(i in 92:nrow(dt)){
  test1[[i]] = malsim_actual_f(dt2,dat_row = i)
  
  bednet_events_scheduled = data.frame(
    timestep = c(0, 3, 6, 6) * 365 + c(0,
                                     dt2$days_after_jan_2017[i],
                                     dt2$net_schedule_2020[i],
                                     dt2$net_delivered_2020[i]),
    name=c("background", 
           "trial_nets",
           "2020_scheduled",
           "2020_delivered")
  )
  
  check1[i,1:3] = c(bednet_events_scheduled$timestep[2:4])
}

par(mfrow=c(3,4))
par(mar=c(3,3,2,1))
for(i in 1:12){
  site_cluster = i
  plot(test1[[site_cluster]]$prev_match ~ test1[[site_cluster]]$timestep,pch="",ylim=c(0,1),
       main = print(dt2$cluster[i]))
  lines(test1[[site_cluster]]$prev_match ~ test1[[site_cluster]]$timestep,col="red")
  points(dt2$Prevalence_baseline_2_10_yrs[site_cluster] ~ c(check1[i,1]),pch=19)
  points(dt2$Prevalence_6m[site_cluster] ~ c(check1[i,1] + 6*30),pch=8)
  points(dt2$Prevalence_12m[site_cluster] ~ c(check1[i,1] + 12*30),pch=8)
  points(dt2$Prevalence_18m[site_cluster] ~ c(check1[i,1] + 18*30),pch=8)
  points(dt2$Prevalence_25m[site_cluster] ~ c(check1[i,1] + 25*30),pch=8)
  abline(v = check1[i,1],lty=2)
  abline(v = check1[i,2],lty=2)
  abline(v = check1[i,3],lty=2,col="blue")
  
}


###########################################
##
## These are the median data for Jaffer et al 2023 (submission year)

check1 = as.data.frame(check1)
colnames(check1) = c("trial net days","scheduled 2020 days","delivered 2020 days")

## store these results in simulations folder
saveRDS(test1, file="simulations/median_outputs_cluster_orderv1.RData")
saveRDS(check1, file="simulations/scheduled_delivered_nets_2020v1.RData")
