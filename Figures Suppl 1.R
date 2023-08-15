
######################################
##
## Supplementary Figures - clusters and data
##
#######################################
par(mfrow = c(2,4))
for(i in 1:nrow(dt2)){
  ## Figure 1c
  ## Validation example of a cluster with the full uncertainty and
  ## points overlaid
  plot(test1[[i]]$prev_match ~ test1[[i]]$timestep, xlim = c(2*365,365*10),
       ylim = c(0, 1), ylab = "Prevalence children 2 to 10 years (%)",
       xlab = "Time",xaxt = "n",bty = "n",pch="",yaxt="n",
       main = print(dt2$Health_Sub_District.1[i]))
  axis(1, at = c(2*365,4*365,6*365,8*365,10*365),
       labels = c("Jan 2016","Jan 2018","Jan 2020","Jan 2022","Jan 2024"))
  axis(2, las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))
  
  for(j in 1:20){
    uncertainty_4 = readRDS(paste0("simulations/uncertainty_outputs_cluster_",print(dt2$cluster[i]),".RData"))
    lines(uncertainty_4[[j]]$prev_scheduled ~ test1[[i]]$timestep,col=adegenet::transp("grey",0.4))
    lines(uncertainty_4[[j]]$prev_match ~ test1[[i]]$timestep,col=adegenet::transp("aquamarine",0.4))
    
  }
  lines(test1[[i]]$prev_scheduled ~ test1[[i]]$timestep, col="darkred",cex=1.5)
  lines(test1[[i]]$prev_match ~ test1[[i]]$timestep, col="aquamarine4",cex=1.5)
  
  date_estimate = c(3*365+dt2$days_after_jan_2017[dt2$cluster == dt2$cluster[i]],
                    3*365+dt2$days_after_jan_2017[dt2$cluster == dt2$cluster[i]]+6*30,
                    3*365+dt2$days_after_jan_2017[dt2$cluster == dt2$cluster[i]]+12*30,
                    3*365+dt2$days_after_jan_2017[dt2$cluster == dt2$cluster[i]]+18*30,
                    3*365+dt2$days_after_jan_2017[dt2$cluster == dt2$cluster[i]]+25*30)
  observed_estimate = c(dt2$Prevalence_baseline_2_10_yrs[dt2$cluster == dt2$cluster[i]],
                        dt2$Prevalence_6m[dt2$cluster == dt2$cluster[i]],
                        dt2$Prevalence_12m[dt2$cluster == dt2$cluster[i]],
                        dt2$Prevalence_18m[dt2$cluster == dt2$cluster[i]],
                        dt2$Prevalence_25m[dt2$cluster == dt2$cluster[i]])
  
  points(observed_estimate ~ date_estimate,col = c("red","grey40","blue","darkblue","black"),
         pch = c(19,8,8,8,8),cex=1.5)
  
  polygon(c(6*365, 9*365, 
            9*365, 6*365),
          c(0.6,0.6,0.7,0.7),border=NA,col = adegenet::transp("blue",0.4))
  
}


