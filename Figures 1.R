# Figures for Jaffer et al 2023 (submission year)

## Figure 1 Demonstrating the modelling process

# par(oma=c(2, 2, 2, 2), mar=c(3, 3, 3, 3))

## read in the median data outputs
test1 = readRDS("simulations/median_outputs_cluster_order.RData")
check1 = readRDS("simulations/scheduled_delivered_nets_2020.RData")

## set up layout
layout(matrix(c(1, 2,  1, 3), ncol=2))

## Figure 1a
plot(test1[[4]]$prev_scheduled ~ test1[[4]]$timestep, xlim = c(2*365,7.5*365),
     ylim = c(0,1),ylab = "Prevalence in 2 to 10 years (%)",
     yaxt = "n", xaxt = "n",
     xlab = "", bty = "n",pch="")
axis(1, at = c(2*365,3*365,4*365,5*365,6*365,7*365),
     labels = c("Jan 2016","Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021"))
axis(2, las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))

text(2.5*365,0.05,"TIME LINE")
polygon(c(2*365+134, 3*365+134, 
          3*365+134, 2*365+134),
        c(0.4,0.4,0.5,0.5),border=NA,col = adegenet::transp("blue",0.4))

for(i in 1:20){
        uncertainty_4 = readRDS(paste0("simulations/uncertainty_outputs_cluster_4.RData"))
        lines(uncertainty_4[[i]]$prev_scheduled ~ test1[[4]]$timestep,col=adegenet::transp("grey",0.4))
        lines(uncertainty_4[[i]]$prev_match ~ test1[[4]]$timestep,col=adegenet::transp("aquamarine",0.4))
        
}
lines(test1[[4]]$prev_scheduled ~ test1[[4]]$timestep, col="darkred",cex=1.5)
lines(test1[[4]]$prev_match ~ test1[[4]]$timestep, 
      col = "aquamarine4", lty=2, lwd = 1.5)


text(3*365+150,0.8,"Mass campaign 2017", col = "grey30")
arrows(x0 = check1[4,1],
       y0 = 0.4,
       x1 = check1[4,1],
       y1 = 0,angle = 30,code = 2,length = 0.1,lwd = 3,col="grey40")


text(check1[4,2],0.8,"Mass campaign 2020", col = "darkred")
text(check1[4,2],0.72,"as scheduled", col = "darkred")
arrows(x0 = check1[4,2],
       y0 = 0.4,
       x1 = check1[4,2],
       y1 = 0,angle = 30,code = 2,length = 0.1,lwd = 3,col="darkred")


text(check1[4,3],0.6,"Mass campaign 2020", col = "darkgreen")
text(check1[4,3],0.52,"as delivered", col = "darkgreen")
arrows(x0 = check1[4,3],
       y0 = 0.35,
       x1 = check1[4,3],
       y1 = 0,angle = 30,code = 2,length = 0.1,lwd = 3,col="darkgreen",lty=4)

date_estimate = c(3*365+dt2$days_after_jan_2017[dt2$cluster == 4],
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+6*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+12*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+18*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+25*30)
observed_estimate = c(dt2$Prevalence_baseline_2_10_yrs[dt2$cluster == 4],
                      dt2$Prevalence_6m[dt2$cluster == 4],
                      dt2$Prevalence_12m[dt2$cluster == 4],
                      dt2$Prevalence_18m[dt2$cluster == 4],
                      dt2$Prevalence_25m[dt2$cluster == 4])

points(observed_estimate ~ date_estimate,col = c("red","grey40","blue","darkblue","black"),
       pch = c(19,8,8,8,8),cex=1.5)


## Figure 1b
## A measured prev versus modelled prevalence
## point colours and shapes to distinguish
## timing (baseline, 6 months, 12 months, 18 months 25 months) (shape)
## and net type (colour)
xx = seq(0,1,0.2)
yy = seq(0,1,0.2)

plot(yy ~ xx, xaxt = "n", yaxt = "n",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled",
     bty = "n", pch = "")
abline(a= 0,b = 1, lty= 2)
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))



## median runs
dt2 = read.csv("raw data/dt2.csv",header = TRUE)
## baseline, 6, 12, 18, 25 months
for(i in 1:nrow(dt2)){
        simulated_estimate = c(test1[[i]]$prev_match[dt2$days_after_jan_2017[i]],
                               test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                               test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                               test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                               test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+25*30])
        observed_estimate = c(dt2$Prevalence_baseline_2_10_yrs[i],
                              dt2$Prevalence_6m[i],dt2$Prevalence_12m[i],dt2$Prevalence_18m[i],dt2$Prevalence_25m[i])
        
        ## generate estimates for each cross sectional survey
        uncertainty_a = readRDS(paste0("simulations/uncertainty_outputs_cluster_",print(dt2$cluster[i]),".RData"))
        baseline_uncertainty = c(uncertainty_a[[1]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[2]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[3]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[4]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[5]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[6]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[7]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[8]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[9]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[10]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[11]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[12]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[13]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[14]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[15]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[16]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[17]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[18]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[19]]$prev_match[dt2$days_after_jan_2017[i]],
                                 uncertainty_a[[20]]$prev_match[dt2$days_after_jan_2017[i]])
        
        at6m_uncertainty = c(uncertainty_a[[1]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[2]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[3]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[4]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[5]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[6]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[7]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[8]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[9]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[10]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[11]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[12]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[13]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[14]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[15]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[16]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[17]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[18]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[19]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                             uncertainty_a[[20]]$prev_match[dt2$days_after_jan_2017[i]+6*30])
        
        a12m_uncertainty = c(uncertainty_a[[1]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[2]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[3]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[4]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[5]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[6]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[7]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[8]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[9]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[10]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[11]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[12]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[13]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[14]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[15]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[16]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[17]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[18]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[19]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                             uncertainty_a[[20]]$prev_match[dt2$days_after_jan_2017[i]+12*30])
        
        a18m_uncertainty = c(uncertainty_a[[1]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[2]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[3]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[4]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[5]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[6]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[7]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[8]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[9]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[10]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[11]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[12]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[13]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[14]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[15]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[16]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[17]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[18]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[19]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                             uncertainty_a[[20]]$prev_match[dt2$days_after_jan_2017[i]+18*30])
        
        a25m_uncertainty = c(uncertainty_a[[1]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[2]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[3]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[4]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[5]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[6]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[7]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[8]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[9]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[10]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[11]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[12]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[13]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[14]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[15]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[16]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[17]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[18]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[19]]$prev_match[dt2$days_after_jan_2017[i]+25*30],
                             uncertainty_a[[20]]$prev_match[dt2$days_after_jan_2017[i]+25*30])
        
        ## segments from uncertainty
        segments(x0 = min(baseline_uncertainty),
                 x1 = max(baseline_uncertainty),
                 y0 = dt2$Prevalence_baseline_2_10_yrs[i],
                 y1 = dt2$Prevalence_baseline_2_10_yrs[i],
                 col = "red")
        
        # segments(x0 = min(at6m_uncertainty),
        #          x1 = max(at6m_uncertainty),
        #          y0 = dt2$Prevalence_6m[i],
        #          y1 = dt2$Prevalence_6m[i],
        #          col = "grey")
        # 
        # segments(x0 = min(a12m_uncertainty),
        #          x1 = max(a12m_uncertainty),
        #          y0 = dt2$Prevalence_12m[i],
        #          y1 = dt2$Prevalence_12m[i],
        #          col = "blue")
        # 
        # segments(x0 = min(a18m_uncertainty),
        #          x1 = max(a18m_uncertainty),
        #          y0 = dt2$Prevalence_18m[i],
        #          y1 = dt2$Prevalence_18m[i],
        #          col = "darkblue")
        # 
        segments(x0 = min(a25m_uncertainty),
                 x1 = max(a25m_uncertainty),
                 y0 = dt2$Prevalence_25m[i],
                 y1 = dt2$Prevalence_25m[i],
                 col = "black")
        
        # points(observed_estimate ~ simulated_estimate,col = c("red","grey","blue","darkblue","black"),
        #        pch = 19)

        points(c(observed_estimate[1],observed_estimate[5]) ~ 
                       c(simulated_estimate[1],simulated_estimate[5]),col = c("red","black"),
               pch = 19)
        
}

legend("topleft",legend = c("baseline (calibration)",
                            # "6 months",
                            # "12 months",
                            # "18 months",
                            "25 months"),
       col = c("red",
               # "grey","blue","darkblue",
               "black"),
       pch = 19,bty = "n")

## Figure 1c
## Validation example of a cluster with the full uncertainty and
## points overlaid
plot(test1[[4]]$prev_match ~ test1[[4]]$timestep, xlim = c(2*365,365*10),
     ylim = c(0, 1), ylab = "Prevalence children 2 to 10 years (%)",
     xlab = "Time",xaxt = "n",bty = "n",pch="",yaxt="n",
     main = "Buliisa Sub-District")
axis(1, at = c(2*365,4*365,6*365,8*365,10*365),
     labels = c("Jan 2016","Jan 2018","Jan 2020","Jan 2022","Jan 2024"))
axis(2, las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))

for(i in 1:20){
        uncertainty_4 = readRDS(paste0("simulations/uncertainty_outputs_cluster_4.RData"))
        lines(uncertainty_4[[i]]$prev_scheduled ~ test1[[4]]$timestep,col=adegenet::transp("grey",0.4))
        lines(uncertainty_4[[i]]$prev_match ~ test1[[4]]$timestep,col=adegenet::transp("aquamarine",0.4))
        
}
lines(test1[[4]]$prev_scheduled ~ test1[[4]]$timestep, col="darkred",cex=1.5)
lines(test1[[4]]$prev_match ~ test1[[4]]$timestep, col="aquamarine4",cex=1.5)

date_estimate = c(3*365+dt2$days_after_jan_2017[dt2$cluster == 4],
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+6*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+12*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+18*30,
                  3*365+dt2$days_after_jan_2017[dt2$cluster == 4]+25*30)
observed_estimate = c(dt2$Prevalence_baseline_2_10_yrs[dt2$cluster == 4],
                      dt2$Prevalence_6m[dt2$cluster == 4],
                      dt2$Prevalence_12m[dt2$cluster == 4],
                      dt2$Prevalence_18m[dt2$cluster == 4],
                      dt2$Prevalence_25m[dt2$cluster == 4])

points(observed_estimate ~ date_estimate,col = c("red","grey40","blue","darkblue","black"),
       pch = c(19,8,8,8,8),cex=1.5)

polygon(c(6*365, 9*365, 
          9*365, 6*365),
        c(0.6,0.6,0.7,0.7),border=NA,col = adegenet::transp("blue",0.4))

