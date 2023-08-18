##########################
##
## Figure 3
##
##########################


## median runs
dt2 = read.csv("raw data/dt2.csv",header = TRUE)

dd = expand.grid(cluster = dt2$cluster)

## baseline, 6, 12, 18, 25 months
for(i in 1:nrow(dt2)){
  simulated_estimate = c(test1[[i]]$prev_match[dt2$days_after_jan_2017[i]],
                         test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+6*30],
                         test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+12*30],
                         test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+18*30],
                         test1[[i]]$prev_match[dt2$days_after_jan_2017[i]+25*30])
  observed_estimate = c(dt2$Prevalence_baseline_2_10_yrs[i],
                        dt2$Prevalence_6m[i],dt2$Prevalence_12m[i],dt2$Prevalence_18m[i],dt2$Prevalence_25m[i])
  dd[i,2:11] = c(observed_estimate,simulated_estimate)
  
}

## generate absolute difference in estimate
for(i in 1:5){
  for(j in 1:nrow(dt2)){
    
    dd[j,11+i] = abs(dd[j,1+i] - dd[j,6+i])
    
  }
  
}
df = merge(dt2, dd, by = "cluster")
df$nets_early = ifelse(df$net_delivered_2020 - df$net_schedule_2020 < 0, 1, 0)


#######################
##
## Estimating the population of people receiving nets early

## Number of clusters with early nets
100 - sum(df$nets_early)

## Number of people receivng nets early
sum(df$HSD.population.at.2017.18.UCC[df$nets_early == 1])

## Number of people receiving nets laet
sum(df$HSD.population.at.2017.18.UCC[df$nets_early == 0])

sum(df$HSD.population.at.2017.18.UCC[df$nets_early == 0])/sum(df$HSD.population.at.2017.18.UCC)

## How early/late were people receiving nets?
df$delivery_early_late_category = ifelse(df$net_delivered_2020 - df$net_schedule_2020 > 90, "3 months or more late",
                                         ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= 90 & df$net_delivered_2020 - df$net_schedule_2020 > 60, "2 to 3 months late",
                                                ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= 60 & df$net_delivered_2020 - df$net_schedule_2020 > 30, "1 to 2 months late",
                                                       ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= 30 & df$net_delivered_2020 - df$net_schedule_2020 > 0, "1 month late",
                                                              ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= 0 & df$net_delivered_2020 - df$net_schedule_2020 > -30, "1 month early",
                                                                     ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= -30 & df$net_delivered_2020 - df$net_schedule_2020 > -60, "1 to 2 months early",
                                                                            ifelse(df$net_delivered_2020 - df$net_schedule_2020 <= -60 & df$net_delivered_2020 - df$net_schedule_2020 > -90, "2 to 3 months early","3 months or more early")))))))

### Figure 3

par(mfrow=c(2,2))
par(mar=c(4,6,5,2))

## a) shows
barplot(c(mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "3 months or more early"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "2 - 3 months early"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "1 to 2 months early"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "1 month early"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "1 month late"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "1 to 2 months late"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "2 to 3 months late"]),
          mean(df$HSD.population.at.2017.18.UCC[df$delivery_early_late_category == "3 month or more late"])),
        ylab="",cex.axis = 1.2,cex.lab = 1.2,
        xlab = "Received nets early (months)    Received nets late (months)",
        xaxt="n",yaxt="n",ylim=c(0,200000))
        # col = c("white","grey80","grey50","grey20",
        #         "lightblue","blue","darkblue","black",NA))
# col = c("white","grey80","grey50","grey20",
# adegenet::transp("darkblue",c(0.3,0.5,0.7,0.9)),NA))
mtext("Mean Population estimate for",
      side=2,line = 3.8)
mtext("Health sub District (thousands)",
      side=2,line = 2.6)

abline(v=4.9,lty=2)

axis(1, at=seq(0.8,9.1,length=8),
     labels=c(">3","2-3","1-2","0-1","0-1","1-2","2-3",">3"))
axis(2,las=2,at=c(0,50000,100000,150000,200000),
     labels = c(0,50,100,150,200))
        
        


## b) to show
## read in the median data outputs
test1 = readRDS("simulations/median_outputs_cluster_orderv1.RData")
check1 = readRDS("simulations/scheduled_delivered_nets_2020.RData")

plot((test1[[4]]$n_clin_match/10000)*dt2$HSD.population.at.2017.18.UCC[4] ~ test1[[4]]$timestep, xlim = c(2*365,365*10),
     ylab = "",
     ylim=c(0,800),
     xlab = "Time",xaxt = "n",bty = "n",pch="",yaxt="n",
     main = "Buliisa Sub-District (Population estimate 124876)")
axis(1, at = c(2*365,4*365,6*365,8*365,10*365),
     labels = c("Jan 2016","Jan 2018","Jan 2020","Jan 2022","Jan 2024"))
axis(2, las = 2, at = seq(0,800,100))
mtext("All-age clinical incidence",
      side=2,line = 3.2)

# for(i in 1:20){
#   uncertainty_4 = readRDS(paste0("simulations/uncertainty_outputs_cluster_4.RData"))
#   lines(uncertainty_4[[i]]$n_clin_scheduled ~ test1[[4]]$timestep,col=adegenet::transp("grey",0.4))
#   lines(uncertainty_4[[i]]$n_clin_match ~ test1[[4]]$timestep,col=adegenet::transp("aquamarine",0.4))
#   
# }
lines((test1[[4]]$n_clin_scheduled/10000)*dt2$HSD.population.at.2017.18.UCC[4] ~ test1[[4]]$timestep, col="darkred",cex=1.5)
lines((test1[[4]]$n_clin_match/10000)*dt2$HSD.population.at.2017.18.UCC[4] ~ test1[[4]]$timestep, col="aquamarine4",cex=1.5)

polygon(c(6*365, 9*365, 
          9*365, 6*365),
        c(1550,1550,1450,1450),border=NA,col = adegenet::transp("blue",0.4))

polygon(c(5*365, 6*365, 
          6*365, 5*365),
        c(800,800,750,750),border=NA,col = adegenet::transp("orange",0.7))
polygon(c(6*365, 7*365, 
          7*365, 6*365),
        c(750,750,700,700),border=NA,col = adegenet::transp("darkblue",0.4))
polygon(c(7*365, 8*365, 
          8*365, 7*365),
        c(700,700,650,650),border=NA,col = adegenet::transp("darkblue",0.4))
polygon(c(8*365, 9*365, 
          9*365, 8*365),
        c(650,650,600,600),border=NA,col = adegenet::transp("darkblue",0.4))


## c)
clin_inc_1000 = expand.grid(year = c(2019,2020,2021,2022))
 

for(i in 1:nrow(df)){

  df$sum_inc_clinical_2019_sched[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_scheduled[(5*365):(6*365)])/10000
  df$sum_inc_clinical_2020_sched[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_scheduled[(6*365):(7*365)])/10000
  df$sum_inc_clinical_2021_sched[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_scheduled[(7*365):(8*365)])/10000
  df$sum_inc_clinical_2022_sched[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_scheduled[(8*365):(9*365)])/10000
  
  df$sum_inc_clinical_2019[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_match[(5*365):(6*365)])/10000
  df$sum_inc_clinical_2020[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_match[(6*365):(7*365)])/10000
  df$sum_inc_clinical_2021[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_match[(7*365):(8*365)])/10000
  df$sum_inc_clinical_2022[i] = df$HSD.population.at.2017.18.UCC[i]*sum(test1[[i]]$n_clin_match[(8*365):(9*365)])/10000

  df$mean_inc_clinical_2019[i] = df$HSD.population.at.2017.18.UCC[i]*mean(test1[[i]]$n_clin_match[(5*365):(6*365)])/10000
  df$mean_inc_clinical_2020[i] = df$HSD.population.at.2017.18.UCC[i]*mean(test1[[i]]$n_clin_match[(6*365):(7*365)])/10000
  df$mean_inc_clinical_2021[i] = df$HSD.population.at.2017.18.UCC[i]*mean(test1[[i]]$n_clin_match[(7*365):(8*365)])/10000
  df$mean_inc_clinical_2022[i] = df$HSD.population.at.2017.18.UCC[i]*mean(test1[[i]]$n_clin_match[(8*365):(9*365)])/10000

  df$sum_inc_clinical_2019_1000[i] = 1000*sum(test1[[i]]$n_clin_match[(5*365):(6*365)])/10000
  df$sum_inc_clinical_2020_1000[i] = 1000*sum(test1[[i]]$n_clin_match[(6*365):(7*365)])/10000
  df$sum_inc_clinical_2021_1000[i] = 1000*sum(test1[[i]]$n_clin_match[(7*365):(8*365)])/10000
  df$sum_inc_clinical_2022_1000[i] = 1000*sum(test1[[i]]$n_clin_match[(8*365):(9*365)])/10000
  
  df$mean_inc_clinical_2019_1000[i] = 1000*mean(test1[[i]]$n_clin_match[(5*365):(6*365)])/10000
  df$mean_inc_clinical_2020_1000[i] = 1000*mean(test1[[i]]$n_clin_match[(6*365):(7*365)])/10000
  df$mean_inc_clinical_2021_1000[i] = 1000*mean(test1[[i]]$n_clin_match[(7*365):(8*365)])/10000
  df$mean_inc_clinical_2022_1000[i] = 1000*mean(test1[[i]]$n_clin_match[(8*365):(9*365)])/10000
  
}

clin_inc_1000$early_delivery = c(mean(df$sum_inc_clinical_2019_1000[df$nets_early == 1]),
                                 mean(df$sum_inc_clinical_2020_1000[df$nets_early == 1]),
                                 mean(df$sum_inc_clinical_2021_1000[df$nets_early == 1]),
                                 mean(df$sum_inc_clinical_2022_1000[df$nets_early == 1]) )

clin_inc_1000$late_delivery = c(mean(df$sum_inc_clinical_2019_1000[df$nets_early == 0]),
                                mean(df$sum_inc_clinical_2020_1000[df$nets_early == 0]),
                                mean(df$sum_inc_clinical_2021_1000[df$nets_early == 0]),
                                mean(df$sum_inc_clinical_2022_1000[df$nets_early == 0]) )

clin_inc_1000

barplot(as.numeric(c(clin_inc_1000[1,2:3],NA,
          clin_inc_1000[2,2:3],NA,
          clin_inc_1000[3,2:3],NA,
          clin_inc_1000[4,2:3])),
        ylab = "",xlab = "Years",
        xaxt="n",cex.axis = 1.2,cex.lab = 1.2,
        ylim=c(0,800),yaxt="n",
        col = rep(c("lightgreen",adegenet::transp("blue",0.7),NA),3))
axis(2,las=2,at=c(0,100,200,300,400))
axis(1,at=c(1.35,4.88,8.45,11.78),
     labels = c(2019,2020,2021,2022))

mtext("Mean all-age clinical cases per ",
      side=2,line = 3.8)
mtext("Health sub District per 1000 people",
      side=2,line = 2.6)

legend("topleft",legend=c("Delivered early","Delivered late"),
       title="ITN mass campaign",bty="n",
       pch=15,col=c("lightgreen","blue"))

points(rnorm(mean=0.8,n = length(df$mean_inc_clinical_2019[df$nets_early == 1]),sd=0.1),
       df$mean_inc_clinical_2019[df$nets_early == 1],col=adegenet::transp("darkgreen",0.5),
       pch=19)
points(rnorm(mean=1.9,n = length(df$mean_inc_clinical_2019[df$nets_early == 0]),sd=0.1),
       df$mean_inc_clinical_2019[df$nets_early == 0],col=adegenet::transp("darkblue",0.5),
       pch=19)
points(rnorm(mean=4.25,n = length(df$mean_inc_clinical_2020[df$nets_early == 1]),sd=0.1),
       df$mean_inc_clinical_2020[df$nets_early == 1],col=adegenet::transp("darkgreen",0.5),
       pch=19)
points(rnorm(mean=5.35,n = length(df$mean_inc_clinical_2020[df$nets_early == 0]),sd=0.1),
       df$mean_inc_clinical_2020[df$nets_early == 0],col=adegenet::transp("darkblue",0.5),
       pch=19)
points(rnorm(mean=7.95,n = length(df$mean_inc_clinical_2021[df$nets_early == 1]),sd=0.1),
       df$mean_inc_clinical_2021[df$nets_early == 1],col=adegenet::transp("darkgreen",0.5),
       pch=19)
points(rnorm(mean=9.05,n = length(df$mean_inc_clinical_2021[df$nets_early == 0]),sd=0.1),
       df$mean_inc_clinical_2021[df$nets_early == 0],col=adegenet::transp("darkblue",0.5),
       pch=19)
points(rnorm(mean=11.65,n = length(df$mean_inc_clinical_2022[df$nets_early == 1]),sd=0.1),
       df$mean_inc_clinical_2022[df$nets_early == 1],col=adegenet::transp("darkgreen",0.5),
       pch=19)
points(rnorm(mean=12.75,n = length(df$mean_inc_clinical_2022[df$nets_early == 0]),sd=0.1),
       df$mean_inc_clinical_2022[df$nets_early == 0],col=adegenet::transp("darkblue",0.5),
       pch=19)



## d)
sum_clin_cases = expand.grid(year = c(2019,2020,2021,2022))
sum_clin_cases$scheduled = c(sum(df$sum_inc_clinical_2019_sched),
                             sum(df$sum_inc_clinical_2020_sched),
                             sum(df$sum_inc_clinical_2021_sched),
                             sum(df$sum_inc_clinical_2022_sched))
                             
sum_clin_cases$delivered = c(sum(df$sum_inc_clinical_2019),
                             sum(df$sum_inc_clinical_2020),
                             sum(df$sum_inc_clinical_2021),
                             sum(df$sum_inc_clinical_2022))


barplot(as.numeric(c(sum_clin_cases[1,2:3],NA,
                     sum_clin_cases[2,2:3],NA,
                     sum_clin_cases[3,2:3],NA,
                     sum_clin_cases[4,2:3])),
        ylim = c(0,6000000),ylab = "",xlab = "Years",
        xaxt="n",cex.axis = 1.2,cex.lab = 1.2,
        yaxt="n",
        col = c(adegenet::transp(c("darkred","aquamarine2"),0.2),NA,rep(c("darkred",adegenet::transp("aquamarine4",0.7),NA),3))
        )
axis(2,las=2,at=seq(0,4000000,500000),labels=seq(0,4000,500))
axis(1,at=c(1.35,4.88,8.45,11.78),
     labels = c(2019,2020,2021,2022))

mtext("Total cases (thousands)",
      side=2,line = 3.2)

legend("topright",legend=c("Scheduled","Delivered"),
       title="ITN mass campaign",bty="n",
       pch=15,col=c("darkred","aquamarine4"))


par(xpd=NA,cex = 1.11)

text(x = -20, y = 18000000,"(A)")
text(x = -2, y = 18000000,"(B)")
text(x = -20, y = 8000000,"(C)")
text(x = -2, y = 8000000,"(D)")
