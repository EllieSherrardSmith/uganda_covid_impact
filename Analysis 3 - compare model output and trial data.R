#########################
##
## Analysis 3
##
##########################

## statistical comparison of estimated prevalence at 
## 6, 12, 18 and 25 months

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

## LMs for respective simulated prevalence versus observed
# summary(lm(dd[,2] ~ dd[,7] + 0)) ## baseline... calibrated
summary(lm(dd[,3] ~ dd[,8] + 0))
summary(lm(dd[,4] ~ dd[,9] + 0))
summary(lm(dd[,5] ~ dd[,10] + 0))
summary(lm(dd[,6] ~ dd[,11] + 0))

length(which(dd[,13] < 0.1))
length(which(dd[,14] < 0.1))
length(which(dd[,15] < 0.1))
length(which(dd[,16] < 0.1))

df = merge(dt2, dd, by = "cluster")
df$nets_early = ifelse(df$net_delivered_2020 - df$net_schedule_2020 < 0, 1, 0)

## Supplementary Figure 2
par(mfrow=c(2,2))

mod6months = glm(df$V3 ~ df$V8 + df$nets_early)
prev_mod = seq(0,max(df$V3,na.rm=TRUE),0.01)
y_nets_early_no = mod6months$coefficients[1] + mod6months$coefficients[2] * prev_mod + 0
y_nets_early_ye = mod6months$coefficients[1] + mod6months$coefficients[2] * prev_mod + mod6months$coefficients[3]

plot(df$V3 ~ df$V8,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 6-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V3[df$nets_early == 0] ~ df$V8[df$nets_early == 0], col = "grey40", pch = 19)
points(df$V3[df$nets_early == 1] ~ df$V8[df$nets_early == 1], col = "grey40")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "grey40")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "grey40")
axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))
##
mod12months = glm(df$V4 ~ df$V9 + df$nets_early)
prev_mod = seq(0,max(df$V4,na.rm=TRUE),0.01)
y_nets_early_no = mod12months$coefficients[1] + mod12months$coefficients[2] * prev_mod + 0
y_nets_early_ye = mod12months$coefficients[1] + mod12months$coefficients[2] * prev_mod + mod12months$coefficients[3]

plot(df$V4 ~ df$V9,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 12-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V4[df$nets_early == 0] ~ df$V9[df$nets_early == 0], col = "blue", pch = 19)
points(df$V4[df$nets_early == 1] ~ df$V9[df$nets_early == 1], col = "blue")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "blue")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "blue")
axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))


##

mod18months = glm(df$V5 ~ df$V10 + df$nets_early)
prev_mod = seq(0,max(df$V5,na.rm=TRUE),0.01)
y_nets_early_no = mod18months$coefficients[1] + mod18months$coefficients[2] * prev_mod + 0
y_nets_early_ye = mod18months$coefficients[1] + mod18months$coefficients[2] * prev_mod + mod18months$coefficients[3]

plot(df$V5 ~ df$V10,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 18-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V5[df$nets_early == 0] ~ df$V10[df$nets_early == 0], col = "darkblue", pch = 19)
points(df$V5[df$nets_early == 1] ~ df$V10[df$nets_early == 1], col = "darkblue")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "darkblue")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "darkblue")

axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))


##
mod25months = glm(df$V6 ~ df$V11 + df$nets_early)
prev_mod = seq(0,max(df$V6,na.rm=TRUE),0.01)
y_nets_early_no = mod25months$coefficients[1] + mod25months$coefficients[2] * prev_mod + 0
y_nets_early_ye = mod25months$coefficients[1] + mod25months$coefficients[2] * prev_mod + mod25months$coefficients[3]

plot(df$V6 ~ df$V11,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 25-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V6[df$nets_early == 0] ~ df$V11[df$nets_early == 0],  pch = 19)
points(df$V6[df$nets_early == 1] ~ df$V11[df$nets_early == 1], )
lines(y_nets_early_no ~ prev_mod, lty = 1)
lines(y_nets_early_ye ~ prev_mod, lty = 2)

axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))


## And try with an interaction term
## Supplementary Figure 2
par(mfrow=c(2,2))

mod6months = glm(df$V3 ~ df$V8 * df$nets_early)
prev_mod = seq(0,max(df$V3,na.rm=TRUE),0.01)
y_nets_early_no = mod6months$coefficients[1] + mod6months$coefficients[2] * prev_mod + 0 + mod6months$coefficients[4] * prev_mod * 0 
y_nets_early_ye = mod6months$coefficients[1] + mod6months$coefficients[2] * prev_mod + mod6months$coefficients[3] + mod6months$coefficients[4] * prev_mod * 1

plot(df$V3 ~ df$V8,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 6-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V3[df$nets_early == 0] ~ df$V8[df$nets_early == 0], col = "grey40", pch = 19)
points(df$V3[df$nets_early == 1] ~ df$V8[df$nets_early == 1], col = "grey40")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "grey40")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "grey40")
axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))

##
mod12months = glm(df$V4 ~ df$V9 * df$nets_early)
prev_mod = seq(0,max(df$V4,na.rm=TRUE),0.01)
y_nets_early_no = mod12months$coefficients[1] + mod12months$coefficients[2] * prev_mod + 0 + mod12months$coefficients[4] * prev_mod * 0 
y_nets_early_ye = mod12months$coefficients[1] + mod12months$coefficients[2] * prev_mod + mod12months$coefficients[3] + mod12months$coefficients[4] * prev_mod * 1

plot(df$V4 ~ df$V9,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 12-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V4[df$nets_early == 0] ~ df$V9[df$nets_early == 0], col = "blue", pch = 19)
points(df$V4[df$nets_early == 1] ~ df$V9[df$nets_early == 1], col = "blue")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "blue")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "blue")
axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))


##

mod18months = glm(df$V5 ~ df$V10 * df$nets_early)
prev_mod = seq(0,max(df$V5,na.rm=TRUE),0.01)
y_nets_early_no = mod18months$coefficients[1] + mod18months$coefficients[2] * prev_mod + 0 + mod18months$coefficients[4] * prev_mod * 0 
y_nets_early_ye = mod18months$coefficients[1] + mod18months$coefficients[2] * prev_mod + mod18months$coefficients[3] + mod18months$coefficients[4] * prev_mod * 1

plot(df$V5 ~ df$V10,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 18-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V5[df$nets_early == 0] ~ df$V10[df$nets_early == 0], col = "darkblue", pch = 19)
points(df$V5[df$nets_early == 1] ~ df$V10[df$nets_early == 1], col = "darkblue")
lines(y_nets_early_no ~ prev_mod, lty = 1, col = "darkblue")
lines(y_nets_early_ye ~ prev_mod, lty = 2, col = "darkblue")

axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))


##
mod25months = glm(df$V6 ~ df$V11 * df$nets_early)
prev_mod = seq(0,max(df$V6,na.rm=TRUE),0.01)
y_nets_early_no = mod25months$coefficients[1] + mod25months$coefficients[2] * prev_mod + 0 + mod25months$coefficients[4] * prev_mod * 0 
y_nets_early_ye = mod25months$coefficients[1] + mod25months$coefficients[2] * prev_mod + mod25months$coefficients[3] + mod25months$coefficients[4] * prev_mod * 1

plot(df$V6 ~ df$V11,pch="",ylim=c(0,1),xlim = c(0,1),
     xaxt="n",yaxt="n",
     main = "At 25-months",
     ylab = "Prevalence 2 to 10 years as measured",
     xlab = "Prevalence 2 to 10 years as modelled")
points(df$V6[df$nets_early == 0] ~ df$V11[df$nets_early == 0],  pch = 19)
points(df$V6[df$nets_early == 1] ~ df$V11[df$nets_early == 1], )
lines(y_nets_early_no ~ prev_mod, lty = 1)
lines(y_nets_early_ye ~ prev_mod, lty = 2)

axis(1, at = seq(0,1,0.2),labels = seq(0,100,20))
axis(2, las = 2, at = seq(0,1,0.2),labels = seq(0,100,20))

