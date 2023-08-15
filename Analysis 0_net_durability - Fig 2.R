###################################################
##
## 1 Site parameter set up

## Work out seasonal pattern relative to timing of net distributions

## Fit function to net use estimates for each of 104 clusters


library(rstan)
library(adegenet)

## Looping through the 104 clusters 
PARMS_USAGE = expand.grid(1:1000)

## Specify the times when data were observed
time_obs = c(0,6/12,12/12,18/12,25/12)

# Coverage 6-months: 71% ## From supplementary figures Staedke et al 2020
# Coverage 12-months: 63%
# Coverage 18-months: 49%
## Specify the data observed (including the uncertainties as this provides more data points for a better fit)
## where more specific data for different clusters are available this can be improved.

net_data = read.csv("raw data/net_use.csv",header=TRUE)
net_data25 = read.csv("raw data/For Ellie_requested info for 25m_final.csv",header=TRUE)
net_dataAll = merge(net_data, net_data25, by.x = "cluster", all.x = TRUE)
net_dataAll$nets_25m = net_dataAll$LLIN_use/100
net_dataAll = net_dataAll[complete.cases(net_dataAll),]

for(i in 1:nrow(net_dataAll)){

  ## ORIGINAL USES NET_DATA ONLY AND 6-18 MONTHS
  standard_net_usage = c(0.95,as.numeric(net_dataAll[i,c(4,6,8,11)]))
  
  # standard_net_usage = c(0.71,0.63,0.49,
  #                        0.62,0.58,0.44, ##these are the estimated uncertainty ranges
  #                        0.72,0.68,0.54) ##these are the estimated uncertainty ranges
  
  ## Transform data for the model fitting
  y_standard_net_usage = log(standard_net_usage)
  
  ## Specify a time series to project outcomes across
  time_m = seq(0,3,0.01)
  
  ## Specify the exponential decay model in Rstan
  stanmodelcode <- "
  data {
  int<lower=0> N;
  int<lower=0> N2; //the size of the new_X matrix
  vector[N] y;
  vector[N] x;
  vector[N2] New_x;
  
  }
  parameters {
  real beta0;
  real beta1;
  real sigma;
  }
  transformed parameters {
  vector[N] m;
  m = beta0 + beta1 * x;
  } 
  model {
  // priors
  beta0 ~ cauchy(0, 10); 
  beta1 ~ cauchy(0, 10); 
  
  // likelihood
  y ~ normal(m, sigma);   
  }
  generated quantities {
  vector[N2] y_pred;
  y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
  }
  "
  stanDso <- stan_model(model_code = stanmodelcode) 
  
  ## Put data into a list
  dat_standard <- list(N = length(time_obs), 
                       N2 = length(seq(0,3,0.01)),
                       y = y_standard_net_usage, 
                       x = time_obs,
                       New_x = seq(0,3,0.01))
  
  ## Run the statistical model
  fit <- sampling(stanDso, 
                  data = dat_standard, ## specify data
                  iter = 4000,         ## speify iterations
                  warmup=2000)         ## specify warm up (default is 4 chains here)
  
  
  ## plotting the posterior distribution for the parameters
  post_beta<-As.mcmc.list(fit,pars="beta0")
  plot(post_beta)
  
  ## gradient is fit to the data for alpha
  ## standard_net_usage ~ exp(-alpha*time_obs)
  b0 <- extract(fit, 'beta0')
  b0<- unlist(b0, use.names=FALSE)
  b1 <- extract(fit, 'beta1')
  b1<- unlist(b1, use.names=FALSE)
  
  ## Back translate output estimates
  y_predicted = mean(b0) + mean(b1)*time_m; 
  y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)
  
  ## and select some uncertainty for the Bayesian estimates
  y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
  y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)
  
  ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
  parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)]) ##
  parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
  
  PARMS_USAGE[,i] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)
  
}



standard_net_usage = c(0.95,as.numeric(net_data[i,c(4,6,8)]))

# standard_net_usage = c(0.71,0.63,0.49,
#                        0.62,0.58,0.44, ##these are the estimated uncertainty ranges
#                        0.72,0.68,0.54) ##these are the estimated uncertainty ranges

## Transform data for the model fitting
y_standard_net_usage = log(standard_net_usage)

## Specify a time series to project outcomes across
time_m = seq(0,3,0.01)

## Specify the exponential decay model in Rstan
stanmodelcode <- "
data {
int<lower=0> N;
int<lower=0> N2; //the size of the new_X matrix
vector[N] y;
vector[N] x;
vector[N2] New_x;

}
parameters {
real beta0;
real beta1;
real sigma;
}
transformed parameters {
vector[N] m;
m = beta0 + beta1 * x;
} 
model {
// priors
beta0 ~ cauchy(0, 10); 
beta1 ~ cauchy(0, 10); 

// likelihood
y ~ normal(m, sigma);   
}
generated quantities {
vector[N2] y_pred;
y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
}
"
stanDso <- stan_model(model_code = stanmodelcode) 

## Put data into a list
dat_standard <- list(N = length(time_obs), 
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01))

## Run the statistical model
fit <- sampling(stanDso, 
                data = dat_standard, ## specify data
                iter = 4000,         ## speify iterations
                warmup=2000)         ## specify warm up (default is 4 chains here)


## plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

## Back translate output estimates
y_predicted = mean(b0) + mean(b1)*time_m; 
y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)

## and select some uncertainty for the Bayesian estimates
y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)


###################################################
##
## Some fits are poorer than others
## Check

## Final output plotted to confirm the fit
par(mfrow = c(1,1))
plot(standard_net_usage[1:5] ~ ## plotting the mean estimates first
       time_obs[1:5],pch="",
     ylab="Households with at least one net per two occupants (%)", ## specify data used
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.4,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.4,cex.axis=1.4)

# polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
# lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)

# ## add the range in the uncertainty from the observed data for the trial if available
# for(i in 1:4){
#   segments(x0=time_obs[i],x1=time_obs[i],
#            y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
# }

for(i in 1:90){
  standard_net_usage = c(0.95,as.numeric(net_dataAll[i,c(4,6,8,11)]))
  points(standard_net_usage[1:5] ~ time_obs[1:5])
  
}

# ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
# parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.2) & b1 <= quantile(b1,0.8)]) ##
# parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

# PARMS_USAGE[,1] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)


for(i in c(1:54,56:90)){
  # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE[,i])
  D_LOW = quantile(PARMS_USAGE[,i],0.01)
  D_UPP = quantile(PARMS_USAGE[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  lines(aa ~ time_m,col="blue")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # 
}

##odd fits
# for(i in c(9,18,30,52,69,89,100,104)){
  for(i in c(55)){
    # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE[,i])
  D_LOW = quantile(PARMS_USAGE[,i],0.01)
  D_UPP = quantile(PARMS_USAGE[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  lines(aa ~ time_m,col="red")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # 
}
oddies = c(9,18,30,52,65,66,76,69,83,89,100,104)
standard_net_usage_outliers = array(dim=c(4,length(oddies)))
for(i in 1:length(oddies)){
  standard_net_usage_outliers[,i] = c(0.95,as.numeric(net_data[oddies[i],c(4,6,8)]))
}


PARMS_USAGE2 = array(dim=c(1000,104))
for(i in 1:104){
  PARMS_USAGE2[,i] = ifelse(PARMS_USAGE[,i] > 15,15,PARMS_USAGE[,i])
  
}

PARMS_USAGE3 = array(dim=c(1000,104))
for(i in 1:104){
  PARMS_USAGE3[,i] = ifelse(PARMS_USAGE2[,i] < 0,1.5,PARMS_USAGE2[,i])
  
}
D = colMeans(PARMS_USAGE3)
hist(D,breaks = 50)

write.csv(PARMS_USAGE3)

######################
## Final output plotted to confirm the fit
par(mfrow = c(1,2))
plot(standard_net_usage[1:4] ~ ## plotting the mean estimates first
       time_obs[1:4],pch="",
     # ylab="Households with >= net per two occupants (%)", ## specify data used
     ylab="Estimated adherrence to net use (%)", ## specify data used
     xlab="Time in years since mass distribution",yaxt="n",ylim=c(0,1),cex.lab=1,cex.axis=1,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1,cex.axis=1)

# polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
# lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)

# ## add the range in the uncertainty from the observed data for the trial if available
# for(i in 1:4){
#   segments(x0=time_obs[i],x1=time_obs[i],
#            y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
# }
net_timing = read.csv("data/net_distribution.csv",header=TRUE)
## Additional vectors
net_timing$net_type_INPUTA = ifelse(net_timing$Net_Type == "Olyset Net", 1,
                                    ifelse(net_timing$Net_Type == "PermaNet 2.0", 2,
                                           ifelse(net_timing$Net_Type == "Olyset Plus", 3, 4)))
net_timing$net_type_INPUTB = ifelse(net_timing$Net_Type == "Olyset Net", 1,
                                    ifelse(net_timing$Net_Type == "PermaNet 2.0", 2,
                                           ifelse(net_timing$Net_Type == "Olyset Plus", 1, 2)))

cols_nets = c("blue","darkred","darkgreen","orange")
for(i in 1:104){
  standard_net_usage = c(0.95,as.numeric(net_data[i,c(4,6,8)]))
  points(standard_net_usage[1:4] ~ time_obs[1:4],
         col=cols_nets[net_timing$net_type_INPUTA[i]])
  
}

# ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
# parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.2) & b1 <= quantile(b1,0.8)]) ##
# parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

# PARMS_USAGE[,1] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)

D = colMeans(PARMS_USAGE3)
colours_nets = c(adegenet::transp("blue",0.2),"darkgreen",adegenet::transp("darkred"),"orange")
as.numeric(net_data$net_type)

colour_match = ifelse(as.numeric(net_data$net_type) == 1,colours_nets[1],
                      ifelse(as.numeric(net_data$net_type) == 2,colours_nets[2],
                             ifelse(as.numeric(net_data$net_type) == 3,colours_nets[3],colours_nets[4])))

for(i in 1:104){
  # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE3[,i])
  D_LOW = quantile(PARMS_USAGE3[,i],0.01)
  D_UPP = quantile(PARMS_USAGE3[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  # lines(aa ~ time_m,col="blue")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # polygon(c(time_m,rev(time_m)),
  #         c(aL,rev(aU)),col=adegenet::transp("blue",0.05),border = NA) 
  lines(aa ~ time_m,col=cols_nets[net_timing$net_type_INPUTA[i]])
  
}
for(i in 1:104){
  standard_net_usage = c(0.95,as.numeric(net_data[i,c(4,6,8)]))
  points(standard_net_usage[1:4] ~ time_obs[1:4],pch=19,col=colour_match[i])
  
}
D = colMeans(PARMS_USAGE3)
# colours_nets = c(adegenet::transp("orange",0.2),"orange",adegenet::transp("darkgreen"),"darkgreen")
as.numeric(net_data$net_type)

colour_match = ifelse(as.numeric(net_data$net_type) == 1,colours_nets[1],
                      ifelse(as.numeric(net_data$net_type) == 2,colours_nets[2],
                             ifelse(as.numeric(net_data$net_type) == 3,colours_nets[3],colours_nets[4])))
hist(D,breaks = 104,main="",yaxt="n",border=F,
     ylab = "Frequency",xlab="Adherence to net use (years)",col=colour_match)


abline(v=median(D),col="blue",lwd=2)
axis(2,las=2,at=c(0:7))
abline(v=quantile(D,0.9),col="blue",lwd=2,lty=2)
abline(v=quantile(D,0.1),col="blue",lwd=2,lty=2)
