
#########################################
#   Mark Recovery Model - YOY NP 2022   #
#       Script by Thornton Ritz         # 
#########################################

## Read in associated packages ##
library(tidyverse)
library(jagsUI)
library(readxl)
library(dplyr)
library(jagshelper)

## Data preparation ##
# YOY NP is my capture history file, for all fish there is an associated site variable
# that is used to identify the groupings to create separate M-ARRAY by grouping 
# Each wetland (4 wetlands) has 4 sites, 2 reference (no restoration) and 2 spawning pools (restoration). 
# Combined reference sites and spawning pools within wetlands 
# 4 days stocking + 42 (no nets) + 18 days (nets)

setwd("/Users/tritz1/Desktop/Mark_Recovery_Wetland")
YOY_NP<- read_excel("YOY_NP_M_Array_Wetland_2.xlsx")
YOY_NP$Site<- as.factor(YOY_NP$Site)

## Environmental Covariate Data Frames ##

DO<- read_excel("YOY_NP_DO_mean_Wetland.xlsx")
DO<-as.matrix(DO)

TEMP<- read_excel("YOY_NP_TEMP_mean_Wetland.xlsx")
TEMP<-as.matrix(TEMP)

## Create grouping variables ##
# Generate a matrix for each grouping that is then used for each MARRAY #

## Chippewa Creek ##
Chip <- YOY_NP |>
  filter(Site %in% c("CHIP_REF3", "CHIP_REF4", "CHIP_SP2")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(Chip)
Chip<- as.numeric(Chip)
dim(Chip) <- x

## Cranberry Creek ##
Cran <- YOY_NP |>
  filter(Site %in% c("CM_REF4" , "CM_REF5", "CM_SP1", "CM_SP2")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(Cran)
Cran<- as.numeric(Cran)
dim(Cran) <- x

## French Creek ##
French<- YOY_NP |>
  filter(Site %in% c("FC_REF6", "FC_REF7", "FC_SP5", "FC_SP6")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(French)
French<- as.numeric(French)
dim(French) <- x

## Point Vivian ##
PV<- YOY_NP |>
  filter(Site %in% c("PV_REF1", "PV_REF2", "PV_SP4", "PV_SP5")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(PV)
PV<- as.numeric(PV)
dim(PV) <- x



## Removal of site 
YOY_NP <- YOY_NP |>
  select(-Site)

YOY_NP<- as.matrix(YOY_NP)
x<- dim(YOY_NP)
dim(YOY_NP) <- x

## Recovery Rate Contributions ##
# Did we have nets in the water? ( Yes - 1 No - 0) (42 days post stocking of unsampled growth period)
# Did we capture a YOY NP in those nets? ( Yes - 1 No - 0) (Net)
# How many of those caught did we keep? ( 0.0 - 1.0   % kept to be recovered) (RR)
# How many of kept fish had a OTC mark indicating recapture (0.0 - 1.0   % marked) (Mark)
# 8x64 matrix for 8 groupings and 64 days ##

Net <- read_excel("YOY_NP_CATCH_REC_Wetland.xlsx")
Net<- Net |>
  select(-Site)
Net<- as.matrix(Net)

RR<- read_excel("YOY_NP_PERCENT_RETAIN_Wetland.xlsx")
RR<- RR |>
  select(-Site)
RR<- as.matrix(RR)

Mark<- read_excel("YOY_NP_Mark_Rate_Wetland_2.xlsx")
Mark<- Mark |>
  select(-Site)
Mark<- as.matrix(Mark)

## First marking day ##
get.first <- function(x) min(which(x!=0))
f <- apply(YOY_NP, 1, get.first)

# 8.3. The mark-recovery model fitted with the multinomial likelihood 
# 8.3.1. Constant parameters
# Define function to create an m-array based for mark-recovery (MR) data
marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Create vector with occasion of marking 
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  # Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  # Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# Replace empty columns with zeros as place holders for MARRAY 
# only one stocking occasion per site, new fish are not stocked every sampling day#

marr1 <- marray.dead(Chip)
marr1 <- marr1 %>% replace(is.na(.), 0)
View(marr)

marr2 <- marray.dead(Cran)
marr2 <- marr2 %>% replace(is.na(.), 0)

marr3 <- marray.dead(French)
marr3 <- marr3 %>% replace(is.na(.), 0)

marr4 <- marray.dead(PV)
marr4 <- marr4 %>% replace(is.na(.), 0)

site<- rep(0,61)

## Specify model in BUGS language 

sink("mr-mnl.jags")
cat("
model {

# Priors and constraints
 
for (s in 1:n.sites) {
for (t in 1:n.occasions){
# Linear predictor 
surv[s,t] <- ilogit(logit_surv_mu[s] +beta_DO * DO[s,t]+
                          beta_TEMP * TEMP[s,t])
}
  
surv.mean[s] <- mean(surv[s,])
surv.int.mean[s] <- surv.mean[s]^64

logit_surv_mu[s] ~ dnorm(grand.surv.logit.mu, grand.surv.tau)


} 
capture ~ dunif(0,1)
int.capture <- 1-(1-capture)^64
beta_DO ~ dnorm(0,0.001)
beta_TEMP ~ dnorm(0,0.001)
grand.surv.logit.mu ~ dnorm(0,0.001)
grand.surv.tau <- 1/grand.surv.sd^2
grand.surv.sd ~ dunif(0,10)

for(s in 1:n.sites) {
for (t in 1:n.occasions){

  rec[s,t] <- (Net[s,t] * capture)

}
  rec.mean[s] <- mean(rec[s,]) # Mean daily recovery rate...will be low because many 0s

}

# Define the multinomial likelihood
for (t in 1:n.occasions){
   marr1[t,1:(n.occasions+1)] ~ dmulti(pr1[t,], rel1[t])
   marr2[t,1:(n.occasions+1)] ~ dmulti(pr2[t,], rel2[t])
   marr3[t,1:(n.occasions+1)] ~ dmulti(pr3[t,], rel3[t])
   marr4[t,1:(n.occasions+1)] ~ dmulti(pr4[t,], rel4[t])

   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:n.occasions){
   pr1[t,t] <- (1-surv[1,t])*rec[1,t]
   pr2[t,t] <- (1-surv[2,t])*rec[2,t]
   pr3[t,t] <- (1-surv[3,t])*rec[3,t]
   pr4[t,t] <- (1-surv[4,t])*rec[4,t]

   # Above main diagonal
   for (j in (t+1):n.occasions){
      pr1[t,j] <- prod(surv[1,t:(j-1)])*(1-surv[1,j])*rec[1,j]
      pr2[t,j] <- prod(surv[2,t:(j-1)])*(1-surv[2,j])*rec[2,j]
      pr3[t,j] <- prod(surv[3,t:(j-1)])*(1-surv[3,j])*rec[3,j]
      pr4[t,j] <- prod(surv[4,t:(j-1)])*(1-surv[4,j])*rec[4,j]

      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr1[t,j] <- 0
      pr2[t,j] <- 0
      pr3[t,j] <- 0
      pr4[t,j] <- 0

      } #j
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr1[t,n.occasions+1] <- 1-sum(pr1[t,1:n.occasions])
   pr2[t,n.occasions+1] <- 1-sum(pr2[t,1:n.occasions])
   pr3[t,n.occasions+1] <- 1-sum(pr3[t,1:n.occasions])
   pr4[t,n.occasions+1] <- 1-sum(pr4[t,1:n.occasions])

   } #t
}
",
fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr1 = marr1, marr2 = marr2, marr3 = marr3, marr4 = marr4, 
                  n.occasions = dim(marr1)[2]-1, 
                  rel1 = rowSums(marr1), rel2 = rowSums(marr2), rel3 = rowSums(marr3), rel4 = rowSums(marr4),
                  Net=Net, DO=DO, TEMP=TEMP, RR=RR, Mark=Mark,
                  n.sites=4)

# Initial values

inits <- function(){list(logit_surv_mu=rnorm(4,0,1),
                         grand.surv.logit.mu=rnorm(1),
                         grand.surv.sd=runif(1),
                      beta_DO=rnorm(1,0,1),
                      beta_TEMP=rnorm(1,0,1),
                      capture=runif(1,0,1)
                         )
                        }  

# Parameters monitored
parameters <- c("surv.mean", "rec.mean", "beta_DO", "beta_TEMP", "logit_surv_mu", "capture",
                "surv.int.mean", "capture", "int.capture",
                "grand.surv.logit.mu", "grand.surv.sd")
# MCMC settings # this are low so I can run it quickly to assess ##
ni <- 50000
nt <- 3
nb <- 25000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr <- jags(jags.data, inits, parameters, "mr-mnl.jags", n.chains = nc, n.thin = nt,
           n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors

print(mr, digits = 5)
plot(mr)

mr.sum<- mr$summary

densityplot(mr)
library(lattice)
xyplot(mr)
write.csv(mr.sum, "model.sum.COV_2.csv")



Surv<- as.data.frame(mr$sims.list$surv.mean)


write.csv(Surv, "Survival_Simulations_2.csv")






TEMP_rng<- range(TEMP[3,])
TEMP_seq<- seq(TEMP_rng[1], TEMP_rng[2], length.out=100)
TEMP_std<- (TEMP_seq-mean(TEMP[3,]))/ sd(TEMP[3,])
head(TEMP_std)

surv_est<- mr$mean$logit_surv_mu[3] + mr$mean$beta_TEMP[3] * TEMP_std
surv_est<- plogis(surv_est)
head(surv_est)  

view(surv_est)
plot(TEMP_seq, surv_est, type = 'l', xlab="TEMP", ylab="Surv")

nsamples<- mr$mcmc.info$n.samples

surv_post<- matrix(NA, nrow = nsamples, ncol=100)
 for(i in 1:100){
   surv_post[,i] <- mr$sims.list$logit_surv_mu[3] + mr$sims.list$beta_TEMP[3] * TEMP_std[i]
 }

surv_post<- plogis(surv_post) 

surv_lower<- apply(surv_post, 2, quantile, 0.025)
surv_upper<- apply(surv_post, 2, quantile, 0.975)


plot.new()
plot(TEMP_seq, surv_est, type='l', xlab="TEMP", ylab="Survival")
polygon(c(TEMP_seq, rev(TEMP_seq)),
        c(surv_lower, rev(surv_upper)),
        col="grey90", border=NA)
lines(TEMP_seq, surv_est)
