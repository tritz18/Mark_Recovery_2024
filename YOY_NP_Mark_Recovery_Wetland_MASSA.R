
#########################################
#   Mark-Recovery Model - MASSA  #
#       Script by Thornton Ritz         # 
#########################################

#### Read in associated packages ####
library(tidyverse)
library(jagsUI)
library(readxl)
library(dplyr)

#### Data preparation ##
## YOY NP is my capture history file used to create M-array, ##
## associated site variable used to split m-arrays up by wetland ####
setwd("/Users/thorn/OneDrive/Desktop/Mark_Recovery_Wetland")
YOY_NP<- read_excel("Massa_MARRAY.xlsx")
YOY_NP$Site<- as.factor(YOY_NP$Site)

#### Split up capture history by wetland to run for Type specific estimates ##
## Generate a matrix for each grouping that is then used for each MARRAY ####

## Non-excavated ##
NON_EXC <- YOY_NP |>
  filter(Site %in% c("Non-excavated")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(NON_EXC)
NON_EXC<- as.numeric(NON_EXC)
dim(NON_EXC) <- x

## Channel ##
CHAN <- YOY_NP |>
  filter(Site %in% c("Channel")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(CHAN)
CHAN<- as.numeric(CHAN)
dim(CHAN) <- x

# SPawning pool ##
SP <- YOY_NP |>
  filter(Site %in% c("Spawning_Pool")) |>
  select(-Site) |>
  as.matrix() 
x<-dim(SP)
SP<- as.numeric(SP)
dim(SP) <- x


#### Removal of site variable, split up by wetland already #### 
YOY_NP <- YOY_NP |>
  select(-Site)

YOY_NP<- as.matrix(YOY_NP)
x<- dim(YOY_NP)
dim(YOY_NP) <- x

#### Recovery Rate Contributions ##
## Did we have nets in the water? ( Yes - 1 No - 0) (42 days post stocking of unsampled growth period) ##
## Did we capture a YOY NP in those nets? ( Yes - 1 No - 0) (Net) ####

Net <- read_excel("CATCH_REC_MASSA.xlsx")
Net<- Net |>
  select(-Site)
Net<- as.matrix(Net)

#### First marking day ####
get.first <- function(x) min(which(x!=0))
f <- apply(YOY_NP, 1, get.first)

#### The mark-recovery model fitted with the multinomial likelihood  ##
## Define function to create an m-array based for mark-recovery (MR) data ####
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

#### Replace empty columns with zeros as place holders for MARRAY ##
## only one stocking occasion per site, new fish are not stocked every sampling day ####

marr1 <- marray.dead(CHAN)
marr1 <- marr1 %>% replace(is.na(.), 0)

marr2 <- marray.dead(NON_EXC)
marr2 <- marr2 %>% replace(is.na(.), 0)

marr3 <- marray.dead(SP)
marr3 <- marr3 %>% replace(is.na(.), 0)

#### Specify model in BUGS language ####

sink("mr-mnl.jags")
cat("
model {

# Priors and constraints
 
for (s in 1:n.sites) {
for (t in 1:n.occasions){
    surv[s,t] <- surv.mean[s]

  
}
  surv.mean[s] ~ dunif(0,1)

# Mean daily surv. by wetland, period surv. calculation   
surv.int.mean[s] <- surv.mean[s]^46
} 

# Capture indicates prob. to capture ind. (does it go into net), vague prior  
capture ~ dunif(0,1)

# Period capture calculation 
int.capture <- 1-(1-capture)^46

for(s in 1:n.sites) {
for (t in 1:n.occasions){

  rec[s,t] <- (Net[s,t] * capture)

}
  rec.mean[s] <- mean(rec[s,]) 

}

# Define the multinomial likelihood
for (t in 1:n.occasions){
   marr1[t,1:(n.occasions+1)] ~ dmulti(pr1[t,], rel1[t])
   marr2[t,1:(n.occasions+1)] ~ dmulti(pr2[t,], rel2[t])
   marr3[t,1:(n.occasions+1)] ~ dmulti(pr3[t,], rel3[t])

   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:n.occasions){
   pr1[t,t] <- (1-surv[1,t])*rec[1,t]
   pr2[t,t] <- (1-surv[2,t])*rec[2,t]
   pr3[t,t] <- (1-surv[3,t])*rec[3,t]

   # Above main diagonal
   for (j in (t+1):n.occasions){
      pr1[t,j] <- prod(surv[1,t:(j-1)])*(1-surv[1,j])*rec[1,j]
      pr2[t,j] <- prod(surv[2,t:(j-1)])*(1-surv[2,j])*rec[2,j]
      pr3[t,j] <- prod(surv[3,t:(j-1)])*(1-surv[3,j])*rec[3,j]

      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr1[t,j] <- 0
      pr2[t,j] <- 0
      pr3[t,j] <- 0

      } #j
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr1[t,n.occasions+1] <- 1-sum(pr1[t,1:n.occasions])
   pr2[t,n.occasions+1] <- 1-sum(pr2[t,1:n.occasions])
   pr3[t,n.occasions+1] <- 1-sum(pr3[t,1:n.occasions])

   } #t
}
",
fill = TRUE)
sink()


#### Bundle data ####
jags.data <- list(marr1 = marr1, marr2 = marr2, marr3 = marr3, 
                  n.occasions = dim(marr1)[2]-1, 
                  rel1 = rowSums(marr1), rel2 = rowSums(marr2), rel3 = rowSums(marr3),
                  Net=Net,n.sites=3)

#### Initial values ####
inits <- function(){list(capture=runif(1,0,1),
                         surv.mean=rep(runif(1,0,1),3))
                        }  

# Parameters monitored
parameters <- c("surv.mean", "rec.mean", "capture", "surv.int.mean",
                 "int.capture")
# MCMC settings # this are low so I can run it quickly to assess ##
ni <- 5000
nt <- 2
nb <- 2500
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
