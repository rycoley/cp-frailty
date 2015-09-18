##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Estimating Effectiveness in HIV Prevention Trials with a Hierarachical Compound Poisson Frailty Model""
##### part of dissertation work at University of Washington


### For primary simulation
#correct priors
#run in parallel with multicore on SGE computing cluster


### WORKFLOW: Get simulation settings, load libraries,  define parameters and hyperparameters for data generation, define sampling algorithm settings, source necessary R scripts, run one replication of simulation, save results in csv files




###get simulation settings

#get task id to set seed
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))

#export environment variables proporiton not at risk (PNR) and event rate (ER) in command line of qsub 
## Set probability *not* at risk
PNR<-Sys.getenv("PNR")
PNR<-as.numeric(PNR)
PNR

##Set event rate
ER<-Sys.getenv("ER")
ER<-as.numeric(ER)
ER


### load libraries
library(coda)
library(survival)


### define parameters, hyperparameters for data generation
rho<--log(PNR) 

#hyperparameters 
#rate and shape parameters for beta prior on probability not at risk 
rho.b<-2
rho.a<-(rho.b*PNR)/(1-PNR) 

#shape and rate parameters for gamma prior on baseline hazard
lam.s<-0.05*50
lam.r<-50

#shape and rate parameters for gamma prior on eta
eta.s<-1
eta.r<-0.5



### define sampler settings
A<-10000 #first gibbs sampler iteration to save
B<- 25000 #total number of sampler iterations
ke<- 1 #number thinned


### sources R scripts
source("cpe-source.R") #log (posterior) likelihoods and sampling functions
source("cpe-to-run.R")	#do.one() function does data generation and sampling algorithm 
	

res<-do.one(seed=SEED)



beta.name<-paste("beta-cpe-cpt-",ER,"-",PNR,"-",SEED,".csv",sep="")
lambda.name<-paste("lam-cpe-cpt-",ER,"-",PNR,"-",SEED,".csv",sep="")
eta.name<-paste("eta-cpe-cpt-",ER,"-",PNR,"-",SEED,".csv",sep="")
rho.name<-paste("rho-cpe-cpt-",ER,"-",PNR,"-",SEED,".csv",sep="")
HD.name<-paste("HD-cpe-cpt-",ER,"-",PNR,"-",SEED,".csv",sep="")

write.csv(res$bet.res, beta.name) 
write.csv(res$lam.res, lambda.name) 
write.csv(res$eta.res, eta.name) 
write.csv(res$rho.res, rho.name) 
write.csv(res$HD, HD.name) 


