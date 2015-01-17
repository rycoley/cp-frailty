##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Estimating Effectiveness in HIV Prevention Trials with a Hierarachical Compound Poisson Frailty Model""
##### part of dissertation work at University of Washington


### For primary simulation
#correct priors
#export environment variables proporiton not at risk (PNR) and event rate (ER) in command line of qsub 
#run in parallel with multicore on SGE computing cluster

#G<-Sys.getenv("G")
#G<-as.numeric(G)
G<-1 #I use this to set the seed if I want to run multiple simulations with same PNR and ER

## Set probability *not* at risk
PNR<-Sys.getenv("PNR")
PNR<-as.numeric(PNR)
PNR

##Set event rate
ER<-Sys.getenv("ER")
ER<-as.numeric(ER)
ER

library(multicore)


rho<--log(PNR) 

#priors 
#rate and shape parameters for beta prior on probability not at risk 
rho.b<-2
rho.a<-(rho.b*PNR)/(1-PNR) 

#shape and rate parameters for gamma prior on baseline hazard
lam.s<-0.05*50
lam.r<-50

#shape and rate parameters for gamma prior on eta
eta.s<-1
eta.r<-0.5

S<- 500 #number of simulations
A<-10000 #first gibbs sampler iteration to save
B<- 25000 #total number of sampler iterations
ke<- 1 #number thinned


source("cpe-source.R") #log (posterior) likelihoods and sampling functions
source("cpe-to-run.R")	#do.one() function does data generation and sampling algorithm 
	
write.csv( "completed", paste("status-",PNR,"-",ER,".csv",sep="")) #set up file to monitor which samplers have finished
	
	
res<-mclapply(1:S,do.one,mc.cores=20) #multicore apply function for 20 cores 


res.beta<-res.lam<-res.eta<-res.rho<-matrix(nrow=S,ncol=4) #matrices for saimpling results
res.HD<-matrix(nrow=S,ncol=3) #matrix for hellinger distance

#fill matrices with results from S simulations
for(s in 1:S){
	res.beta[s,]<-res[[s]]$bet.res
	res.lam[s,]<-res[[s]]$lam.res
	res.rho[s,]<-res[[s]]$rho.res
	res.eta[s,]<-res[[s]]$eta.res
	res.HD[s,]<-res[[s]]$HD}


#write results to csv file
beta.name<-paste("beta-cpe-cpt-",ER,"-",PNR,".csv",sep="")
lambda.name<-paste("lam-cpe-cpt-",ER,"-",PNR,".csv",sep="")
eta.name<-paste("eta-cpe-cpt-",ER,"-",PNR,".csv",sep="")
rho.name<-paste("rho-cpe-cpt-",ER,"-",PNR,".csv",sep="")
HD.name<-paste("HD-cpe-cpt-",ER,"-",PNR,".csv",sep="")

write.csv(res.beta, beta.name) 
write.csv(res.lam, lambda.name) 
write.csv(res.eta, eta.name) 
write.csv(res.rho, rho.name) 
write.csv(res.HD, HD.name) 





