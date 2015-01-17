##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Estimating Effectiveness in HIV Prevention Trials with a Hierarachical Compound Poisson Frailty Model""
##### part of dissertation work at University of Washington

### Challenge simulations
#Non-informative priors
#export environment variables proporiton not at risk (PNR) and event rate (ER) in command line of qsub 
#run in parallel with multicore on SGE computing cluster

#code similar to cpe-cpt.R, see for annotated code

#G<-Sys.getenv("G")
#G<-as.numeric(G)
G<-1


PNR<-Sys.getenv("PNR")
PNR<-as.numeric(PNR)
PNR

ER<-Sys.getenv("ER")
ER<-as.numeric(ER)
ER

library(multicore)

rho<--log(PNR) 

rho.b<-1
rho.a<-1

S<- 500
A<-10000
B<- 25000 
ke<- 1

source("cpe-source3.R")
source("cpe-to-run3.R")	


res<-mclapply(1:S,do.one,mc.cores=20)

res.beta<-res.lam<-res.eta<-res.rho<-matrix(nrow=S,ncol=4)
res.HD<-matrix(nrow=S,ncol=3)

for(s in 1:S){
	res.beta[s,]<-res[[s]]$bet.res
	res.lam[s,]<-res[[s]]$lam.res
	res.rho[s,]<-res[[s]]$rho.res
	res.eta[s,]<-res[[s]]$eta.res
	res.HD[s,]<-res[[s]]$HD}


beta.name<-paste("beta-cpe-cpt3-",ER,"-",PNR,".csv",sep="")
lambda.name<-paste("lam-cpe-cpt3-",ER,"-",PNR,".csv",sep="")
eta.name<-paste("eta-cpe-cpt3-",ER,"-",PNR,".csv",sep="")
rho.name<-paste("rho-cpe-cpt3-",ER,"-",PNR,".csv",sep="")
HD.name<-paste("HD-cpe-cpt3-",ER,"-",PNR,".csv",sep="")

write.csv(res.beta, beta.name) 
write.csv(res.lam, lambda.name) 
write.csv(res.eta, eta.name) 
write.csv(res.rho, rho.name) 
write.csv(res.HD, HD.name) 



