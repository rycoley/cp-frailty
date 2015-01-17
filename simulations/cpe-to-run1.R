##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Estimating Effectiveness in HIV Prevention Trials with a Hierarachical Compound Poisson Frailty Model""
##### part of dissertation work at University of Washington

#For estimation with informative prior, incorrect on PNR
#Accompanies cpe-cpt1.R and cpe-source1.R

#This file runs a single simulation, including data generation based on random seed and runs gibbs sampler to get parameter estimate

#function is similar to that defined in cpe-to-run.R, see there for annotated code

library(compiler)
enableJIT(1)

library(survival)

do.one<-function(seed){

print(c(seed,"assigned"))

true.rho<-rho

true.lam<-lambda<-0.05
true.bet<-beta<-log(0.5)

true.eta<-eta<-2
nu<-rho*eta

set.seed(seed)

N<-rep(0,n) 
for(i in 1:n){N[i]<-rpois(1,rho)}

Z<-rep(0,n) 
for(i in 1:n){Z[i]<-rgamma(1, shape=N[i]*eta, rate=nu)}

Y<-rep(NA,n)
for(i in 1:n){if(!N[i]==0){
		Y[i]<-rexp(1, rate=(lambda*exp(X[i]*beta)*Z[i]))	}
		else {Y[i]<-1000}	}

Tf<-quantile(Y,ER) 
delta<-rep(0,n)
delta[Y<=Tf]<-1		
Y[Y>Tf]<-Tf

Tf.mo<-c(1:length(months))[months==min(months[Tf<months])]

pdf.true<-get.surv.pdf(cdf=get.surv.mo(Y=get.ymo(Y, Tf.mo=Tf.mo),delta=delta, Tf.mo=Tf.mo), Tf.mo=Tf.mo)


### GIBB'S SAMPLER
keep<-seq(A,B,ke)
K<-length(keep)

rho.mat<-eta.mat<-lambdas<-betas<-vector(length=K)
rhohat<- runif(1,0.5,1) 
etahat<- runif(1,1,3) 
lam<- runif(1,0.025,0.075) 
bet<- 0 

rho<-eta<-lambda<-beta<-NULL

Nhat<-Zhat<-vector(length=n)
N<-Z<-NULL

pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)

	for(i in 1:n){NiZi<-get.NiZi(deltai=delta[i], X=X[i], pk0=pk0, pk1=pk1, rhoi=rhohat, etai=etahat, lexbi=lam*exp(X[i]*bet), Yi=Y[i], N0=1)
		Nhat[i]<-NiZi$Nc
		Zhat[i]<-NiZi$Zc}
HD<-vector(length=K)

for(j in 2:B){
		
	rhohat<-slice.rho(rho0=rhohat, eta=etahat, Nhat=Nhat, Zhat=Zhat, rho.a=rho.a, rho.b=rho.b)

	etahat<-slice.eta(eta0=etahat, rho=rhohat, Nhat=Nhat, Zhat=Zhat, eta.s=eta.s, eta.r=eta.r)
	
	pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
	pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)
	
	for(i in 1:n){NiZi<-get.NiZi(deltai=delta[i], X=X[i], pk0=pk0, pk1=pk1, rhoi=rhohat, etai=etahat, lexbi=lam*exp(X[i]*bet), Yi=Y[i], N0=Nhat[i])
		Nhat[i]<-NiZi$Nc
		Zhat[i]<-NiZi$Zc}
		
	lam<-get.lambda(delta=delta, Zs=Zhat, bet=bet, Y=Y, lam.s=lam.s, lam.r=lam.r)
	bet<-slice.beta(beta0=bet, lam=lam, Zhat=Zhat, delta=delta, Y=Y)


		if(j %in%keep){#print(j)
		k.num<-c(1:K)[j==seq(A,B,ke)]
		rho.mat[k.num]<-rhohat
		eta.mat[k.num]<-etahat
		lambdas[k.num]<-lam
		betas[k.num]<-bet

		HD[k.num]<-get.yrep.surv(bet=bet, lam=lam, rhohat=rhohat, etahat=etahat,  pdf.true=pdf.true, Tf=Tf, Tf.mo=Tf.mo)
		}}

bet.res<-c(median(betas), cover(x=betas, true=true.bet))
lam.res<-c(median(lambdas), cover(x=lambdas, true=true.lam))
rho.res<-c(median(rho.mat), cover(x=rho.mat, true=true.rho))
eta.res<-c(median(eta.mat), cover(x=eta.mat, true=true.eta))

HD<-quantile(HD,p=c(0.5, 0.025, 0.975))


results<-list(bet.res=bet.res, lam.res=lam.res, rho.res=rho.res, eta.res=eta.res, HD=HD, seed=seed) 

print(c(seed,"completed"))
#write.csv(c(read.csv(paste("status-",PNR,"-",ER,".csv",sep="")),seed), paste("status-",PNR,"-",ER,".csv",sep=""))

return(results)}



