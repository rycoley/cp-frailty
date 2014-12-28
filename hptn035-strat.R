##### rycoley@gmail.com
#####This is 2 of 2 code files for re-creating the application analysis in Coley and Brown (2014) ``Estimating effectiveness in HIV prevention trials with a hierarchical compound Poisson frailty model"
##### part of dissertation work at University of Washington


#site-stratified analysis, proportional shift in hazard before 6 months
#prior for baseline hazard (lambda[j]) dependent on expected number of exposure processes (rho[j]) for each site j. 
#assume similar baseline hazard for those at-risk across sites


#unfortunately, this trial data is not yet publicly available. This code should show how to do a similar analysis with interval-censored event times
data<-read.csv("data035-cp.csv")
X<-data$tx #binary in this application

n0<-sum(X==0) #number in each treatment arm
n1<-sum(X==1)
n<-n0+n1

sn<-data$sn #factor variable for each site
G<-length(unique(sn)) #number of sites

n.site<-vector(length=G) #number of participants within each site
for(g in 1:G){n.site[g]<-sum(data$sn==g)}

V.site<-matrix(nrow=n, ncol=G) #design matrix for sites
for(g in 1:G){V.site[,g]<-sn==g}

##we want to compare observed and predcited survival up to 3 years
##for quarters
Tf<-12 # 12 quarters=3 years
quarters<-seq(1/4,Tf/4,1/4)

##in months
Tf.mo<-36 #
months<-seq(1/12,Tf.mo/12,1/12) #36
months<-round(months,4)


source("hptn035-strat-source.R") #accompanying code with functions for log posteriors and sampling 

###
#Y<-data$time
delta<-rep(0,n) 
delta[data$pos==1]<-1 #indicator of event

neg.time<-data$neg.time/365 #final follow-up time for those without event observed

#event times for those with event observed
pos.min.time<-data$pos.min.time/365  #last negative test
pos.max.time<-data$pos.max.time/365 #first positive test
pos.int.time<-data$pos.int.time/365


##here, I get kaplan meier survival for observed data. will compare to data simulated under the posterior
library(survival)
Y.mid<-rep(NA,n)
Y.mid[delta==0]<-neg.time[delta==0]
Y.mid[delta==1]<-pos.int.time[delta==1]

#quarters
Y.q<-get.yq(Y=Y.mid)
P.obs<-get.surv.q(Y.q,delta)
P.tx.obs<-get.surv.tx.q(Y.q,delta)
pdf.obs<-get.surv.pdf(cdf=P.obs)
pdf.placebo<-get.surv.pdf(cdf=P.tx.obs[,1])
pdf.tx<-get.surv.pdf(cdf=P.tx.obs[,2])
pdf.true<-as.matrix(cbind(pdf.obs, pdf.placebo, pdf.tx))


Y<-Y.mid<-Y.q<-Y.mo<-NULL
pdf.obs<-pdf.placebo<-pdf.tx<-NULL


##function that performs Gibbs sampling algorithm and compares predicted to observed survival
#inputs are B- total number of iterations, A- number burn-in, ke- number to thin, priors, starting seed for replication

do.one<-function(A, B, ke, theta.s, theta.r, tau.a, tau.b, rho.a, rho.b, eta.s, eta.r, lam0.s, lam0.r, seed){ 
		
set.seed(seed)

#set-up vectors/matrices for storing sampled values
keep<-seq(A,B,ke)
K<-length(keep)

rho.mat<-lambdas<-eta.mat<-matrix(nrow=K, ncol=G)
betas<-thetas<-taus<-lambda0s<-vector(length=K)

N.mat<-Z.mat<-matrix(nrow=K, ncol=n) #can save time by not saving these 
Nhat<-Zhat<-vector(length=n)

surv.all<-surv.placebo<-surv.pro<-matrix(nrow=K, ncol=Tf.mo)
HD<-matrix(nrow=K, ncol=3)


#get initial values
rhohat<- runif(G,0.5,1) #<-rho.start
etahat<- runif(G,1,3) #<-nu.start
lam<- runif(G,0.05,0.15) #<-lambda.start
bet<- 0 #<-beta.start
theta<-runif(1,0.15, 0.25)
tau<-runif(1,8,16)
lam0<-1

#n-length vector giving the rho/eta/lambda for each subj based on site
rho.vec<-as.vector(V.site%*%rhohat)
eta.vec<-as.vector(V.site%*%etahat)
lam.vec<-as.vector(V.site%*%lam)
		
Y.mat<-matrix(nrow=K, ncol=n)
Y<-rep(NA,n)
Y[delta==0]<-neg.time[delta==0]
Y[delta==1]<-pos.int.time[delta==1] #use midpoint for initial event times
Y.mat[1,]<-round(Y,4)

ch.vec<-get.cuml.haz(Y=Y, lam.vec=lam.vec, lam0=lam0)

for(i in 1:n){
	if(delta[i]==0){ #only those without event observed
		Nc<-rpois(1, woah(rhoi=rho.vec[i], etai=eta.vec[i], chexbi=ch.vec[i]*exp(X[i]*bet)))
		if(Nc==0){
			Nhat[i]<-0
			Zhat[i]<-0	}}
	if(delta[i]==1){
		N0<-rpois(1,1) 
		while(N0==0){N0<-rpois(1,1)} #only time truncated poisson is used
		mh.res<- MH.Ni(N0=N0, rhoi=rho.vec[i], etai=eta.vec[i], chexbi=ch.vec[i]*exp(X[i]*bet))		
		Nc<- mh.res[1]
		}
	if(Nc>0){
		Zc<- rgamma(1, shape=(Nc*eta.vec[i] + delta[i]), rate=((rho.vec[i]*eta.vec[i]) + ch.vec[i]*exp(X[i]*bet)))		 
		Nhat[i]<-Nc
		Zhat[i]<-Zc}	}
		N.mat[1,]<-Nhat
		Z.mat[1,]<-round(Zhat,4)

##GIBB'S SAMPLER
for(j in 2:B){

	if(j %in% keep){k<-c(1:K)[j==seq(A,B,ke)]}
	
	#sample event times from observed interval
	for(i in 1:n){if(delta[i]==1){
	Y[i]<-slice.Yi(Y0=Y[i], min.y=pos.min.time[i], max.y=pos.max.time[i], Zi=Zhat[i], lami=lam.vec[i], lam0=lam0, exbi=exp(X[i]*bet))}}	
	if(j %in% keep){Y.mat[k,]<-round(Y,4)}

	ch.vec<-get.cuml.haz(Y=Y, lam.vec=lam.vec, lam0=lam0)
	
	#sample mean number of exposure processes for each site
	rhohat<-slice.rho(rho0=rhohat, eta=etahat, Nhat=Nhat, Zhat=Zhat, rho.a=rho.a, rho.b=rho.b, theta=theta, tau=tau, lam=lam)
	if(j %in% keep){rho.mat[k,]<-rhohat}
	rho.vec<-as.vector(V.site%*%rhohat)

	#sample shape parameter for gamma-distributed risk associated with each exposure process
	etahat<-slice.eta(eta0=etahat, rho=rhohat, Nhat=Nhat, Zhat=Zhat, eta.s=eta.s, eta.r=eta.r)
	if(j %in% keep){eta.mat[k,]<-etahat}
	eta.vec<-as.vector(V.site%*%etahat)

	#sample number of exposure processes and total frailty for each individual, as described in Methods section
for(i in 1:n){
	#Nc<-Zc<-NULL
	if(delta[i]==0){ #only those without event observed
		Nc<-rpois(1, woah(rhoi=rho.vec[i], etai=eta.vec[i], chexbi=ch.vec[i]*exp(X[i]*bet)))	
		if(Nc==0){
			Nhat[i]<-0
			Zhat[i]<-0	}}
	if(delta[i]==1){
		#N0<-mh.res<-NULL
		mh.res<- MH.Ni(N0=Nhat[i], rhoi=rho.vec[i], etai=eta.vec[i], chexbi=ch.vec[i]*exp(X[i]*bet))		
		Nc<- mh.res[1]
#		AP.mat[j,i]<- mh.res[2]
		}
	if(Nc>0){
		Zc<- rgamma(1, shape=(Nc*eta.vec[i] + delta[i]), rate=((rho.vec[i]*eta.vec[i]) + ch.vec[i]*exp(X[i]*bet)))		 
		Nhat[i]<-Nc
		Zhat[i]<-Zc}	}

	if(j%in%keep){	
		N.mat[k,]<-Nhat
		Z.mat[k,]<-round(Zhat,4)}	

	#sample expected baseline hazard among exposured individuals
	theta<-slice.theta(theta=theta, rhoh=rhohat, tau=tau, lam=lam, theta.s=theta.s, theta.r=theta.r)
	if(j %in% keep){thetas[k]<-theta}
	
	#sample rate parameter for prior on site-specific baseline hazards
	tau<-slice.tau(tau0=tau, theta=theta, rhoh=rhohat, lam=lam, tau.a=tau.a, tau.b=tau.b)
	if(j %in% keep){taus[k]<-tau}

	#sample site-specific baseline hazard
	for(g in 1:G){lam[g]<-get.lambda(deltas=delta[sn==g], Xs=X[sn==g], Zs=Zhat[sn==g], bet=bet, Ys=Y[sn==g], lam0=lam0, theta=theta, rho=rhohat[g], tau=tau)}
	if(j %in% keep){lambdas[k,]<-lam}
	lam.vec<-as.vector(V.site%*%lam)
	
	#sample shared proporitonal adjustment to constant baseline hazard for first 6 months 
	lam0<-get.lambda0(deltas=delta, Xs=X, Zs=Zhat, bet=bet, Ys=Y, lam.vec=lam.vec, lam0.s=lam0.s, lam0.r=lam0.r)
	if(j %in% keep){lambda0s[k]<-lam0}
	
	ch.vec<-get.cuml.haz(Y=Y, lam.vec=lam.vec, lam0=lam0)
	
	#sample intervention effectiveness
	bet<-slice.beta(beta0=bet, ch.vec=ch.vec, Zhat=Zhat, delta=delta)
	if(j %in% keep){betas[k]<-bet}
	
	if(j %in% keep){
		#simulate date under posterior using observed censoring patterns
		prob.cens<-get.prob.cens(delta=delta,Y=Y)		
		prob.cens[Tf,]<-prob.cens[Tf,]+prob.cens[(Tf+1),]
		prob.cens[(Tf+1),]<-rep(0,G)

		surv.rep<-get.yrep.surv(bet=bet, lam=lam, lam0=lam0, rhohat=rhohat, etahat=etahat, prob.cens=prob.cens,pdf.true=pdf.true)

		surv.all[k,]<-surv.rep$surv.all
		surv.placebo[k,]<-surv.rep$surv.tx[,1]
		surv.pro[k,]<-surv.rep$surv.tx[,2]
	
		HD[k,]<-surv.rep$HD

		print(j)}
}

res<-list(bet=betas, lam=lambdas, lam0=lambda0s, rho=rho.mat, eta=eta.mat, theta=thetas, tau=taus, N.mat=N.mat, Z.mat=Z.mat, Y.mat=Y.mat, surv.all=surv.all, surv.placebo=surv.placebo, surv.pro=surv.pro, HD=HD)
return(res)}


#example for implementing and saving results

res<-do.one(A=100000, B=300000, ke=10, theta.s=10, theta.r=28, tau.a=2, tau.b=20, rho.a=c(7.039, 9.107, 7.279, 5.080, 7.170, 8.645, 7.360), rho.b=rep(2,G), eta.s=rep(2,G), eta.r=rep(0.5,G), lam0.s=3, lam0.r=4, seed=SEED)

write.csv(list(bet=res$bet, lam=res$lam, lam0=res$lam0, eta=res$eta, rho=res$rho, tau=res$tau, theta=res$theta), paste("res-035-strat-parameters-",SEED,".csv",sep=""))

write.csv(list(surv.all=res$surv.all, surv.placebo=res$surv.placebo, surv.pro=res$surv.pro), paste("res-035-strat-survival-",SEED,".csv",sep=""))

write.csv(res$HD,paste("res-035-strat-HD-",SEED,".csv",sep=""))

write.csv(res$N.mat,"res-035-strats-Ns-1.csv")
write.csv(res$Z.mat,"res-035-strats-Zs-1.csv")
write.csv(res$Y.mat,"res-035-strats-Ys-1.csv")

res<-NULL
