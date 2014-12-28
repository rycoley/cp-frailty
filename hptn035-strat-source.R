
##### rycoley@gmail.com
#####This is 1 of 2 code files for re-creating the application analysis in Coley and Brown (2014) ``Estimating effectiveness in HIV prevention trials with a hierarchical compound Poisson frailty model"
##### part of dissertation work at University of Washington


#This file contains the posterior log likelihoods and sampling functions for all model parameters.  

## This is a site-stratified analysis that assumes a piecewise constant baseline hazard (pre- and post- 6 months) for each site and compound Poisson frailty (one frailty distribution/site)

#log-posterior likelihood for parameter in proportional hazards function (beta)
#in application, hazard ratio associated with intervention
#prior(beta)~N(0,10) 
lik.beta<-function(bet, ch.vec, Zs, delta){#
	XB<-X[Zs>0]*bet
	return(sum(delta[Zs>0]*XB) - sum(Zs[Zs>0]*exp(XB)*ch.vec[Zs>0]) + log(dnorm(bet,0,10)) )}

#slice sampler for beta
slice.beta<-function(beta0, ch.vec, Zhat, delta){#delta, rhohat, nuhat, nor, 
	z<-J<-K<-L<-R<-beta.star<-NULL
	z<- lik.beta(bet=beta0, ch.vec=ch.vec, Zs=Zhat, delta=delta) - rexp(1,1) #want to make a slice with all betas that have a log-likelihood above this
	w<- 0.05 #width of each step
	m<- 10 #number of total steps out
	L<- beta0 - (w*runif(1,0,1)) #starting lower bound of slice
	R<- L + w #starting upper bound of slice
	J<- floor(m*runif(1,0,1)) #keep track of number of steps
	K<- (m-1)-J
	while(lik.beta(bet=L,  ch.vec=ch.vec, Zs=Zhat, delta=delta) > z & J>0){L<- L-w
		J<- J-1} #extend slice to contain all betas below starting beta with log likelihood above z
	while(lik.beta(bet=R,  ch.vec=ch.vec, Zs=Zhat, delta=delta) > z & K>0){R<- R+w
		K<- K-1}
	beta.star<- runif(1,L,R) #sample new beta from this slice 
	while(lik.beta(bet=beta.star, ch.vec=ch.vec, Zs=Zhat, delta=delta) < z){ 
		if(beta.star<beta0){L<-beta.star}
		if(beta.star>beta0){R<-beta.star}
		beta.star<- runif(1,L,R)}#repat sampling to obtain sample of beta that has log-likelihood above z
	return(beta.star)}

##see paper for definition of baseline hazard
#sample constant baseline hazard for time>6 months for each site, lambda[j]
#posterior is gamma when prior(lambda)~gamma
get.lambda<-function(deltas, Xs, Zs, bet, Ys, lam0, theta, rho, tau){
	shape<- sum(deltas) + (tau*theta*(1-exp(-rho)))
	rate<- sum(Zs[Ys<=0.5] * exp(as.vector(Xs[Ys<=0.5]*bet)) * lam0 * Ys[Ys<=0.5]) + sum(Zs[Ys>0.5] * exp(as.vector(Xs[Ys>0.5]*bet)) * (0.5*lam0 + (Ys[Ys>0.5] - 0.5))) + tau
	return(rgamma(1,shape=shape, rate=rate))}

#sample shared proportional adjustment for baseline hazard before 6 months 
#here, lambda0. in paper, lambda_{6mo}
get.lambda0<-function(deltas,Xs,Zs,bet,Ys,lam.vec,lam0.s,lam0.r){
	shape<- sum(deltas[Ys<=0.5]) + lam0.s
	rate<- sum(Zs[Ys<=0.5] * exp(as.vector(Xs[Ys<=0.5]*bet)) * lam.vec[Ys<=0.5] * Ys[Ys<=0.5]) + sum(Zs[Ys>0.5] * exp(as.vector(Xs[Ys>0.5]*bet)) * lam.vec[Ys>0.5] * 0.5) + lam0.r
	return(rgamma(1, shape=shape, rate=rate))}


#log prior density for number of exposure processes for individuals in site j, rho[j]
#prior(1-exp(-rho[j]))) ~ beta(A,B)
#here, rho.a=A, rho.b=B  
ldrho<-function(rhoh, rho.a,rho.b){ #log likelihood of prior #this is for proprotion not at risk
	return(-rho.a*rhoh + (rho.b-1)*log(1-exp(-rhoh)))}

#log-posterior likelihood for each rho[j]
lik.rho<-function(rhoh, g, eta, Ns, Zs, rho.a, rho.b, theta, tau, lam.g){
	nu<-rhoh*eta
	A<-tau*theta*(1-exp(-rhoh))
	return(-(n.site[g]*rhoh)+ sum(Ns[Ns>0]*log(rhoh)) + sum(Ns[Ns>0]*eta*log(nu))   -sum(Zs[Zs>0]*nu) + ldrho(rhoh=rhoh, rho.a, rho.b) + A*log(tau) - lgamma(A) + A*log(lam.g) )}

#slice sampler for each rho[j]
#comparable to slice.beta. see notes there for more details on function
slice.rho<-function(rho0, eta, Nhat, Zhat, rho.a, rho.b, theta, tau, lam){#delta, nuhat, lexbT, 
	for(g in 1:G){
	z<-J<-K<-L<-R<-rho.star<-NULL
	rho.star<-rho0
	z<- lik.rho(rhoh=rho0[g], g=g, eta=eta[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], rho.a=rho.a[g], rho.b=rho.b[g], theta=theta, tau=tau, lam.g=lam[g]) -rexp(1,1) 	
	w<- 0.1
	m<- 10
	L<- rho0[g] - (w*runif(1,0,1))
	R<- L + w
	L<- max(L,1e-10)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-10 & lik.rho(rhoh=L, g=g, eta=eta[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], rho.a=rho.a[g], rho.b=rho.b[g], theta=theta, tau=tau, lam.g=lam[g]) > z & J>0){
		L<- max(L-w,1e-10)
		J<- J-1}
		L<-max(L,1e-10)
	while(lik.rho(rhoh=R, g=g, eta=eta[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], rho.a=rho.a[g], rho.b=rho.b[g], theta=theta, tau=tau, lam.g=lam[g]) > z & R<10 & K>0){
		R<- R+w
		K<- K-1}
		R<-min(R,10)
	rho.star[g]<- runif(1,L,R)
	while(lik.rho(rhoh=rho.star[g], g=g, eta=eta[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], rho.a=rho.a[g], rho.b=rho.b[g], theta=theta, tau=tau, lam.g=lam[g]) < z){
		if(rho.star[g] < rho0[g]){L<- rho.star[g]}
		if(rho.star[g] > rho0[g]){R<- rho.star[g]}
		rho.star[g]<- runif(1,L,R)}
		rho0<-rho.star}
	return(rho.star)	}

#log-posterior likelihood for eta[j], the site-specific shape parameter of the gamma random variable for amount of risk associated with each exposure process 
#prior(eta[j])~Gamma(O_{eta[j]}, T_{eta[j]}) = Gamma(eta.s, eta.r)
lik.eta<-function(etah, rho, Ns, Zs, eta.s, eta.r){
	nu<-rho*etah
	return(sum(Ns[Ns>0]*etah*log(nu)) - sum(lgamma(Ns[Ns>0]*etah)) + sum(Ns[Zs>0]*etah*log(Zs[Zs>0])) -sum(Zs[Zs>0]*nu) + log(dgamma(etah, shape=eta.s, rate=eta.r)) )}

#slice sampler for each eta[j]
#comparable to slice.beta. see notes there for more details on function
slice.eta<-function(eta0, rho, Nhat, Zhat, eta.s, eta.r){#rhohat, delta,  lexbT, 
	for(g in 1:G){
	z<-J<-K<-L<-R<-eta.star<-NULL
	eta.star<-eta0
	z<- lik.eta(etah=eta0[g], rho=rho[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], eta.s=eta.s[g], eta.r=eta.r[g]) -rexp(1,1) 	
	w<- 0.1
	m<- 10
	L<- eta0[g] - (w*runif(1,0,1))
	R<- L + w
	L<-max(L,1e-20)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-20 & lik.eta(etah=L, rho=rho[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], eta.s=eta.s[g], eta.r=eta.r[g]) > z  & J>0){L<- max(L-w,1e-20)
		J<- J-1}
	while(lik.eta(etah=R, rho=rho[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], eta.s=eta.s[g], eta.r=eta.r[g]) > z & K>0){R<- R+w
		K<- K-1}
	eta.star[g]<- runif(1,L,R)
	while(lik.eta(etah=eta.star[g], rho=rho[g], Ns=Nhat[sn==g], Zs=Zhat[sn==g], eta.s=eta.s[g], eta.r=eta.r[g]) < z){
		if(eta.star[g] < eta0[g]){L<- eta.star[g]}
		if(eta.star[g] > eta0[g]){R<- eta.star[g]}
		eta.star[g]<- runif(1,L,R)}
		eta0<-eta.star}
	return(eta.star)	}

#log-posterior likelihood for the baseline hazard among those at-risk, as explained in Application section of paper
#prior(theta)~Gamma(O_{theta}, T_{theta}) = Gamma(theta.s, theta.r)
lik.theta<-function(theta, rhoh, tau, lam, theta.s, theta.r){
	A<-tau*theta*(1-exp(-rhoh))
	return( sum(A*log(tau)) - sum(lgamma(A)) + sum(A*log(lam)) + log(dgamma(theta, shape=theta.s, rate=theta.r)))}

#slice sampling function for theta
#comparable to slice.beta. see notes there for more details on function
slice.theta<-function(theta, rhoh, tau, lam, theta.s, theta.r){#rhohat, delta,  lexbT, 
	z<-J<-K<-L<-R<-theta.star<-NULL
	z<- lik.theta(theta=theta, rhoh=rhoh, tau=tau, lam=lam, theta.s=theta.s, theta.r=theta.r) -rexp(1,1) 
	w<- 0.01
	m<- 5
	L<- theta - (w*runif(1,0,1))
	R<- L + w
	L<- max(L,1e-20)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-20 & lik.theta(theta=L, rhoh=rhoh, tau=tau, lam=lam, theta.s=theta.s, theta.r=theta.r) > z & J>0){L<-max(L-w,1e-20)
		J<- J-1}
		L<-max(L,1e-20)
	while(lik.theta(theta=R, rhoh=rhoh,tau=tau, lam=lam, theta.s=theta.s, theta.r=theta.r) > z & K>0){R<- R+w
		K<- K-1}
	theta.star<- runif(1,L,R)
	while(lik.theta(theta=theta.star, rhoh=rhoh, tau=tau, lam=lam, theta.s=theta.s, theta.r=theta.r) < z){
		if(theta.star < theta){L<- theta.star}
		if(theta.star > theta){R<- theta.star}
		theta.star<- runif(1,L,R)}
	return(theta.star)	}	

#log-posterior likelihood for the rate parameter for the gamma prior on all lambda[j], as explained in Application section of paper
#prior(tau)~Inverse Gamma(tau.a, tau.b)
lik.tau<-function(tau, theta, rhoh, lam, tau.a, tau.b){
	A<-tau*theta*(1-exp(-rhoh))
	return(sum(A*log(tau)) -sum(lgamma(A)) + sum(A*log(lam)) - sum(tau*lam) + (-tau.a-1)*log(tau) - (tau.b/tau))}

#slice sampler for tau
#comparable to slice.beta. see notes there for more details on function
slice.tau<-function(tau0, theta, rhoh, lam, tau.a, tau.b){#rhohat, delta,  lexbT, 
	z<-J<-K<-L<-R<-tau.star<-NULL
	z<- lik.tau(tau=tau0, theta=theta, rhoh=rhoh, lam=lam, tau.a=tau.a, tau.b=tau.b) -rexp(1,1) 
	w<- 0.2
	m<- 10
	L<- tau0 - (w*runif(1,0,1))
	R<- L + w
	L<- max(L,1e-20)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-20 & lik.tau(tau=L, theta=theta, rhoh=rhoh, lam=lam, tau.a=tau.a, tau.b=tau.b) > z & J>0){L<-max(L-w,1e-20)
		J<- J-1}
		L<-max(L,1e-20)
	while(lik.tau(tau=R, theta=theta, rhoh=rhoh, lam=lam, tau.a=tau.a, tau.b=tau.b) > z & K>0){R<- R+w
		K<- K-1}
	tau.star<- runif(1,L,R)
	while(lik.tau(tau=tau.star, theta=theta, rhoh=rhoh, lam=lam, tau.a=tau.a, tau.b=tau.b) < z){
		if(tau.star < tau0){L<- tau.star}
		if(tau.star > tau0){R<- tau.star}
		tau.star<- runif(1,L,R)}
	return(tau.star)	}	

#log-transformed density for the number of exposure processes for participant i, marginalized over Z[i]
#follows eqn (8) in Methods section
lik.Ni<-function(Ni,rho,eta,chexbi){ #for deltai=1
	nu<-rho*eta
	return(Ni*log(rho) -lgamma(Ni+1) + log(Ni*eta) + Ni*eta*log(nu) - ((Ni*eta)+1)*log(nu + chexbi)  ) 	}

#proposal probabilities for metropolis-hastings algorithm
pN<-function(given){
	if (given==1){return(1/2)}
	if (given>1){return(1/3)}}

#metropolis-hastings algorithm to sample N[i], the number of exposure processes for participant i
MH.Ni<-function(N0, rhoi, etai, chexbi){
	U<-runif(1,0,1) #sample candidate, conditional on obesrved event
	if(N0>1){ 
		if(U < 1/3 ){Nc<- N0-1 }
		if(U > 1/3 & U < 2/3){Nc<- N0
			return(c(Nc,1))}
		if(U > 2/3){Nc<- N0+1}}
	if(N0==1){
		if(U < 1/2){Nc<- N0
			return(c(Nc,1))}
		if(U > 1/2){Nc<- N0+1}}
	
	accept<- rbinom(1,1, min( (exp(lik.Ni(Ni=Nc,rho=rhoi, eta=etai, chexbi=chexbi)-lik.Ni(Ni=N0, rho=rhoi, eta=etai, chexbi=chexbi) )) * (pN(Nc)/pN(N0)) ,1) )
	if(accept==1){return(c(Nc,1))}
	if(accept==0){return(c(N0,0)) }}

#mean parameter for Poisson dist on N[i] when no event observed, delta[i]=0
#follows eqn (8) in Methods section
woah<-function(rhoi, etai, chexbi){ #
	nu<-rhoi*etai
	return(rhoi*(nu/(nu+chexbi))^etai)} 

#log-likelihood for event time 
lik.Yi<-function(y, Zi, lami, lam0, exbi){
	if(y<=0.5){return(-Zi*exbi*lami*lam0*y)}
	else{return(-Zi*exbi*lami*(lam0*0.5 + (y-0.5)))}}

#slice sampler for interval-censored event time.
slice.Yi<-function(Y0, min.y, max.y, Zi, lami, lam0, exbi){
	z<-J<-K<-L<-R<-y.star<-NULL
	z<- lik.Yi(y=Y0, Zi=Zi, lami=lami, lam0=lam0, exbi=exbi) - rexp(1,1)
	w<- 0.05
	m<- 10
	L<- Y0 - (w*runif(1,0,1))
	R<- L + w
	L<-max(L,min.y)
	R<-min(R,max.y)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>min.y & lik.Yi(y=L, Zi=Zi, lami=lami, lam0=lam0, exbi=exbi) > z & J>0){L<- L-w
		J<- J-1}
		L<-max(L, min.y)
	while(R<max.y & lik.Yi(y=R, Zi=Zi, lami=lami, lam0=lam0, exbi=exbi) > z & K>0){R<- R+w
		K<- K-1}
		R<-min(R, max.y)
	y.star<- runif(1,L,R)
	while(lik.Yi(y=y.star, Zi=Zi, lami=lami, lam0=lam0, exbi=exbi) < z){
		if(y.star<Y0){L<-y.star}
		if(y.star>Y0){R<-y.star}
		y.star<- runif(1,L,R)}
	return(y.star)}

#function to obtain cumulative hazard, Lambda(t)= intergral_0^t lambda(t)* Z[i]* exp(X[i]*beta)
get.cuml.haz<-function(Y,lam.vec,lam0){
	ch.vec<-vector(length=n)
	ch.vec[Y<=0.5]<-lam.vec[Y<=0.5]*lam0*Y[Y<=0.5]
	ch.vec[Y>0.5]<-lam.vec[Y>0.5]*(lam0*0.5 + (Y[Y>0.5]-0.5))
	return(ch.vec)}



#the following are functions to simulate survival under the posterior based on current parameter samples and observed censoring patterns

##for analysis in paper, simulated event times are discretized to quarters. next, survival functions are discretized to months for calculating hellinger distance.

#round event time to months within a year
get.ymo<-function(Y){
	Y.mo<-rep((1/12),n)
	for(j in 2:Tf.mo){
		Y.mo[Y>months[(j-1)] & Y<=months[j]]<-months[j]}
	Y.mo[Y>months[Tf.mo]]<-months[Tf.mo]
	return(Y.mo)}

#kaplan meier survival estimate for event times (in months)
get.surv.mo<-function(Y.mo,delta){
	so<-Surv(time=Y.mo, event=delta)
	km<-survfit(formula=so~1)
	if(length(km$surv)==Tf.mo){return(km$surv)}
	else{km.surv<-vector(length=Tf.mo)
		for(j in 1:Tf.mo){if(sum(km$time==months[j])>0){km.surv[j]<-km$surv[km$time==months[j]]}
			else{if(j==1){km.surv[j]<-1}
				else{km.surv[j]<-km.surv[(j-1)]}}}
				return(km.surv)}}

#kaplan meier survival for event times (in months) within treatment arms		
get.surv.tx.mo<-function(Y.mo,delta){
	so<-Surv(time=Y.mo, event=delta)
	km<-survfit(formula=so~X)
	if(length(km$time)==(Tf.mo*2)){return(cbind(km$surv[1:Tf.mo], km$surv[(Tf.mo+1):(Tf.mo*2)]))} 
	else{ so<-km<-NULL
		surv.mat<-matrix(nrow=Tf.mo,ncol=2)
		surv.mat[,1]<-get.surv.mo(Y.mo=Y.mo[X==0], delta=delta[X==0])
		surv.mat[,2]<-get.surv.mo(Y.mo=Y.mo[X==1], delta=delta[X==1])
		return(surv.mat)} }

##code for quarters
#event times in quarters, for calculating Hellinger distance
get.yq<-function(Y){
	Y.q<-rep((1/4),n)
	for(j in 2:Tf){
		Y.q[Y>quarters[(j-1)] & Y<=quarters[j]]<-quarters[j]}
	Y.q[Y>quarters[Tf]]<-quarters[Tf]
	return(Y.q)}

#kaplan meier survival for event times (in quarters)
get.surv.q<-function(Y.q,delta){ 
	so<-Surv(time=Y.q, event=delta)
	km<-survfit(formula=so~1)
	if(length(km$surv)==Tf){return(km$surv)}
	else{km.surv<-vector(length=Tf)
		for(j in 1:Tf){if(sum(km$time==quarters[j])>0){km.surv[j]<-km$surv[km$time==quarters[j]]}
			else{if(j==1){km.surv[j]<-1}
				else{km.surv[j]<-km.surv[(j-1)]}}}
				return(km.surv)}}
#kaplan meier survival for event times (in quarters) within treatment arm
get.surv.tx.q<-function(Y.q,delta){
	so<-Surv(time=Y.q, event=delta)
	km<-survfit(formula=so~X)
	if(length(km$time)==(Tf*2)){return(cbind(km$surv[1:Tf], km$surv[(Tf+1):(Tf*2)]))} 
	else{ so<-km<-NULL
		surv.mat<-matrix(nrow=Tf,ncol=2)
		surv.mat[,1]<-get.surv.q(Y.q[X==0], delta[X==0])
		surv.mat[,2]<-get.surv.q(Y.q[X==1], delta[X==1])
		return(surv.mat)} }

#calculate density for discrete time from survival
get.surv.pdf<-function(cdf){
	pdf<-rep(0,(Tf+1))
	pdf[1]<-1-cdf[1]
	for(t in 2:Tf){
		pdf[t]<-cdf[(t-1)]-cdf[t]} 
	pdf[(Tf+1)]<-cdf[Tf]
	return(pdf)}

#time<-seq(0,(Tf/12),(1/12))
time<-seq(0,(Tf/4),(1/4))

#obtain observed censoring patterns	
get.time.cens<-function(Y, Nt){
	time.cens<-matrix(nrow=n, ncol=(Nt-1))
	for(s in 1:(Nt-1)){
		time.cens[,s]<-time[s]<= Y & Y<time[(s+1)] & delta==0}
		return(time.cens)	}

get.prob.cens<-function(delta, Y){
	Nt<-length(time)
	t.cens<-get.time.cens(Y=Y, Nt=Nt)

	prob.cens<-matrix(nrow=Nt, ncol=G)
	prob.cens[1:(Nt-1),]<-t(t.cens)%*%V.site #R.mat

	for(g in 1:G){
		prob.cens[Nt,g]<-sum(Y>=(Tf/4) & sn==g & delta==0)
		prob.cens[,g]<-prob.cens[,g]/sum(sn==g & delta==0)}
	return(prob.cens)}	

#calculate hellinger distance
get.HD<-function(p,q){ #hellinger distance
	hd<-rep(0,ncol(p))
	for(j in 1:ncol(p)){
		hd[j]<-(1/sqrt(2)) * sqrt( sum( (sqrt(p[,j])-sqrt(q[,j]))^2 ) )}
		return(hd)}


#simulate survival times given current parameter samples, censoring
#calculated hellinger distance between observerd survival and simulated  	
get.yrep.surv<-function(bet, lam, lam0, rhohat, etahat, prob.cens, pdf.true){
	
	set.seed(1)
		
	N.rep<-Z.rep<-rep(0,n)
	for(i in 1:n){
		N.rep[i]<-rpois(1,rhohat[sn[i]])
		if(N.rep[i]>0){Z.rep[i]<-rgamma(1,shape=N.rep[i]*etahat[sn[i]], rate=(etahat[sn[i]]*rhohat[sn[i]]))}}

		Y.rep<-rep(NA,n)
		for(i in 1:n){if(Z.rep[i]>1e-20){
			Y.rep[i]<-rexp(1, rate=lam[sn[i]]*lam0*exp(X[i]*bet)*Z.rep[i])
			if(Y.rep[i]>0.5){
				Y.rep[i]<-0.5+rexp(1,rate=lam[sn[i]]*Z.rep[i]*exp(X[i]*bet))}}
			else{Y[i]<-5000}}

	delta.rep<-rep(0,n)

	for(i in 1:n){
		T.cens<-NULL
		T.cens<-quarters[rmultinom(1,1,prob.cens[,sn[i]])==1]
		if(is.na(Y.rep[i])){Y.rep[i]<-T.cens}
		else{delta.rep[i]<-Y.rep[i]<T.cens
			Y.rep[i]<-min(Y.rep[i],T.cens)}}

	Y.q<-get.yq(Y=Y.rep)
#	Y.q<-round(Y.q)

	P.rep<-get.surv.q(Y.q,delta.rep)
	P.tx.rep<-get.surv.tx.q(Y.q,delta.rep)

	pdf.rep<-as.matrix(cbind(get.surv.pdf(cdf=P.rep), get.surv.pdf(cdf=P.tx.rep[,1]), get.surv.pdf(cdf=P.tx.rep[,2])))
	HD<-get.HD(p=pdf.true, q=pdf.rep)

	Y.mo<-get.ymo(Y=Y.rep)
	surv.all<-round(get.surv.mo(Y.mo=Y.mo, delta=delta.rep),5)
	surv.tx<-round(get.surv.tx.mo(Y.mo=Y.mo, delta=delta.rep),5)

return(list(HD=HD, surv.all=surv.all, surv.tx=surv.tx))}





