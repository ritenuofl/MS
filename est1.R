################## (updated) ######################
########## load all the needed packages #############
library(mstate)
library(seqinr)
### incorporating (fixed) covariates into markov chain
## and using mstate and lm to estimate parameters

#RM: this file is a storehouse of functions that we need
#later
#RM Just now: 2 functions one for generation, 
#one for msprep
set.seed(123)

## we have N observations,pi0 would be N*k matrix; Q would be a list of transition
## rate matrix (k*k) for each pt depending on their Z (e.g., Z~N(0,1))
#k<-3
#N<-10
#Z<-rbinom(N,1,0.5)
#print

gen<-function(N=200,k=3,pr=.5){
  beta=matrix(.5,k,k); diag(beta)=0
  #RM pr- chance of binary
#cov_ctmc<-function(N,k,pi0,q0,Z) {
  Z<-rbinom(N,1,pr)
pi0=matrix(rep(c(0.4,0.5,0.1),N),byrow=TRUE,ncol=k)
q0<-matrix(runif(k*k,0,10),ncol=k)  ### fixed param;baseline Q matrix

sum_q<-c()
for (i in 1:k){
  sum_q[i]<-sum(q0[i,][(row(q0)!=col(q0))[i,]])
}
diag(q0)<--sum_q
#beta=0.5 ##RM: already in function argument
for (i in 1:N){ ## i=1
   #RM: not wrong at all but maybe  b_fun may not 
  # be put in a separate func. Since you do not need to
  # use it for  many different occassions repeatedly.

    #loop
  #bfun<-function(h,l,Z) {
  Q =q0*exp(beta*Z[i])
   #RM:Just returns a matrix at one go instead of h,l piecewise
  #RM: Q will vary for each i, but u may not have to store every Q
  diag(Q)=0
  diag(Q)=-apply(Q,1,sum)
  #RM: just do negative row sum
    
#RM1: again,overall good. Since we are not looking at all these different
#Q matrices we do not maybe even need to store them 
# and keep a  list of them.  Yes all Qs are different and this way 
#coding captures this variation among Qs, but for generation
#as well as estimation, we would be  looking only at betas,we can get 
#individual Qs whenever we want by the mat= step
########### finish simulate Q matrices ##############  
################### start simulate Y and T ###############
YT<-list()
## t and y are vectors
#for (i in 1:N) {## i=1;i=2
  #RM: Just using one for loop, that is,in this setup
  #generating q s and data at one go.
  t0<-0
  tmax<-1
  t<-y<-c()
  y0<-sample(k,1,prob=pi0[i,])              ## initial state
  y[1]<-y0
  #RM: The Gillespie chunk below is perfect
  while (t0<tmax) {
    tstar<-rexp(1,-Q[y0,y0])
    t0<-t0+tstar
    weights=Q[y0,]
    weights[y0]=0
    y0<-sample(k,1,prob=weights/sum(weights))
    y<-c(y,y0)
    t<-c(t,t0)
    ## str(t);str(y)
    yt<-cbind(y=y[-length(y)],t)
  }
  YT[[i]]<-yt
}
return(YT)
}
#check generation code works
YT1= gen()
#str(YT1)    
#plot(YT1[[1]][,"t"],YT1[[1]][,"y"]) 
#plot(YT1[[2]][,"t"],YT1[[2]][,"y"]) 

