#Author: Joseph Feldman
#Date: January 14th, 2021
#Description: Code for the MCMC algorithm used for estimation of a Gaussian Copula using the extended rank-probit likelihood
## Specifically, the code produces posterior samples of Gaussian copula correlation matrices as well as intercepts, which are then 
## used in 2_Synthesizer_demo.R for construction of synthetic data sets through posterior predictive sampling.
#---------------------------------------------------------------------------------------------------------------------------------
#####Inputs:
# Y: (N x p) dataframe to be synthesized with binary and categorical variables encoded as factors. In this exercise, we 
##create a simulated data set with one categorical variable, one binary variable, one continuous variable, and a count variable 
##that is excluded for targeted synthesis. The resultant dataframe that is used to fit thecopula is Y_cop
# nsamp: number of iterations for the Markov chain.
# odens: output density- number of iterations between saved samples
# plugin.threshold: if the number of unique values of a variable exceeds
## this integer, then plug-in the empirical distribution as the marginal in MCMC.
#####Outputs:
# C.psamp: Posterior samples of gaussian copula correlation matrix
# alpha.psamp: Posterior samples of mean vector from gaussian copula

########## R libraries required for model

library(truncnorm)
library(MASS)
source('sampleMGP.R') # for sampling multiplicative gamma process parameters within mcmc

########## create simulated data set

# create categorical variable
cat = sample(c("red", "orange","yellow","green"),1000, replace = T, prob = c(.3,.1,.4,.2))

# create variable variable

cond.bin<-function(cat){
  if(cat == "red"){
    return(rbinom(1,1,.8))
  }
  else if(cat == "orange"){
    return(rbinom(1,1,.2))
  }
  else if(cat == "yellow"){
    return(rbinom(1,1,.5))
  }
  else{
    return(rbinom(1,1,.9))
  }
}
bin<- sapply(cat,cond.bin)

#continuous
cont <- ifelse(bin ==1, rnorm(sum(bin ==1),1,2),rnorm(sum(bin == 0),0,.5))

#create Y_cop dataframe
Y_cop<- data.frame(cat = as.factor(cat), bin = as.factor(bin), cont = cont)

#create count variable as a (complex) function of the previous variables 
Y_onehot<- model.matrix(~  cat + bin + cont + bin*cont, data = Y_cop)

lambdas = round(apply(Y_onehot,1,function(x) sum(x * c(310,4,3,2,3,.5,2))) + cont^2)

count  = rpois(1000, lambdas)

Y<- data.frame(Y_cop, count)

########## Data Wrangling

#### Only proceed to this chunk if the factor columns aren't organized in descending order by number of levels
#### with numeric columns coming after
#Pull out columns that are factors
factor.cols<- which(sapply(Y_cop,is.factor))

#Arrange Factors in order of number of Factors
if(length(factor.cols)>1){
  inc_f = Y_cop[,names(sort(sapply(Y_cop[,sapply(Y_cop, is.factor)], nlevels), decreasing = TRUE))]
  Y_cop<- data.frame(inc_f, Y_cop[,-factor.cols])
  factor.cols <- which(sapply(Y_cop,is.factor) == TRUE)
}
####Skip if ordering complete

#Create One hot encoded data frame for rank-probit likelihood
cat<- NULL
cat_col = 1 
n_cat = n_bin = 0
is_bin = is_cat = NULL 


for(i in factor.cols){
  if(nlevels(Y_cop[,i]) >2){
    one.hot.cat<- model.matrix(~Y_cop[,i] - 1)
    colnames(one.hot.cat)<- levels(Y_cop[,i])
    cat<- cbind(cat,one.hot.cat)
    cat_col = cat_col + nlevels(Y_cop[,i])
    n_cat = n_cat +nlevels((Y_cop[,i]))
    is_cat <- c(is_cat, rep(1,nlevels(Y_cop[,i])))
    is_bin<- c(is_bin, rep(0,nlevels(Y_cop[,i])))
    
  }
  else{
    one.hot.bin<- model.matrix(~Y_cop[,i])[,2]
    cat<-cbind(cat,one.hot.bin)
    colnames(cat)[cat_col] = colnames(Y_cop)[i]
    cat_col = cat_col + 1
    is_cat <- c(is_cat, rep(0,1))
    is_bin<- c(is_bin, rep(1,1))
  }
  
}


one.hot.cols<- ncol(cat) # number of one hot encoded categories
Y_mod<- cbind(cat,Y_cop[,-factor.cols]) # data frame used for sampler
n<-dim(Y_mod)[1] #N
p<-dim(Y_mod)[2] #p
is_cat_bin<- c(rep(1,one.hot.cols), rep(0,p - one.hot.cols)) #For extended rank-probit likelihood, determines re-sampling step for column j
vnames<-colnames(Y_mod)
Y_mod<-as.matrix(Y_mod)
colnames(Y_mod)<-vnames

########## Initialize Values 
# Initialize Z
set.seed(1)

R<-NULL
for(j in 1:p) { R<-cbind(R, match(Y_mod[,j],sort(unique(Y_mod[,j])))) }
Rlevels<-apply(R,2,max,na.rm=TRUE)
Ranks<- apply(Y_mod,2,rank,ties.method="max",na.last="keep")
N<-apply(!is.na(Ranks),2,sum)
U<- t( t(Ranks)/(N+1))
Z<-qnorm(U)
Zfill<-matrix(rnorm(n*p),n,p)
Z[is.na(Y_mod)]<-Zfill[is.na(Y_mod)]


# Randomly initialize categorical latent columns

Z[,1:one.hot.cols] = rnorm(n*one.hot.cols)

 

#Initialize mean vector
alpha =c(qnorm(colMeans(Y_mod[,1:one.hot.cols])), rep(0,p - one.hot.cols))



# initialize covariance matrix

k.star =min(which(summary(prcomp(Z))$importance[3,] > 0.95)) #min(p, ceiling(log(p))) # 

# Initialize Factor Model Parameters: these are 

# Local df:
nu = 3; 

# Global shrinkage:
a1 = 2; a2 = 3;

# Error variance:
a.sigma = 1; b.sigma = 0.3


# Use SVD for simplicity
svd0 = svd(Z); # plot(1-cumsum(svd0$d^2)/sum(svd0$d^2)); #plot(svd0$d[-1]/svd0$d[-min(p,n)])

# Factor loadings (p x k.star):
Lambda = Lambda.s = as.matrix(svd0$v[,1:k.star])

# Factors (n x k.star):
eta = svd0$u[,1:k.star]%*%diag(svd0$d[1:k.star], k.star) #y%*%Lambda

# Residuals (n x p):
eps = Z - (tcrossprod(eta, Lambda))

# Diagonal error variance (p-dimensional vector):
Sigma.diag = apply(eps, 2, var) 


# Local and global precision parameters:
phi.jh = 1/Lambda^2; tau.h = delta.h = rep(1, k.star)

# Use of Empirical CDF Threshold
plugin.threshold = 100
plugin.marginal = (apply(Y_mod, 2, function(x) {
  length(unique(x))
}) > plugin.threshold)

##########Storage
nsamp = 10000
odens = 1
C.psamp<- array(0,c(nsamp/odens,p,p))
alpha.psamp<- array(0,c(nsamp/odens,p))


########## Start MCMC

for(ns in 1:nsamp){
  
  
  cp.eta = crossprod(eta)
  for(j in 1:p){
    chQj = chol(diag(phi.jh[j,]*tau.h, k.star) + cp.eta/Sigma.diag[j])
    lj = crossprod(eta, Z[,j] - alpha[j])/Sigma.diag[j]
    Lambda[j,] = backsolve(chQj,forwardsolve(t(chQj), lj) + rnorm(k.star))
  }
  
  
  # Step 2: sample the error variances 
  eps = Z - sweep(tcrossprod(eta, Lambda),2,alpha,'+')
  
  Sigma.diag = apply(eps, 2, function(x) 1/rgamma(n = 1, shape = a.sigma + n/2, 
                                                  rate = b.sigma + 1/2*sum(x^2))) 
  
  # Step 3: sample the factors
  chQeta = chol(diag(k.star) + crossprod(Lambda, diag(1/Sigma.diag))%*%Lambda)
  leta = tcrossprod(crossprod(Lambda, diag(1/Sigma.diag)), sweep(Z,2,alpha,'-'))
  eta = t(backsolve(chQeta,forwardsolve(t(chQeta), leta) + rnorm(n*k.star))) #for(i in 1:n) eta[i,]= backsolve(chQeta,forwardsolve(t(chQeta), leta[,i]) + rnorm(k.star))
  
  # Step 4: sample phi.jh
  phi.jh = matrix(rgamma(n = p*k.star, shape = (nu + 1)/2,
                         rate = (nu + Lambda^2*matrix(rep(tau.h, each = p), nr = p))/2), nr = p) #for(h in 1:k.star){for(j in 1:p) phi.jh[j,h] = rgamma(n = 1, shape = (nu + 1)/2, rate = (nu + Lambda[j,h]^2*tau.h[h])/2)
  
  # Step 5: sample tau.h via delta.h
  delta.h = sampleMGP(theta.jh = sqrt(phi.jh)*Lambda, delta.h = delta.h, a1 = a1, a2 = a2)
  tau.h = cumprod(delta.h)
  
  # Step 6: Sample alpha

  eps = Z - tcrossprod(eta, Lambda)
  for(j in 1:one.hot.cols){
    lalpha = sum(eps[,j])/Sigma.diag[j]
    Qalpha = (n/Sigma.diag[j] + 1)^-1
    alpha[j] = rnorm(1,Qalpha*lalpha,sqrt(Qalpha))
  }

  
  # Step 7: Re-sample Z
  for(j in sample(p)){ 
    #marginal moments
    alphaj = alpha[j] +Lambda[j,]%*%t(eta)
    sdj = sqrt(Sigma.diag[j])
    if(is_cat_bin[j] == 0){ #if j is numeric
      if(!plugin.marginal[j]){
        for(r in 1:Rlevels[j]){
          
          ir<-(1:n)[R[,j] == r & !is.na(R[,j])]
          lb<- suppressWarnings(max(Z[R[,j] == r-1,j], na.rm = T))
          ub <- suppressWarnings(min(Z[R[,j] == r+1,j], na.rm = T))
          Z[ir,j] = qnorm(runif(length(ir), pnorm(lb, alphaj[ir],sdj),pnorm(ub,alphaj[ir],sdj)),alphaj[ir],sdj)}
      }
    }
    
    else{ #if j is a categorical sub-level
      
      Z[which(Y_mod[,j] == 1),j] =   rtruncnorm(n = 1,a =0 ,mean = alphaj[which(Y_mod[,j] == 1)],sd = sdj)
      Z[which(Y_mod[,j] == 0), j] = rtruncnorm(n = 1,b = 0,mean = alphaj[which(Y_mod[,j] == 0)],sd = sdj)
      
    }
    
  }
  
  
  
  
  ### Save Results
  if(ns%%odens == 0){
    Omega = tcrossprod(Lambda) + diag(Sigma.diag)
    C<-Omega/(sqrt(diag(Omega)) %*% t(sqrt(diag(Omega)))) #copula correlation matrix
    C.psamp[ns/odens,,] <- C
    alpha.psamp[ns/odens,]= diag(diag(Omega^-.5))%*%alpha #scale intercept vector
    if ((ns%%100) == 0) {
      cat(round(100 * ns/nsamp), "percent done ", date(), 
          "\n")
    }
  }
  
  if(any(is.infinite(Z)) | any(is.na(Z))){
    print('oo')
    break
  }
  
}
