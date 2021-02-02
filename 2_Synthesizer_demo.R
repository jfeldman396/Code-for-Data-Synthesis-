
#Author: Joseph Feldman
#Date: January 14th, 2021
#Description: Algorithm for producing synthetic observations of variables modeled with the copula in Sampler.R. Specifically, 
## the algorithm uses posterior samples of model parameters for posterior predictive sampling, with special marginal synthesis for
## categorical variables, as described in the paper. This algorithm produces synthetic observations of Y_cop,
## which are then used as synthetic covariates for simulation of targeted response variables in 3_NonLinearSynthesizer_demo.R
#------------------------------------------------------------------------------------------------------------------------------------
#### Script for Generating Synthetic Data North Carolina Data
##Inputs
# Y: (N x p) dataframe of mixed types with binary and categorical variables encoded as factors. For this demo,
## we use exclude the targeted variable, count, from Y and synthesize Y_cop from 1_Sampler_demo.R.
# Y_mod: (N x p*) dataframe with column expansions for categorical variables in Y_cop, used in MCMC.
# one.hot.cols: number of one-hot-encoded columns in Y_mod
# C.psamp: (Nsamp x p* x p*) array of posterior samples of copula correlation matrix
# alpha.psamp: (Nsamp x p*) array of posterior samples of mean vector
# Nsims: number of synthetic samples to create
# plugin.marginal: a logical of length p. Determines which variables to use marginal distributions
# and which to use density estimates for synthesis
# Burn: Burn-in samples to discard
##Outputs
# syn: Synthetic data set of exposure profiles (Nsims x p)
##########################################################################

######### R libraries required for algorithm
library(TruncatedNormal) #note v 2.2 or later is required for this algorithm
library(MASS)
library(BMS)

Burn = 5000

######### Get union over all categorical variables and empirical frequencies.

# Pull out levels of categorical variables with > 2 levels and create the multivariate categorical variable
cat.vbls = which(sapply(Y_cop[,factor.cols], nlevels)>2) #categorical columns
cat.cols = nlevels(Y_cop[,cat.vbls]) #length of categorical columns
one.hot.cols = cat.cols + length(which(sapply(Y_cop, nlevels)==2))
lvls = vector("list", length(cat.vbls)) # names of each level among categorical variables
for(v in 1:length(cat.vbls)){
  lvls[[v]] = t(sapply(Y_cop[,cat.vbls[v]],levels))[1,]
}

if(length(cat.vbls != 0)){ 
  dim_sim = NULL
  for(i in 1:length(cat.vbls)){
    dim_sim = c(dim_sim, nlevels(Y_cop[,cat.vbls[i]])) #for simulating categorical realizations
  }
  
  if(length(cat.vbls) == 1){ #if there is only one  categorical variable in the data set
    cat.mat<- model.matrix(~Y_cop[,cat.vbls] - 1)
    probs = colMeans(cat.mat) # Empirical proportions for simulation
  }
  else{#if there is more than one categorical variable
    intrct = interaction(Y_cop[,cat.vbls[1]]) #iterate over interactions to get the multivariate union
    for(i in 2:length(cat.vbls)){
      to_intrct = Y_cop[,cat.vbls[i]]
      intrct = interaction(intrct,to_intrct)
    }
    
    cat.mat = model.matrix(~intrct -1) #binarize the categorical variable
    probs = colMeans(cat.mat)  #get the empirical probabilities of each category
    
  }
  
}

########## Simulate data by sampling from posterior predictive distribution

sim.Z = NULL #realizations of latent variables
cat.vec = NULL#storage for categorical variables
Nsims = dim(Y)[1]
sim.cats<- t(rmultinom(Nsims,1,probs)) #Fix categories

for(i in 1:Nsims){
  if(length(cat.vbls)>0){
    #Which category the observation belongs to
    cat.obs<- array(sim.cats[i,], dim_sim) 
    hot<- which(cat.obs ==1, arr.ind = T) 
    
    c<- rep(0,1)
    ub = NULL
    lb = NULL
    for(j in 1:length(lvls)){ # get truncation for simulation of z.cat and fixed categoriess for RPL
      u = rep(0, dim_sim[j])
      l = rep(-Inf, dim_sim[j])
      l[hot[j]] = 0
      u[hot[j]] = Inf
      ub = c(ub,u)
      lb = c(lb, l)
      c[j] = lvls[[j]][hot[j]]
      
    }
    cat.vec = rbind(cat.vec,c)
    
    # Simulation of z_cat
    C.i = C.psamp[i+Burn,,]
    alpha.i = alpha.psamp[i+Burn,]
    z.cat<- rtmvnorm(1, mu = as.vector(alpha.i[1:cat.cols]), sigma =C.i[1:cat.cols,1:cat.cols],lb = lb, ub = ub)
    
    # conditional moments
    C.s = C.i[-c(1:cat.cols),c(1:cat.cols)]%*%solve(C.i[c(1:cat.cols),c(1:cat.cols)])
    C.c = C.i[-c(1:cat.cols),-c(1:cat.cols)] - C.s%*%C.i[1:cat.cols,-c(1:cat.cols)]
    
    alpha.c = alpha.i[(cat.cols + 1):p] + C.s%*%t(t(z.cat) - alpha.i[1:cat.cols])
    
    #simulation of z_noncat | z_cat
    z.con<-mvrnorm(n =1 , alpha.c, Sigma = C.c)
    
    samp = c(z.cat,z.con)
    sim.Z = rbind(sim.Z,samp)
    
    
  }
  if (i%%100 == 0) {
    cat(round(100 * i/Nsims), "percent done ", date(), 
        "\n")
  }
}



######## Transform latent realizations to correct scale

#Factorize categorical variables for synthetic data
colnames(cat.vec) = colnames(Y_cop[,cat.vbls])
cat.vec = data.frame(cat.vec)

for(i in 1:ncol(cat.vec)){
  
  cat.vec[,i] = factor(cat.vec[,i], levels = lvls[[i]])
  
}

#Binary Variables
bin.cols = which(apply(Y_mod[,-c(1:cat.cols)],2,function(x)length(unique(x))) == 2)
bin<-ifelse(sim.Z[,cat.cols + bin.cols]>0,1,0)
bin<-factor(bin, levels = c("0","1"))

# Numerical Variables
num<- array(0,c(nrow(cat.vec),(p- one.hot.cols)))

quant<- (sapply(sim.Z[,-c(1:one.hot.cols)],pnorm))
col = one.hot.cols + 1

for(i in 1:(p-one.hot.cols)){
  if(plugin.marginal[col]){
    #ensure data is on correct scale: if there are no negative values in the observed data, 
    #there should be no negative observations in the synthetic data
    if(min(Y_mod[,col])>=0){
      sim_col<- exp(BMS::quantile.density(density(log(Y_mod[,col])), probs = quant))  
    }
    else{
      sim_col<- BMS::quantile.density(density(Y_mod[,col]), probs = quant)
    }
    col = col +1
    num[,i] =sim_col
  }
  else{
    
    sim_col<- quantile(Y_mod[,col], probs = quant[,i],type = 1)
    col = col +1
    num[,i] =sim_col
    
  }
}




#join factors, binary and  numerical  data 
syn = cbind.data.frame(cbind(cat.vec,bin))
syn = data.frame(syn, num)
colnames(syn) = c('cat','bin', 'cont')


