#Author: Joseph Feldman
#Date: January 14th, 2021
#Description: Algorithm for producing synthetic observations of targeted response variables - in this case the "count" variable 
## from the simulated data set given synthetic Y_cop
#-------------------------------------------------------------------------------------------------------------------------------------------------------

##Inputs
# Y: full data set, including targeted count variable
# syn: Synthetic Y_cop data produced from 2_Synthesizer_demo.R
## Outputs:  
# syn: synthetic data with newly appended outcomes of interest

###############################################################################################

########## R libraries required for algorithm
library(dbarts)

########## Arrange Data into training and testing, with synthetic data being test data
Y_full <- Y

Y_cols<- c('count') #Response Variables of interest
Y_train = t(t(Y_full[,colnames(Y_full) %in% Y_cols]))
X_train =Y_full[,!(colnames(Y_full) %in% Y_cols)]

X_test = syn[,colnames(X_train)]

#ensure binary variables are numeric
X_test[,c("bin")] = sapply(X_test[,c("bin")],as.numeric)
X_train[,c("bin")] = sapply(X_train[,c("bin")],as.numeric)

########## Initialize RL
R<-NULL
p = dim(Y_train)[2]
n = dim(Y_train)[1]
for(j in 1:p) { R<-cbind(R, match(Y_train[,j],sort(unique(Y_train[,j])))) }
Rlevels<-apply(R,2,max,na.rm=TRUE)
Ranks<- apply(Y_train,2,rank,ties.method="max",na.last="keep")
N<-apply(!is.na(Ranks),2,sum)
U<- t( t(Ranks)/(N+1))
Z<-qnorm(U)



########## Create dataframe for dbarts
dat<- data.frame(Z = Z[,1], X_train)
  


## Hyperparameters
sigest = apply(Z, 2, sd)
sigdf = 3;sigquant = 0.95; k = 2.0; power =2.0;base = 0.95
n.trees = 200

# Initialize samplers with data:
control = dbartsControl(n.chains = 1, n.burn = 0, n.samples = 1,
                        n.trees = n.trees)



sampler = dbarts(Z ~ .,data = dat, test = X_test,
                  control = control,
                  tree.prior = cgm(power, base),
                  node.prior = normal(k),
                  resid.prior = chisq(sigdf, sigquant),
                  sigma = sigest[1])


Burn = 100
Nsims = 1000 + Burn

########## Storage
f.samp = array(0,c(Nsims -Burn, dim(X_test)[1])) #tree samples
sigmas =  array(0,c(Nsims - Burn,1)) #error sd 

########## Sampler: switch  to dat2 and sampler2 upon synthesizing reading scores
dat$Z = Z[,1]
sampler$setResponse(dat$Z)
samp = sampler$run(updateState = TRUE) #initialize f and sigma
for(ns in 1:Nsims){
  muj = samp$train
  sdj = samp$sigma
  #resample Z
  for(r in 1:Rlevels[1]){
    ir<-(1:n)[R[,1] == r & !is.na(R[,1])]
    lb<- suppressWarnings(max(Z[R[,1] == r-1,1], na.rm = T))
    ub <- suppressWarnings(min(Z[R[,1] == r+1,1], na.rm = T))
    Z[ir,1] = qnorm(runif(length(ir), pnorm(lb, muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
  }
  #sample trees
  dat$Z = Z[,1]
  sampler$setResponse(dat$Z)
  samp = sampler$run(updateState = TRUE)
  if(ns > Burn){
    
    #save test tree estimates
    f.samp[ns - Burn,] =  samp$test
    sigmas[ns - Burn,] = samp$sigma
    if (ns%%100 == 0) {
      cat(round(100 * ns/Nsims), "percent done ", date(), 
          "\n")
    }
  }
  
}

########## Copula link function to get data back on scale

count = quantile(Y$count, probs = pnorm(rnorm(dim(syn)[1],mean = colMeans(f.samp), sd = mean(sigmas)) , mean = 0, sd = 1), type=1)
syn[,"count"] = count

########## Basic checks on synthetic data utility

#univariate
summary(syn$count)
summary(Y$count)
plot(density(syn$cont), main = "Synthetic (black) and Observed (red) Continuous Vbl. Density",
     ylim = range(0,.4),xlim = range(Y$cont, syn$cont))
lines(density(Y$cont), col = 2)
table(syn$cat)
table(Y$cat)

#bivariate
cor(syn$cont,syn$count)
cor(Y$cont,Y$count)
xtabs(~syn$cat + syn$bin)/dim(syn)[1]
xtabs(~Y$cat + Y$bin)/dim(Y)[1]


