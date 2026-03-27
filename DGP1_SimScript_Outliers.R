### Simulation for DGP1 section - Outlier data

library(zCompositions)
library(coda.base)
library(pROC)
library(caret)
library(doParallel)
library(doRNG)


rm(list=ls())

# Calculate the number of cores
no_cores <- 18
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# FUNCTIONS ##################################################################

## Generate alr-coord time series from symbolised Markov-Switching time series
## (for 2 states + 1 alr coord)

create.data <- function(n.series , l.series, prob00 = 0.9, prob11 = 0.7, alpha = 1.5, sigma = 1){
  
  nse <- n.series # Number of states
  capt <- l.series # Length of the series
  pphi <- 1
  T <- capt + pphi
  ns <- 2 # Number of primitive states
  nk <- pphi + 1 # First obs for evaluation
  captst <- capt - pphi # captst is the effective sample size
  
  alpha <- alpha # overall slope of the series, difference between means
  sig <- sigma # time series variability
  p00k <- prob00 # prob inertia state 0 - expansion
  p11k <- prob11 # prob inertia state 1 - recession
  
  x2 <- runif(T)   # MS trend
  
  ff <- rep(0,T)  
  for (i in 2:T){
    if (ff[i-1] < 0.5){ # Expansion
      if (x2[i] < p00k){ff[i] <- 0}
      else {ff[i] <- 1}
    }
    if (ff[i-1] >= 0.5){ # Recession 
      if (x2[i] < p11k){ff[i] <- 1}
      else {ff[i] <- 0}
    }
  }
  
  for (it in 2:T-1){
    if((ff[it-1]==0) && (ff[it]==1) && (ff[it+1]==0)){ff[it]=0}
    if((ff[it-1]==1) && (ff[it]==0) && (ff[it+1]==1)){ff[it]=1}
  }
  
  y <- matrix(0, nrow = 1, ncol = nse)
  for (it in 2:T){
    if ((ff[it]==0)){
      yit=alpha+y[it-1]+sqrt(sig)*matrix(rnorm(nse), nrow = 1)
    }
    if ((ff[it]==1)){
      yit=-alpha+y[it-1]+sqrt(sig)*matrix(rnorm(nse), nrow = 1)
    }
    y <- rbind(y, yit)
  }
  
  y <- y[2:nrow(y),]
  colnames(y) <- paste("S",1:ncol(y),sep="")
  ff <- ff[2:(length(ff))]
  
  return(data.frame(State = ff, y))
  
}

make.outlier <- function(dat, p){
  
  # dat: including first State column
  # p: p% of time series including outlying values
  
  outs <- function(x){
    b <- boxplot(x, plot = FALSE)$stats[c(1,5)] 
    if (runif(1,-1,1) > 0) {return(b[2] + b[2]*sample(c(0.2,0.4,0.6),1))}
    else {return(b[1] - b[1]*sample(c(0.2,0.4,0.6),1))}
  }
  
  if (p == 0){return(dat)} # To include de no outliers case
  
  State <- dat[,"State"] # Exclude state column
  dat <- dat[,-1]
  
  dat <- t(dat) # Series by rows for this
  rows <- sample(1:nrow(dat),round(p*nrow(dat),0)) # random selection of p% series (rows) with outliers
  perc.cols <- runif(length(rows),0.1,0.3) # Between 10-30% elements of a series containing outliers for each row
  num.cols <- round(perc.cols*ncol(dat),0)
  for (i in 1:length(rows)){
    cols <- sample(ncol(dat),num.cols[i]) # random selection of cols with imposed outliers applied to each row
    dat[rows[i], cols] <- apply(dat[, cols], 2, outs)
  }
  
  
  return(data.frame(State=State,t(dat))) # Back to series by columns
}

symbolicTScoords_1ALR <- function(X, m = 3){ # Generate symbolic time series + alr coord representation
  
  # X: time series by columns, including first State column
  # m: length of the m-histories
  
  symbolize <- function(x, m){
    
    # x: time series
    # m: m-history length
    
    s<-embed(x, m)
    nas <- !complete.cases(s)
    s <- t(apply(s, 1, order)) # orders
    s <- do.call("paste", c(as.data.frame(s), sep = ""))
    s[nas] <- NA
    return(as.factor(s))
  }
  
  perm <- function(v) {
    n <- length(v)
    if (n == 1) v
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
      X
    }
  }
  S <- perm(c(1,2,3)) # All possible symbols
  S <- do.call("paste", c(as.data.frame(S), sep = ""))
  
  sym.freqs <- function(x, S){
    
    # x: Symbols across series at a time point
    # S: set of all possible symbols
    
    levels(x) <- c(levels(x),setdiff(S,levels(x)))
    x <- factor(x, levels = S)
    freq <- table(x[x%in%S])
    freq
  }
  
  # SYMBOLIZE TIME SERIES (BASED ON 3-HISTORIES) + COMPOSITION
  
  X.sym <- apply(X[,-1],2,function(x) symbolize(x,m=m))
  X.sym.freq <- t(apply(X.sym,1,function(x) sym.freqs(x,S)))
  X.sym.prop <- t(apply(X.sym.freq, 1, function(x) x/sum(x)))
  X.sym.prop.sub <- cbind(Other = rowSums(X.sym.prop[,!(colnames(X.sym.prop)%in%c("321"))]),
                          "321"=X.sym.prop[,"321"])
  
  X.sym.prop.sub <- data.frame(X.sym.prop.sub)
  
  # Deal with zeros if any
  if (any(X.sym.prop.sub==0)){

    mins <- apply(X.sym.prop.sub,2,function(x) min(x[x!=min(x)])) # Thresholds for replacement

    X.sym.prop.sub.imp <- multLN(X.sym.prop.sub, label = 0, dl = mins)
  }
  if (any(X.sym.prop.sub==0)==FALSE){
    X.sym.prop.sub.imp <- X.sym.prop.sub
  }
  
  # alr coordinates
  X.sym.prop.sub.imp.coord <- coordinates(X.sym.prop.sub.imp, basis = 'alr')
  
  return(cbind(State=X$State[m:length(X$State)],X.sym.prop.sub.imp.coord))
  
}

## MS model fit to alr coords
# (State 1 is recession, greyed area in plot)
fitSYM_1ALR_TS_2STATES <- function(X){
  
  # X: matrix [states, alr.coord]
  
  # Required functions
  
  ofn <- function(th2,y){
    th<-trans(th2)
    captst<-nrow(y)
    pt0<-1/2
    pt1<-1/2
    vari<-th[3]
    p00 <- th[4]
    p10<- 1-p00
    p11 <- th[5]
    p01 <- 1-p11
    nume<-(2*pi)*(vari^0.5)
    ft<-0
    for (it in 1:captst){
      e0t=y[it,]-th[1]
      e1t=y[it,]-th[2]
      ft0 <- (1/nume)*exp(-0.5*e0t^2/vari)
      ft1 <- (1/nume)*exp(-0.5*e1t^2/vari)
      ftotal <- pt0*ft0+pt1*ft1
      
      ptt0 <- pt0*ft0/ftotal
      ptt1 <- pt1*ft1/ftotal
      
      prob0[it]<<-ptt0
      prob1[it]<<-ptt1
      
      ft <- ft + log(ftotal)
      
      pt0 <- p00*ptt0+p01*ptt1
      pt1 <- p10*ptt0+p11*ptt1
    }
    return(-ft)
  }
  trans <- function(teta){
    th <-teta
    th[3] <- teta[3]^2
    th[4]<- teta[4]^2/(1+teta[4]^2)
    th[5]<- teta[5]^2/(1+teta[5]^2)
    return(th)
    
  }
  
  BC<-X$State # vector of states
  AY <- cbind(alr1=X[,-1]) # column of alr coords
  T<-length(AY)
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  
  NBER<-BC
  state0<-which(NBER==0)
  state1<-which(NBER==1)
  
  mu0<-mean(AY[state0]) # Initial conditions based on actual states
  mu1<-mean(AY[state1])
  
  v<-1
  p_00 <- 1.5;p_11<-1.5
  startval<-c(mu0,mu1,v,p_00,p_11)
  
  # Optimisation
  
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  opti <- optim(startval,ofn,gr="BFGS",AY)
  theta <- opti$par
  thn <- trans(theta)
  
  res <- as.data.frame(cbind(State=BC,prob0,prob1)) # Probability of states over time
  
  # Identify prob linked to State 1
  # Trick to match up state with prob column
  BS1 <- mean((res$State-res[,"prob0"])^2, na.rm = T)
  BS2 <- mean((res$State-res[,"prob1"])^2, na.rm = T)
  refprob <- ifelse(BS1 < BS2, "prob0", "prob1")
  
  return(list(res,refprob))
  
}

## Performance
distortion <- function(fit.result){
  
  refprob <- fit.result[[2]]
  fit.result <- fit.result[[1]]
  
  # Performance assessment
  AUROC <- roc(fit.result$State,fit.result[,refprob])$auc
  BS <- mean((fit.result$State-fit.result[,refprob])^2, na.rm = T)
  
  obs <- as.factor(fit.result$State) # Need to be factors to use caret
  pred <- as.factor(ifelse(fit.result[,refprob] >= 0.5,1,0))
  pred <- factor(pred,levels=levels(obs)) # Just in case a level is missing or different order
  
  cm <- confusionMatrix(pred, obs, mode = "everything", positive="1") # Pred vs Obs
  cm$table <- round(prop.table(cm$table),4)
  metrics <- data.frame("AUC"=AUROC,"Brier"=BS,matrix(cm$byClass,nrow=1,dimnames=list("",names(cm$byClass))))
  
  return(metrics)
}

#############################
# Performing the simulations
#############################

## Parameters definition

nsim <- 500 # number of simulations

# N and D are data matrix no. rows and columns

N <- c(250) 
D <- c(250)

p11 <- 0.9
p22 <- c(0.3,0.7) # Inertia of the second state
sigma <- c(1,2) # variance of the error term of the MSM
alpha <- c(0.5,1.5) # business cycle regimes, distance between slopes at each state {1,3} -> {0.5,1.5} because around zero
p <- c(0,0.2,0.4,0.6,0.8) # prop of series affected
method <- c("Symbolic")

Parameters.temp <- expand.grid(N,D,p22,sigma,alpha,p,method)
Parameters.temp <- Parameters.temp[order(Parameters.temp$Var1,Parameters.temp$Var2,
                                         Parameters.temp$Var3,Parameters.temp$Var4,
                                         Parameters.temp$Var5,Parameters.temp$Var6,Parameters.temp$Var7),]

set.seed(534534) # Generate random number to generate the same data set for different methods in the otherwise same scenario
Parameters.temp$Seed <- rep(round(rnorm(nrow(Parameters.temp)/length(method),500,200)),each=length(method))  

rownames(Parameters.temp) <- NULL
colnames(Parameters.temp)<-c("N","D","p22","sigma","alpha","%","Method","Seed")

Parameters <- list()
ncases <- nrow(Parameters.temp)

for(case in 1:ncases){
  # Configuration parameters
  N <- Parameters.temp[case,1]
  D <- Parameters.temp[case,2]
  p22 <- Parameters.temp[case,3]
  sigma <- Parameters.temp[case,4]
  alpha <- Parameters.temp[case,5]
  p <- Parameters.temp[case,6]
  method <- Parameters.temp[case,7]
  seed <- Parameters.temp[case,8]
  
  Parameters[[case]]<-list(N=N,D=D,p22=p22,sigma=sigma,alpha=alpha,p=p,method=method,seed=seed)
}

## Performing all the simulations

set.seed(12213123)

results <- list()

for (case in 1:ncases){
  
  N <- Parameters[[case]]$N
  D <- Parameters[[case]]$D
  p22 <- Parameters[[case]]$p22
  sigma <- Parameters[[case]]$sigma
  alpha <- Parameters[[case]]$alpha
  p <- Parameters[[case]]$p
  method <- Parameters[[case]]$method
  seed <- Parameters[[case]]$seed
  
  print(paste("case:",case))
  
  set.seed(seed)
  
  TAB <- foreach(k = 1:nsim,
                 .packages = c("zCompositions","coda.base","pROC","caret"),
                 .combine=rbind) %dorng%
    {
      
      #simulation for the current case
      
      res <- list()
      for (j in 1:100){ # Try up to 100 times, just in case of numerical error
        res[[j]] <- try(create.data(n.series = N, l.series = D, prob00 = p11, prob11 = p22,
                                    alpha = alpha, sigma = sigma), silent = TRUE)
        if ((!is.character(res[[j]])) && (length(table(res[[j]][,"State"]))==2)){ # Check that State includes two states
          data <- res[[j]]
          break
        }
      }
      # Impose outliers
      mdata <- make.outlier(dat = data, p = p) # data set with outliers imposed
      # Fit MSM
      if (method == "Symbolic"){
        symcoordata <- symbolicTScoords_1ALR(mdata, m = 3)
        res <- list()
        for (j in 1:100){ # Try up to 50 times, just in case of numerical error
          res[[j]] <- try(fitSYM_1ALR_TS_2STATES(symcoordata), silent = TRUE)
          if (!is.character(res[[j]])){
            fitmodel <- res[[j]]
            break
          }
        }
      }
      if (method == "Standard"){
        symcoordata <- symbolicTScoords_1ALR(mdata, m = 3) # Exclude actual state column
        res <- list()
        for (j in 1:100){ # Try up to 100 times, just in case of numerical error
          res[[j]] <- try(fitSYM_1ALR_TS_2STATES(symcoordata), silent = TRUE)
          if (!is.character(res[[j]])){
            fitmodel <- res[[j]]
            break
          }
        }
      }
      
      # Distortion measures
      metrics <- distortion(fitmodel)
      
      return(metrics)
    }
  
  # Create row with details of scenario and results
  results[[case]] <- data.frame("Case"=case,"n"=N,"D"=D,"p22"=p22,"Sigma"=sigma,
                                "Alpha"=alpha,"p"=p,
                                "Method"=as.character(method),
                                "AUC"=colMeans(TAB, na.rm = T)[1],
                                "Brier"=colMeans(TAB, na.rm = T)[2],
                                "Sensitivity"=colMeans(TAB, na.rm = T)[3],
                                "Specificity"=colMeans(TAB, na.rm = T)[4],
                                "Pos Pred Value"=colMeans(TAB, na.rm = T)[5],
                                "Neg Pred Value"=colMeans(TAB, na.rm = T)[6],
                                "Precision"=colMeans(TAB, na.rm = T)[7],
                                "Recall"=colMeans(TAB, na.rm = T)[8],
                                "F1"=colMeans(TAB, na.rm = T)[9],
                                "Prevalence"=colMeans(TAB, na.rm = T)[10],
                                "Detection Rate"=colMeans(TAB, na.rm = T)[11],
                                "Detection Prevalence"=colMeans(TAB, na.rm = T)[12],
                                "Balanced Accuracy"=colMeans(TAB, na.rm = T)[13])
  
}

results <- do.call(rbind,results)
rownames(results) <- NULL

stopCluster(cl)

