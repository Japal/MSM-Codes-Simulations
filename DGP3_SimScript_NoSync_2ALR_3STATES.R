### Simulation for DGP3 section - asynchronous time series

# Some series enter recession earlier and some (the same or others) leave recession later than the majority

library(zCompositions)
library(coda.base)
library(pROC)
library(caret)
library(doParallel)
library(doRNG)

rm(list=ls())

# Calculate the number of cores
no_cores <- detectCores() - 2
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# FUNCTIONS ##################################################################

## Generate alr-coord time series from symbolised Markov-Switching time series

create.data.async.v3 <- function(n.series , l.series, async.prop = 0, async.gap = 1,
                                 prob00 = 0.9, prob11 = 0.7, alpha = 1, sigma = 1){
  
  # async.prop: prop of asynchronous time series
  # async.gap: number of lagged time points of asynchronous time series
  #            (it is up to the given number, e.g. 3 means lags of 1, 2 or 3 time points)
  
  # Compute lagged (delayed) time series
  lag.ts <- function(x,lag) {embed(c(rep(NA,lag), x), lag+1)[,lag+1]}
  
  # n.series: number of states
  # l.series: length of series
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
  
  #theta <- c(mu0k,mu1k,phi0k,phi1k,sig0k,sig1k,p00k,p11k)
  x2 <- runif(T)   # MS trend
  
  ff <- rep(0,T)  # ff=0 in expansions and ff=1 in recessions
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
  
  # Identify start and end of change of actual state
  change.id <- diff(ff)
  change.id <- c(0,change.id) # To give same length as ff
  change.entry <- which(change.id==1) # To mark exactly first and last time points of state 1
  change.exit <- which(change.id==-1)
  change.id[change.exit] <- 0
  change.id[change.exit - 1] <- -1
  change.exit <- which(change.id==-1) # Update for the following loop
  
  # Generate % asynced series by expanding to both left and right the recession state up to as determined by async.gap
  
  if (async.prop > 0){
    
    y.async <- matrix(0, nrow = T, ncol = round(nse*async.prop,0))
    gap <- sample(async.gap:async.gap, size = ncol(y.async), replace = T) # gap for each async series.
    
    for (j in 1:ncol(y.async)){
      
      ff.async <- ff
      
      for (i in 1:length(change.entry)){ # For each recession start
        if ((change.entry[i] - gap[j]) > 0){ # Do not apply if below the start of the series
          ff.async[seq(change.entry[i] - gap[j], change.entry[i])] <- 1
        }
      }
      for (i in 1:length(change.exit)){ # For each recession end
        if ((change.exit[i] + gap[j]) <= length(ff)){ # Do not apply if beyond the end of the series
          ff.async[seq(change.exit[i], change.exit[i] + gap[j])] <- 1
        }
      }
      
      tmp <- matrix(0, ncol = 1)
      
      for (it in 2:T){
        if ((ff.async[it]==0)){ # Expansion
          yit <- alpha + tmp[it-1]+sqrt(sig)*matrix(rnorm(1), nrow = 1)
        }
        if ((ff.async[it]==1)){ # Recession
          yit <- -alpha + tmp[it-1]+sqrt(sig)*matrix(rnorm(1), nrow = 1)
        }
        tmp <- rbind(tmp, yit)
      }
      y.async[,j] <- tmp
    }
    
    y.async <- y.async[2:nrow(y.async),]
    
  }
  
  # Generate synced series
  
  if (async.prop > 0){  
    y <- matrix(0, nrow = 1, ncol = nse - ncol(y.async))
  }
  else{
    y <- matrix(0, nrow = 1, ncol = nse)
  }
  
  for (it in 2:T){
    if ((ff[it]==0)){ # Expansion
      yit = alpha + y[it-1]+sqrt(sig)*matrix(rnorm(ncol(y)), nrow = 1)
    }
    if ((ff[it]==1)){ # Recession
      yit= -alpha + y[it-1]+sqrt(sig)*matrix(rnorm(ncol(y)), nrow = 1)
    }
    y <- rbind(y, yit)
  }
  
  y <- y[2:nrow(y),]
  
  # Asynced and synced together
  if (async.prop > 0){
    y <- cbind(y.async,y)
  }
  
  colnames(y) <- paste("S",1:ncol(y),sep="")
  ff <- ff[2:(length(ff))]
  
  return(data.frame(State = ff, y))
  
}
symbolicTScoords_2ALR <- function(X, m = 3){ # Generate symbolic time series + alr coord representation
  
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
  # Reduce to (123,321,Other)-subcomposition - 321 is expansion
  X.sym.prop.sub <- cbind("123"=X.sym.prop[,"123"],"321"=X.sym.prop[,"321"],
                           Other = rowSums(X.sym.prop[,!(colnames(X.sym.prop)%in%c("123","321"))]))
  
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
  
  return(cbind(State=X$State[m:length(X$State)], X.sym.prop.sub.imp.coord))
  
}

## MS model fit to alr coords
# Includes identification of probability vector (out of the 3 corresponding to the fitted states) 
# associated to the original recession state (greyed area in plots, state = 1)
fitSYM_2ALR_TS_3STATES <- function(X){
  
  # X: matrix [states, alr.coord]
  
  # Required functions
  
  ofn<-function(th2,y){
    th<-trans(th2)
    captst<-nrow(y)
    pt0<-1/3
    pt1<-1/3
    pt2<-1/3
    vari<-matrix(NaN,2,2)
    vari[1,1]<-th[8]
    vari[2,2]<-th[9]
    vari[1,2]<-th[7]
    vari[2,1]<-th[7]
    p00 <- th[10]
    p10 <- th[13]
    p20 <- 1-p00-p10
    p01 <- th[14]
    p11 <- th[11]
    p21 <- 1-th[14]-th[11]
    p02 <- th[15]
    p12 <- 1-th[12]-th[15]
    p22 <- th[12]
    nume<-(2*pi)*(det(vari)^0.5)
    ft<-0
    for (it in 1:captst){
      e0t=y[it,]-cbind(th[1],th[2])
      e1t=y[it,]-cbind(th[3],th[4])
      e2t=y[it,]-cbind(th[5],th[6])     
      ft0 <- (1/nume)*exp(-0.5*e0t%*%solve(vari)%*%t(e0t))
      ft1 <- (1/nume)*exp(-0.5*e1t%*%solve(vari)%*%t(e1t))
      ft2 <- (1/nume)*exp(-0.5*e2t%*%solve(vari)%*%t(e2t))
      ftotal <- pt0*ft0+pt1*ft1+pt2*ft2
      
      ptt0 <- pt0*ft0/ftotal
      ptt1 <- pt1*ft1/ftotal
      ptt2 <- pt2*ft2/ftotal
      
      prob0[it]<<-ptt0
      prob1[it]<<-ptt1
      prob2[it]<<-ptt2
      
      ft <- ft + log(ftotal)
      
      pt0 <- p00*ptt0+p01*ptt1+p02*ptt2
      pt1 <- p10*ptt0+p11*ptt1+p12*ptt2
      pt2 <- p20*ptt0+p21*ptt1+p22*ptt2
    }
    return(-ft)
  }
  trans <- function(teta){
    th <-teta
    th[8] <- teta[8]^2
    th[9] <- teta[9]^2
    th[10]<- teta[10]^2/(1+teta[10]^2+teta[13]^2)
    th[13]<- teta[13]^2/(1+teta[10]^2+teta[13]^2)
    th[11] <- teta[11]^2/(1+teta[11]^2+teta[14]^2) 
    th[14] <- teta[14]^2/(1+teta[11]^2+teta[14]^2) 
    th[15] <- teta[15]^2/(1+teta[15]^2+teta[12]^2) 
    th[12] <- teta[12]^2/(1+teta[15]^2+teta[12]^2) 
    return(th)
    
  }
  trans2<- function(teta){
    th <-teta
    th[8] <- teta[8]^2
    th[9] <- teta[9]^2
    th[10]<- exp(teta[10])/(1+exp(teta[10]))
    th[11] <- exp(teta[11])/(1+exp(teta[11])+exp(teta[12])) 
    th[12] <- exp(teta[12])/(1+exp(teta[11])+exp(teta[12])) 
    th[13] <- exp(teta[13])/(1+exp(teta[13])) 
    return(th)
    
  }
  
  je <- 0           # je=0 to skip hessian calculations 
  
  BC<-X$State # vector of states
  AY <-cbind(X$alr1,X$alr2)
  T<-length(AY)
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  prob2<-rep(0,captst)
  
  NBER<-BC
  state0<-which(NBER==1)
  state1<-which(NBER==1)
  state2<-which(NBER==0)
  mu0<-colMeans(AY[state0,])
  mu1<-colMeans(AY[state1,])
  mu2<-colMeans(AY[state2,])
  mu01<--1.5;mu02<-1 
  mu11<-0;mu12<--0.5
  mu21<-3;mu22<--1
  cova<-0.5;v2<-1;v3<-1
  p_00 <-1.5;p_11<-1.5;p_22<-1.5;
  p_10<-1.5;p_01<-0.5;p_02<-0.5;
  
  startval<-c(mu0[1],mu0[2],mu1[1],mu1[2],mu2[1],mu2[2],cova,v2,v3,p_00,p_11,p_22,p_10,p_01,p_02)
  nth<-length(startval)
  
  # Optimisation
  
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  prob2<-rep(0,captst)
  opti <- optim(startval,ofn,gr="BFGS",AY)
  theta <- opti$par
  thn <- trans(theta)
  
  # Standard deviations

  if (je!=0){
    options=list(fnscale=1
                 ,parscale=rep.int(1, nth)
                 ,ndeps=rep.int(1e-3, nth)
    )
    h<-optimHess(theta,ofn,AY,gr=NULL,control = options)
    hi<-solve(h);
    stdor<-diag(hi)^.5;
    gr<-gradient(trans,theta,centered = FALSE, pert = 1e-8)
    Hfin<-gr%*%hi%*%t(gr)
    std=diag(Hfin)^.5
  }
  
  res <- as.data.frame(cbind(State=BC,prob0,prob1,prob2))
  
  # Identify prob linked to State 1
  # Trick to match up state with prob column
  BS1 <- mean((res$State-res[,"prob0"])^2, na.rm = T)
  BS2 <- mean((res$State-res[,"prob1"])^2, na.rm = T)
  BS3 <- mean((res$State-res[,"prob2"])^2, na.rm = T)
  refprob <- names(res)[which.min(c(BS1,BS2,BS3)) + 1]
  
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

nsim <- 500

N <- c(250) 
D <- c(250) 
p11 <- 0.97 
p22 <- 0.95
# NOTE: p11 and p22 set high to avoid excessive transitions that can mixed up with the increasing async gap
sigma <- 1 # c(1,2) # variance of the error term of the MSM
alpha <- 2.5 # distance between slopes at each state {1,3} -> {0.5,1.5} because around zero
p <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # prop of series out of sync
gap <- c(1,5,10,15) # misalignment gap
method <- c("Symbolic")

Parameters.temp <- expand.grid(N,D,p22,sigma,alpha,p,gap,method)
Parameters.temp <- Parameters.temp[order(Parameters.temp$Var1,Parameters.temp$Var2,
                                         Parameters.temp$Var3,Parameters.temp$Var4,
                                         Parameters.temp$Var5,Parameters.temp$Var6,Parameters.temp$Var7,
                                         Parameters.temp$Var8),]

set.seed(534534) # Generate random number to generate the same data set for different methods in the otherwise same scenario
Parameters.temp$Seed <- rep(round(rnorm(nrow(Parameters.temp)/length(method),500,200)),each=length(method))  

rownames(Parameters.temp) <- NULL
colnames(Parameters.temp)<-c("N","D","p22","sigma","alpha","%","Gap","Method","Seed")

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
  gap <- Parameters.temp[case,7]
  method <- Parameters.temp[case,8]
  seed <- Parameters.temp[case,9]
  
  Parameters[[case]]<-list(N=N,D=D,p22=p22,sigma=sigma,alpha=alpha,p=p,gap=gap,method=method,seed=seed)
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
  gap <- Parameters[[case]]$gap
  method <- Parameters[[case]]$method
  seed <- Parameters[[case]]$seed
  
  print(paste("case:",case))
  
  set.seed(seed) # To generate the same data sets for different methods in the otherwise same scenario
  
  TAB <- foreach(k = 1:nsim,
                 .packages = c("zCompositions","coda.base","pROC","caret"),
                 .combine=rbind) %dorng%
    {
      
      #simulation for the current case
      
      res <- list()
      for (j in 1:100){ # Try up to 100 times, just in case of numerical error
        res[[j]] <- try(create.data.async.v3(n.series = N, l.series = D, async.prop = p, async.gap = gap,
                                          prob00 = p11, prob11 = p22,
                                          alpha = alpha, sigma = sigma), silent = TRUE)
        if (!is.character(res[[j]])){
          data <- res[[j]]
          break
        }
      }
      # Fit MSM
      if (method == "Symbolic"){
        symcoordata <- symbolicTScoords_2ALR(data, m = 3)
        res <- list()
        for (j in 1:100){ # Try up to 50 times, just in case of numerical error
          res[[j]] <- try(fitSYM_2ALR_TS_3STATES(symcoordata), silent = TRUE)
          if (!is.character(res[[j]])){
            fitmodel <- res[[j]]
            break
          }
        }
      }
      if (method == "Standard"){
        symcoordata <- symbolicTScoords_2ALR(mdata, m = 3) # Exclude actual state column
        res <- list()
        for (j in 1:100){ # Try up to 100 times, just in case of numerical error
          res[[j]] <- try(fitSYM_2ALR_TS_3STATES(symcoordata), silent = TRUE)
          if (!is.character(res[[j]])){
            fitmodel <- res[[j]]
            break
          }
        }
      }
      
      # Distortion measures

      metrics <- tryCatch(distortion(fitmodel),
                          error=function(e) rep(NA,13))
      return(metrics)
    }
  
  # Create row with details of scenario and results
  results[[case]] <- data.frame("Case"=case,"n"=N,"D"=D,"p22"=p22,"Sigma"=sigma,
                                "Alpha"=alpha,"p"=p,"Gap"=gap,
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
