### Simulation for DGP2 section - Breaks

library(zCompositions)
library(coda.base)
library(foreach)
library(doParallel)

rm(list=ls())

n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

# FUNCTIONS ##################################################################

## Generate alr-coord time series from symbolised Markov-Switching time series

create.data.breaks <- function(n.series , l.series, alpha = 1, sigma = 1){
  
  nse <- n.series # Number of series
  capt <- l.series # Length of the series
  pphi <- 1
  T <- capt + pphi
  ns <- 4 # Number of primitive states
  nk <- pphi + 1 # First obs for evaluation
  captst <- capt - pphi # captst is the effective sample size
  mu1k <- 1
  mu2k <- mu1k + alpha*3
  mu3k <- mu1k - alpha*1
  mu4k <- mu1k - alpha*2
  phik <- 0.3
  sigk <- sigma
  
  y <- matrix(0, nrow = 1, ncol = nse)
  ff<-0
  
  for (it in nk:floor(T/4)){
    yit <- mu1k + phik*(y[it-1]-mu1k) + sqrt(sigk)*matrix(rnorm(nse), nrow = 1)
    y <- rbind(y, yit)
    ff<- rbind(ff, 1)
  }
  for (it in (floor(T/4)+1):floor(T/2)){
    yit <- mu2k + phik*(y[it-1]-mu2k) + sqrt(sigk)*matrix(rnorm(nse), nrow = 1)
    y <- rbind(y, yit)
    ff<- rbind(ff, 2)
  }
  for (it in (floor(T/2)+1):floor(3*T/4)){
    yit <- mu3k + phik*(y[it-1]-mu3k) + sqrt(sigk)*matrix(rnorm(nse), nrow = 1)
    y <- rbind(y, yit)
    ff<- rbind(ff, 3)
  }
  for (it in (floor(3*T/4)+1):T){
    yit <- mu4k + phik*(y[it-1]-mu4k) + sqrt(sigk)*matrix(rnorm(nse), nrow = 1)
    y <- rbind(y, yit)
    ff<- rbind(ff, 4)
  }
  
  y <- y[nk:nrow(y),]
  colnames(y) <- paste("S",1:ncol(y),sep="")
  ff <- ff[nk:length(ff)]
  
  return(data.frame(State = ff, y))
  
}

symbolicTScoords_2ALR <- function(X, m = 3){ # Generate symbolic time series + alr coord representation
  
  # X: time series by columns, including first State column
  # m: length of the m-histories
  
  symbolize <- function(x, m){
    
    # x: time series
    # m: m-history length
    
    #s <- t(apply(embed(x, m), 1, rev)) # m-histories
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
fitSYM_2ALR_TS_3STATES.breaks <- function(X, width = 0){
  
  # X: matrix [states, alr.coord]
  
  # Required functions
  
  # width: width of the surrounding of the break point to consider when checking obs - pred matchup
  
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
    p10 <- 1-p00
    p20 <- 0
    p01 <- th[11]
    p11 <- th[12]
    p21 <- 1-p01-p11
    p02 <- 0
    p12 <- 1-th[13]
    p22 <- th[13]
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
    th[10]<- teta[10]^2/(1+teta[10]^2)
    th[11] <- teta[11]^2/(1+teta[11]^2+teta[12]^2) 
    th[12] <- teta[12]^2/(1+teta[11]^2+teta[12]^2) 
    th[13] <- teta[13]^2/(1+teta[13]^2) 
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
  
  BC<-X$State # vector of states
  AY <-cbind(X$alr1,X$alr2)
  T<-length(AY)
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  prob2<-rep(0,captst)
  
  mu_b1_s0<--3;mu_b1_s1<--1;mub1_s2<-0; 
  mu_b2_s0<-1;mu_b2_s1<--1;mub2_s2<--2;
  cova<-0.5;v2<-1;v3<-1
  p_00 <- 0.5;p_01<-0.5;p_11<-3;p_22<-1
  startval<-c(mu_b1_s0,mu_b2_s0,mu_b1_s1,mu_b2_s1,mub1_s2,mub2_s2,cova,v2,v3,p_00,p_01,p_11,p_22)
  
  # Optimisation
  
  captst<-nrow(AY)
  prob0<-rep(0,captst)
  prob1<-rep(0,captst)
  prob2<-rep(0,captst)
  opti <- optim(startval,ofn,gr="BFGS",AY)
  theta <- opti$par
  thn <- trans(theta)
  
  res <- as.data.frame(cbind(State=BC,prob0,prob1,prob2))
  
  # Split obs state into steady/break (0 vs 1) ("status" used to refer to this binary partition)
  obs.status <- c(0,diff(res$State)) # Mark start of new original state
  # Split pred state (up, steady, down) into steady/break
  pred.status <- apply(res[,-1],1, which.max) # Identify highest probability for each observation
  pred.status <- ifelse(pred.status == names(which.max(table(pred.status))), 0, 1) # Identify steady as the most frequent and assign steady/break labels
  
  break.at <- which(obs.status==1)
  matches <- rep(0,3)
  
  for(i in 1:3){ # Define matching-up neighbourhood around the three break points
    if(any(pred.status[(break.at[i]-width):(break.at[i]+width)] == 1)){matches[i] <- 1}
  }
  
  return(matches)
}

#############################
# Performing the simulations
#############################

## Parameters definition

nsim <- 500

# N and D are data matrix no. rows and columns
N <- c(250) 
D <- c(250) 
sigma <- c(1,2,3,4) # variance of the error term of the MSM
alpha <- c(0.5,1,2,3,4,5) # separation of means
method <- c("Symbolic")
width <- c(1,2,3) # Width of neighbourhood around the break point to check match-up

Parameters.temp <- expand.grid(N,D,sigma,alpha,method,width)
Parameters.temp <- Parameters.temp[order(Parameters.temp$Var1,Parameters.temp$Var2,
                                         Parameters.temp$Var3,Parameters.temp$Var4,
                                         Parameters.temp$Var5,Parameters.temp$Var6),]

set.seed(534534) # Generate random number to generate the same data set for different methods in the otherwise same scenario
Parameters.temp$Seed <- rep(round(rnorm(nrow(Parameters.temp)/length(method),500,200)),each=length(method))  

rownames(Parameters.temp) <- NULL
colnames(Parameters.temp)<-c("N","D","sigma","alpha","Method","Width","Seed")

Parameters <- list()
ncases <- nrow(Parameters.temp)

for(case in 1:ncases){
  # Configuration parameters
  N <- Parameters.temp[case,1]
  D <- Parameters.temp[case,2]
  sigma <- Parameters.temp[case,3]
  alpha <- Parameters.temp[case,4]
  method <- Parameters.temp[case,5]
  width <- Parameters.temp[case,6]
  seed <- Parameters.temp[case,7]
  
  Parameters[[case]]<-list(N=N,D=D,sigma=sigma,alpha=alpha,method=method,width=width,seed=seed)
}

## Performing all the simulations

set.seed(12213123)

results <- list()

for (case in 1:ncases){

  N <- Parameters[[case]]$N
  D <- Parameters[[case]]$D
  sigma <- Parameters[[case]]$sigma
  alpha <- Parameters[[case]]$alpha
  method <- Parameters[[case]]$method
  width <- Parameters[[case]]$width
  seed <- Parameters[[case]]$seed
  
  print(paste("case:",case))
  
  set.seed(seed) # To generate the same data sets for different methods in the otherwise same scenario
  
  pred <- foreach(k = 1:nsim,
                  .combine = 'c',
                  .packages = c("zCompositions","coda.base")) %dopar% {
  
        #simulation for the current case
        data <- create.data.breaks(n.series = N, l.series = D,
                                   alpha = alpha, sigma = sigma)
        # Fit MSM
          symcoordata <- symbolicTScoords_2ALR(data, m = 3)
          res <- list()
          for (j in 1:100){ # Try up to 50 times, just in case of numerical error
            res[[j]] <- try(fitSYM_2ALR_TS_3STATES.breaks(symcoordata, width = width), silent = TRUE)
            if (!is.character(res[[j]])){
              fitmodel <- res[[j]]
              break
            }
          }
    return(fitmodel)      
  }
  
  results[[case]] <- data.frame("Case"=case,"n"=N,"D"=D,"Sigma"=sigma,
                                "Alpha"=alpha,
                                "Method"=as.character(method),
                                "Width"=width,
                                "Accuracy"=mean(pred, na.rm = T))
}

results <- do.call(rbind,results)
rownames(results) <- NULL

