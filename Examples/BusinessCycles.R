### Simulate Markov switching time series
### Generate associated compositional symbolic time series in coordinates (CSC)



#install.packages("mixtools")
#library(mixtools)
library(zCompositions)
library(coda.base)
library(pROC)
library("rootSolve")
remove(list=ls())   # Remove global environment
cat("\f")           # Clear the screen
graphics.off()      # Close the current graphical device
set.seed(2)
require("readxl")
require("writexl")
#set.seed(2344)
#############################################################
#                   REQUIRED FUNCTIONS                      #
#############################################################

#-----------------------------------------------------------#
#         Extract symbols from time series                  #
#-----------------------------------------------------------#
# Extract symbols from time series (based on m-histories)
symbolize <- function(x, m){
  
  # x: time series
  # m: m-history length
  
  #s <- t(apply(embed(x, m), 1, rev)) # m-histories
  s<-embed(x, m)
  s <- t(apply(s, 1, order)) # orders
  return(do.call("paste", c(as.data.frame(s), sep = "")))
  
}

#-----------------------------------------------------------#
#            Compute all possible symbols in S              #
#-----------------------------------------------------------#
## Space of all symbols 3-histories
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

#-----------------------------------------------------------#
#    Compute frequencies of (all possible) symbols in S     #
#-----------------------------------------------------------#
# Compute frequencies of (all possible) symbols across time series at each time point
sym.freqs <- function(x, S){
  
  # x: Symbols across series at a time point
  # S: set of all possible symbols
  
  levels(x) <- c(levels(x),setdiff(S,levels(x)))
  x <- factor(x, levels = S)
  freq <- table(x[x%in%S])
  freq
}

#-----------------------------------------------------------#
#      Procedures to estimate Markov-switching models       #
#-----------------------------------------------------------#

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

#############################################################
#                                                           #
#                   The code start here                     #
#                                                           #
#############################################################
je<-0           # je=0 to skip hessian calculations 

#############################################################
#                      Loading data                         #
#############################################################

mydata2<-read_xlsx("coincident_indicators.xlsx",sheet=1)
mydata <- mydata2[4:nrow(mydata2),3:ncol(mydata2)]
fecha<-data.matrix(mydata2[4:nrow(mydata2):(nrow(mydata2)-2),1])   # Recall m=2
NBER<-data.matrix(mydata2[4:nrow(mydata2),2])
AY<-data.matrix(mydata)

# nAY<-nrow(AY)
# AY2<-100*(AY[2:nAY,]-AY[1:(nAY-1),])/AY[1:(nAY-1),]
# NBER<-NBER[2:nrow(NBER),]
# AY<-AY2

X <- cbind(NBER,AY)

indica<-ts(X[,-1],start=c(1979,04),frequency=12)
states2<-ts(NBER,start=c(1979,04),frequency=12)

plot(states2,type='h',col="grey75",lwd=2.0,
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(indica,type='l',ylab="",xlab="",xaxt="n")


#############################################################
#      Symbolise time series (based on 3-histories)         #
#############################################################
## 2. Symbolise time series (based on 3-histories)
X.sym <- apply(X[,-1],2,function(x) symbolize(x,m=3))

## 3. Symbol frequencies across series at each time point (composition)
X.sym.freq <- t(apply(X.sym,1,function(x) sym.freqs(x,S)))
# In proportions
X.sym.prop <- t(apply(X.sym.freq, 1, function(x) x/sum(x)))
  # Example distribution for first row of X.sym.prop
barplot(X.sym.prop[1,])

## 4. Work with subcomposition (123,321,Other)
X.sym.prop.sub <- cbind("123"=X.sym.prop[,"123"],"321"=X.sym.prop[,"321"],
                        Other = rowSums(X.sym.prop[,!(colnames(X.sym.prop)%in%c("123","321"))]))

Pi_123<-X.sym.prop[,"123"]
Pi_321<-X.sym.prop[,"321"]
Pi_oth<-rowSums(X.sym.prop[,!(colnames(X.sym.prop)%in%c("123","321"))])
Pi_213<-X.sym.prop[,"213"]
Pi_132<-X.sym.prop[,"132"]
Pi_231<-X.sym.prop[,"231"]
Pi_312<-X.sym.prop[,"312"]


NBER<-NBER[3:length(NBER)] # Recall m=3
#states2<-ts(NBER,start=c(1979,06),frequency=12)

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 123",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_123,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 321",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_321,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of others",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 213",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_213,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

# plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 213 on onthers",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)   
# matplot(Pi_213/Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 132",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_132,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

# plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 132 on onthers",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)   
# matplot(Pi_132/Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 231",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(Pi_231,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

# plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 231 on onthers",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)   
# matplot(Pi_231/Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 312",
     ylab="",xlab="",yaxt="n")
par(new = TRUE)   
matplot(Pi_312,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

# plot(states2,type='h',col="grey75",lwd=2.0,main="Proportion of 312 on onthers",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)   
# matplot(Pi_312/Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))


#############################################################
#               Imputation of zeroes                        #
#############################################################
X.sym.prop.sub <- data.frame(X.sym.prop.sub)
# Deal with zeros if any
# if (any(X.sym.prop.sub==0)){
# zPatterns(X.sym.prop.sub, label = 0)
# mins <- apply(X.sym.prop.sub,2,function(x) min(x[x!=min(x)]))    # Thresholds for replacement
# X.sym.prop.sub.imp <- lrEM(X.sym.prop.sub, label = 0, dl = mins) # Impute zeros
# }
if (any(X.sym.prop.sub==0)){
  zPatterns(X.sym.prop.sub, label = 0,plot = FALSE)
  mins <- apply(X.sym.prop.sub,2,function(x) min(x[x!=min(x)]))    # Thresholds for replacement
  X.sym.prop.sub.imp <- multLN(X.sym.prop.sub, label = 0, dl = mins) # Impute zeros
}

if (any(X.sym.prop.sub==0)==FALSE){
  X.sym.prop.sub.imp<-X.sym.prop.sub
}

#############################################################
#       ALR coordinates of symbols (123,321,Other)          #
#############################################################

# Sequential binary partition to define basis for custom balances
sbp <- sbp_basis(list(b1 = Other ~ X123 + X321, 
                 b2 = X123 ~ X321),
                 data = X.sym.prop.sub)

X.sym.prop.sub.imp.coord <- coordinates(X.sym.prop.sub.imp, basis = 'alr')
# b2 is balance coordinate of symbol 123 vs. 321 (relative weight of one vs. the other)

plot(states2,type='h',col="grey75",lwd=2.0,main="Balance 123 vs. others",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(X.sym.prop.sub.imp.coord$alr1,type='l',ylab="",xlab="",xaxt="n",lwd=2)

plot(states2,type='h',col="grey75",lwd=2.0,main="Balance 321 vs. others",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(X.sym.prop.sub.imp.coord$alr2,type='l',ylab="",xlab="",xaxt="n",lwd=2)


#############################################################
#                                                           #
#                 Markov-switching in symbols               #
#                                                           #
#############################################################
AY <-cbind(X.sym.prop.sub.imp.coord$alr1,X.sym.prop.sub.imp.coord$alr2)
#AY <-cbind(X.sym.prop.sub[,1],X.sym.prop.sub[,2])
T<-length(AY)
captst<-nrow(AY)
prob0<-rep(0,captst)
prob1<-rep(0,captst)
prob2<-rep(0,captst)

#############################################################
#             INITIAL PARAMS' VALUES                        #
#############################################################
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

#startval<-c(mu01,mu02,mu11,mu12,mu21,mu22,cova,v2,v3,p_00,p_01,p_11,p_22)
startval<-c(mu0[1],mu0[2],mu1[1],mu1[2],mu2[1],mu2[2],cova,v2,v3,p_00,p_11,p_22,p_10,p_01,p_02)
nth<-length(startval)

#############################################################
#             OPTIMIZATION PROCEDURES                       #
#############################################################
captst<-nrow(AY)
prob0<-rep(0,captst)
prob1<-rep(0,captst)
prob2<-rep(0,captst)
opti <- optim(startval,ofn,gr="BFGS",AY)
theta <- opti$par
thn <- trans(theta)

#############################################################
#                   Standard deviations                     #
#############################################################
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

#############################################################
#                       GRAPHS                              #
#############################################################
plot(states2,type='h',col="grey75",lwd=2.0,main="Prob of state 0",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(prob0,type='l',ylab="",xlab="",xaxt="n",lwd=2)

plot(states2,type='h',col="grey75",lwd=2.0,main="Prob of state 1",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(prob1,type='l',ylab="",xlab="",xaxt="n",lwd=2)

plot(states2,type='h',col="grey75",lwd=2.0,main="Prob of state 2",
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)   
matplot(prob2,type='l',ylab="",xlab="",xaxt="n",lwd=2)

#############################################################
#                       AUROC                              #
#############################################################
roc_fun <- roc(NBER,prob0)
AUROC<-roc_fun$auc
BS<-mean((NBER-prob0)^2)
cat("\f")
print("Parameter estimates")
thn
if (je!=0){
  print("Standard deviations")
  std
}
print("AUROC")
AUROC
print("Brier Score")
BS


