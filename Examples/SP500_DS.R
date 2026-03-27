### Simulate Markov switching time series
### Generate associated compositional symbolic time series in coordinates (CSC)



#install.packages("openxlsx")
#library(mixtools)
library(zCompositions)
library(coda.base)
library(pROC)
library("rootSolve")
library(openxlsx)
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

#############################################################
#                                                           #
#                   The code start here                     #
#                                                           #
#############################################################
je<-0           # je=0 to skip hessian calculations 

#############################################################
#                      Loading data                         #
#############################################################
mydata0<-read.xlsx("SP500_DS.xlsx",sheet=1)  # 1973.01-2024.01
VIX<-read.xlsx("SP500_DS.xlsx",sheet=2)      # 1973.01-2024.01
mydata0<-mydata0[2:nrow(mydata0),]  # 1973.01 is empty 
VIX<-data.matrix(VIX[2:nrow(VIX),2])
fecha<-data.matrix(mydata0[,1])
NBER<-data.matrix(mydata0[,2])
mydatakk <- mydata0[,3:ncol(mydata0)]
mydata2 <- mydatakk
mydata2[2:nrow(mydata2), ] <- (mydatakk[2:nrow(mydatakk), ] - mydatakk[1:(nrow(mydatakk) - 1), ]) / mydatakk[1:(nrow(mydatakk) - 1), ]
mydata2<-mydata2[2:nrow(mydata2),]
fecha<-fecha[2:nrow(fecha),]
NBER<-NBER[2:nrow(NBER),]
VIX<-VIX[2:nrow(VIX),]
column_means <- apply(mydata2, 2, function(x) mean(x, na.rm = TRUE))
for (i in nrow(mydata2)) {
  if (!any(is.na(mydata2[i,]))) {
    mydata2[i,] <- mydata2[i,] - column_means
  }
}
mydata<-mydata2^2

AY<-as.matrix(mydata)
X<-data.matrix(AY)

indica<-ts(X,start=c(1973,02),frequency=12)

plot(indica[,1],type='l',col="grey75",lwd=0.1,
     ylab="",xlab="",ylim = NULL,yaxt="n")
par(new = TRUE)
matplot(indica,type='l',ylab="",xlab="",ylim = c(0,0.5),xaxt="n" )

#############################################################
#      Symbolise time series (based on 3-histories)         #
#############################################################
## 2. Symbolise time series (based on 3-histories)
X.sym <- apply(X,2,function(x) symbolize(x,m=3))
X.sym2 <- X.sym

for(j in 1:ncol(X.sym)){
  pos<-which(is.na(X[,j]))
  X.sym2[pos,j]<-NA
}

## 3. Symbol frequencies across series at each time point (composition)
X.sym.freq <- t(apply(X.sym2,1,function(x) sym.freqs(na.omit(x),S)))
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


matplot(Pi_123,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1),main="123")

matplot(Pi_321,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1),main="321")

matplot(Pi_oth,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1),main="others")

#matplot(Pi_213,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

#matplot(Pi_132,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

#matplot(Pi_231,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

#matplot(Pi_312,type='l',ylab="",xlab="",xaxt="n",lwd=2,ylim=c(0,1))

matplot(X.sym.prop.sub,type="l", ylab="",main="Proportion of symbols")
legend("topleft", legend=c("123", "321", "oth"), col=c(1:3), lty=1:3, inset=0.01,cex=0.5)

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
# sbp <- sbp_basis(list(b1 = Other ~ X123 + X321, 
#                 b2 = X123 ~ X321),
#                 data = X.sym.prop.sub)

X.sym.prop.sub.imp.coord <- coordinates(X.sym.prop.sub.imp, basis = 'alr')
# b2 is balance coordinate of symbols 123 and 321 (relative to others)

alr1t<-ts(X.sym.prop.sub.imp.coord$alr1,start=c(1973,03),frequency=12)
plot(alr1t,type='l',ylab="",xlab="",lwd=2,main="ARL1")

alr2t<-ts(X.sym.prop.sub.imp.coord$alr2,start=c(1973,03),frequency=12)
plot(alr2t,type='l',ylab="",xlab="",lwd=2,main="ARL2")

matplot(cbind(X.sym.prop.sub.imp.coord$alr1,X.sym.prop.sub.imp.coord$alr2),type='l',ylab="",xlab="",xaxt="n",lwd=2)
legend("bottomleft", legend=c("123/other", "321/other"), col=c(1:2), lty=1:2, inset=0.01,cex=0.8)


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
# state0<-which(NBER==1)
# state1<-which(NBER==1)
# state2<-which(NBER==0)
# mu0<-colMeans(AY[state0,])
# mu1<-colMeans(AY[state1,])
# mu2<-colMeans(AY[state2,])

mu01<--3;mu02<--0.6    # 321 (up)
mu11<--0.5;mu12<--2.75    # 123 (down)
mu21<--1.5;mu22<--1.5  # Others

#mu21<--2;mu22<--2  # Others

cova<-0.5;v2<-1;v3<-1
p_00 <-1.5;p_11<-1.5;p_22<-1.5;    
p_10<-0.5;p_01<-0.5;p_02<-0.5;

startval<-c(mu01,mu02,mu11,mu12,mu21,mu22,cova,v2,v3,p_00,p_01,p_11,p_22)
#startval<-c(mu0[1],mu0[2],mu1[1],mu1[2],mu2[1],mu2[2],cova,v2,v3,p_00,p_11,p_22,p_10,p_01,p_02)
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
VIXt<-ts(VIX,start=c(1973,03),frequency=12)
states2<-ts(NBER,start=c(1973,03),frequency=12)
plot(states2,type='h',col="grey75",lwd=2.0,main="Prob S=0",
     ylab="",xlab="",ylim = NULL,yaxt="n")
lines(VIXt, col = "red", lty = 2)
par(new = TRUE)   
prob0t<-ts(prob0,start=c(1973,03),frequency=12)
matplot(prob0t,type='l',ylab="",xlab="",xaxt="n",lwd=2)
par(new = TRUE)
plot(VIXt, col = "red", lty = 2,type="l")


# plot(states2,type='h',col="grey75",lwd=2.0,main="Prob S=1",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)
# prob1t<-ts(prob1,start=c(1973,03),frequency=12)
# matplot(prob1t,type='l',ylab="",xlab="",xaxt="n",lwd=2)
# 
# plot(states2,type='h',col="grey75",lwd=2.0,main="Prob S=2",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)
# prob2t<-ts(prob2,start=c(1973,03),frequency=12)
# matplot(prob2t,type='l',ylab="",xlab="",xaxt="n",lwd=2)
# 
# plot(states2,type='h',col="grey75",lwd=2.0,main="",
#      ylab="",xlab="",ylim = NULL,yaxt="n")
# par(new = TRUE)
# matplot(cbind(alr1t,alr2t),type='l',ylab="",xlab="",xaxt="n",lwd=2)
# legend("bottomleft", legend=c("arl1", "arl2"), col=c(1:2), lty=1:2, inset=0.01,cex=0.6)



