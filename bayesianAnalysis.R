library(rjags)
library(mnormt)
library(MASS)

setwd("~/Downloads/PartyNationalization/")


## Define useful functions

#Mean center function
stdize <- function(x){
  (x-mean(x,na.rm=TRUE))
}


#Calculate multinomial Cook's Distances
cooks.multinom <- function(obj){
  O <- dim(obj$model)[1]
  V <- vcov(obj)
  theta.hat <- c(t(coef(obj)))
  upd.model <- function(n){
    temp <- update(natModelM,data=natModelM$model[-n,])
    theta.hat.i <- c(t(coef(temp)))
    theta.diff <- theta.hat - theta.hat.i
    return(t(theta.diff)%*%V%*%theta.diff)
    }
  c.distances <- apply(t(1:O),2,upd.model)
}


## Read data in and preprare it for analysis
natData <- read.table("fullData.csv",header=TRUE,sep=",")
natData$Experience <-  natData$Experience + 1

natData$NatParB <- natData$NatCoalPar/100
natData$Homogeneity <- natData$PartyHomogeneityCab

scale <- 1e8  #Make scale hundreds of millions of national currency
natData$Transfers2.t1 <- ceiling(natData$Transfers2.t1/scale)
natData$PublicGoods.t1 <- ceiling(natData$PublicGoods.t1/scale)
natData$TotExpenditures2.t1 <- ceiling(natData$TotExpenditures2.t1/scale)
natData$Others.t1 <- ceiling((natData$TotExpenditures2.t1)-
                             (natData$Transfers2.t1 + natData$PublicGoods.t1)) #Round up.

#Construct outcome variable
outcome <- cbind(natData$Others.t1,
                   natData$PublicGoods.t1, #TARGETABLE!! (2 is targetable)
                   natData$Transfers2.t1) #NON-TARGETABLE!! (3 is nontargetable)
stdData <- apply(natData[,-1],2,stdize)

## SEM requested by Reviewer 2 
## WARNING!!
## THIS IS THE WRONG MODELING STRATEGY. CONTACT AUTHORS FOR BETTER EXPLANATION
##

#firstStage <- lm(stdize(NatParB)~stdize(Homogeneity)
#                     +stdize(Age)
#                     +stdize(Unemploy)
#                     +stdize(AveMag)
#                     +stdize(ETH_FRAC)
#                     +stdize(FISC_DECEN)
#                     +stdize(Leg_ENP)
#                    ,data=natData)
#effects <- array(NA,c(100,2))
#meanPreds <- predict(firstStage)
#sigmaPreds <- diag(summary(firstStage)$sigma^2,length(meanPreds))
#for(i in 1:100){
#firstStageSamp <- mvrnorm(1,mu=meanPreds,Sigma=sigmaPreds)
#natModelSecStage <- multinom(outcome~NatParB
#                     +stdize(Age)
#                     +stdize(Unemploy)
#                     +stdize(AveMag)
#                     +stdize(ETH_FRAC)
#                     +stdize(FISC_DECEN)
#                     +stdize(Leg_ENP)
#                     ,maxit=1000
#                     ,data=natData)
#effects[i,] <- c(summary(natModelSecStage)$coefficients[1,2]
#                 ,summary(natModelSecStage)$coefficients[2,2])
#}


#And continuing with our thing...

#Prepare data to pass to JAGS
n.beta <- 10
sim.nat <- seq(min(stdize(natData$NatParB),na.rm=TRUE),max(stdize(natData$NatParB),na.rm=TRUE),length.out=100)

nat.data <- list(n.obs = dim(natData)[1],
                 n.beta = n.beta,
                 mu.beta2 = rep(0,n.beta),
                 mu.beta3 = rep(0,n.beta),
                 sim.nat = sim.nat,
                 sigma.beta = 0.001*diag(n.beta),
                 y=outcome,
                 X=cbind(1
                   ,stdize(natData$NatParB)
                   ,stdize(natData$Homogeneity)
                   ,stdize(natData$NatParB)*stdize(natData$Homogeneity)
                   ,stdize(natData$Age)
                   ,stdize(natData$Unemploy)
                   ,stdize(natData$AveMag)
                   ,stdize(natData$ETH_FRAC)
                   ,stdize(natData$FISC_DECEN)
                   ,stdize(natData$Leg_ENP)       
                  ),
                   n=natData$TotExpenditures2.t1)

nat.inits <- function(){
  list(beta = cbind(0,matrix(rmnorm(2,varcov=1*diag(n.beta)),nrow=n.beta,ncol=2)))
}


## Estimate Null model
nat.dataNull <- list(n.obs = dim(natData)[1],
                     mu.beta = rep(0,2),
                     sigma.beta = 0.001*diag(2),
                     y=outcome,
                     n=natData$TotExpenditures2.t1)

nat.initsNull <- function(){
  list(beta = c(0,matrix(rmnorm(1,varcov=1*diag(2)),nrow=1,ncol=2)))
}

load.module("dic")  #To generate model comparison 


## Estimate reported model
natModel <- jags.model(file="natModel.jag",
                       data=nat.data,
                       n.chains=2,
                       n.adapt=1.5e6,
                       inits=nat.inits)

results <- coda.samples(natModel,
                        n.iter=2e6,
                        thin=100,
                        variable.names=c("beta","deviance","mareffSim"))

natModelNull <- jags.model(file="natModelNull.jag",
                       data=nat.dataNull,
                       n.chains=2,
                       n.adapt=1.5e6,
                       inits=nat.initsNull)

resultsNull <- coda.samples(natModelNull,
                        n.iter=2e6,
                        thin=100,
                        variable.names=c("deviance"))

#load(file="bayesianModel.RData")
save(results,resultsNull, file="bayesianModel.RData")



#### Prepare data for producing predicted shares plots
probs <- c(0.05,0.5,0.95) # Set level of credible intervals
simNat <- seq(0,1,length.out=100)[19:80]  #Generate simulated nationalization values
predict.x.max <- array(NA,c(1e4,n.beta,100))
predict.x.mean <- array(NA,c(1e4,n.beta,100))
predict.x.min <- array(NA,c(1e4,n.beta,100))

post.samples <- apply(results[[1]],2,sample,size=1e4)

index <- 1:2e4
huh.vector <- stdize(seq(1/100,1,length.out=100))


## Generate predicted values for maximal homogeneity levels
for(i in 1:100){
  predict.x.max[,,i] <- matrix(cbind(1
                                     ,huh.vector[i]
                                     ,max(stdize(natData$Homogeneity),na.rm=TRUE)
                                     ,huh.vector[i]*max(stdize(natData$Homogeneity),na.rm=TRUE)
                                     ,mean(stdize(natData$Age),na.rm=TRUE)
                                     ,mean(stdize(natData$Unemploy),na.rm=TRUE)
                                     ,mean(stdize(natData$AveMag),na.rm=TRUE)
                                     ,mean(stdize(natData$ETH_FRAC),na.rm=TRUE)
                                     ,mean(stdize(natData$FISC_DECEN),na.rm=TRUE)
                                     ,mean(stdize(natData$Leg_ENP),na.rm=TRUE)                                     
                                     )
                           ,nrow=1e4,ncol=n.beta,byrow=TRUE)
}



pred.expeta.max <- array(NA,c(1e4,3,100))
for(i in 1:100){
  for(j in 1:1e4){
    pred.expeta.max[j,,i] <- exp(c(predict.x.max[j,,i]%*%post.samples[j,1:n.beta],predict.x.max[j,,i]%*%post.samples[j,(n.beta+1):(n.beta*2)],predict.x.max[j,,i]%*%post.samples[j,((n.beta*2)+1):(n.beta*3)]))
  }
}
pred.shares.max <- apply(pred.expeta.max,3,apply,1,function(x){c(x[2]/sum(x),x[3]/sum(x))})
ps.max.target <- pred.shares.max[index%%2!=0,]
ps.max.target <- apply(ps.max.target,2,quantile,probs=probs)
ps.max.nontarget <- pred.shares.max[index%%2==0,]
ps.max.nontarget <- apply(ps.max.nontarget,2,quantile,probs=probs)
  

## Generate predicted values for mean homogeneity levels

for(i in 1:100){
  predict.x.mean[,,i] <- matrix(cbind(1
                                      ,huh.vector[i]
                                      ,mean(stdize(natData$Homogeneity),na.rm=TRUE)
                                      ,huh.vector[i]*mean(stdize(natData$Homogeneity),na.rm=TRUE)
                                      ,mean(stdize(natData$Age),na.rm=TRUE)
                                      ,mean(stdize(natData$Unemploy),na.rm=TRUE)
                                      ,mean(stdize(natData$AveMag),na.rm=TRUE)
                                      ,mean(stdize(natData$ETH_FRAC),na.rm=TRUE)
                                      ,mean(stdize(natData$FISC_DECEN),na.rm=TRUE)
                                      ,mean(stdize(natData$Leg_ENP),na.rm=TRUE)                                     
                                     )
                           ,nrow=1e4,ncol=n.beta,byrow=TRUE)
}
pred.expeta.mean <- array(NA,c(1e4,3,100))
for(i in 1:100){
  for(j in 1:1e4){
    pred.expeta.mean[j,,i] <- exp(c(predict.x.mean[j,,i]%*%post.samples[j,1:n.beta],predict.x.mean[j,,i]%*%post.samples[j,(n.beta+1):(n.beta*2)],predict.x.mean[j,,i]%*%post.samples[j,((n.beta*2)+1):(n.beta*3)]))
  }
}
pred.shares.mean <- apply(pred.expeta.mean,3,apply,1,function(x){c(x[2]/sum(x),x[3]/sum(x))})
ps.mean.target <- pred.shares.mean[index%%2!=0,]
ps.mean.target <- apply(ps.mean.target,2,quantile,probs=probs)
ps.mean.nontarget <- pred.shares.mean[index%%2==0,]
ps.mean.nontarget <- apply(ps.mean.nontarget,2,quantile,probs=probs)


## Generate predicted values for minimal homogeneity levels

for(i in 1:100){
  predict.x.min[,,i] <- matrix(cbind(1
                                     ,huh.vector[i]
                                     ,min(stdize(natData$Homogeneity),na.rm=TRUE)
                                     ,huh.vector[i]*min(stdize(natData$Homogeneity),na.rm=TRUE)
                                     ,mean(stdize(natData$Age),na.rm=TRUE)
                                     ,mean(stdize(natData$Unemploy),na.rm=TRUE)
                                     ,mean(stdize(natData$AveMag),na.rm=TRUE)
                                     ,mean(stdize(natData$ETH_FRAC),na.rm=TRUE)
                                      ,mean(stdize(natData$FISC_DECEN),na.rm=TRUE)
                                      ,mean(stdize(natData$Leg_ENP),na.rm=TRUE)                                     
                                     )
                           ,nrow=1e4,ncol=n.beta,byrow=TRUE)
}
pred.expeta.min <- array(NA,c(1e4,3,100))
for(i in 1:100){
  for(j in 1:1e4){
    pred.expeta.min[j,,i] <- exp(c(predict.x.min[j,,i]%*%post.samples[j,1:n.beta],predict.x.min[j,,i]%*%post.samples[j,(n.beta+1):(n.beta*2)],predict.x.min[j,,i]%*%post.samples[j,((n.beta*2)+1):(n.beta*3)]))
  }
}
pred.shares.min <- apply(pred.expeta.min,3,apply,1,function(x){c(x[2]/sum(x),x[3]/sum(x))})
ps.min.target <- pred.shares.min[index%%2!=0,]
ps.min.target <- apply(ps.min.target,2,quantile,probs=probs)
ps.min.nontarget <- pred.shares.min[index%%2==0,]
ps.min.nontarget <- apply(ps.min.nontarget,2,quantile,probs=probs)



#### Plot predicted values with credible bands

pdf("Graphs/PredictedShares.pdf",height=6,width=14)
par(mfrow=c(1,2))
plot(ps.max.target[2,19:80]~simNat,
     type='n',
     ylim=c(0,1),
     ylab="% of Targetable Expenditures",
     xlab="Government Nationalization",
     main="Targetable Expenditures",
     axes=FALSE)

polygon(c(simNat,rev(simNat)),
c(ps.max.target[3,19:80],rev(ps.max.target[1,19:80])),
col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

polygon(c(simNat,rev(simNat)),
c(ps.mean.target[3,19:80],rev(ps.mean.target[1,19:80])),
col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

polygon(c(simNat,rev(simNat)),
c(ps.min.target[3,19:80],rev(ps.min.target[1,19:80])),
col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)


axis(1)
axis(2)

lines(ps.max.target[2,19:80]~simNat,
     type='l',
     lwd=1.7,
     col="white")
lines(ps.mean.target[2,19:80]~simNat,
      type='l',
      lwd=1.7,
      lty=2,
      col="white")
lines(ps.min.target[2,19:80]~simNat,
      type='l',
      lwd=1.7,
      lty=3,
      col="white")

legend("topleft",
       bty="n",
       lty=c(1,2,3),
       lwd=c(2,2,2),
       col=c("black","black","black"),
       legend=c("High","Average","Low")
       ,title="Constituency Similarity")

plot(ps.max.nontarget[1,19:80]~simNat,
     type='n',
     ylim=c(0,1),
     ylab="% of Non-Targetable Expenditures",
     xlab="Government Nationalization",
     main="Non-Targetable Expenditures",
     axes=FALSE)
polygon(c(simNat,rev(simNat)),
c(ps.max.nontarget[3,19:80],rev(ps.max.nontarget[1,19:80])),
col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

polygon(c(simNat,rev(simNat)),
c(ps.min.nontarget[3,19:80],rev(ps.min.nontarget[1,19:80]))
,col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

polygon(c(simNat,rev(simNat)),
c(ps.mean.nontarget[3,19:80],rev(ps.mean.nontarget[1,19:80])),
col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

axis(1)
axis(2)
lines(ps.max.nontarget[2,19:80]~simNat,
     type='l',
     lwd=1.7,
     col="white")
lines(ps.mean.nontarget[2,19:80]~simNat,
      type='l',
      lwd=1.7,
      lty=2,
      col="white")
lines(ps.min.nontarget[2,19:80]~simNat,
      type='l',
      lwd=1.7,
      lty=3,
      col="white")
legend("bottomright",
       bty="n",
       lty=c(1,2,3),
       lwd=c(2,2,2),
       col=c("black","black","black"),
       legend=c("High","Average","Low")
       ,title="Constituency Similarity")
dev.off()




### Plot marginal effects of homogeneity

#sim.nat <- seq(min(natData$NatParB), max(natData$NatParB),length.out=100)
#theNames <- rownames(summary(results,quantiles=probs)$quantiles)
#TIndex <- grepl("mareffSim\\[.*,2\\]",theNames)
#NTIndex <-grepl("mareffSim\\[.*,3\\]",theNames)                     
#mareffT<-summary(results,quantiles=probs)$quantiles[TIndex,]
#mareffNT<-summary(results,quantiles=probs)$quantiles[NTIndex,]

#pdf("Graphs/MarEffSim(PartHomCabCoalNatPar80).pdf",height=6,width=6)
#plot(mareffNT[,2]~sim.nat,
#     type='n',
#     ylim=c(min(mareffNT[,1]),max(mareffNT[,3])),
#     ylab="Effect of Constituency Similarity on Non-Targetable Expenditures",
#     xlab="Government Nationalization",
#     main="",
#     axes=FALSE)
#polygon(c(sim.nat,rev(sim.nat)),
#c(mareffNT[,3],rev(mareffNT[,1])),
#col=rgb(red=0,green=0,blue=0,alpha=150,max=255),border=NA)

#abline(h=0, lty=2)

#axis(1)
#axis(2)
#lines(mareffNT[,2]~sim.nat,
#     type='l',
#     lwd=1.7,
#     col="white")
#dev.off()
