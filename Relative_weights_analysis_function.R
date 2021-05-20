rwa <- function(thedata,iter = 5000){
  
multRegress<-function(mydata){
  numVar<<-NCOL(mydata)
  Variables<<- names(mydata)[2:numVar]
  
  mydata<-cor(mydata, use="complete.obs")
  RXX<-mydata[2:numVar,2:numVar]
  RXY<-mydata[2:numVar,1]
  
  RXX.eigen<-eigen(RXX)
  D<-diag(RXX.eigen$val)
  delta<-sqrt(D)
  
  lambda<-RXX.eigen$vec%*%delta%*%t(RXX.eigen$vec)
  lambdasq<-lambda^2
  beta<-solve(lambda)%*%RXY
  rsquare<<-sum(beta^2)
  
  RawWgt<-lambdasq%*%beta^2
  import<-(RawWgt/rsquare)*100
  
  result<<-data.frame(Variables, Raw.RelWeight=RawWgt, Rescaled.RelWeight=import)
}


multBootstrap<-function(mydata, indices){
  mydata<-mydata[indices,]
  multWeights<-multRegress(mydata)
  return(multWeights$Raw.RelWeight)
}

multBootrand<-function(mydata, indices){
  mydata<-mydata[indices,]
  multRWeights<-multRegress(mydata)
  multReps<-multRWeights$Raw.RelWeight
  randWeight<-multReps[length(multReps)]
  randStat<-multReps[-(length(multReps))]-randWeight
  return(randStat)
}

multBootcomp<-function(mydata, indices){
  mydata<-mydata[indices,]
  multCWeights<-multRegress(mydata)
  multCeps<-multCWeights$Raw.RelWeight
  comp2Stat<-multCeps-multCeps[1]
  comp2Stat<-comp2Stat[-1]
  Labels2<<-Labels[-1]
  return(comp2Stat)
}

#bootstrapping
# install.packages("boot")
library(boot)

mybootci<-function(x){
  boot.ci(multBoot,conf=0.95, type="bca", index=x)
}

runBoot<-function(num){
  INDEX<-1:num
  test<-lapply(INDEX, FUN=mybootci)
  test2<-t(sapply(test,'[[',i=4)) #extracts confidence interval
  CIresult<<-data.frame(Variables, CI.Lower.Bound=test2[,4],CI.Upper.Bound=test2[,5])
}
myRbootci<-function(x){
  boot.ci(multRBoot,conf=0.95,type="bca",index=x)
}

runRBoot<-function(num){
  INDEX<-1:num
  test<-lapply(INDEX,FUN=myRbootci)
  test2<-t(sapply(test,'[[',i=4))
  CIresult<<-data.frame(Labels,CI.Lower.Bound=test2[,4],CI.Upper.Bound=test2[,5])
}

myCbootci<-function(x){
  boot.ci(multC2Boot,conf=0.95,type="bca",index=x)
}

runCBoot<-function(num){
  INDEX<-1:num
  test<-lapply(INDEX,FUN=myCbootci)
  test2<-t(sapply(test,'[[',i=4))
  CIresult<<-data.frame(Labels2,CI.Lower.Bound=test2[,4],CI.Upper.Bound=test2[,5])
}

myGbootci<-function(x){
  boot.ci(groupBoot,conf=0.95,type="bca",index=x)
}

runGBoot<-function(num){
  INDEX<-1:num
  test<-lapply(INDEX,FUN=myGbootci)
  test2<-t(sapply(test,'[[',i=4))
  CIresult<<-data.frame(Labels,CI.Lower.Bound=test2[,4],CI.Upper.Bound=test2[,5])
}


Labels<-names(thedata)[2:length(thedata)]
multRegress(thedata)
RW.Results<-result #The Raw and Rescaled Weights
RSQ.Results<-rsquare #R-squared For the Model
#Bootstrapped Confidence interval around the individual relative weights, 5000 bootstrapped resamples
multBoot<-boot(thedata, multBootstrap, iter)
multci<-boot.ci(multBoot,conf=0.95, type="bca")
runBoot(length(thedata[,2:numVar]))
CI.Results<-CIresult  #BCa Confidence Intervals around the raw weights
#Bootstrapped Confidence interval tests of Significance, 5000 bootstrapped resamples
randVar<-rnorm(length(thedata[,1]),0,1)
randData<-cbind(thedata,randVar)
multRBoot<-boot(randData,multBootrand, iter)
multRci<-boot.ci(multRBoot,conf=0.95, type="bca")
runRBoot(length(randData[,2:(numVar-1)]))
CI.Significance<-CIresult #BCa Confidence Interval Tests of significance; If Zero is not included, Weight is Significant

return(
  list(
    #R-squared For the Model
RSQ.Results,
    #The Raw and Rescaled Weights
RW.Results,
    #BCa Confidence Intervals around the raw weights
CI.Results,
    #BCa Confidence Interval Tests of significance
    #If Zero is not included, Weight is Significant
CI.Significance)
  )
# )
}

