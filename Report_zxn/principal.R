

setwd("F:/商业数据挖掘_光华/homework/2")

#library(sas7bdat)
############ read data files
#zip1=read.sas7bdat("zipdemo1.sas7bdat")
#save(zip1,file="zip1.rda")

##### load the data files
load("zip1.rda")
summary(zip1)
sum(!complete.cases(zip1))

##### conduct the multiple imputation
#require(mi)
#zip1_com = mi(zip1[,-1])
#comp_zip1 = mi.data.frame(zip1_com,m=1) 
#zip1_comp = cbind(zip1[,1],comp_zip1)
#save(zip1_comp,file="zip1_comp.rda")
load("zip1_comp.rda")
summary(zip1_comp)

#### Principal Component Analysis
dd=zip1_comp[,-1]
zip_prin=princomp(dd,cor=T,scores=T)
screeplot(zip_prin,type="lines")
summary(zip_prin)
loadings(zip_prin)

#### Factor Analysis
dd_std=scale(dd)
zip_fac=factanal(dd_std,factor=19,rotation="varimax",n.obs=nrow(dd),control=list(trace=T))
zip_fac


zip_fac$STATISTIC

zip_fac$PVAL

zip_fac$criteria

loadings(zip_fac)

apply(zip_fac$loadings,2,function(x) return(sum(x^2)))

zip_fac=factanal(dd_std,factor=19,rotation="varimax",n.obs=nrow(dd),control=list(trace=T))

flag=1; i=15
while(flag==1)
{
  zip_fac=factanal(dd_std,factor=i,rotation="varimax",n.obs=nrow(dd))
  fac_var=apply(zip_fac$loadings,2,function(x) return(sum(x^2)))
  cat("nFactor:",i,"\n","Factor variance:",fac_var,"\n","\n")
  if (all(fac_var>=1))
    flag=0
  i=i-1
}

# Determine Number of Factors to Extract
library(nFactors)
ev <- eigen(cor(dd_std)) # get eigenvalues
ap <- parallel(subject=nrow(dd_std),var=ncol(dd_std),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

nBartlett(x=ev$values, N= nrow(dd_std), alpha=0.05, details=TRUE)
nBentler(x=ev$values, N= nrow(dd_std), alpha=0.05, details=TRUE)
nCng(x=ev$values,  model="factors", details=TRUE)

dd_std=scale(dd)
zip_fac=factanal(dd_std,factor=6,rotation="varimax",n.obs=nrow(dd))
loadings(zip_fac)
zip_fac=factanal(dd_std,factor=6,rotation="promax",n.obs=nrow(dd))
loadings(zip_fac)


#=========================================================================
### other packages and methods.

library(psych)
library(GPArotation)
fit <- fa(dd_std, nfactors=7, rotate="varimax" ,max.iter = 30,fm="pa")
factor.stats(dd_std,fit,n.obs=nrow(dd))
fit # print results

#using the Harman 24 mental tests, compare a principal factor with a principal components solution
pc <- principal(Harman74.cor$cov,4,rotate="varimax")
pa <- fa(Harman74.cor$cov,4,fm="pa" ,rotate="varimax")  #principal axis 
uls <- fa(Harman74.cor$cov,4,rotate="varimax")          #unweighted least squares is minres
wls <- fa(Harman74.cor$cov,4,fm="wls")       #weighted least squares

#to show the loadings sorted by absolute value
print(uls,sort=TRUE)

#then compare with a maximum likelihood solution using factanal
mle <- factanal(covmat=Harman74.cor$cov,factors=4)
factor.congruence(list(mle,pa,pc,uls,wls))
#note that the order of factors and the sign of some of factors may differ 

#finally, compare the unrotated factor, ml, uls, and  wls solutions
wls <- fa(Harman74.cor$cov,4,rotate="none",fm="wls")
pa <- fa(Harman74.cor$cov,4,rotate="none",fm="pa")
minres <-  factanal(factors=4,covmat=Harman74.cor$cov,rotation="none")
mle <- fa(Harman74.cor$cov,4,rotate="none",fm="mle")
uls <- fa(Harman74.cor$cov,4,rotate="none",fm="uls")
factor.congruence(list(minres,mle,pa,wls,uls))
#in particular, note the similarity of the mle and min res solutions
#note that the order of factors and the sign of some of factors may differ 



#an example of where the ML and PA and MR models differ is found in Thurstone.33.
#compare the first two factors with the 3 factor solution 
Thurstone.33 <- as.matrix(Thurstone.33)
mle2 <- fa(Thurstone.33,2,rotate="none",fm="mle")
mle3 <- fa(Thurstone.33,3 ,rotate="none",fm="mle")
pa2 <- fa(Thurstone.33,2,rotate="none",fm="pa")
pa3 <- fa(Thurstone.33,3,rotate="none",fm="pa")
mr2 <- fa(Thurstone.33,2,rotate="none")
mr3 <- fa(Thurstone.33,3,rotate="none")
factor.congruence(list(mle2,mr2,pa2,mle3,pa3,mr3))

v9 <- sim.hierarchical()
f3 <- fa(v9,3)
factor.stats(v9,f3,n.obs=500)
f3o <- fa(v9,3,fm="pa",rotate="Promax")
factor.stats(v9,f3o,n.obs=500)




