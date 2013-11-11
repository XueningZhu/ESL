mm=matrix(c(0,0.9,0.8,0,0,0.7,0,0,0),ncol=3,byrow=T)
mat=mm+t(mm)+diag(1,3)
X=mvrnorm(n = 1000, mu=c(1,2,4), Sigma=mat)
data=X

### missing at random
num=1:1000
missing_ind=NULL
for (i in 1:3)
{
  missing_ind=rbind(missing_ind,cbind(sample(num,100),i))
  num=setdiff(num,missing_ind[i,])
}
data[missing_ind]=NA


#http://stackoverflow.com/questions/18395368/sampling-from-an-inverse-gamma-distribution-in-r
MI.mcmc<-function(data)
{
  v0=1; s0_2=1
  a0=v0/2; b0=1/2*v0*s0_2
  mu0=matrix(0,nrow=ncol(data)-1);Gamma0=diag(1,ncol(data)-1)
  data_com=data[complete.cases(data),]
  n=nrow(data_com)
  var_na=apply(data,2,function(x) {
    return(which(is.na(x)))
  })
  
  mu=rep(list(mu=mu0),ncol(data_com))
  Gamma=rep(list(Gamma=Gamma0),ncol(data))
  a=rep(list(a=a0),ncol(data_com))
  b=rep(list(b=b0),ncol(data_com))
  y=rep(list(NULL ),ncol(data_com))
  iter=1
  while(iter<=100)
  {
    show(iter)
    for (i in 1:ncol(data))
    {
      mu0=matrix(mu[[i]][,ncol(mu[[i]])],ncol=1)
      Gamma0=Gamma[[i]][(nrow(Gamma[[i]])-ncol(data)+2):nrow(Gamma[[i]]),]
      a0=a[[i]][length(a[[i]])]
      b0=b[[i]][length(b[[i]])]
      
      mu1=solve(t(data_com[,-i])%*%data_com[,-i]+Gamma0)%*%(Gamma0%*%mu0+t(data_com[,-i])%*%matrix(data_com[,i],ncol=1))
      Gamma1=(t(data_com[,-i])%*%data_com[,-i]+Gamma0)
      a1=a0+n/2
      b1=b0+1/2*(sum(data_com[,i]*data_com[,i])+t(mu0)%*%Gamma0%*%mu0-t(mu0)%*%Gamma0%*%mu0)
      
      sigma1 = sqrt(1/rgamma(1,a1,b1))
      beta1=mvrnorm(n=1, mu1, sigma1^2*solve(Gamma1))
      y1=data[var_na[,i],-i]%*%beta1+rnorm(length(var_na[,i]),mean=0,sd=sigma1)
      data[var_na[,i],i]=y1
      
      mu[[i]]=cbind(mu[[i]],mu1)
      Gamma[[i]]=rbind(Gamma[[i]],Gamma1)
      a[[i]]=c(a[[i]],a1)
      b[[i]]=c(b[[i]],b1)
      y[[i]]=cbind(y[[i]],y1)
    }
    iter=iter+1
  }
  return(list(Data=data,Y=y,Mu=mu,Gamma=Gamma,A=a,B=b))
}
res=MI.mcmc(data)

data1=data
data1[missing_ind]=do.call(c,lapply(res$Y,function(x)
{
  zz=apply(x[,(ncol(x)-10):ncol(x)],1,mean,na.rm=T)
  return(zz)
}))

library(ggplot2)
XX=data.frame(X,type="origin")
dd=data.frame(data,type="missing")

m <- ggplot(XX,  aes(x=X3))  ### X2 / X3 is also ok
m <- m + geom_histogram(data=XX,color="red",fill="red",alpha=0.3, binwidth = 0.5)
m <- m + geom_histogram(data=dd,color="black",fill="black",alpha=0.1, binwidth = 0.5)
m + geom_histogram(data=dd1[is.na(dd$X1),],color="blue",fill="blue",alpha=0.3, binwidth = 0.5)

