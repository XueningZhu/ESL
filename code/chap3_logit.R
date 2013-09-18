
require(fields) 

hmap.func <- function(a, f, xlab, ylab) 
{
  image.plot(a, a, outer(a, a, f), zlim = c(0, 1), xlab = xlab, ylab = ylab)
}

g <- function(x, y) 1 / (1 + 2^(x^3+ y+x*y))

d <- data.frame(x1 = rnorm(1000), x2 = rnorm(1000)
                ,x3 = rnorm(1000), x4 = rnorm(1000))
d$y = with(d, ifelse(runif(1000) < g(x1, x2), 1, 0))


# plot signal (left hand plot below)
a = seq(-2, 2, len = 100)
hmap.func(a, g, "x1", "x2")


## logistic regression: binary classification

get.p<-function(beta,X)
{
  xx=exp(X%*%matrix(beta,ncol=1))
  p=xx/(1+xx)
  return(p)
}
logit.beta<-function(beta0,X,Y,eps)
{
  P=get.p(beta0,X)
  W=diag(as.vector(P*(1-P)))
  beta1=beta0+solve(t(X)%*%W%*%X)%*%t(X)%*%(Y-P)
  
  while(sum((beta1-beta0)^2)/sum(beta1^2)>eps^2)
  {
    cat("beta1:",beta1,"\n")
    beta0=beta1
    P=get.p(beta0,X)
    W=diag(as.vector(P*(1-P)))
    beta1=beta0+solve(t(X)%*%W%*%X)%*%t(X)%*%(Y-P)
  }
  return(beta1)
}

logit.predict<-function(X,Y,eps=10^(-1))
{
  X=as.matrix(cbind(1,X))
  Y=matrix(Y,ncol=1)
  beta0=sample(c(-1, 1), ncol(X), replace=TRUE)
  beta=logit.beta(beta0,X,Y,eps=eps)
  pred=get.p(beta,X)
  return(list(Beta=beta,Pred=pred))
}

logit.res=logit.predict(d[,-1],d$y,eps=10^(-2))

f<-function(x,y,beta)
{
  zero=matrix(0,nrow=length(x),ncol=nrow(beta)-3)
  X=as.matrix(cbind(1,x,y,zero))
  p=get.p(beta,X)
  return(p)
}
image.plot(a, a, outer(a, a, f,beta=logit.res$Beta), zlim = c(0, 1), xlab = "x1", ylab = "x2")

aa=outer(a, a, f,beta=logit.res$Beta)