
#### for chap 4 LDA & QDA (Linear Discriminant Analysis)

### generate the dataset: different mean/ same covariance
x1=rnorm(1000)
x2=rnorm(1000)

d1=data.frame(x=x1+x2,y=x1-x2,group=1)
d2=data.frame(x=2+x1+x2,y=2+x1-x2,group=2)
#d2=data.frame(x=8+x1-2*x2,y=-1+2*x1+x2,group=3)
data=rbind(d1,d2)

### estimation with no bias
pai=summary(factor(data$group))/nrow(data)
l_data=split(data,factor(data$group))
mu=sapply(l_data,function(x) return(apply(x,2,mean)))
l_mat=lapply(1:length(l_data),function(i) {y=t(as.matrix(l_data[[i]][,1:2]))-mu[1:2,i]; return(y%*%t(y))})
Sigma=do.call("+",l_mat)/(nrow(data)-length(pai))

### Linear discriminant function (covariance is identical)
LFun <- function(x,y,mu,sigma,pai)
{
  mu=matrix(mu,ncol=1)
  return(cbind(x,y)%*%solve(sigma)%*%mu-as.numeric(t(mu)%*%solve(sigma)%*%mu*1/2)+log(pai))
}
LDA <- function(x,y,mu,sigma,pai)
{
  z=sapply(1:length(x),function(i) {return(max(LFun(x[i],y[i],mu[1:2,1],sigma,pai[1]),
                                             LFun(x[i],y[i],mu[1:2,2],sigma,pai[2])))})
  return(z)
}

##color: http://research.stowers-institute.org/efg/R/Color/Chart/ColorChart.pdf
LDA.col <- function(x,y,mu,sigma,pai)
{
  z=sapply(1:length(x),function(i) {return(c("#FF7256","#BCEE68")[which.max(c(LFun(x[i],y[i],mu[1:2,1],sigma,pai[1]),
                                               LFun(x[i],y[i],mu[1:2,2],sigma,pai[2])))])})
  return(z)
}

### for the 3d plot (discriminant hyperplane)
library(rgl)
a = 5
x = seq(-a, a, 0.05)
y = seq(-a, a, 0.05)
z = outer(x,y,LDA, mu, sigma=Sigma,pai=pai)
z.col = as.character(matrix(outer(x,y,LDA.col, mu, sigma=Sigma,pai=pai),nrow=1))
persp3d(x, y, z,col=z.col)

### the LDA Classifiers
  ### solve the implicit function for y: uniroot is employed here
f <- function (x,sigma,mu,pai) 
{
  g <- function(y) t(rbind(x,y))%*%solve(sigma)%*%(mu[1:2,2]-mu[1:2,1])-
    1/2*t(mu[1:2,2])%*%solve(sigma)%*%mu[1:2,2]+
    1/2*t(mu[1:2,1])%*%solve(sigma)%*%mu[1:2,1]-
    log(pai[1])+log(pai[2])
  y=uniroot ( g , c ( -100 , 100 ) )
  return(y)
}

l_y=sapply(data$x,function(x)
  {
  return(f(x,sigma=Sigma,mu,pai)$root)
})

data$l_y=l_y

### for the LDA visualization: 2d
library(ggplot2)
p <- ggplot(data, aes(x, y))
p <- p + geom_point(aes(colour=factor(group)))
p + geom_line(aes(x=x,y=l_y),col="red",size=1.5)

### ==============================================================
### QDA: Quadratic Discriminant Analysis
### generate the dataset: different mean & covariance
x1=rnorm(1000)
x2=rnorm(1000)

d1=data.frame(x=x1+x2,y=x1-x2,group=1)
d2=data.frame(x=3+x1+x2,y=4+x1-x2,group=2)
d3=data.frame(x=8+x1-2*x2,y=-1+2*x1+x2,group=3)
data=rbind(d1,d2,d3)
p <- ggplot(data, aes(x, y))
p + geom_point(aes(colour=factor(group)))

### estimation with no bias (some test for mu and sigma further)
pai=summary(factor(data$group))/nrow(data)
l_data=split(data,factor(data$group))
mu=sapply(l_data,function(x) return(apply(x,2,mean)))
l_mat=lapply(1:length(l_data),function(i) {y=t(as.matrix(l_data[[i]][,1:2]))-mu[1:2,i]; return(y%*%t(y))})
Sigma=Reduce("+",l_mat)/(nrow(data)-length(pai))

### ==========
l_y12=sapply(data$x,function(x)
{
  return(f(x,sigma=Sigma,mu[,1:2],pai[1:2])$root)
})

data$l_y12=l_y12

l_y23=sapply(data$x,function(x)
{
  return(f(x,sigma=Sigma,mu[,2:3],pai[2:3])$root)
})

data$l_y23 = l_y23

l_y13=sapply(data$x,function(x)
{
  return(f(x,sigma=Sigma,mu[,c(1,3)],pai[c(1,3)])$root)
})

data$l_y13 = l_y13


### ====================

### group1 & group3: counter intuitive (actually low accuracy)
p <- ggplot(data, aes(x, y))
p <- p + geom_point(aes(colour=factor(group)))
p <- p + geom_line(aes(x=x,y=l_y12),col="red",size=1.5)
p <- p + geom_line(aes(x=x,y=l_y23),col="black",size=1.5)
p <- p + geom_line(aes(x=x,y=l_y13),col="blue",size=1.5)
p

### QDA
Sigma_d=lapply(1:length(l_mat),function(i) return(l_mat[[i]]/(nrow(l_data[[i]])-1)))

QFun <- function(x,y,mu,sigma,pai)
{
  res1=sapply(1:length(x),function(i) {
    X=rbind(x[i],y[i])-mu
    res=-1/2*log(abs(det(sigma)))-1/2*t(X)%*%solve(sigma)%*%X+log(pai)
  })
  return(res1)
}
QDA <- function(x,y,mu,sigma,pai)
{
  N=length(pai)
  z=sapply(1:length(x),function(i) {
    res=sapply(1:N,function(j) {return(QFun(x[i],y[i],mu[1:2,j],sigma[[j]],pai[j]))})
    return(max(res))
  })
  return(z)
}

##color: http://research.stowers-institute.org/efg/R/Color/Chart/ColorChart.pdf
QDA.col <- function(x,y,mu,sigma,pai)
{
  N=length(pai)
  z=sapply(1:length(x),function(i) 
    {return(c("#FF7256","#BCEE68","#6495ED" )[which.max(sapply(1:N,function(j) 
      {return(QFun(x[i],y[i],mu[1:2,j],sigma[[j]],pai[j]))}))])})
  return(z)
}

### for the 3d plot (discriminant hyperplane)
library(rgl)
x = seq(-5, 7, 0.1)
y = seq(-5, 7, 0.1)


z = outer(x,y,QDA, mu, sigma=Sigma_d,pai=pai)
z.col = as.character(matrix(outer(x,y,QDA.col, mu, sigma=Sigma_d,pai=pai),nrow=1))
persp3d(x, y, z,col=z.col)

z2 = outer(x,y,QFun, mu=mu[1:2,2], sigma=Sigma_d[[2]],pai=pai[2])
QFun(1,3,mu=mu[1:2,1], sigma=Sigma_d[[1]],pai=pai[1])
persp3d(x, y, z2,col=z.col)

f_QDA <- function (x,sigma,mu,pai) 
{
  N=length(pai)
  comb = combn(N,2)
  res=apply(comb,2,
        function(com)
          {
          tryCatch({
            g <- function(y) {
              QFun(x,y,mu[1:2,com[1]],sigma[[com[1]]],pai[com[1]])-
                QFun(x,y,mu[1:2,com[2]],sigma[[com[2]]],pai[com[2]])
            }
            if(all(com==c(1,3)))
              interval=c(-5,0)
            else
              interval=c(0,5)
            y=uniroot ( g , interval )$root
            
          }, error = function(err) {
            # warning handler picks up where error was generated
            print(paste("MY_ERROR:  ",err))
            y=NA
            return(y)
          }, finally = {
            print(y)
          }) # END tryCatch
          
        })
  return(res)
}
#res=f_QDA(1,Sigma_d,mu,pai)

Q_y=t(sapply(data$x,function(x)
{
  return(f_QDA(x,sigma=Sigma_d,mu,pai))
}))
colnames(Q_y)=paste("Q_y",apply(combn(3,2),2,paste,collapse=""),sep="")
data1=data.frame(data,Q_y)

p <- ggplot(data1, aes(x, y))
p <- p + geom_point(aes(colour=factor(group)))
p <- p + geom_line(aes(x=x,y=Q_y12),col="red",size=1.5)
p <- p + geom_line(aes(x=x,y=Q_y23),col="black",size=1.5)
p <- p + geom_line(aes(x=x,y=Q_y13),col="blue",size=1.5)
p








