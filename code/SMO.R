require(fields) # for heatmap plot with legend

# heatmap wrapper, plotting func(x, y) over range a by a
hmap.func <- function(a, f, xlab, ylab) 
{
  image.plot(a, a, outer(a, a, f), zlim = c(0, 1), xlab = xlab, ylab = ylab)
}

# define class signal
f <- function(x, y) 1 / (1 + 2^(x^3+ y+x*y))

predict.SMO<-function(x,y,alpha,data,Y,b,ker=Gauss.ker)
{
  X=c(x,y,0,0)
  kk=apply(data,1,ker,z=X)
  res=sum(Y*alpha*kk)+b
  return(res)
}

predict.SMO(x=a[1],y=a[3],alpha=res$Alpha,data=d[1:100,1:4],Y=d$y[1:100],b=res$B,ker=Gauss.Ker)

aa=sapply(a,function(x)
  {
  rr=sapply(a,predict.SMO,y=x,alpha=res$Alpha,
         data=d[1:100,1:4],Y=d$y[1:100],b=res$B,ker=Gauss.Ker)
  return(rr)
})

image.plot(a, a, aa, zlim = c(min(aa), max(aa)))

bb=outer(a, a, f)
logit_aa=exp(aa)/(1+exp(aa))
image.plot(a, a, bb, zlim = c(min(bb), max(bb)))
image.plot(a, a, logit_aa, zlim = c(min(logit_aa), max(logit_aa)))
#g <- function(x, y) 1 / (1 + 2^(cos(x^3)-x*y+sin(y)))

# create training data
d <- data.frame(x1 = rnorm(1000), x2 = rnorm(1000)
                ,x3 = rnorm(1000), x4 = rnorm(1000))
d$y = with(d, ifelse(runif(1000) < f(x1, x2), 1, -1))

# plot signal (left hand plot below)
a = seq(-2, 2, len = 100)
hmap.func(a, f, "x1", "x2")

####### the svm classifier
Gauss.Ker<-function(x,z,sigma=1) ### x z is vectors
{
  res=exp(-sum((x-z)^2)/(2*sigma^2))
  return(res)
}
Ker.mat<-function(Ker,data)
{
  res=apply(data,1,function(x)
    {
    return(apply(data,1,Ker,z=x))
  })
#  rr=do.call(cbind,res)
  return(res)
}

#aa=Ker.mat(Gauss.Ker,data=d[,1:4])


is.support<-function(alpha,C,eps)
{
  is_support=sapply(alpha,function(x) {
    return(x>eps&x<C-eps)
  })
  return(is_support)
}
is.KKT<-function(y,alpha,eps,C,gg)
{
  ss=sum(y*alpha)
  if (abs(ss)>eps) return(F)
  if (any(-eps>alpha|alpha>C+eps)) return(F)
  yg=y*gg
  
  is_KKT=sapply(1:length(alpha),function(i) {
    if (abs(alpha[i])<=eps&yg[i]>=1) return(T)
    if (alpha[i]>0&alpha[i]<C&abs(yg[i]-1)<eps) return(T)
    if (abs(alpha[i]-C)<eps&yg[i]<=1) return(T)
    return(F)
  })
  return(is_KKT)
}

first.alpha<-function(y,alpha,eps,C,gg,is_support) ## C is the cost parameter
{
  yg=y*gg
  
  if(any(is_support)) 
  {
    aa=alpha[is_support]
    ygs=yg[is_support]
    is_KKT=sapply(1:length(aa),function(i) {
      #   if (alpha[i]==0&yg[i]>=1) return(T)
      if (aa[i]>0&aa[i]<C&abs(ygs[i]-1)<eps) return(T)
      return(F)
      #  if (abs(alpha[i]-C)<eps&yg[i]<=1) return(T)
    })
    if (!all(is_KKT))
    {
      return(sample(which(is_support)[!is_KKT],1))
    }
  }
  
  inds=!is_support
  aa=alpha[inds]
  isyg=yg[inds]
  
  is_KKT=sapply(1:length(aa),function(i) {
    if (abs(aa[i])<=eps&isyg[i]>=1) return(T)
    #  if (alpha[i]>0&alpha[i]<C&abs(yg[i]-1)<eps) return(T)
    if (abs(aa[i]-C)<eps&isyg[i]<=1) return(T)
    return(F)
  })
  return(sample(which(inds)[!is_KKT],1))
}

second.alpha<-function(E,alpha1_ind)
{
  inds=(1:length(E))[-alpha1_ind]
  ee=E[-alpha1_ind]
  if (E[alpha1_ind]>0) 
  {
    mm=min(ee)
    alpha2_ind=inds[sample(which(ee==mm),1)]
  }
  else
  {
    mm=max(ee)
    alpha2_ind=inds[sample(which(ee==mm),1)]
  }
  return(alpha2_ind)
}

Update.alphas<-function(alpha1_ind,alpha2_ind,y,alpha,C,E,ker_mat)
{
  yy=y[c(alpha1_ind,alpha2_ind)]
  ee=E[c(alpha1_ind,alpha2_ind)]
  if (yy[1]!=yy[2])
  {
    L=max(0,alpha[alpha2_ind]-alpha[alpha1_ind])
    H=min(C,C+alpha[alpha2_ind]-alpha[alpha1_ind])
  }
  else
  {
    L=max(0,alpha[alpha2_ind]+alpha[alpha1_ind]-C)
    H=min(C,alpha[alpha2_ind]+alpha[alpha1_ind])
  }
  kk=ker_mat[alpha1_ind,alpha1_ind]+ker_mat[alpha2_ind,alpha2_ind]-
    2*ker_mat[alpha1_ind,alpha2_ind]
  a2=alpha[alpha2_ind]+yy[2]*(ee[1]-ee[2])/kk
  
  alpha2=a2
  if (a2>H)
  {
    alpha2=H
  }
  if (a2<L)
  {
    alpha2=L
  }
  alpha1=alpha[alpha1_ind]+yy[1]*yy[2]*(alpha[alpha2_ind]-alpha2)
  return(c(alpha1,alpha2))
}
update.b<-function(E,y,b,ker_mat,alpha1_ind,alpha2_ind,alpha_old,alpha_new,eps,C)
{
  ee=E[c(alpha1_ind,alpha2_ind)]
  yy=y[c(alpha1_ind,alpha2_ind)]
  k11=ker_mat[alpha1_ind,alpha1_ind]
  k22=ker_mat[alpha2_ind,alpha2_ind]
  k12=ker_mat[alpha1_ind,alpha2_ind]
  
  b1=-ee[1]-yy[1]*k11*(alpha_new[1]-alpha_old[1])-yy[2]*k12*(alpha_new[2]-alpha_old[2])+b
  b2=-ee[2]-yy[1]*k12*(alpha_new[1]-alpha_old[1])-yy[2]*k22*(alpha_new[2]-alpha_old[2])+b
  if(all(abs(alpha_new)<eps)|all(abs(alpha_new-C)<eps))
  {
    b=mean(c(b1,b2))
  }
  else
  {
    ss=sapply(1:2,function(i) return(alpha_new[i]>0&alpha_new[i]<C))
    
    if (!all(ss))
    {
      return(b)
    }
    b=c(b1,b2)[ss][1]
    if (all(ss))
      cat("b1-b2:",b1-b2,"\n")
  }
  return(b)
}

g<-function(alpha,y,ker_mat,col,b) ### col is the x colind
{
  res=sum(alpha*y*ker_mat[,col])+b
  return(res)
}


update.E<-function(E,is_support,y,alpha,ker_mat,b)
{
  gg=sapply(1:length(alpha),g,alpha=alpha,y=y,
            ker_mat=ker_mat,b=b)
  E=gg-y
  return(E)
}


SMO<-function(data,y,C,eps=10^(-4)) ###data is the features; y is the respond variable
{
  alpha=rep(0,nrow(data))
  b=0
  ker_mat=Ker.mat(Gauss.Ker,data=data)
  gg=sapply(1:length(alpha),g,alpha=alpha,y=y,
            ker_mat=ker_mat,b=b)
  E=gg-y
  is_kkt=is.KKT(y,alpha,eps,C,gg) 
  iter=1
  while(!all(is_kkt))
  {
    iter=iter+1
    is_support=is.support(alpha,C,eps)
    E=update.E(E,is_support,y,alpha,ker_mat,b)
    
    alpha1_ind=first.alpha(y,alpha,eps,C,gg,is_support)
 #   nkkt_ind=which(!is_kkt)
    alpha2_ind=second.alpha(E,alpha1_ind)
    
    alpha_new=Update.alphas(alpha1_ind,alpha2_ind,y,alpha,C,E,ker_mat)
    b=update.b(E,y,b,ker_mat,alpha1_ind,alpha2_ind,
             alpha_old=alpha[c(alpha1_ind,alpha2_ind)],alpha_new,eps,C)
    
    alpha[c(alpha1_ind,alpha2_ind)]=alpha_new
    gg=sapply(1:length(alpha),g,alpha=alpha,y=y,
              ker_mat=ker_mat,b=b)

    is_kkt=is.KKT(y,alpha,eps,C,gg)
    cat("un kkt",sum(!is_kkt),"\n")
  }
  return(list(Alpha=alpha,B=b,Fit=gg,Ker_mat=ker_mat))
}
res=SMO(data=d[1:100,1:4],y=d$y[1:100],C=0.5,eps=10^(-4))










