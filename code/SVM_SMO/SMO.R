require(fields) # for heatmap plot with legend

# heatmap wrapper, plotting func(x, y) over range a by a
hmap.func <- function(a, f, xlab, ylab) 
{
  image.plot(a, a, outer(a, a, f), zlim = c(0, 1), xlab = xlab, ylab = ylab)
}

# define class signal
f <- function(x, y) 1 / (1 + 2^(x^3+ y+x*y))

# create training data
d <- data.frame(x1 = rnorm(200), x2 = rnorm(200)
                ,x3 = rnorm(200), x4 = rnorm(200))
d$y = with(d, ifelse(runif(200) < f(x1, x2), 1, -1))

# plot signal (left hand plot below)

hmap.func(seq(-2, 2, len = 50), f, "x1", "x2")


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
tt=Ker.mat(Gauss.Ker,data=d[,1:4])


is.support<-function(alpha,C,eps)
{
  is_support=sapply(alpha,function(x) {
    return(x>eps&x<C-eps)
  })
  return(which(is_support))
}
is.KKT<-function(y,alpha,eps,C,gg)
{
  ss=sum(y*alpha)
  if (abs(ss)>eps) return(F)
  if (any(-eps>alpha|alpha>C+eps)) return(F)
  yg=y*gg
  
  kkt_stats=rep(0,length(alpha))
  is_KKT=sapply(1:length(alpha),function(i) {
      if (abs(alpha[i])<=eps) { if(yg[i]>=1) return(c(T,kkt_stats[i])); kkt_stats[i]=1}
      if (alpha[i]>eps&alpha[i]<C-eps) {if(abs(yg[i]-1)<eps) return(c(T,kkt_stats[i])); kkt_stats[i]=2}
      if (abs(alpha[i]-C)<eps) {if(yg[i]<=1) return(c(T,kkt_stats[i])); kkt_stats[i]=3}
      return(c(F,kkt_stats[i]))
 #  if (abs(alpha[i])<=eps&yg[i]>=1) { return(T)}
 #  if (alpha[i]>0&alpha[i]<C&abs(yg[i]-1)<eps) return(T)
 #  if (abs(alpha[i]-C)<eps&yg[i]<=1) return(T)
 #   return(F)
      })
  cat(names(table(is_KKT[2,])),"||| kkt_stats:",table(is_KKT[2,]),"\n")
  cat("un kkt",sum(!is_KKT[1,]),"\n")
  return(is_KKT[1,])
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

Update<-function(y,alpha,C,E,ker_mat,eps,b)
{
#  flag=T
  alpha_changed=1
  examineAll=T
  
  while(examineAll|alpha_changed>0)
  {
    if (examineAll)
    {
      alpha1_inds=1:length(alpha)
    }
    else
    {
      alpha1_inds=is.support(alpha,C,eps)
    }
    alpha_changed=0
    ######################################
    for (i in 1:length(alpha1_inds))
    {
      alpha1_ind=alpha1_inds[i]
      alpha2_ind=second.alpha(E,alpha1_ind)
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
      if (abs(ee[1]-ee[2])<eps|abs(L-H)<eps)
      {
        next
      }
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
      
      alpha_new=c(alpha1,alpha2)
      b=update.b(E,y,b,ker_mat,alpha1_ind,alpha2_ind,
                 alpha_old=alpha[c(alpha1_ind,alpha2_ind)],alpha_new,eps,C)
      
      alpha_old=alpha
      alpha[c(alpha1_ind,alpha2_ind)]=alpha_new
      gg=sapply(1:length(alpha),g,alpha=alpha,y=y,
                ker_mat=ker_mat,b=b)
      E=update.E(E,y,alpha,ker_mat,b)
      cat("delta alpha:",sqrt(mean(sum((alpha_old-alpha)^2))),"\n")
      alpha_changed=alpha_changed+1
      cat("alpha changed:",alpha_changed,"\n")
      cat("\n")
    }
    
    if (examineAll)
      examineAll = F
    else if (alpha_changed == 0)
      examineAll = T
  }
  
  return(list(Alpha=c(alpha1,alpha2),B=b,E=E,GG=gg))
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
#    if (all(ss))
#      cat("b1-b2:",b1-b2,"\n")
  }
  return(b)
}

g<-function(alpha,y,ker_mat,col,b) ### col is the x colind
{
  res=sum(alpha*y*ker_mat[,col])+b
  return(res)
}


update.E<-function(E,y,alpha,ker_mat,b)
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
  res=Update(y,alpha,C,E,ker_mat,eps,b)
#  is_kkt=is.KKT(y,alpha,eps,C,gg) 
  return(res)
}
res=SMO(data=d[,1:4],y=d$y,C=1,eps=10^(-3))










