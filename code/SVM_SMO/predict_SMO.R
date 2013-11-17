
predict.SMO<-function(x,y,alpha,data,Y,b,ker=Gauss.ker)
{
  X=c(x,y,0,0)
  kk=apply(data,1,ker,z=X)
  res=sum(Y*alpha*kk)+b
  return(res)
}
a = seq(-2, 2, len = 50)
predict.SMO(x=a[1],y=a[3],alpha=res$Alpha,data=d[,1:4],Y=d$y,b=res$B,ker=Poly.Ker)
aa_poly=sapply(a,function(x)
{
  rr=sapply(a,predict.SMO,y=x,alpha=res_poly$Alpha,
            data=d[,1:4],Y=d$y,b=res$B,ker=Poly.Ker)
  return(rr)
})

logit_aa_poly=exp(aa_poly)/(1+exp(aa_poly))
image.plot(a, a, logit_aa_poly, zlim = c(min(logit_aa_poly), max(logit_aa_poly)))

aa_gauss=sapply(a,function(x)
{
  rr=sapply(a,predict.SMO,y=x,alpha=res_gauss$Alpha,
            data=d[,1:4],Y=d$y,b=res$B,ker=Gauss.Ker)
  return(rr)
})

#bb=outer(a, a, f)
#image.plot(a, a, bb, zlim = c(min(bb), max(bb)))

logit_aa_gauss=exp(aa_gauss)/(1+exp(aa_gauss))
image.plot(a, a, logit_aa_gauss, zlim = c(min(logit_aa_gauss), max(logit_aa_gauss)))
#g <- function(x, y) 1 / (1 + 2^(cos(x^3)-x*y+sin(y)))

library(AUC)
## to calculate the auc
res=res_gauss ### can also be res_poly
fit=(res$Fit-min(res$Fit))/(max(res$Fit)-min(res$Fit))
yy=d$y
yy[yy==-1]=0
plot(sensitivity(fit,factor(yy)))
roc_smo=roc(fit,factor(yy))
plot(roc_smo)
auc(roc_smo)


