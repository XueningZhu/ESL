

### max entropy problem
### sort the dataset with increasing y
### f is a feature matrix (N*length(table(y)))
### w is the same dimension matrix as f

cond.prob<-function(w,f)  
{
  ind=Y==y
  zw=sapply(1:ncol(f),function(i)
    {return(exp(sum(w[,i]*f[,i])))})
  res=zw/sum(zw)
  return(res)
}
gradient<-function(px,pw,ff,pxy)
{
  
}










