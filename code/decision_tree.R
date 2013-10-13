
loan=data.frame(
#  id=1:15,
  y=c(0,0,1,1,0,0,0,1,1,1,1,1,1,1,0),
  age=rep(1:3,each=5),
  work=c(0,0,1,1,0,0,0,1,0,0,0,0,1,1,0),
  house=c(0,0,0,1,0,0,0,1,1,1,1,1,0,0,0),
  credit=c(1,2,2,1,1,1,2,2,3,3,3,1,1,3,1)
  )


#### the decision tree
entropy<-function(x)
{
  xx=table(x)/length(x)
  res=-sum(xx*log2(xx))
  return(res)
}

entropy.cond<-function(cate,cond)
{
  H=tapply(cate,factor(cond),entropy)
  res=sum(table(cond)/length(cond)*H)
  return(res)
}

info.gain<-function(cate,cond)
{
  res=entropy(cate)-entropy.cond(cate,cond)
  return(res)
}

decision.tree<-function(data,y,rows,vars,eps,fit,route)
{
#  dd=data[rows,cols]
#  yy=y[rows]
  t_y=table(y)
  if (length(t_y)==1)
  {
    fit[rows]=y
    res=list(Vars=t_y,Fit=fit)
    return(res)
  }
  info_gain=apply(data,2,function(x)
    {
    res=info.gain(y,x)
    return(res)
  })
  max_i=which.max(info_gain)
  if (info_gain[max_i]<=eps)
  {
    fit[rows]=as.numeric(names(t_y)[which.max(t_y)])
    res=list(Vars=list(fit[rows[1]]),Fit=fit)
    return(res)
  }
#  cols=setdiff(cols,max_i)
#  eval(parse(text=paste("vars",route,"=list(colnames(data)[max_i])",sep="")))

  ind=split(1:nrow(data),factor(data[,max_i]))
  for (i in 1:length(ind))
  {
    rou=paste(route,"$",colnames(data)[max_i],"_",names(ind)[i],sep="")
    res=decision.tree(data[ind[[i]],-max_i],y[ind[[i]]],rows=rows[ind[[i]]],vars,eps,fit,rou="")
    
    eval(parse(text=paste("vars",rou,"=res$Vars",sep="")))  
    fit=res$Fit
  }
#  eval(parse(text=paste("names(vars)","=c('split_var',names(ind))",sep="")))
  res=list(Vars=vars,Fit=fit)
  return(res)
}

res=decision.tree(data=loan[,-1],y=loan$y,rows=1:nrow(loan),vars=list(),
                  eps=0,fit=rep(0,nrow(loan)),route="")



