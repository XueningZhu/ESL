
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
  
  ind=split(1:nrow(data),factor(data[,max_i]))
  for (i in 1:length(ind))
  {
    rou=paste(route,"$",colnames(data)[max_i],"_",names(ind)[i],sep="")
    res=decision.tree(data[ind[[i]],-max_i],y[ind[[i]]],rows=rows[ind[[i]]],vars,eps,fit,rou="")
    
    eval(parse(text=paste("vars",rou,"=res$Vars",sep="")))  
    fit=res$Fit
  }
  res=list(Vars=vars,Fit=fit)
  return(res)
}

res=decision.tree(data=loan[,-1],y=loan$y,rows=1:nrow(loan),vars=list(),
                  eps=0,fit=rep(0,nrow(loan)),route="")

#### for pruning a decision tree


Loss<-function(vec,N_leaf,alpha)
{
  loss=entropy(vec)+alpha*N_leaf
  return(loss)
}

pruning<-function(father_loss,sons_loss,point,y,rows,fit)
{
  if (father_loss-sum(sons_loss)>0)
  {
    attr(point,"loss")=sum(sons_loss)
    return(list(Point=point,Fit=fit))
  }
  else
  {
    res=table(y[rows])
    attr(res,"loss")=father_loss
    max_ind=which.max(res)
    fit[rows]=rep(as.numeric(names(res)[max_ind]),length(rows))
    return(list(Point=res,Fit=fit))
  }
}

pruning.tree<-function(point,data,y,rows,fit,alpha)
{
  is_leaf=!sapply(point,is.list)
  spl_var=gsub("_[0-9]+","",names(point))[1]
  father_loss=Loss(y[rows],1,alpha)
  spl_groups=split(rows,factor(data[rows,spl_var]))
  if (all(is_leaf))
  {
    sons_loss=sapply(spl_groups,function(x)
      {return(Loss(y[x],N_leaf=1,alpha))})    
    res=pruning(father_loss,sons_loss,point,y,rows,fit)
    return(res)
  }
  sons_loss=NULL
  for (i in 1:length(is_leaf))
  {
    if (is_leaf[i])
    {
      tmp=Loss(y[spl_groups[[i]]],N_leaf=1,alpha)
      sons_loss=c(sons_loss,tmp)
    }
    else
    {
      tmp=pruning.tree(point[[i]],data,y,spl_groups[[i]],fit,alpha)
      point[[i]]=tmp$Point
      sons_loss=c(sons_loss,attr(point[[i]],"loss"))
      fit=tmp$Fit
    }
  }
  res=pruning(father_loss,sons_loss,point,y,rows,fit)
  return(res)
}

res=decision.tree(data=loan[,-1],y=loan$y,rows=1:nrow(loan),vars=list(),
                  eps=0,fit=rep(0,nrow(loan)),route="")

res_pru=pruning.tree(point=res$Vars,data=loan[,-1],y=loan$y,rows=1:nrow(loan),fit=res$Fit,alpha=1)











