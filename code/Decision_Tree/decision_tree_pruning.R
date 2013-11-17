alpha.calculate<-function(father_entro,sons_entro)
{
  alpha=(father_entro-sum(sons_entro))/(length(sons_entro)-1)
  return(alpha)
}

alpha.compare<-function(father_alpha,sons_alpha,point)
{
  alpha=c(father_alpha,sons_alpha)
  min_ind=which.min(alpha)
  attr(point,"alpha")=alpha[min_ind]
  if (min_ind==1)
  {
    attr(point,"min.point")=0
  }
  else
  {
    attr(point,"min.point")=min_ind[1]-1
  }
  return(point)
}

Up.tree<-function(point,data,y,rows)
{
  is_leaf=!sapply(point,is.list)
  spl_var=gsub("_[0-9]+","",names(point))[1]
  father_entro=entropy(y[rows])
  spl_groups=split(rows,factor(data[rows,spl_var]))
  sons_entro=sapply(spl_groups,function(x)
  {return(entropy(y[x]))})  
  father_alpha=alpha.calculate(father_entro,sons_entro)
  
  if (all(is_leaf))
  {
    
    attr(point,"alpha")=father_alpha
    attr(point,"min.point")=0
    attr(point,"entropy")=sum(sons_entro)
   # res=pruning(father_loss,sons_loss,point,y,rows,fit)
    return(point)
  }
  sons_alpha=NULL
  sons_entro1=NULL
  for (i in 1:length(is_leaf))
  {
    if (is_leaf[i])
    {
#      tmp=entropy(y[spl_groups[[i]]])
      sons_alpha=c(sons_alpha,100000)
      sons_entro1=c(sons_entro1,sons_entro[i])
    }
    else
    {
      tmp=Up.tree(point[[i]],data,y,spl_groups[[i]])
      point[[i]]=tmp
      #sons_loss=c(sons_loss,attr(point[[i]],"loss"))
      sons_entro1=c(sons_entro1,attr(point[[i]],"entropy"))
      sons_alpha=c(sons_alpha,attr(point[[i]],"alpha"))
    }
  }
  father_alpha=alpha.calculate(father_entro,sons_entro1)
  attr(point,"alpha")=father_alpha
  attr(point,"entropy")=sum(sons_entro1)
  res=alpha.compare(father_alpha,sons_alpha,point)
  return(res)
}

pruning<-function(point,y,rows,fit)
{
  res=table(y[rows])
  max_ind=which.max(res)
  fit[rows]=rep(as.numeric(names(res)[max_ind]),length(rows))
  return(list(Point=res,Fit=fit))
}


Down.tree<-function(point,y,rows,fit)
{
  min_point=attr(point,"min.point")
  if (min_point==0)
  {
    res=pruning(point,y,rows,fit)
   # point=res$Point
  #  fit=res$Fit
    return(res)
  }
  res=Down.tree(point[[min_point]],y,rows=attr(point,"spl_groups")[[min_point]],fit)
  point[[min_point]]=res$Point
  fit=res$Fit
  return(list(Point=point,Fit=fit))
}

tree_alpha=Up.tree(point=res[[1]],data=loan[,-1],y=loan$y,rows=1:nrow(loan))
aa=Down.tree(point=tree_alpha,y=loan$y,rows=1:nrow(loan),fit=res$Fit)
debug(pruning.tree)


tree.sequence<-function(point=point,data=data,y=y,Fit)
{
  tree_seq=list(point)
  alpha=10000
  fit=matrix(Fit,ncol=1)
  i=2
  while(is.list(point))
  {
    tree_alpha=Up.tree(point=point,data=data,y=y,rows=1:nrow(data))
    alpha=c(alpha,attr(tree_alpha,"alpha"))
    res=Down.tree(point=tree_alpha,y=y,rows=1:nrow(data),fit=Fit)
    
    point=res$Point
    tree_seq[[i]]=point
    Fit=res$Fit
    fit=cbind(fit,Fit)
    i=i+1
  }
  return(list(Tree_seq=tree_seq,Fit=fit,Alpha=c(alpha,attr(point,"alpha"))))
}

tree_seq=tree.sequence(point=res$Vars,data=loan[,-1],y=loan$y,Fit=res$Fit)













