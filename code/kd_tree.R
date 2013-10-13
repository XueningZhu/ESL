
### the experimental data points
points=matrix(c(2,3,5,4,9,6,4,7,8,1,7,2,10,12,13,6),byrow=T,ncol=2)

### construct kd-tree

### construct kd-tree
### split--split dimension; rows -- rows used in kd-tree Res -- for the recursive design
kd.tree<-function(points,split,rows,Res)
{
  vec=points[rows,split]
  if (length(vec)%%2==0)
    root_split=median(vec[-1])
  else
    root_split=median(vec)
  
  root=rows[vec==root_split][1]
  left=rows[vec<root_split]
  right=setdiff(rows,c(left,root))
  if (length(left)==0) left=0
  if (length(right)==0) right=0
  
  if (length(left)<=1&length(right)<=1)
  {
    res=rbind(Res,c(root,left,right))
    return(res)
  }
  if (length(left)<=1&length(right)>1)
  {
    tmp=kd.tree(points,split%%ncol(points)+1,right,Res)
    res=rbind(c(root,left,tmp[nrow(tmp),1]),tmp)
    return(res)
  }
  if (length(left)>1&length(right)<=1)
  {
    tmp=kd.tree(points,split%%ncol(points)+1,left,Res)
    res=rbind(tmp,c(root,tmp[nrow(tmp),1],right))
    return(res)
  }
  
  tmp1=kd.tree(points,split%%ncol(points)+1,left,Res)
  tmp2=kd.tree(points,split%%ncol(points)+1,right,tmp1)
  
  res=rbind(tmp2,c(root,tmp1[nrow(tmp1),1],tmp2[nrow(tmp2),1]))
  
  show(res)
  return(res)
}


kd.tree.main<-function(points)
{
  dim_sd=apply(points,2,sd)
  spl_ord=order(dim_sd,decreasing=T)
  points_new=points[,spl_ord]
  kd_tree=kd.tree(points_new,split=1,rows=1:nrow(points),Res=NULL)
  return(list(Spl_ord=spl_ord,Kd_tree=kd_tree))
}


### for searching a kd-tree
### searching the leaf in the point region from top to bottom (using the rules of kd-tree construction)


is.leaf<-function(num_point,leaves)
{
  return(is.element(num_point,leaves))
}

search.leaf.main<-function(points,KdTree,point)
{
  spl_ord=KdTree$Spl_ord
  kd_tree=KdTree$Kd_tree
  
  point=point[spl_ord]
  points=points[,spl_ord]
  leaves=setdiff(c(kd_tree[,2],kd_tree[,3]),c(kd_tree[,1],0))
  
  leaf = search.leaf.sub(1,kd_tree,point,points,kd_line=kd_tree[nrow(kd_tree),],leaves)
  return(leaf)
}


search.leaf.sub<-function(ord,kd_tree,point,points,kd_line,leaves)
{
  roots=kd_line
  Next=ifelse(point[ord]<points[roots[1],ord],roots[2],roots[3])
  if(is.leaf(Next,leaves))
    return(Next)
  line=which(Next==kd_tree[,1])
  show(line)
  
  new_line=kd_tree[line,]
  kd_tree=kd_tree[1:line,]
  return(search.leaf.sub(ord%%ncol(points)+1, kd_tree, point, points ,kd_line=new_line, leaves))
}



### to trace back the tree until arriving at the root


Euler.dist<-function(x,y)
{
  return(sqrt(sum((x-y)^2)))
}

find.father<-function(leaf,kd_tree)
{
  father_line=which(is.element(kd_tree[,2],leaf)|is.element(kd_tree[,3],leaf))
  return(father_line)
}

is.dup<-function(father_point,ord,point,min_dist)
{
  res=father_point[ord]>(point[ord]-min_dist)&father_point[ord]<(point[ord]+min_dist)
  return(res)
}

is.root<-function(father,kd_tree)
{
  res=father==kd_tree[nrow(kd_tree),1]
  return(res)
}

kd.traceback.main<-function(points,leaf,point,kd_tree)
{
  
  
  min_dist=Euler.dist(leaf,point)
  min_point=leaf
  first_ord=(nrow(kd_tree)-1)%%ncol(points)
  if (first_ord==0) first_ord=ncol(points)

  res=kd.traceback.sub(leaf,point,points,ord=first_ord,min_dist,min_point,kd_tree)
  return(res)
}

kd.traceback.sub<-function(leaf,point,points,ord,min_dist,min_point,kd_tree)
{
  leaf_line=find.father(leaf,kd_tree)
  father=kd_tree[leaf_line,1]
  fat_dist=Euler.dist(points[father,],point)
  
  cat("fat_point:",father,"\n")
  cat("fat_dist:",fat_dist,"\n")
  
  if (fat_dist<min_dist)
  {
    min_point=father
    min_dist=fat_dist
  }
  if (is.dup(father_point=points[father,],ord,point,min_dist))
  {
    show("dup!")
    son=setdiff(kd_tree[leaf_line,2:3],leaf)
    son_dist=Euler.dist(points[son,],point)
    
    cat("son_point:",son,"\n")
    cat("son_dist:",son_dist,"\n")
    
    son_line=which(kd_tree[,1]==son)
    son_leaf=search.leaf.main(points,KdTree=kd_tree[1:son_line,],point)
    cat("son_leaf:",son_leaf,"\n")
    
    son_res=kd.traceback.main(points,leaf=son_leaf,point,kd_tree=kd_tree[1:son_line,])
    min_point=son_res$Min_point
    min_dist=son_res$Min_dist
    
  }
  cat("min_point:",min_point,"\n")
  cat("min_dist:",min_dist,"\n")
  
  if (is.root(father,kd_tree))
    return(list(Min_point=min_point,Min_dist=min_dist))
  
  kd_tree=kd_tree[leaf_line:nrow(kd_tree),]
  
  ord_next=(ord-1)%%nrow(points)
  if (ord_next==0) ord_next=1
  
  res=kd.traceback.sub(leaf=father,point,points,ord_next,min_dist,min_point,kd_tree)
  return(res)
}


main.search<-function(points,kd_tree,point)
{
  spl_ord=kd_tree[[1]]
  points=points[,spl_ord]
  point=point[spl_ord]
  leaf=search.leaf.main(points,KdTree=kd_tree,point)
  cat("leaf:",leaf,"\n")
  
  res=kd.traceback.main(points,leaf,point,kd_tree=kd_tree[[2]])
  return(list(First_leaf=leaf,Kd_res=res))
}


kd_tree=kd.tree.main(points)
res=main.search(points,kd_tree,point=c(3,4.5))
kd_res=res[[2]]
leaf=res[[1]]



