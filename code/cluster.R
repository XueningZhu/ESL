



setwd("F:/商业数据挖掘_光华/homework/3")
firm = read.table("Firm.txt")

rownames(firm)=firm[,9]
firm=firm[,1:8]

library(lattice)

pseudo.F<-function(data,max_cluster)
{
  N=nrow(firm)
  pseudo_F=rep(0,max_cluster-1)
  for(k in 2:max_cluster)
  {
    kk=kmeans(firm,centers=k)
    wss=sum(kk$withinss)
    pseudo_F[k-1]=((kk$totss-wss)/(k-1))/(wss/(N-k))
  }
  return(pseudo_F)
}

pf=NULL
for (i in 1:4)
{
  pf=c(pf,pseudo.F(data=firm,max_cluster=15))
}
gg=factor(rep(1:4,each=14),levels=1:4,labels=paste("experiment",1:4))
xyplot(pf~rep(2:15,4)|gg,type="b",xlab="Cluster Number",ylab="Pseudo F")

# Determine number of clusters

pseudo.F<-function(data,max_cluster,measure,plot=T)
{
  N=nrow(firm)
  pseudo_F=rep(0,max_cluster-1)
  hist_data=NULL
  for(k in 2:max_cluster)
  {
    kk=kmeans(firm,centers=k)
    wss=sum(kk$withinss)
    tt=NULL
    for (i in 1:100)
    {
      kk=kmeans(firm,centers=k)
      wss=sum(kk$withinss)
      tt=c(tt,((kk$totss-wss)/(k-1))/(wss/(N-k)))
    }
    
    if (plot)
    {
      hist_data=c(hist_data,tt)
    }
    pseudo_F[k-1]=measure(tt)
  }
  if (plot) return(list(pseudo_F=pseudo_F,hist_data=hist_data))
  return(pseudo_F)
}

pseudo_res=pseudo.F(data=firm,max_cluster=15,measure=mean)
pseudo_F=pseudo_res[[1]]
pseudo_hist=pseudo_res[[2]]
group=factor(rep(2:15,each=100),levels=2:15,
             labels=paste(2:15,"cluster"))
histogram(~pseudo_hist | group, xlab="Pseudo F Statistics for 100 times")

plot(2:15, pseudo_F, type="b", xlab="Number of Clusters",
     ylab="Pseudo F")

pseudo_F=pseudo.F(data=firm,max_cluster=15,measure=median,plot=F) ### without the hist plot
plot(2:15, pseudo_F, type="b", xlab="Number of Clusters",
     ylab="Pseudo F")


wss <- (nrow(firm)-1)*sum(apply(firm,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(firm, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

kk=kmeans(firm,centers=3)
split(rownames(firm),kk$cluster)

#### to check the stability of the clustering
Kmeans.repeat<-function(data,cluster_num,expriment_times=2,repeat_times=100)
{
  MAT=list()
  for (j in 1:expriment_times)
  {
    mat=matrix(1,nrow(data),nrow(data))
    rownames(mat)=rownames(data)
    colnames(mat)=rownames(data)
    for (i in 1:repeat_times)
    {
      kk=kmeans(data,centers=cluster_num)
      mat[t(combn(which(kk$cluster==1),2))]=mat[t(combn(which(kk$cluster==1),2))]+1
    }
    MAT[[j]]=mat
  }
  return(MAT)
}
MAT=Kmeans.repeat(data=firm,cluster_num=3,expriment_times=2,repeat_times=100)

MAT[[1]]-MAT[[2]]


mat=MAT[[1]]
dis=1/t(mat)
diag(dis)=0
dis=as.dist(dis)
fit <- hclust(dis, method="average") 
plot(fit) # display dendogram
# draw dendogram with red borders around the 5 clusters 
clust1=rect.hclust(fit, k=3, border="red")
clust1



###### for scaled firm data
firm_scale=scale(firm)
pseudo_F=pseudo.F(data=firm_scale,max_cluster=15,measure=median,plot=F)
plot(2:15, pseudo_F, type="b", xlab="Number of Clusters",
     ylab="Pseudo F")

wss <- (nrow(firm_scale)-1)*sum(apply(firm_scale,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(firm_scale, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


mat=Kmeans.repeat(data=firm_scale,cluster_num=3,expriment_times=1,repeat_times=100)[[1]]
dis=1/t(mat)
diag(dis)=0
dis=as.dist(dis)
fit <- hclust(dis, method="average")
plot(fit) # display dendogram
# draw dendogram with red borders around the 5 clusters 
clust2=rect.hclust(fit, k=3, border="red")
clust2

get.cluster_id<-function(data,clust)
{
  cluster_id=rep(0,nrow(data))
  cluster_id[unlist(clust)]=rep(1:3,sapply(clust,length))
  cluster_id=as.factor(cluster_id)
  return(cluster_id)
}
library("MASS")
library("ggplot2")
dist=dist(firm_scale)

firm_mds = isoMDS(dist)
x = firm_mds$points[,1]
y = firm_mds$points[,2]
cluster_id=get.cluster_id(data=firm,clust=clust2)
g=ggplot(data.frame(x,y,cluster_id),aes(x,y,label = rownames(firm),color=cluster_id))
g+geom_point(shape=16,size=3)+
  geom_text(hjust=-0.1,vjust=0.5,angle=7,size=6)

#计算距离矩阵
d=dist(firm_scale)
#利用经典多维标度方法建模
firm_cmd=cmdscale(d,k=5,eig=T)
#计算前两个维度的特征值的比，超过0.9，意味着用两个维度即可表示出整个数据结构。
sum((firm_cmd$eig[1:2])^2)/sum((firm_cmd$eig)^2)
x <- firm_cmd$points[,1]
y <- firm_cmd$points[,2]
#绘制两维坐标
cluster_id=get.cluster_id(data=firm,clust=clust2)
g=ggplot(data.frame(x,y,cluster_id),aes(x,y,label = rownames(firm),color=cluster_id))
g+geom_point(shape=16,size=3)+
  geom_text(hjust=-0.1,vjust=0.5,angle=7,size=6)

# Ward Hierarchical Clustering
d <- dist(firm, method = "euclidean") # distance matrix

fit <- hclust(d, method="single") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

fit <- hclust(d, method="average") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

fit <- hclust(d, method="complete") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")





