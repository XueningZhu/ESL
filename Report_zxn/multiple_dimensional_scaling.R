library("MASS")
library("ggplot2")
load("zip2.rda")
dist=1-cor(zip2[,-1])

zip2_mds = isoMDS(dist)
x = zip2_mds$points[,1]
y = zip2_mds$points[,2]
g=ggplot(data.frame(x,y),aes(x,y,label = colnames(zip2)[-1]))
g+geom_point(shape=16,size=3,colour='red')+
  geom_text(hjust=-0.1,vjust=0.5,angle=7,alpha=0.5)



