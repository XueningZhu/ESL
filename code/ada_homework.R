hire=data.frame(health=c(0,0,1,1,1,0,1,1,1,0),work=c(1,3,2,1,2,1,1,1,3,2),
                potential=c(3,1,2,3,3,2,2,1,1,1),level=c(-1,-1,-1,-1,-1,-1,1,1,-1,-1))

hire=data.frame(health=factor(c(0,0,1,1,1,0,1,1,1,0)),work=factor(c(1,3,2,1,2,1,1,1,3,2)),
                potential=factor(c(3,1,2,3,3,2,2,1,1,1)),level=factor(c(-1,-1,-1,-1,-1,-1,1,1,-1,-1)))
hire$level=as.factor(hire$level)



hire=data.frame(health=c(0,0,1,1,1,0,1,1,1,0),work=c(1,3,2,1,2,1,1,1,3,2),
                potential=c(3,1,2,3,3,2,2,1,1,1),level=c(-1,-1,-1,-1,-1,-1,1,1,-1,-1))

library(ada)


x1=rnorm(100)
x2=rnorm(100)
y=x1^2+x2-sin(x1*x2)
y = as.numeric(y>1)
y=y*2-1
hire = data.frame(level=y,x1=x1,x2=x2)

hire_ada<-ada(level~.,data=hire,iter=10,loss="ada",nu=0,bag.frac=0,type="discrete",verbose=T)
predict(hire_ada,hire,type="prob")
varplot(hire_ada)

plot(hire_ada)

## fit discrete ada boost to a simple example
data(iris)
##drop setosa
iris[iris$Species!="setosa",]->iris
##set up testing and training data (60% for training)
n<-dim(iris)[1]
trind<-sample(1:n,floor(.6*n),FALSE)
teind<-setdiff(1:n,trind)
iris[,5]<- as.factor((levels(iris[,5])[2:3])[as.numeric(iris[,5])-1])
##fit 8-split trees
gdis<-ada(Species~.,data=iris[trind,],iter=20,nu=1,type="discrete")
##add testing data set
gdis=addtest(gdis,iris[teind,-5],iris[teind,5])
##plot gdis
plot(gdis,TRUE,TRUE)
##variable selection plot
varplot(gdis)
##pairwise plot
pairs(gdis,iris[trind,-5],maxvar=2)