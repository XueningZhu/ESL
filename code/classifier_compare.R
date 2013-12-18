
library(AUC)
library(ada)
library(e1071)
library(randomForest)
library(kknn)
library(gbm)
library(klaR)
library(ggplot2)

data.generation<-function(n,test_ratio=1/5)
{
  ##### data preparation
  p<-10
  f<-function(x,a,b,d){
    return( a*(x-b)^2+d )
  }
  x1<-runif(n/2,0,4)
  y1<-f(x1,-1,2,1.7)+runif(n/2,-1,1)
  x2<-runif(n/2,2,6)
  y2<-f(x2,1,4,-1.7)+runif(n/2,-1,1)
  y<-c(rep(0,n/2),rep(1,n/2))
  dat<-data.frame(y=factor(y),x1=c(x1,x2),x2=c(y1,y2),
                  matrix(rnorm(n*8),ncol=8))
  names(dat)<-c("y",paste("x",1:10,sep=""))
  plot(dat$x1,dat$x2,pch=c(1:2)[y+1], col=c(1,8)[y+1],
       xlab=names(dat)[2],ylab=names(dat)[3])
  ggplot(dat,aes(x1,x2))+
    geom_point(aes(colour=factor(y)),size=5)
  
  indtest<-sample(1:n,floor(n*test_ratio),FALSE)
  
  train<-dat[-indtest,]
  test<-dat[indtest,]
  return(list(train=train,test=test))
}

a=data.generation(1000)


#### cross validation
gety=function(formula,data)#get response from a formula, code copied from "lm"
{
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  response<-model.response(mf)
  return(unname(response))
}

MSE<-function(y,fit)
{
  if(is.factor(y))
    y=as.numeric(y)-1
  return(mean((y-fit)^2))
}

CV.sub<-function(ind,Train,Test,Train_y,Test_y,formula,pred_alone=T,func,...)
{
  Trainset=Train[-ind,]
  Testset=Train[ind,]
  if (!pred_alone)
  {
    Model=func(formula,Trainset,Testset,Test,...)
    pred_fit=Model$pred_fit
    pred=Model$pred
    Test_pred=Model$test_pred
  }
  else
  {
    Model=func(formula,Trainset,...)
    pred_fit=predict(Model,Trainset,type="prob")
    pred=predict(Model,Testset,type="prob")
    Test_pred=predict(Model,Test,type="prob")
    if (is.matrix(pred))
    {
      pred_fit=pred_fit[,ncol(pred_fit)]
      pred=pred[,ncol(pred)]
      Test_pred=Test_pred[,ncol(Test_pred)]
    }
  }
  fit_auc=auc(roc(pred_fit,as.factor(Train_y[-ind])))
  fit_mse=MSE(Train_y[-ind],pred_fit)
  
  cv_auc=auc(roc(pred,as.factor(Train_y[ind])))
  cv_mse=MSE(Train_y[ind],pred)
  
  Test_auc=auc(roc(Test_pred,as.factor(Test_y)))
  Test_mse=MSE(Test_y,Test_pred)
  return(list(fit_auc=fit_auc,cv_auc=cv_auc,Test_auc=Test_auc,
              fit_mse=fit_mse,cv_mse=cv_mse,Test_mse=Test_mse,
              Test_pred=Test_pred))
}


CrossValid<-function(Train,Test,cv=5,formula,func,...)
{
  Train_y=gety(formula,Train)
  Test_y=gety(formula,Test)
  group=split(1:nrow(Train),rep(1:cv,length.out=nrow(Train)))
  Train_auc=NULL
  Test_auc=NULL
  Test_pred=NULL
  
  res=lapply(group,CV.sub,Train,Test,Train_y,Test_y,formula,func,...)
  fit_auc=sapply(res,function(x) return(x$fit_auc))
  cv_auc=sapply(res,function(x) return(x$cv_auc))
  Test_auc=sapply(res,function(x) return(x$Test_auc))
  
  fit_mse=sapply(res,function(x) return(x$fit_mse))
  cv_mse=sapply(res,function(x) return(x$cv_mse))
  Test_mse=sapply(res,function(x) return(x$Test_mse))
  
  Test_pred=sapply(res,function(x) return(x$Test_pred))
  return(list(fit_auc=fit_auc,cv_auc=cv_auc,Test_auc=Test_auc,
              fit_mse=fit_mse,cv_mse=cv_mse,Test_mse=Test_mse,
              Test_pred=Test_pred))
}


#ada_sim=CrossValid(train,test,cv=5,y~.,ada,iter=50,loss="e",type="discrete")


Mykknn<-function(formula,Trainset,Testset,Test,...)
{
  y=gety(formula,Trainset)
  model=kknn(formula,train=Trainset,
             test=data.frame(rbind(Trainset,Testset,Test)),...)
  
  prob=model$prob[,2]
  pred_fit=prob[1:nrow(Trainset)]
  prob=prob[-(1:nrow(Trainset))]
  
  pred=prob[1:nrow(Testset)]
  test_pred=prob[-(1:nrow(Testset))]
  
  return(list(pred_fit=pred_fit,pred=pred,test_pred=test_pred))
}

logit<-function(formula,Trainset,Testset,Test,...)
{
  model=glm(formula,family=binomial,data=Trainset,...)
  pred_fit=predict(model,Trainset,type="response")
  pred=predict(model,Testset,type="response")
  test_pred=predict(model,Test,type="response")
  return(list(pred_fit=pred_fit,pred=pred,test_pred=test_pred))
}

logit.trans<-function(x)
{
  res=exp(x)/(1+exp(x))
  return(res)
}

Mygbm<-function(formula,Trainset,Testset,Test,...)
{
  col=as.character(formula)[2]
  y1=gety(formula,Trainset)
  Trainset[,col]=as.numeric(y1)-1
  model=gbm(formula,distribution = "bernoulli",data=Trainset,...)
  pred_fit = logit.trans(predict(model,Trainset,...))
  pred = logit.trans(predict(model,Testset,...))
  test_pred = logit.trans(predict(model,Test,...))
  return(list(pred_fit=pred_fit,pred=pred,test_pred=test_pred))
}


MyNB<-function(formula,Trainset,Testset,Test,...)
{
  model=NaiveBayes(formula,data=Trainset,...)
  pred_fit = predict(model,Trainset)$posterior[,2]
  pred = predict(model,Testset)$posterior[,2]
  test_pred = predict(model,Test)$posterior[,2]
  return(list(pred_fit=pred_fit,pred=pred,test_pred=test_pred))
}
Mysvm<-function(formula,Trainset,Testset,Test,...)
{
  model=svm(formula,data=Trainset,probability=T,...)
  
  pred_fit = attr(predict(model,Trainset,probability=T),"probabilities")[,2]
  pred = attr(predict(model,Testset,probability=T),"probabilities")[,2]
  test_pred = attr(predict(model,Test,probability=T),"probabilities")[,2]
  return(list(pred_fit=pred_fit,pred=pred,test_pred=test_pred))
}

#### model benchmark
data=data.generation(1000)
train=data$train
test=data$test
logit_sim=CrossValid(Train=train,Test=test,
                     cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,logit)
knn_sim=CrossValid(Train=train,Test=test,cv=ceiling(nrow(train)/nrow(test))+1,factor(y)~.,pred_alone=F,Mykknn)
NB_sim=CrossValid(Train=train,Test=test,
                  cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,MyNB, usekernel=T,kernel="optcosine")
svm_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,Mysvm,probabiliy=T)
ada_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=T,ada,iter=50,loss="e",type="discrete") 
gbm_sim=CrossValid(Train=train,Test=test,
                   cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,Mygbm,n.trees=500)
rf_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=T,randomForest,ntree=300)
res=list(logit_sim,knn_sim,NB_sim,svm_sim,ada_sim,gbm_sim,rf_sim)

score=c("fit_auc","cv_auc","Test_auc","fit_mse","cv_mse","Test_mse")
model_name=c("logit","kknn","NB","svm","ada","gbm","rf")
rr=do.call(c,lapply(res,function(x){
  y=sapply(1:(length(x)-1),function(i){
    return(mean(x[[i]]))
  })
  names(y)=score
  return(y)
}))
auc_ind=grep("auc",names(rr))
rr_auc=data.frame(model=factor(rep(model_name,each=3)),index=names(rr[auc_ind]),score=rr[auc_ind])
rr_mse=data.frame(model=factor(rep(model_name,each=3)),index=names(rr[-auc_ind]),score=rr[-auc_ind])

ggplot(rr_auc, aes(model, score, fill = index)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text = element_text(size=rel(1.5),colour="black"),
        axis.title.y = element_text(size = rel(1.5), angle = 90,face="bold"),
        axis.title.x = element_text(size = rel(1.5),face="bold"))+
  guides(fill = guide_legend(title="Index",
                           title.theme=element_text(size=15,angle=0,face="bold"),
                           label.theme = element_text(size=12,face="bold",angle = 0)))
                           
ggplot(rr_mse, aes(model, score, fill = index)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text = element_text(size=rel(1.5),colour="black"),
        axis.title.y = element_text(size = rel(1.5), angle = 90,face="bold"),
        axis.title.x = element_text(size = rel(1.5),face="bold"))+
  guides(fill = guide_legend(title="Index",
                             title.theme=element_text(size=15,angle=0,face="bold"),
                             label.theme = element_text(size=12,face="bold",angle = 0)))


#### model comparision

compare<-function(Ns,n_sim=3)
{
  score=c("fit_auc","cv_auc","Test_auc","fit_mse","cv_mse","Test_mse")
  NN=lapply(Ns,function(n){
    tmp=lapply(1:n_sim,function(i){
      cat("n:",n,"sim:",i,"\n")
      data=data.generation(n)
      train=data$train
      test=data$test
      logit_sim=CrossValid(Train=train,Test=test,
                           cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,logit)
      knn_sim=CrossValid(Train=train,Test=test,cv=ceiling(nrow(train)/nrow(test))+1,factor(y)~.,pred_alone=F,Mykknn)
      NB_sim=CrossValid(Train=train,Test=test,
                        cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,MyNB, usekernel=T,kernel="optcosine")
      svm_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,Mysvm,probabiliy=T)
      ada_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=T,ada,iter=50,loss="e",type="discrete") 
      gbm_sim=CrossValid(Train=train,Test=test,
                         cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=F,Mygbm,n.trees=500)
      rf_sim=CrossValid(train,test,cv=ceiling(nrow(train)/nrow(test))+1,y~.,pred_alone=T,randomForest,ntree=300)
      model_name=paste(c("logit","kknn","NB","svm","ada","gbm","rf"),"_size_",n,"_sim_",i,sep="")
      
      
      res=list(logit_sim,knn_sim,NB_sim,svm_sim,ada_sim,gbm_sim,rf_sim)
      res1=lapply(1:6,function(j){
        assign(score[j],
               unlist(lapply(1:length(res),function(i){
                     x=res[[i]][[j]]; names(x)=rep(model_name[i],length(x))
                    return(x)})))
        eval(parse(text=score[j]))
      })
    })
    rr=lapply(1:length(score),function(ii){
      assign(score[ii],
             unlist(lapply(1:length(tmp),function(jj){
                         return(tmp[[jj]][[ii]])})))
      return(eval(parse(text=score[ii])))
    })
    return(rr)
  })
  res=lapply(1:length(score),function(ii){
    assign(score[ii],
           unlist(lapply(1:length(NN),function(jj){
             return(NN[[jj]][[ii]])})))
    return(eval(parse(text=score[ii])))
  })
  return(res)
}

debug(compare)
sim=compare(Ns=c(500,800,1000,1500,2000,5000),n_sim=3)



plot(sim[[2]],sim[[5]])
ll=lm(sim[[2]]~sim[[5]])
fit_draw=data.frame(model=factor(gsub("_\\w+","",names(sim[[1]]))),auc=sim[[1]],mse=sim[[4]])
cv_draw=data.frame(model=factor(gsub("_\\w+","",names(sim[[2]]))),auc=sim[[2]],mse=sim[[5]])
pred_draw=data.frame(model=factor(gsub("_\\w+","",names(sim[[3]]))),auc=sim[[3]],mse=sim[[6]])

ggplot(fit_draw,aes(auc,mse))+
  geom_point(aes(colour=model),size=5)+
  theme(axis.text = element_text(size=rel(1.5),colour="black"),
      axis.title.y = element_text(size = rel(1.5), angle = 90,face="bold"),
      axis.title.x = element_text(size = rel(1.5),face="bold"))+
  guides(fill = guide_legend(title="Index",
                             title.theme=element_text(size=15,angle=0,face="bold"),
                             label.theme = element_text(size=12,face="bold",angle = 0)))

ggplot(cv_draw,aes(auc,mse))+
  geom_point(aes(colour=model),size=5)+
  theme(axis.text = element_text(size=rel(1.5),colour="black"),
        axis.title.y = element_text(size = rel(1.5), angle = 90,face="bold"),
        axis.title.x = element_text(size = rel(1.5),face="bold"))+
  guides(fill = guide_legend(title="Index",
                             title.theme=element_text(size=15,angle=0,face="bold"),
                             label.theme = element_text(size=12,face="bold",angle = 0)))

ggplot(pred_draw,aes(auc,mse))+
  geom_point(aes(colour=model),size=5)+
  theme(axis.text = element_text(size=rel(1.5),colour="black"),
        axis.title.y = element_text(size = rel(1.5), angle = 90,face="bold"),
        axis.title.x = element_text(size = rel(1.5),face="bold"))+
  guides(fill = guide_legend(title="Index",
                             title.theme=element_text(size=15,angle=0,face="bold"),
                             label.theme = element_text(size=12,face="bold",angle = 0)))
