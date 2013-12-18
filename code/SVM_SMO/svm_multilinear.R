

library(e1071)

f <- function(x, y) 1 / (1 + 2^(x^3+ y+x*y))
d <- data.frame(x1 = rnorm(1000), x2 = rnorm(1000))
d$y = with(d, ifelse(runif(1000) < f(x1, x2), 1, -1))
dd=d

dd$x5=d$x1+100*d$x2  ### can also be tried with dd$x5=d$x1+100000*d$x2

model1=svm(y~.,data=d,kernel='radial')
pred1=predict(model1)
model2=svm(y~.,data=dd,kernel='radial')
pred2=predict(model2)

library(AUC)

roc=roc(pred1,factor(d$y))
plot(roc)
auc(roc)

roc=roc(pred2,factor(d$y))
plot(roc)
auc(roc)

coef1=rep(0,nrow(d))
coef2=rep(0,nrow(d))
coef1[model1$index]=model1$coefs
coef2[model2$index]=model2$coefs

mean((coef1-coef2)^2)


