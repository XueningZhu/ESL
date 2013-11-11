

setwd("F:/商业数据挖掘_光华/homework/2")
require(knitr)
######
library(sas7bdat)

############ read data files
zip1=read.sas7bdat("zipdemo1.sas7bdat")
zip2=read.sas7bdat("zipdemo2.sas7bdat")
shopping=read.table("Shopping.txt",header=T,sep="\t")


save(zip1,file="zip1.rda")
save(zip2,file="zip2.rda")
save(shopping,file="shopping.rda")

###########################################

load("shopping.rda")

for (i in 1:10)
{
  shopping[,i]=as.logical(shopping[,i])
}
sum_shop=summary(shopping)
write.csv(sum_shop,file="sum_shop.csv")
############# for the association rules mining
shopping[,1:10]=apply(shopping[,1:10],2,function(x)
  {
  x=as.logical(x)
  y=table(x)/length(x)
  if (y[2]<=0.2)
    x[!x]=NA
  return(x)
})
summary(shopping)


#### 1 without negative transactions
#### unsupervised (11-15 are descriptive variables)
#for (i in 1:10)
#{
#  shopping[,i]=as.logical(shopping[,i])
#  shopping[!shopping[,i],i]=NA 
#}

library(arules)
shop_trans=as(shopping, "transactions")
summary(shop_trans)


#image(shop_trans[1:10])
  ### find the frequent terms
itemFrequencyPlot(shop_trans, support = 0.2, cex.names = 0.8)

itemFreq=itemFrequency(shop_trans, type = "relative")
quantile(itemFreq)
#apriori
rules <- apriori(shop_trans, parameter = list(support = 0.1,confidence = 0.7))
rule_lif=subset(rules, subset = lift > 1.2)
inspect(head(sort(rule_lif, by = "confidence"), n = 20))
#fsets <- eclat(shop_trans, parameter = list(support = 0.05), control = list(verbose = FALSE))

summary(rule_lif)

### convert to dataframe
Rules=as(rule_lif,"data.frame")
#Rules=Rules[Rules$lift>1,]
Rules$rules=as.character(Rules$rules)
rule_spl=do.call(rbind,sapply(Rules$rules,strsplit,"=>"))
rownames(rule_spl)=NULL
colnames(rule_spl)=c("lhs","rhs")
Rules=data.frame(Rules,rule_spl)
Rules=Rules[order(Rules$support,Rules$confidence,decreasing=T),]

supervise=Rules[grepl("(TRUE)|(FALSE)",Rules$rhs,perl=T),]
             ,]
supervise[1:20,1:4]

cate_supvise=lapply(colnames(shopping[1:10]),function(x)
  {
  y=Rules[grepl(x,Rules$rhs,perl=T),1:4]
  cat(x,"\n")
  show(y[1:min(nrow(y),10),])
  return(y)
})

string=paste("(",paste(colnames(shopping[11:15]),collapse=")|("),")",sep="")
people=Rules[grepl(string,Rules$lhs)
             &!grepl("(TRUE)|(FALSE)",Rules$lhs)
             &grepl("(TRUE)|(FALSE)",Rules$rhs)
             &!grepl(string,Rules$rhs),]
people[1:20,1:4]

item=Rules[!grepl(string,Rules$lhs)
             &grepl("(TRUE)|(FALSE)",Rules$lhs)
             &grepl("(TRUE)|(FALSE)",Rules$rhs)
             &!grepl(string,Rules$rhs),]
item[1:20,1:4]

pp=Rules[grepl(string,Rules$lhs)
           &!grepl("(TRUE)|(FALSE)",Rules$lhs)
           &!grepl("(TRUE)|(FALSE)",Rules$rhs)
           &grepl(string,Rules$rhs),]
pp
gender=Rules[grepl("GENDER",Rtules$lhs)&
               !grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)&
               grepl("TRUE",Rules$rhs)
             ,]
age=Rules[grepl("Age",Rules$lhs)&
            !grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)&
            grepl("TRUE",Rules$rhs)
          ,]

marital=Rules[grepl("MARITAL",Rules$lhs)
             &!grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)
           # &grepl("TRUE",Rules$rhs)
          ,]

children=Rules[grepl("CHILDREN",Rules$lhs)
              &!grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)
              # &grepl("TRUE",Rules$rhs)
              ,]

children=Rules[grepl("CHILDREN",Rules$lhs)
               &!grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)
               # &grepl("TRUE",Rules$rhs)
               ,]

working=Rules[grepl("WORKING",Rules$lhs)
               &!grepl("(TRUE)|(FALSE)",Rules$lhs,perl=T)
               # &grepl("TRUE",Rules$rhs)
               ,]

### the lift > 1.2
rule_lif=subset(rules, subset = lift > 1.2)
inspect(rules[1:100])
### sorted rules


### positive & negative rules all included



summary(rules)
summary(fsets)

colnames(shopping)
rule=as(rules,"data.frame")

### mining the unfrequency terms
rules <- apriori(shop_trans, parameter = list(support = 0.01,confidence = 0.6))
rule=subset(rules, subset = lift > 1.2 )
summary(rule)

Rules=as(rule,"data.frame")
Rules$rules=as.character(Rules$rules)
rule_spl=do.call(rbind,sapply(Rules$rules,strsplit,"=>"))
rownames(rule_spl)=NULL
colnames(rule_spl)=c("lhs","rhs")
Rules=data.frame(Rules,rule_spl)
Rules=Rules[order(Rules$support,Rules$confidence,decreasing=T),]

infreq_item=itemFreq[itemFreq<quantile(itemFreq,0.25)]#&grepl("(TRUE|FALSE)",itemFreq,perl=T)]
items=sapply(names(infreq_item),function(x)
{
  res=grepl(x,Rules$rhs)
  return(res)
})

min_freq=apply(items,1,function(x) {
  return(min(infreq_item[x]))
})

min_freq[!is.finite(min_freq)]=NA
save(min_freq,file="min_freq.rda")
load("min_freq.rda")
Rules[which(!is.na(min_freq))[1:20],]

infreq_item=itemFreq[itemFreq<quantile(itemFreq,0.25)&grepl("(TRUE|FALSE)",names(itemFreq),perl=T)]
items=sapply(names(infreq_item),function(x)
{
  res=grepl(x,Rules$rhs)
  return(res)
})
min_freq_item=apply(items,1,function(x) {
  return(min(infreq_item[x]))
})
min_freq_item[!is.finite(min_freq_item)]=NA
save(min_freq_item,file="min_freq_item.rda")
load("min_freq_item.rda")

Rules[which(!is.na(min_freq_item)),]










