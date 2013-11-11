Principal Component Analysis ,Exploratory Factor Analysis and Multidimensinal scaling with zip data
--------------------------------------------------------------------
We report the results in this article to conduct Principal Component Analysis and Factor Analysis with the zip dataset. The main parts are as following:

1.  Data description and missing values imputation
2.  Principal component analysis
3.  Factor analysis
4.  Multidimensional scaling

### Data Description
First, we convert the data format to rda in R with the package 'sas7bdat' and save it as 'zip1.rda' in the working path.
Then we can see the details of the description of the data set.



```r
setwd("F:/商业数据挖掘_光华/homework/2")

# library(sas7bdat) read data files
# zip1=read.sas7bdat('zipdemo1.sas7bdat') save(zip1,file='zip1.rda')
library(xtable)
##### load the data files
load("zip1.rda")
print(xtable(head(summary(zip1)[, 1:5])), type = "html")
```

<!-- html table generated in R 2.15.2 by xtable 1.7-1 package -->
<!-- Mon Nov 11 11:42:13 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH>    ZIPCODE </TH> <TH>    DMAWLTHT </TH> <TH>    INCMINDX </TH> <TH>    WEALTHRT </TH> <TH>    PRCWHTE </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> Min.   : 1001   </TD> <TD> Min.   :0.00   </TD> <TD> Min.   :  3.0   </TD> <TD> Min.   :0.00   </TD> <TD> Min.   :  0.0   </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> 1st Qu.:26071   </TD> <TD> 1st Qu.:2.00   </TD> <TD> 1st Qu.: 69.0   </TD> <TD> 1st Qu.:2.00   </TD> <TD> 1st Qu.: 86.0   </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> Median :49052   </TD> <TD> Median :4.00   </TD> <TD> Median : 84.0   </TD> <TD> Median :3.00   </TD> <TD> Median : 97.0   </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> Mean   :49135   </TD> <TD> Mean   :4.01   </TD> <TD> Mean   : 90.3   </TD> <TD> Mean   :3.63   </TD> <TD> Mean   : 87.9   </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> 3rd Qu.:71292   </TD> <TD> 3rd Qu.:6.00   </TD> <TD> 3rd Qu.:104.0   </TD> <TD> 3rd Qu.:5.00   </TD> <TD> 3rd Qu.: 99.0   </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> Max.   :99929   </TD> <TD> Max.   :9.00   </TD> <TD> Max.   :409.0   </TD> <TD> Max.   :9.00   </TD> <TD> Max.   :100.0   </TD> </TR>
   </TABLE>

```r
sum(!complete.cases(zip1))
```

[1] 908


Some NAs are included in some variables. So we simply conduct the multiple imputation with the 'mi' package in R. Since the it takes a long time to do such process, we just save the imputated dataset before. Here we also deplay the codes.

```r
##### conduct the multiple imputation require(mi) zip1_com = mi(zip1[,-1])
##### comp_zip1 = mi.data.frame(zip1_com,m=1) zip1_comp =
##### cbind(zip1[,1],comp_zip1) save(zip1_comp,file='zip1_comp.rda')
load("zip1_comp.rda")
```

```
## Warning: cannot open compressed file 'zip1_comp.rda', probable reason 'No
## such file or directory'
```

```
## Error: cannot open the connection
```

```r
summary(zip1_comp)
```

```
## Error: object 'zip1_comp' not found
```


### Principal Component Analysis
We conduct the principal component analysis using the function princomp. The correlation matrix is used for the estimation.


```r
#### Principal Component Analysis
dd = zip1_comp[, -1]
```

```
## Error: object 'zip1_comp' not found
```

```r
zip_prin = princomp(dd, cor = T, scores = T)
```

```
## Error: object 'dd' not found
```


Then we draw the scree plot to select the appropriate component number. There's a turning point at Comp.5 in the picture, so we simply choose 5 as the component number, which explains
74.98 percent of the variance.


```r
screeplot(zip_prin, type = "lines")
```

```
## Error: object 'zip_prin' not found
```


we summary the results for details of the components. Besides, if we follow the Kaiser Principle, we'd better to choose 6 as the component number, which explains 78.35 percent variance.


```r
summary(zip_prin)
```

```
## Error: object 'zip_prin' not found
```


In addition, we look into the loading matrix to achieve more information.

The first component is about the population age distribution against the wealth indicators. 

The second mainly describes the races distribution and wealth ratings against some other indicators like ages and NCDB. 

The third contains the comparison about the races structure and wealth indexes, ages, salary structures etc. 

Other components can also be analyzed in this way. 

To summary, the variables can be clustered into aspects about wealth, salary structures, social positions, human races structures, age structures etc.


```r
loadings(zip_prin)
```

```
## Error: object 'zip_prin' not found
```


### Factor Analysis
We conduct the Factor Analysis with the MLE method. To avoid the scaling problem, we center and scale the variables with zero mean and 1 sd.

Besides, we use the varimax principle to conduct the rotation.

```r
#### Factor Analysis
dd_std = scale(dd)
```

```
## Error: object 'dd' not found
```

```r
zip_fac = factanal(dd_std, factor = 10, rotation = "varimax", n.obs = nrow(dd), 
    control = list(trace = T))
```

```
## Error: object 'dd_std' not found
```

```r
zip_fac
```

```
## Error: object 'zip_fac' not found
```


#### Several ways to decide Number of Factors
We use the principle that every factor's variance's ratio should be larger than 1/p. Since the variance of the X (all variance) is p (32), we just need to guarantee that the factor's variance is larger than 1.


```r
flag = 1
i = 15
while (flag == 1) {
    zip_fac = factanal(dd_std, factor = i, rotation = "varimax", n.obs = nrow(dd))
    fac_var = apply(zip_fac$loadings, 2, function(x) return(sum(x^2)))
    cat("nFactor:", i, "\n", "Factor variance:", fac_var, "\n", "\n")
    if (all(fac_var >= 1)) 
        flag = 0
    i = i - 1
}
```

```
## Error: object 'dd_std' not found
```


Here 8 is the number of the factors that subject to the constraint.

Besides that, we can also use the eigenvalues to decide the number. According to the principal component analysis before, we can choose 6 as our number of factors, which accounts for 78.34 percent of the sum of eigenvalues. We can also make other choices due to the demanded cumulative proportion in the summary(zip_prin) before.

The parallel analysis is also a way to decide the number of factors, which is developed by Raiche, Riopel, and Blais. We the R package 'nFactors' fot the analysis.


```r
# Determine Number of Factors to Extract
library(nFactors)
ev <- eigen(cor(dd_std))  # get eigenvalues
```

```
## Error: object 'dd_std' not found
```

```r
ap <- parallel(subject = nrow(dd_std), var = ncol(dd_std), rep = 100, cent = 0.05)
```

```
## Error: object 'dd_std' not found
```

```r
nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea)
```

```
## Error: object 'ev' not found
```

```r
plotnScree(nS)
```

```
## Error: object 'nS' not found
```


The result gives 6 as the best number.

There are other methods provided in this package for the number decision process.

```r
nBartlett(x = ev$values, N = nrow(dd_std), alpha = 0.05, details = TRUE)
```

```
## Error: object 'ev' not found
```

```r
nBentler(x = ev$values, N = nrow(dd_std), alpha = 0.05, details = TRUE)
```

```
## Error: object 'ev' not found
```

```r
nCng(x = ev$values, model = "factors", details = TRUE)
```

```
## Error: object 'ev' not found
```


#### Factor Explainations
We choose 6 as the factor numbers. 

We then use the varimax and promax rotation method to achieve the loading matrix.


```r
#### Factor Analysis
dd_std = scale(dd)
```

```
## Error: object 'dd' not found
```

```r
zip_fac = factanal(dd_std, factor = 6, rotation = "varimax", n.obs = nrow(dd))
```

```
## Error: object 'dd_std' not found
```

```r
loadings(zip_fac)
```

```
## Error: object 'zip_fac' not found
```

```r
zip_fac = factanal(dd_std, factor = 6, rotation = "promax", n.obs = nrow(dd))
```

```
## Error: object 'dd_std' not found
```

```r
loadings(zip_fac)
```

```
## Error: object 'zip_fac' not found
```


From the loading matrix we can see the promax rotation are closely related to the varimax rotation in the top3 factors. We now look into more details of the loading matrix of the promax rotation.

The first factor has a high loading value at WEALTHRT, DMAWLTHT, INCMINDX and CEMI, they are mainly about the income index. And it has relative negative loading on PRCBLCK, PRCNCD1 and PRCRENT compared to the first aspect, which are mainly  related to the human races percent.

The second factor are mostly the PRCNCD3, PRCNCD10 aginst PRCNCD1 and PRCOWNO. They are mainly the descriptions of the %NCDB HH.

The third has a high loading on PRC500K, PRC200K, PRC100K and OOMEDHVL. They are mainly about the OOH Value.

The fourth factor has the highest loading on PRC65P, HHMEDAGE, PRC55P and a relative negative value on PRC3544, PRC4554. They are mostly descriptions of the age structure.

The fifth factor is highest loaded on PRCUN18, PRCTHRE and PRCHHFM. They have a closed relationship with the %HH value.

#### Multidimensional Scaling

We use the package MASS to do the multidimensional scaling with data zip2. For details, we use 1-correlation as the distance measure.


```r
library("MASS")
library(ggplot2)
load("zip2.rda")
```

```
## Error: cannot open the connection
```

```r
dist = 1 - cor(zip2[, -1])
```

```
## Error: object 'zip2' not found
```

```r

zip2_mds = isoMDS(dist)
```

```
## Error: default method not implemented for type 'closure'
```

```r
x = zip2_mds$points[, 1]
```

```
## Error: object 'zip2_mds' not found
```

```r
y = zip2_mds$points[, 2]
```

```
## Error: object 'zip2_mds' not found
```

```r
g = ggplot(data.frame(x, y), aes(x, y, label = colnames(zip2)[-1]))
```

```
## Error: object 'x' not found
```

```r
g + geom_point(shape = 16, size = 3, colour = "red") + geom_text(hjust = -0.1, 
    vjust = 0.5, alpha = 0.5, angle = 7)
```

```
## Error: object 'g' not found
```


From the figure, we can see some features are clustered with each other.

The CEMI, PRC25BA, OOHVI, PRCOOHV, ISPSA, WEALTHRT, DMAWLTHT, MEDSCHYR are very closed to each other. That means the wealth ,the education level and the social position are high correlated.

The PRC65P, HHMEDAGE, PRC55P are very closed. They are all descriptions of ages. We may infer that the householder's median age are upon 55.

The PRC200K, PRC100K, OOHVI are very closed to each other. 

There isn't a significant clustering sign of other features.





