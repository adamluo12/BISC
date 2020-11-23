# BISC:Transcriptional-Bursting-Inference-using-Single-Cell-Transcriptoic-data

## Author: Xizhi Luo, Fei Qin, Feifei Xiao, Guoshuai Cai

## Description
Gene expression in mammalian cells is inherently stochastic that mRNA molecules are synthesized in discrete bursts. Although the advent of single cell transcriptomic technology provides great opportunities to explore the phenomenon of transcriptional bursting, current Beta-Poisson framework greatly suffers from substantial technical noise that leads to false estimations and conclusions. Statistical methods that account for complex nature of single cell transcriptomic data are needed to accurately reveal the heterogeneity of gene expression and infer kinetics parameters.
Here, we develop a Bayesian hierarchical framework, BISC, to study the stochastic gene expression kinetics. Gamma-Poisson model was employed to accommodate overdispersion of read counts, which was dynamically characterized by fitting the mean-variance relationship. The reliable estimation of the dispersion parameter is essential for capturing the expression heterogeneity, which improve the estimation of kinetics parameter. Also, we proposed a differential bursting analysis framework to identify heterogeneous bursting kinetics under different studying conditions.

## Required Packages
```r
library(rstan)
library(mgcv)
library(edgeR)
library(splatter)
```

## Running BISC
BISC can take gene expression read count data, to minimize unwanted heterogeneity, we recommended data generated from cells of the same population. BISC has three mdoelling options, 1) Poisson-Beta model (Bursting); 2) Poisson-Beta model with modification to enforce a mean-variance trend (Bursting + BCV); 3) Poisson-Beta model with modifications to enforce a mean-variance trend and include dropout events (Bursting+BCV+Dropout). 

## Examples
(1) Load read count data
```r
load("read_count.RData")
```
(2) Estimate library size factor
```r
sum=colSums(count)
count=count[,sum!=0]
lib.size<-colSums(count)
lib.size=lib.size/mean(lib.size)
N=dim(count)[1]
K=dim(count)[2]
```
(3) Estimate Dropout probability related parameters
```r
splat=splatEstimate(count)
bcv.df=getParam(splat,"bcv.df")
drop.x0=getParam(splat,"dropout.mid")
drop.tau=getParam(splat,"dropout.shape")
```
(4) Estimate log2CPM and raw dispersion, and fit the mean-variance trend
```r
disps <- edgeR::estimateDisp(count)
logcpm=edgeR::aveLogCPM(count)
data_gam=data.frame(logcpm,disps)
formula <- gam(disps~s(logcpm),data=data_gam)

```
(5) Caliberated BCV estimates
```r
bcv=matrix(rep(1,ncol(count)*nrow(count)),ncol=ncol(count))
for (c in 1:ncol(count)) {
  bcv[,c] <- predict(formula,edgeR::cpm(count,log=T,prior.counts=1)[,c])
}
if(bcv.df==Inf){
  bcv=bcv
}else{
  bcv <- bcv*sqrt(bcv.df / rchisq(dim(count)[1], df = bcv.df))
}
```
(6) Running rstan
### Bursting 
```r
fit=stan(file="bhm.stan",data=list(N=N,K=K,y=count,l=as.numeric(lib.size)),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","s"),save_warmup=FALSE)
```
### Bursting + BCV
```r
fit1=stan(file="bhm1.stan",data=list(N=N,K=K,y=count,l=as.numeric(lib.size),bcv=bcv),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","s"),save_warmup=FALSE)
```
### Bursting + BCV + Dropout
```r
fit2=stan(file="bhm2.stan",data=list(N=N,K=K,y=count,l=as.numeric(lib.size),bcv,tau=drop.tau,x0=drop.x0),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","p","s"),save_warmup=FALSE)
```









