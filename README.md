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
BISC can take gene expression read count data, to minimize unwanted heterogeneity, we recommended data generated from cells of the same population. BISC has three mdoelling options, 1) Poisson-Beta model (PB); 2) Poisson-Beta model with modification to enforce a mean-variance trend (PB-trend); 3) Poisson-Beta model with modifications to enforce a mean-variance trend and include dropout events (ZIPB-trend). 

## Install BISC
```r
install_github("adamluo12/BISC/BISC")
```
## Running BISC
BISC takes the read count input values.

## Examples
(1) Load read count data
```r
> count[1:10,1:5]
    NA19098.r3.D02 NA19101.r1.C10 NA19101.r3.B02 NA19239.r3.B02 NA19098.r2.G06
G1               1              0              0              0              0
G2               0              0              0              0              0
G3               1              0              4              3              0
G4               0              1              0              0              0
G5               0              0              0              0              0
G6               0              0              0              0              0
G7               0              0              0              0              0
G8               2              3              0              0              2
G9               4              4              2              3              4
G10              0              0              0              0              0
N=dim(count)[1]
K=dim(count)[2]
```
(2) BISC estimation with 3 models
### PB model
```r
result=BISC_estimate(data=count,model="PB",iter=4000)
```
### PB-trend model
```r
result=BISC_estimate(data=count,model="PB-trend",iter=4000)
```
### ZIPB-trend model
```r
result=BISC_estimate(data=count,model="ZIPB-trend",iter=4000)
```
(3) BISC output: we used the mean values of the posterior estimates of kon, koff and s.
```r
final=apply(result,2,mean)
kon_est=final[1:N]
> kon_est[1:10]
  kon[1]   kon[2]   kon[3]   kon[4]   kon[5]   kon[6]   kon[7]   kon[8]   kon[9]  kon[10] 
1.326592 2.376357 1.876852 1.688580 2.454198 3.058243 1.896000 3.509753 1.847244 2.269638 
koff_est=final[N+1:2*N]
> koff_est[1:10]
 koff[1]  koff[2]  koff[3]  koff[4]  koff[5]  koff[6]  koff[7]  koff[8]  koff[9] koff[10] 
2.959785 4.043145 3.044111 3.532146 2.557586 3.190786 4.101236 2.378761 2.911859 2.846453 
s_est=final[(2*N+1):3*N]
> s_est[1:10]
    s[1]     s[2]     s[3]     s[4]     s[5]     s[6]     s[7]     s[8]     s[9]    s[10] 
39.46621 34.29013 36.84154 41.46866 36.94047 36.57172 38.43477 38.48573 37.84727 39.22281 
```
### Differential Bursting analysis 
To compare bursting kinetics between two studying groups A and B for each gene, BISC formulates a statistical testing framework based on the posterior MCMC samples of bursting parameters (BISC estimation step above).
```r
DB_result=DB(data1=BISC_A,data2=BISC_B,frequency=0,size=0,log2=T)
```
### HeatMap, MA plot and Volcano plot
```r
library(edgeR)
DB_plot(file1,file2,file3,DBdata=DB_result,main="DB frequency Group A vs. B")
```










