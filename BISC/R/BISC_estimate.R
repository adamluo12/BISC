#' Transcriptional bursting estimation with Stan
#'
#' @export
#' @param data read count input data matrix (row: genes; column: cells).
#' @param model estimation model, defaults to "PB-trend". BISC includes three model options: PB model, PB-trend model and ZIPB-trend model
#' @param iter number of MCMC interations, defaults to 4000
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#'

BISC_estimate=function(data,model="PB-trend",iter=4000){
  sum=colSums(count)
  count=count[,sum!=0]
  lib.size<-colSums(count)
  lib.size=lib.size/mean(lib.size)
  N=dim(count)[1]
  K=dim(count)[2]
  splat=splatEstimate(count)
  bcv.df=getParam(splat,"bcv.df")
  drop.x0=getParam(splat,"dropout.mid")
  drop.tau=getParam(splat,"dropout.shape")
  disps <- edgeR::estimateDisp(count)
  logcpm=edgeR::aveLogCPM(count)
  data_gam=data.frame(logcpm,disps)
  formula <- gam(disps~s(logcpm),data=data_gam)
  bcv=matrix(rep(1,ncol(count)*nrow(count)),ncol=ncol(count))
  for (c in 1:ncol(count)) {
    bcv[,c] <- predict(formula,edgeR::cpm(count,log=T,prior.counts=1)[,c])
  }
  if(bcv.df==Inf){
    bcv=bcv
  }else{
    bcv <- bcv*sqrt(bcv.df / rchisq(dim(count)[1], df = bcv.df))
  }
  if(model=="PB"){
    fit=rstan::sampling(stanmodels$bhm,data=list(N=N,K=K,y=count,l=as.numeric(lib.size)),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","s"),save_warmup=FALSE)
  }
  if(model=="PB-trend"){
    fit=rstan::sampling(stanmodels$bhm1,data=list(N=N,K=K,y=count,l=as.numeric(lib.size),bcv=bcv),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","s"),save_warmup=FALSE)
  }
  if(model=="ZIPB-trend"){
    fit=rstan::sampling(stanmodels$bhm2,data=list(N=N,K=K,y=count,l=as.numeric(lib.size),bcv,tau=drop.tau,x0=drop.x0),chains=1,iter =5000,control = list(adapt_delta = 0.99),  pars=c("kon","koff","p","s"),save_warmup=FALSE)
  }
  result=data.matrix(fit)
  return(result)
}