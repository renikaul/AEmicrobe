
#useful libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(gplots)
library(bbmle)

library(stats)
library(rootSolve)
library(mvtnorm)
library(numDeriv)
##Functions used in analysis
##1) a function describing model, 2-parameter Weibull function
weib=function(lam,k,x){
  p<-1-exp(-(x/lam)**k)
  return(p)
}

#a function that calculates the negative log-likelihood based on a binomial distribution 
#using a single x,y value
LL1=function(lam,k,x,y){
  p<-1-exp(-(x/lam)^k)
  loglik <- -sum(dbinom(y,1,p, log=TRUE))
  return(loglik)
}
#using multiple x,y values
LL=function(param,data){
  param[1]->lam
  param[2]->k
  data[[1]]->x
  data[[2]]->y 
  p<-1-exp(-(x/lam)^k)
  loglik <- -sum(dbinom(y,1,p, log=TRUE))
  return(loglik)
}
#Plotting Contours
#image cmd, taking in a matrix



#Estimate the max. likelihood of binomial distribution
mle.estimate=function(data){
  initial<-list(lam=1,k=1)#initial parameters for MLE
  #Fit MLE, just change model name
  model<-mle2(minuslogl=LL1,start=initial, data=data,method='SANN', trace=TRUE)
  #Pull out point estimates for plotting
  l<-lam.est<-coef(model)[1]
  k<-k.est<-coef(model)[2]
  param<-c(lam.est,k.est)

  return(model)
}

#Produce contour plot of parameter estimates with 95% CI around estimate
con.plot=function(window,res,data,model,n,title, ring=1, surface=FALSE, same=TRUE, ylim=c(window[3],window[4]),xlim=c(window[1],window[2])){
  #name param
  min.lam<-window[1]
  max.lam<-window[2]
  min.k<-window[3]
  max.k<-window[4]
  lam.res<-res[1]
  k.res<-res[2]
  #grid size
  lam.seq<-seq(min.lam,max.lam, length.out=lam.res)
  k.seq<-seq(min.k,max.k, length.out=k.res)
  #make grid
  grid<-matrix(NA, nrow=lam.res, ncol=k.res)
  
  for(i in 1:lam.res){
    for(j in 1:k.res){
      grid[i,j]=LL(c(lam.seq[i],k.seq[j]),data)
    }
  }
  #k:y-axis, col
  #lam: x-axis, row
  
  #Calculate where to put 95% CI
  CI<--logLik(model)+1.96
  
if (surface==TRUE) {
  #change col to white to remove surface levels. 
  image(z=grid, x=lam.seq,y=k.seq, col=heat.colors(n))
  #Plot contour map with 95% CI
  #col= changes CI ring
  contour(x=lam.seq,y=k.seq,z=grid,levels=CI,add=TRUE,labcex=1.5, drawlabels = FALSE, nlevels=1, col=ring)
  title(main=title) 
} else {
  contour(x=lam.seq,y=k.seq,z=grid,levels=CI,add=same,ylim=ylim,xlim=xlim,labcex=1.5, drawlabels = FALSE, nlevels=1, col=ring)
}
  invisible(grid)
}

#Calculate range of inflection point
inflection.range=function(model, data){
  hessian(LL, coef(model),data=data)->hess
  solve(hess,diag(c(1,1)))->A
  rmvnorm(1000000,coef(model),A)->sample
  # the second derivative of the two parameter Weibull function is:
  inflect<- function(v){v[1]*((v[2]-1)/v[2])**(1/v[2])}
  points<- apply(sample,1,inflect)
  points <- points[!is.na(points)]
  qnorm(c(0.025,.975),mean=mean(points),sd=sd(points)) -> boundary
  return(boundary)
}

par.range=function(model, data){
  hessian(LL, coef(model),data=data)->hess
  solve(hess,diag(c(1,1)))->A
  rmvnorm(1000000,coef(model),A)->sample
  # return lambda and k values associated with the range
  qnorm(c(0.025,0.975), mean=mean(sample[,2]), sd=sd(sample[,2]))->k.boundry
  qnorm(c(0.025,0.975), mean=mean(sample[,1]), sd=sd(sample[,1]))->l.boundry
  range<-c(k.boundry, l.boundry)
  names(range)<-c('k.min','k.max','lam.min','lam.max')
  return(range)
}


persist.plot=function(y,x,z,model, conplot.out,title){
  #z range of initial population sizes in experiment
  l<- coef(model)[1] #fitted lambda value
  k<-coef(model)[2] #fitted k value
  min<-conplot.out[[2]] #vector of inflection points when k is at 1.001 
  max<- conplot.out[[3]] #vector of inflection points when k is at max value within 95CI
  
  prediction<-1-exp(-(z/l)^k)
  plot(y~jitter(x),col=1, xlim=c(0,10), ylim=c(0,1), main=title)
  points(z,y=prediction,col=2, type='l')
  
  for (i in 1:length(min)){
    line<-min[i]
    abline(v=line, lty=2)
  }
  
  for (i in 1:length(max)){
    line<-max[i]
    abline(v=line, lty=3)
  }
}

modelplot=function(start,end, x,y,model){
  #Plot data with fitted curve, shaded 95% CI of inflection point 
  plot(NA, xlim=c(-2,8), ylim=c(0,1),yaxt='n', xaxt='n', xlab='', ylab='')
  
  
  #shade range of inflection point
  s<-log(start)
  e<-log(end)
  shade<-seq(s,e, by=0.025)
  for(i in 1:length(shade)){
    shade[i]->draw
    abline(v=draw, col='lightgrey')
  }
  
  log.x<-log(x)+1
  cells<-c(1,10,100,1000)
  tick<-log(cells)
  
  points(y~jitter(log.x),col=1)
  axis(2, at=c(0,1),las=1)
  mtext(2, text='Establishment', line=2.5)
  axis(1, at=tick, label=cells)
  mtext(1, text='Initial Population Size (# cells)', line=2.5)
  
  coef(model)[1]->l
  coef(model)[2]->k
  #plot fitted line
  z<-seq(0,3000, by=0.1)
  prediction<-1-exp(-(z/l)^k)
  lines(y=prediction,x=log(z),col=2, type='l',lwd=2)
}





