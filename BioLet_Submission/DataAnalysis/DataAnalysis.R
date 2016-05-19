#This document contains all code needed to preform analysis reported in: 'Experimental demonstration of an Allee effect in microbial populations'
#The document consists of 3 parts:
#1. All functions and packages
#2. Analysis with figures
#3. Analysis with figures for sensitivity analysis (in supplemental material)

###################
## Part 1 
###################
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
con.plot=function(window,res,data,model,n,title, ring=1, plot=TRUE, surface=FALSE, same=TRUE, ylim=c(window[3],window[4]),xlim=c(window[1],window[2])){
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
if (plot==TRUE){
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
  }
  invisible(list(lam=lam.seq, k=k.seq,LLgrid=grid))
}

#Calculate range of inflection point based on hessian
inflection.range.hess=function(model, data){
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


#Determine max/min of k and lambda based on the output of con.plot (test)
paramCI=function(test, model){
CItest<-c()
#go through each cell and ask if within CI
  for (i in 1:length(test$k)){
    temp3<-c()
    for (j in 1:length(test$lam)){
      ifelse(-logLik(model)+1.96 > test$LLgrid[j,i],1,0)-> temp #1 falls within interval
      c(k=test$k[i],lam=test$lam[j],interval=temp)->temp2
      rbind(temp3,temp2)->temp3
    }
    rbind(CItest, temp3)->CItest
  }
  #Only look at those within CI 
  CItest[which(CItest[,3]>0),]->CItest2
  c(k=coef(model)[2],k.min=min(CItest2[,1]), k.max=max(CItest2[,1]),lam=coef(model)[1], lam.min=min(CItest2[,2]), lam.max=max(CItest2[,2]))->out
return(list(CItest2,out))
}

#Calculate range of inflection point based on param values within CI, data=output from con.plot[[1]]
inflection.range=function(data){
  # the second derivative of the two parameter Weibull function is: v1:k, v2:lam
  inflect<- function(v){v[2]*((v[1]-1)/v[1])**(1/v[1])}
  points<- apply(data,1,inflect)
  points <- points[!is.na(points)]
  boundary<-c(min(points),max(points))
  return(boundary)
}


################
## Part 2
###############


#libraries not in scripts
library(reshape2) #for data summary
library(gplots) #for plots with CI
library(Hmisc) # for binomal CI calculation
library(waffle) # for dot plot 

##treat: (1) high carbon, (2) low carbon, (3)low carbon + predation
#load in all data
read.csv("completedata.csv", header=T)->d

#seperate data into 3 treatments
d[d$treat==1,c(1,2)]->highC
d[d$treat==2,c(1,2)]->lowC
d[d$treat==3,c(1,2)]->lowCP


#create summary table (# established per initial population by treatment)
##########################

#prep data by grouping in experimental treatments of "Initial" population size
highC.group<- group_by(highC, x)
lowC.group<- group_by(lowC, x)
lowCP.group<- group_by(lowCP, x)

#summary x=initial pop, N= num replicates, est= number established, average, lower and upper binomial confidence interval based on number of replicates 
highC.stat<-summarise(highC.group, N=length(y), est=sum(y), bin.ave=binconf(x=sum(y), n=36)[1],bin.lci=binconf(x=sum(y), n=36)[2],bin.uci=binconf(x=sum(y), n=36)[3])
lowC.stat<-summarise(lowC.group, N=length(y), est=sum(y), bin.ave=binconf(x=sum(y), n=12)[1],bin.lci=binconf(x=sum(y), n=12)[2],bin.uci=binconf(x=sum(y), n=12)[3])
lowCP.stat<-summarise(lowCP.group, N=length(y), est=sum(y), bin.ave=binconf(x=sum(y), n=12)[1],bin.lci=binconf(x=sum(y), n=12)[2],bin.uci=binconf(x=sum(y), n=12)[3])
###### 

#Fit data to Weibul
##########################

#calculate MLE surface
highC.fit<-mle.estimate(data=highC)
lowC.fit<-mle.estimate(data=lowC)
lowCP.fit<-mle.estimate(data=lowCP)

#get intervals on k and lambda
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=highC, model=highC.fit, title='',plot=FALSE, same=FALSE, ring='black',ylim=c(0,12),xlim=c(0,10))->highC.LL
con.plot(window=c(0,10,0,25), res=c(500,1250), n=10, data=lowC, model=lowC.fit, title='',plot=FALSE, same=FALSE, ring='blue')->lowC.LL
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=lowCP, model=lowCP.fit, title='',plot=FALSE, same=FALSE, ring='red')->lowCP.LL

paramCI(highC.LL, highC.fit)->highC.param
paramCI(lowC.LL, lowC.fit)->lowC.param
paramCI(lowCP.LL, lowCP.fit)->lowCP.param

signif(highC.param[[2]], digits=3)->highC.CI
signif(lowC.param[[2]], digits=3)->lowC.CI
signif(lowCP.param[[2]], digits=3)->lowCP.CI
#####
#FIGURES
########



#Create plot of fitted Weibul curves with data
#################
#pdf('fig1.pdf', width=7, height=5)
postscript("fig1.eps",title='Fig1',horizontal = FALSE, height=5,width=7)
par(oma=c(0,0,0,0), mar=c(5,3,1,1))
x<-seq(0,30000, by=0.05)
plot(NA, ylim=c(-0.2,1.2), xlim=c(0,10.5), xaxt='n', yaxt='n',xlab=expression(Inoculum ~ Size ~ (cells~mL^{-1})), ylab='')
tick3<-c(1,10,100,1000,10000)
axis(1, at=log(tick3), labels=tick3)
axis(2, at=c(0,0.5,1), labels=c('0','0.5','1.0'), hadj=0.5,padj=0.5, las=1)
mtext('Probability of Establishment', side=2, line=1.75)
#Add fitted line with inflection point CI
######
#calculate ranges on inflection points
inflection.range(highC.param[[1]])->highC.range
inflection.range(lowC.param[[1]])->lowC.range
inflection.range(lowCP.param[[1]])->lowCP.range

#shade range of inflection point
shade<-seq(lowCP.range[1],lowCP.range[2], by=0.01)
for(i in 1:length(shade)){
  shade[i]->draw
  abline(v=draw, col='lightpink')
}

abline(v=lowCP.range[1], lty=3, col='lightpink', lwd=2)
abline(v=lowCP.range[2], lty=3, col='lightpink', lwd=2)

shade<-seq(lowC.range[1],lowC.range[2], by=0.01)
for(i in 1:length(shade)){
  shade[i]->draw
  abline(v=draw, col='lightblue1')
}
abline(v=lowC.range[1], lty=3, col='lightblue',lwd=2)
abline(v=lowC.range[2], lty=3, col='lightblue',lwd=2)

shade<-seq(highC.range[1],highC.range[2], by=0.01)
for(i in 1:length(shade)){
  shade[i]->draw
  abline(v=draw, col='gray85')
}
abline(v=highC.range[1], lty=3, col='grey', lwd=2)
abline(v=highC.range[2], lty=3, col='grey', lwd=2)

##HighC
#plot weibull line with given parameters
k<- coef(highC.fit)[2]
l<- coef(highC.fit)[1]
lgy<-1-exp(-(log(x)/l)^k)
lines(lgy~log(x)) 
#add inflection point
l*(((k-1)/k)^(1/k))->xi
yi<-1-exp(-((xi)/l)^k)
points(yi~xi, pch=17)

#lowC
#plot weibull line with given parameters
k<- coef(lowC.fit)[2]
l<- coef(lowC.fit)[1]
lgy<-1-exp(-(log(x)/l)^k)
lines(lgy~log(x), col='blue') 
#add inflection point
l*(((k-1)/k)^(1/k))->zi
yi<-1-exp(-((zi)/l)^k)
points(yi~zi, pch=15, col='blue')

#lowCP
#plot weibull line with given parameters
k<- coef(lowCP.fit)[2]
l<- coef(lowCP.fit)[1]
lgy<-1-exp(-(log(x)/l)^k)
lines(lgy~log(x), col='red') 
#add inflection point
l*(((k-1)/k)^(1/k))->zpi
yi<-1-exp(-((zpi)/l)^k)
points(yi~zpi, pch=19, col='red')
#####

#Add Experimental data 
#########################################

#Add ave data with error
plotCI(x=highC.stat$x, y=highC.stat$bin.ave, ui=highC.stat$bin.uci,gap=0, li=highC.stat$bin.lci, col='black',pch=2, add=TRUE)
plotCI(x=lowC.stat$x, y=lowC.stat$bin.ave, ui=lowC.stat$bin.uci,gap=0, li=lowC.stat$bin.lci, col='blue',pch=0, add=TRUE)
plotCI(x=(lowCP.stat$x+0.1), y=lowCP.stat$bin.ave, ui=lowCP.stat$bin.uci,gap=0, li=lowCP.stat$bin.lci, col='red',pch=1, add=TRUE)

#Add 0/1 data
#points should be spaced 0.02 apart on y-axis

#try adding numbers so less cluttered
#text(x=highC.stat$x, y=1.035, labels=highC.stat$est)
#text(x=highC.stat$x, y=-0.035, labels=(36-highC.stat$est))
# points(x=highC.stat$x, y=rep(1.035, length(highC.stat$x)),pch=2, cex=3)
# points(x=highC.stat$x, y=rep(-0.035, length(highC.stat$x)),pch=2, cex=3)
# 
# text(x=lowC.stat$x+0.04, y=1.035, labels=lowC.stat$est)
# text(x=lowC.stat$x+0.04, y=-0.035, labels=(12-lowC.stat$est))
# points(x=lowC.stat$x+0.04, y=rep(1.038, length(lowC.stat$x)),pch=0, cex=3)
# points(x=lowC.stat$x+0.04, y=rep(-0.035, length(lowC.stat$x)),pch=0, cex=3)
# 
# 
# text(x=lowCP.stat$x+0.04, y=1.09, labels=lowCP.stat$est)
# text(x=lowCP.stat$x+0.04, y=-0.09, labels=(12-lowCP.stat$est))
# points(x=lowCP.stat$x+0.04, y=rep(1.09, length(lowCP.stat$x)),pch=1, cex=3)
# points(x=lowCP.stat$x+0.04, y=rep(-0.09, length(lowCP.stat$x)),pch=1, cex=3)

#brute force waffle plot

xspace<-0.072
yspace<-0.02

#highC
########
#established
xseq<-unique(highC.stat$x)
eseq<-highC.stat$est
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<-1.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y+(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=17, col='black')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=17, col='black')
  }
}

#failed
xseq<-unique(highC.stat$x)
eseq<-(highC.stat$N-highC.stat$est)
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<- -.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y-(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=17, col='black')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=17, col='black')
  }
}
#####

#lowC
#######
xseq<-unique(lowC.stat$x)
eseq<-lowC.stat$est
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<-1.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y+(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=15, col='blue')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=15, col='blue')
  }
  
}

#failed
xseq<-unique(lowC.stat$x)
eseq<-(lowC.stat$N-lowC.stat$est)
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<- -.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y-(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=15, col='blue')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=15, col='blue')
  }
}



#lowCP
#######
xseq<-unique(lowCP.stat$x)+.15
eseq<-lowCP.stat$est
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<-1.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y+(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=19, col='red')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=19, col='red')
  }
  
}

#failed
xseq<-unique(lowCP.stat$x)+0.15
eseq<-(lowCP.stat$N-lowCP.stat$est)
for (i in 1:length(xseq)){
  #starting corner
  x<-xseq[i]
  y<- -.035
  #build grid
  xnumcol<-c(1:4)
  ynumrow<-c(1:10)
  xcol<-x+(xnumcol*xspace)
  yrow<- y-(ynumrow*yspace)
  dpoint<-cbind(c(rep(xcol[1],10),rep(xcol[2],10),rep(xcol[3],10),rep(xcol[4],10)), c(rep(yrow,4)))
  #trim to number established
  if (eseq[i]==1) {
    dpoint[1,]->wpoint
    points(x=wpoint[1], y=wpoint[2], cex=0.4, pch=19, col='red')
  }
  if (eseq[i]>1){
    dpoint[1:eseq[i],]->wpoint
    points(wpoint, cex=0.4, pch=19, col='red')
  }
}

#add legend
legend(x=6, y=0.45, legend=c('High', 'Low', 'Low+Predation'),title='Carbon Resource', col=c(1,4,2), pch=c(17,15,19), lty=1, bty='n')

box()
dev.off()

###############

#Create plot of likelihood surface
###########################
pdf('fig2.pdf', width=7, height=5)
postscript("fig2.eps",title='Fig2',height=5,width=7, horizontal = FALSE)

par(oma=c(0,0,0,0), mar=c(5,4,1,1))
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=highC, model=highC.fit, title='', same=FALSE, ring='black',ylim=c(0,12),xlim=c(0,10))
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=lowC, model=lowC.fit, title='', same=TRUE, ring='blue')
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=lowCP, model=lowCP.fit, title='', same=TRUE, ring='red')
mtext('k',2, line=2.75)
#mtext(expression(paste( lambda, ~~(log[10]~(cells~mL^-1)))), 1, line=3)
mtext(expression(lambda), 1, line=3)

#mtext('Parameter Estimate', 3, line=1)
#add points
points(coef(highC.fit)[2]~coef(highC.fit)[1], pch=17, cex=0.5)
points(coef(lowC.fit)[2]~coef(lowC.fit)[1], pch=15, cex=0.5, col='blue')
points(coef(lowCP.fit)[2]~coef(lowCP.fit)[1], pch=19, cex=0.5, col='red')
#add legend
#legend(x=8,y=8, legend=c('High', 'Low', 'Low+Predation'),title='Carbon Resouce', col=c(1,4,2), pch=c(17,15,19), lty=1, bty='n')
#add hypothesis test
abline(h=1, lty=2)
abline(v=0)
text(x=.2, y=6, labels=expression(paste("Critical Population\n Threshold = 1 cell",~mL^-1, sep='')), srt=90)
text(x=4.5, y=1.25, labels="Allee effect detected")

#Add labels
text(x=1,y=3,labels='High\nCarbon')
text(x=2.45,y=10, labels='Low\nCarbon')
text(x=6,y=6, labels='Low Carbon\n with Predation')


#this should be updated manually if change in data or analysis
#Add parameter estimates
text(x=5.5,y=11.55, labels=expression(hat(k)), cex=1)
text(x=8.25,y=11.55, labels=expression(hat(lambda)), cex=1)
legend(x=4.25,y=11.5, legend=c('2.57 (1.72,3.65)     1.92 (1.64,2.12)',
                               '11.7 (5.42,24.2)     3.19 (3.01,3.43)',
                               ' 5.8  (3.78,8.75)     4.64 (4.33,4.99)'),
       col=c(1,4,2), pch=c(17,15,19), lty=1, bty='n') 
#####
dev.off()


################
# Part 3:Sensitivity analysis: Drop one out
######################################
# Scheme: sequentially drop one replicate/line from data and fit parameters. The estimated parameter space will then be plotted over the likelihood 
#surface. This should result in a cloud of likelihood intervals along with the estimate from the entire data set. This analysis will identifiy any single case that has undo influence.


#load model functions
source("script/mlefunctions2.R")
##treat: (1) high carbon, (2) low carbon, (3)low carbon + predation
#load in all data
read.csv("completedata.csv", header=T)->d

tiff('pointwise.tiff', width=7, height=5, units="in", compression='none', res=100)
par(oma=c(0,0,0,0), mar=c(5,4,1,1))

#Pull out treatment type
#Treatment 1
###############
#high Carbon
d[d$treat==1,]->highC

goodfit<-mle.estimate(data=highC)
con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=highC, model=goodfit, title='', same=FALSE, ring='black',ylim=c(0,12),xlim=c(0,10))
mtext('k',2, line=2.75)
#mtext(expression(paste( lambda, ~~(log[10]~(cells~mL^-1)))), 1, line=3)
mtext(expression(lambda), 1, line=3)

#add legend
#legend(x=8,y=8, legend=c('High', 'Low', 'Low+Predation'),title='Carbon Resouce', col=c(1,4,2), pch=c(17,15,19), lty=1, bty='n')
#add hypothesis test
abline(h=1, lty=2)
abline(v=0)
text(x=.2, y=6, labels=expression(paste("Critical Population\n Threshold = 1 cell",~mL^-1, sep='')), srt=90)
text(x=4.5, y=1.25, labels="Allee effect detected")

#Add labels
text(x=1,y=3,labels='High\nCarbon')
text(x=2.45,y=10, labels='Low\nCarbon')
text(x=6,y=6, labels='Low Carbon\n with Predation')

#Remove single case
for ( j in 1: nrow(highC)){
  highC[-j,]->temp
  #Fit model to temp data
  tempfit<-mle.estimate(data=temp)
  con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=temp, model=tempfit, title='', same=TRUE, ring='grey')
  points(coef(tempfit)[2]~coef(tempfit)[1], pch=17, cex=0.5, col='grey')
  print(j)
}

con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=highC, model=goodfit, title='', same=TRUE, ring='black')
points(coef(goodfit)[2]~coef(goodfit)[1], pch=17, cex=0.5)

#Treatment 2
###############
#low Carbon
d[d$treat==2,]->lowC

poorfit<-mle.estimate(data=lowC)
#Remove single case
for ( j in 1: nrow(lowC)){
  lowC[-j,]->temp
  #Fit model to temp data
  tempfit<-mle.estimate(data=temp)
  con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=temp, model=tempfit, title='', same=TRUE, ring='lightblue1')
  points(coef(tempfit)[2]~coef(tempfit)[1], pch=15, cex=0.5, col='lightblue1')
  print(j)
}

con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=lowC, model=poorfit, title='', same=TRUE, ring='blue')
points(coef(poorfit)[2]~coef(poorfit)[1], pch=15, cex=0.5, col='blue')

#Treatment 3
###############
#low Carbon+Pred
d[d$treat==3,]->lowCP

predfit<-mle.estimate(data=lowCP)
#Remove single case
for ( j in 1: nrow(lowCP)){
  lowCP[-j,]->temp
  #Fit model to temp data 
  tempfit<-mle.estimate(data=temp)
  con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=temp, model=tempfit, title='', same=TRUE, ring='mistyrose1')
  points(coef(tempfit)[2]~coef(tempfit)[1], pch=19, cex=0.5, col='mistyrose1')
  print(j)
}

con.plot(window=c(0,10,0,13), res=c(500,500), n=10, data=lowCP, model=predfit, title='', same=TRUE, ring='red')
points(coef(predfit)[2]~coef(predfit)[1], pch=19, cex=0.5, col='red')

dev.off()