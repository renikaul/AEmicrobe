library(adaptivetau)

#function to describe the transition rates under exponential growth
exp.transrates = function(x, par, t){
  #par=c(lambda, mu)
  #x=prey
  rates=c(
    #birth
    (par[1]*x),
    #death
    (par[2]*x)
  )
  return(rates)
}

#function to describe the transition rates without predation
demo.transrates = function(x, par, t){
  #par=c(lambda, theta, mu)
  #x=prey
  rates=c(
    #birth
    ((par[1]*x^2)/(par[2]+x)),
    #death
    (par[3]*x)
  )
  return(rates)
}

#function to describe the transition rates with predation
predator.transrates = function(x, par, t){
  #par=c(lambda, theta, mu, alpha, handling, predator)
  #x=prey
  rates=c(
    #birth
    (par[1]*x*(x/(par[2]+x))),
    #death
    ((par[3]*x)+(par[6]*par[4]*x)/(1+par[4]*par[5]*x))
  )
  return(rates)
}

#function to run an "experiment" with the same initial population sizes as done before
experiment<-function(reps, transrates, par, tmax){
  
  #transition matrix
  M=matrix(c(1, #host birth
             -1 #host natural death
  ),ncol=2)
  
  sim.out<-c()
  for (j in 1:reps)
  { j-> sim
    
    #Initial conditions of state variables
    Y<-c(1,2,4,8,16,32,64,128,256,512,1024,2048)
    max.out<-c()
    for (i in 1:length(Y)){
      Y[i]-> Y0
      
      stoch<-ssa.adaptivetau(Y0, M , transrates, par, tmax)
      cbind(Y0,theta=par[2], max=max(stoch[,2]), min=min(stoch[,2]))->new.max
      rbind(max.out,new.max)->max.out
    }
    cbind(sim, max.out)->sim.new
    rbind(sim.out, sim.new)->sim.out
  }
  
  ifelse(sim.out[,5]==0,0,1)->persist
  cbind(sim.out, persist)-> sim.out
  return(sim.out)
}


#function to run many experiments and fit to weibul
fit.demo.allee.simulation=function(lambda=2, mu=1, theta.series=3, reps=12, transrates){
#lambda= 2#per capita birth
#theta.series=seq(from=1,to=100, by=5)#Allee strength
#mu= 1#per capita death
sim.data<-c()
est.data<-c()

#loop through different values of theta (simulation and fitting)
for (i in 1:length(theta.series)){
  theta.series[i]->theta
  par<-c(lambda, theta, mu)
#Stochastic simulation of experiment
  #transition matrix
  M=matrix(c(1, #host birth
             -1 #host natural death 
             ),ncol=2)
  
  sim.out<-c()
  for (j in 1:reps){ 
    j-> sim
    #Initial conditions of state variables
    Y<-c(1,2,4,8,16,32,64,128,256,512,1024,2048)
    max.out<-c()
    for (i in 1:length(Y)){
      Y[i]-> Y0
      stoch<-ssa.adaptivetau(Y0, M , transrates, par, 200)
      cbind(Y0,theta=par[2], max=max(stoch[,2]), min=min(stoch[,2]))->new.max
      rbind(max.out,new.max)->max.out
    }
    cbind(sim, max.out)->sim.new
    rbind(sim.out, sim.new)->sim.out  
}
  ifelse(sim.out[,5]==0,0,1)->persist
#return stochastic simulation
  cbind(sim.out, persist)-> sim.temp
#save stochastic data
  rbind(sim.data, sim.temp)->sim.data
  #fit model to current theta value
  sim.temp[,6]->y
  log(sim.temp[,2]+1)->x 
  temp.data<-list(x=x,y=y)
  #calculate MLE surface
  fitted<-mle.estimate(data=temp.data, n=1000, res=c(500,500),window=c(0,8,0,8), plot.surface = FALSE)
  par.range(fitted, temp.data)->par.interval
  
  #save fitted param
  c(theta=theta,coef(fitted)[1],coef(fitted)[2], LL=logLik(fitted), par.interval)-> temp.est
  rbind(est.data, temp.est)->est.data
}
return(list(sim.data, est.data))
}

fit.pred.allee.simulation=function(lambda=2, mu=1, attack=1, handling=2.35, predator=10,theta.series=3, reps=12, transrates){
  #param=c(lambda, mu, attack, handling,predator)

sim.pred<-c()
est.pred<-c()

for (i in 1:length(theta.series)){
  theta.series[i]->theta
  pred.par=c(lambda, theta, mu, attack, handling,predator)
  #Stochastic simulation of experiment
  #transition matrix
  M=matrix(c(1, #host birth
             -1 #host natural death 
  ),ncol=2)
  
  sim.out<-c()
      for (j in 1:reps){ 
        j-> sim
        #Initial conditions of state variables
        Y<-c(1,2,4,8,16,32,64,128,256,512,1024,2048)
        max.out<-c()
            for (i in 1:length(Y)){
              Y[i]-> Y0
              stoch<-ssa.adaptivetau(Y0, M , transrates, pred.par, 200)
              cbind(Y0,theta=par[2], max=max(stoch[,2]), min=min(stoch[,2]))->new.max
              rbind(max.out,new.max)->max.out
            }
        cbind(sim, max.out)->sim.new
        rbind(sim.out, sim.new)->sim.out  
      }
  ifelse(sim.out[,5]==0,0,1)->persist
  #return stochastic simulation
  cbind(sim.out, persist)-> sim.temp
  #save data
  rbind(sim.pred, sim.temp)->sim.pred
  
  #fit model to current theta value
  sim.pred[,6]->y
  log(sim.pred[,2]+1)->x 
  temp.data<-list(x=x,y=y)
  #calculate MLE surface
  fitted<-mle.estimate(data=temp.data, n=1000, res=c(500,500),window=c(0,8,0,8))
  par.range(fitted, temp.data)->par.interval
  
  #save fitted param
  c(theta=theta,coef(fitted)[1],coef(fitted)[2], LL=logLik(fitted), par.interval)-> temp.est
  rbind(est.pred, temp.est)->est.pred
}

return(list(sim.pred, est.pred))
}

fit.exp.simualtion=function(simulations=100, exp.par=c(2,1),replicates=12, rates=exp.transrates){
  #lambda= 2#per capita birth
  #mu= 1#per capita death
  #exp.par<-c(lambda, mu)
  
  sim.exp<-c()
  est.exp<-c()
  for (i in 1:simulations){
    i-> sim
    experiment(reps=replicates, transrates=rates, par=exp.par, tmax=200)->sim.temp
    #save data
    rbind(sim.exp, sim.temp)->sim.exp
    #fit model to current theta value
    sim.temp[,6]->y
    log(sim.temp[,2]+1)->x 
    temp.data<-list(x=x,y=y)
    #calculate MLE surface
    fitted<-mle.estimate(data=temp.data, n=1000, res=c(500,500),window=c(0,8,0,8), plot.surface=FALSE)
    par.range(fitted, temp.data)->par.interval
    
    #save fitted param
    c(sim=i,coef(fitted)[1],coef(fitted)[2], LL=logLik(fitted), par.interval)-> temp.est
    rbind(est.exp, temp.est)->est.exp
  }
  return(list(sim.exp, est.exp))

}