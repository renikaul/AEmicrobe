#This script will take raw H1 outputs, reshape the data and fit a Weibull model using MLE
#Required data files:
#
#
#
#
#

#Required packages


#Require Scripts
source('script/pr.reformat.R') #converts H1 data into easy to use formate
source('script/mlefunctions2.R') #MLE fitting functions
source('script/simfunctions.R') #simulates data with known AE strength

#Files produced:
#data files:

#plots

######################
#Reshape data files##
#####################

#Demographic data is split across 3 96 well plates
#Demographic trial data
read.csv('data/demo_pl1.csv',header=T)->d1
pr.reformat('data/demo_pl2.txt')->d2
pr.reformat('data/demo_pl3.txt')-> d3

###Need to remove uness columns
#add sampling info to data not recorded by plate reader
plate<-c(1)
cbind(plate,d1)->d1
plate<-c(2)
cbind(plate,d2)->d2
plate<-c(3)
cbind(plate, d3)->d3

rbind(d1[,-c(3,4,5)],d2[,-c(3,4)],d3[,-c(3,4)])->data
#reshape data
#md<-melt(d, id=c("plate","Measurment"), variable.name="Well", value.name="Growth")

#read in plate maps
read.csv('data/demoplatemap.csv', header=T)-> p3
#####################################################
#write script to tell me time to threshold
#####################################################
out<-c()


for ( i in 1:3){
  i-> pl
  data[which(data$plate==i),]-> temp2
  na.omit(temp2)->temp2
  temp2[,-c(1:2)]-> temp3
  #establishment threshold is set at OD>0.05 of blank correct read
  apply(temp3,2,function(V){return(which(V>0.05)[1])})->est
  
  #  apply(mask, 2, sum)-> est
  rbind(out, est)-> out
}

plate<-c(1:3)
rownames(out)<-NULL
sample.estab<-cbind(plate,out)
as.data.frame(sample.estab)->sample.estab
#current units is in num time points above threshold, need to convert to time 

#groups wells by treatment. 
melt(sample.estab, id=c("plate"), variable.name='Well', value.name='est')-> est.sum
(est.sum$est)*(5/60)->hour #convert from time point to hours
round(hour,digits=2)->hour
cbind(est.sum,hour)->est.sum

#write script to tell me max value in each population
#######################################################################
yield<-c()

for ( i in 1:3){
  data[which(data$plate==i),]-> temp2
  na.omit(temp2)->temp2
  apply(temp2[,-c(1:3)], 2, max)-> y
  rbind(yield,y)-> yield
}


#Make it pretty

plate<-c(1:3)
rownames(yield)<- NULL
as.data.frame(yield)-> yield
cbind(plate, yield)->max.out
#Summary statistics
#groups wells by treatment. 
melt(max.out, id=c("plate"), variable.name='Well', value.name='yield')-> max.sum

##########################################################################################
#Write script to tell me which populations established, based on max value
##########################################################################################

ifelse(max.sum$yield>0.05, 1,0)->persist
cbind(max.sum, persist)->max.sum
max.sum$yield*max.sum$persist->new.yield
cbind(max.sum, new.yield)->max.sum

###########################################
#Compile all of the data together and calculate summary statistics
######################################################################
inner_join(max.sum, est.sum, by=c('Well','plate'))-> demo.sum
inner_join(demo.sum,p3, by='Well')->demo.sum

demo.group<- group_by(demo.sum, Initial)

# Create a function that calculates 95% confidence intervals for the given
# data vector using a t-distribution
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}

#summarise can't handle NAs, so need to write external function. 
length.na<- function(data){
  new.data<-na.omit(data)
  len<-length(new.data)
  return(len)
}

mean.na<- function(data){
  new.data<-na.omit(data)
  ave<-mean(new.data)
  return(ave)
}

sd.na<- function(data){
  sdeviation<-sd(data, na.rm=TRUE)
  return(sdeviation)
}

demo.stat<-summarise(demo.group, p.N=length(persist), p.ave=mean(persist), p.sd=sd(persist), est.N=length.na(hour), est.ave=mean.na(hour), est.sd=sd.na(hour), y.N=length(new.yield), y.ave=mean(new.yield), y.sd=sd(new.yield))

demo.stat[which(demo.stat$Initial>0),]->demo.stat
(demo.stat$Initial/200)*75->adj.initial #density of initial population (cells/mL)
cbind(adj.initial, demo.stat)->demo.stat

