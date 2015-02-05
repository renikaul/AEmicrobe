##Pre-processing functions

#####################################################
#write script to calculate summary statistics
#####################################################
ttothreshold=function(data,cat.type, n.cat, thres,convert){
#cat.type; how to sort out different types of data (ie. plate or read type)
#n.cat; number of levels within category
#thres: value to use to determine if a population is extablished
#factor to multiple read number to convert to hours.
  tto<-c()
  yield<-c()
  
for ( i in 1:n.cat){
  data[which(cat.type==i),]-> temp2
  na.omit(temp2)->temp2
  temp2[,-c(1:2)]-> temp3#look only at data columns
  
  #establishment threshold is set at thres
  apply(temp3,2,function(V){return(which(V>thres)[1])})->est
  #max yield 
  apply(temp3, 2, max)-> y
  
  #save info for n.cat
  rbind(tto, est)-> tto
  rbind(yield,y)-> yield
}

#time to establishment prep for melt function
plate<-c(1:n.cat)
rownames(tto)<-NULL
sample.estab<-cbind(plate,tto)
as.data.frame(sample.estab)->sample.estab
#groups wells by treatment. 
melt(sample.estab, id=c("plate"), variable.name='Well', value.name='est')-> est.sum
#current units is in num time points above threshold, need to convert to time 
(est.sum$est)*convert->hour #convert from time point to hours
round(hour,digits=2)->hour
cbind(est.sum,hour)->est.sum #final out put for time to threshold

#max yield prep for melt function
plate<-c(1:n.cat)
rownames(yield)<- NULL
as.data.frame(yield)-> yield
cbind(plate, yield)->max.out
#groups wells by treatment. 
melt(max.out, id=c("plate"), variable.name='Well', value.name='yield')-> max.sum
#score persistance based on max.sum
ifelse(max.sum$yield>thres, 1,0)->persist
cbind(max.sum, persist)->max.sum
max.sum$yield*max.sum$persist->new.yield #values of 0 put in for populations that do not establish (easy col to sort downstream)
cbind(max.sum, new.yield)->max.sum #final output for yield 

#combine the threshold and yield into one output
inner_join(max.sum, est.sum, by=c('Well','plate'))-> demo.sum

return(demo.sum)
}

##################################################
#stat calculations within summarise function
###################################################
#Write set of functions to use within the summarise function, summarise can handle very basic statistics, but not data containing NAs.
#These functions are col with NAs. 
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
