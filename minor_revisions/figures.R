#Reply to reviewers of minor edits; most of the comments require knowing the dectection threshold for the 96well and 384 well plate. 
#A serial dilution was made, added to 96 and 384 well plate. Read in plate reader than cells were quantified over at flowcyt

#read in plate reader info
setwd("/Users/renikaul/Dropbox/Vibrio/MicrobialAlleeEffect/AEmicrobe/ReplytoReviewers_JMDedits/minor_revisions/")

source("/Users/renikaul/Dropbox/Rscripts/pr.reformat.R")
source("/Users/renikaul/Dropbox/Rscripts/pr.broken.R")
#source("/Users/renikaul/Dropbox/Rscripts/comparePlot.R")

pr.broken("MSReviewer_response.txt")->W96
pr.broken("MSReviewer_response384.txt")->W384

#define replicates
W96s<-list(c("A1","B1"),
           c("A2","B2"),
           c("A3","B3"),
          
           c("A5","B5"),
           c("A6","B6"),
           c("A7","B7"),
           c("A8","B8"),
           c("A9","B9"),
           c("A10","B10"),
           c("A11","B11"),
           c("A12","B12")
           )
W96.od<-c()
for(i in 1:length(W96s)){
  set=W96[,which(colnames(W96) %in% W96s[[i]])]
  temp<-c(dilution=i,mean=apply(set,1,mean), sd=apply(set,1,sd))
  W96.od<-rbind(W96.od,temp)
}

W384s<-list(c("E9","F9","G9"),
           c("E10","F10","G10"),
           c("E11","F11","G11"),
         
           c("E13","F13","G13"),
           c("E14","F14","G14"),
           c("E15","F15","G15"),
           c("E16","F16","G16"),
           c("E17","F17","G17"),
           c("E18","F18","G18"),
           c("E19","F19","G19"),
           c("E8","F8","G8")
           )

W384.od<-c()
for(i in 1:length(W384s)){
  set=W384[,which(colnames(W384) %in% W384s[[i]])]
  temp<-c(dilution=i,mean=apply(set,1,mean), sd=apply(set,1,sd))
  W384.od<-rbind(W384.od,temp)
}

#remove empty well (Dilution #4)

library(gplots)
plotCI(y=W96.od[,2],x=W384.od[,2], uiw=W96.od[,3],ylab="96", xlab="384", ylim=c(0,0.5),xlim=c(0,1), gap=0)
plotCI(y=W96.od[,2],x=W384.od[,2], uiw=W384.od[,3],ylab="96", xlab="384", ylim=c(0,0.5),xlim=c(0,1), gap=0, add=TRUE, err=x)



#determine detection limit of 96 vs 384

#flow data
read.csv("countdata.csv", header=T)->cnt
#only interested in first 12
cnt[1:12,]->cnt
#add dilution factor
dil<-c(4000,8000,16000,32000,64000,128000)
cnt[,11]*dil->den.uL
#density of overnight culture
mean(den.uL*1000)->den.ml

#combine OD and cell count

#cell cnt of dilution series
dil2<-c(1,10,100,2000,4000,8000,16000,32000,64000,128000)
cell<-den.ml/dil2

#plot OD vs CFU for each plate
png("detection.png", width=480, height=480, units= "px")
plot((W96.od[1:10,2]-W96.od[11,2])~log10(cell), col="blue", ylim=c(0,.35), ylab="OD620", lwd=2, cex=1.25)
lines((W96.od[1:10,2]-W96.od[11,2])~log10(cell), col="blue", lwd=2)
lines((W384.od[1:10,2]-W384.od[11,2])~log10(cell), lwd=2)
points((W384.od[1:10,2]-W384.od[11,2])~log10(cell), pch=2, cex=1.25, lwd=2)
abline(h=0.25, lty=2)
text(y=0.25,x=6, label="Establishment\nThreshold" )
legend('topleft', bty='n', legend=c("96 well plate","384 well plate"), col=c("blue","black"), pch=c(21,2), lwd=2)
dev.off()

##Generation time/ Time to detection
dbtime=function(Nt, No, t){log(2)/(log(Nt/No)/t)}

#from one cell 
dbtime(Nt=10^7, No=1, t=96)
#from HC
dbtime(Nt=10^7, No=5, t=96)
dbtime(Nt=10^7, No=13, t=96)

#In HC
HC_No<-c(5,10,20,40,80,160,320)
HC_t<-c(32,32,32,31,30,29,27)

dbtime(Nt=10^7, No=HC_No, t=HC_t)->HC_db

#In LC
LC_No<-c(13,26,53,106,213,426,853,1706,3413,6826,13653,27306)
LC_t<-c(95,84,63,58,52,48,35,37,35,31,29,26)
dbtime(Nt=10^7, No=LC_No, t=LC_t)->LC_db

#In LCP
LCP_No<-c(13,26,53,106,213,426,853,1706,3413,6826,13653,27306)
LCP_t<-c(86,94,89,87,69,69,37,35,34,31,27,26)
dbtime(Nt=10^7, No=LCP_No, t=LCP_t)->LCP_db

#plot it out
plot(LC_db, type='b', ylim=c(0,5))
points(LCP_db, type='b', col=2)
points(HC_db, type='b', col=3)
