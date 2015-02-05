pr.reformat <- function(file){

  #read the input file
  lines = readLines(file)

  lines[length(lines)+1] <- ''						#add blank line to end
  lines <- sapply(lines,sub,pattern='[[:space:]]+$',replacement='') 	#trim tariling whitespace
  if(length(grep("^Read", lines[1]))==0){
    match=lines[1]
    lines <- sapply(lines, sub, pattern=paste('^',match,sep=''), replacement=paste('Read',match,sep=' '))
  }
  names(lines)<-NULL
  splits <- strsplit(lines, split="\t")					#Split each line by tabs

  out <- vector()

  if(lines[[2]] == ''){
    # Input is in Format 1

    start <- which(sapply(lines,grep,pattern='Read')==1)
    stop <- c(start[-1]-1, length(lines))

    type=vector()
    for(i in start){
      type=c(type, strsplit(lines[i],'[ :]')[[1]][2])
    }

    names <- c('Measurment', 'Time', 'Type', splits[[3]][2:length(splits[[3]])])

    for(i in length(start):1){
      for(j in (start[i]+3):stop[i]){
        #Add a row of data
        if(grep("^[[:digit:]]{1}:[[:digit:]]{2}:[[:digit:]]{2}",lines[[j]]) && any(splits[[j]][-1] != '')){
          tmp <- c((j-start[i]-2), splits[[j]][1], paste('Read',type[i]), splits[[j]][2:length(splits[[j]])])
          if(length(tmp)==length(names)){
            out=rbind(out,tmp)
          }else{
            out=rbind(out,c(tmp,rep('',times=(length(names)-length(tmp)))))
          }
        }
      }
    }

  }else{
    # Input is in Format 2

    i <- 1
    while(i <= length(lines)){
      
      if(length(grep("^Read", lines[i]))!=0){
        type=paste('Read',strsplit(lines[i],'[ :]')[[1]][2])
        i=i+1
        measure=strsplit(lines[i],' ')[[1]][3]
        time=strsplit(lines[i],'[()]')[[1]][2]
        i=i+1
        cn=splits[[i]][-1]
        i=i+1
        tmp=vector()
        rn=vector()
        while(lines[i] != ''){
          rn=c(rn,splits[[i]][1])
          tmp=rbind(tmp,splits[[i]][2:(length(splits[[i]])-1)])
          i=i+1
        }
        dim(tmp)<-NULL
        if(any(tmp != '')){
          out=rbind(out,c(measure, time, type, tmp))
        }
      }
      
      i <- i+1
    }

    names=c('Measurment', 'Time', 'Type')
    for(i in 1:length(cn)){
      for(j in 1:length(rn)){
        names=c(names, paste(rn[j],cn[i],sep=''))
      }
    }

  }

  l=list()
  l[[1]]=as.numeric(out[,1])
  l[[2]]=as.character(out[,2])
  l[[3]]=as.character(out[,3])
  for(i in 4:dim(out)[2]){
    l[[i]]=as.numeric(out[,i])
  }
  ret=as.data.frame(l,stringsAsFactors=FALSE)
  names(ret)<-names
  ret=ret[order(ret$Measurment),]
  
  return(ret)
}
