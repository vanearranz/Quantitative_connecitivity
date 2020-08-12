LoadSpecies<- function(nameSpecies,dataset,simulation, simulationIterations){
  
  data<-read.csv(dataset,header=TRUE)
  
  #Add a column with the name of the species
  data <-cbind(data,Species=nameSpecies)
  
  #Get number of rows of the file
  nRows <-nrow(data)
  
  #Get minimum and maximum
  
  nMin <-min(data$StartPoint)
  nMax <-max(data$EndPoint)
  
  
  #THESE ARE THE PARAMETERS YOU HAVE TO CHANGE
  #set number of Test and iterations per Test
  #nSimulation = How many simulations you want to do
  #nSimulationIterations = How many iterations per simulation
  
  nSimulation <- simulation
  nSimulationIterations <- simulationIterations
  
  
  #Create a matrix to hold the simulation outputs and the random numbers we use to testing so we can manually check the output
  aSimulationOutput<-matrix(99999,nSimulation,nSimulationIterations)
  aSimulationNumbers<-matrix(99999,nSimulation,nSimulationIterations)
  
  for (x in 1:nSimulation){
    aSegmentBucket<-c()
    aRandomNumbers<-c()
    y<-1
    while(y <= nSimulationIterations){
      #Get random number in between min and max
      nRandom <-runif(1,nMin,nMax)
      aSegmentFound<-which(data$EndPoint<=nRandom)
      
      
      #if length of aSegmentFound is Zero means the random number falls into the 1st segment
      if (length(aSegmentFound)>0){
        nSegment <-max(aSegmentFound) + 1
      }
      else{
        nSegment<-1
      }
      
      #Check that the segment does not exists in the array aSegmentBucket, this will ensure to have ALWAYS different segments
      if (!(nSegment %in% aSegmentBucket)){
        aSegmentBucket <-c(aSegmentBucket,nSegment)
        aRandomNumbers <-c(aRandomNumbers,nRandom)
        y<-y+1
      }
      
    }
    
    
    aSimulationOutput[x,]<-c(aSegmentBucket)
    aSimulationNumbers[x,]<-c(aRandomNumbers)
    
  }
  
  
  ##OUTPUT THE RANDOM RUMBERS FOR SIMULATION
  write.csv(aSimulationOutput,file=paste0(getwd(),"/output/",nameSpecies,"SegmentOutput.csv"))
  write.csv(aSimulationNumbers,file=paste0(getwd(),"/output/",nameSpecies,"RandomNumbers.csv"))
  
  
  Segments<-read.csv("Species/Segments.csv",header = TRUE)
  oSegments<-Segments
  
  
  for (x in 1:nSimulation){
    aResults<-c()
    
    for (y in 1:nSimulationIterations){
      
      oSegmentData<-data[which(data$Segment==aSimulationOutput[x,y]),]
      aResults$a[between(oSegments$EndPoint,oSegmentData$StartPoint,oSegmentData$EndPoint)==TRUE]<-1
      aResults$a[is.na(aResults$a)]<-0
      names(aResults)[names(aResults)=="a"]<-paste0("Simulation_",nSimulationIterations)
      
    }
    
    dfSimulation<-as.data.frame(Reduce("+",aResults))
    colnames(dfSimulation)<-paste0("Simulation_",x)
    oSegments<-cbind(oSegments,dfSimulation)
    
  }
  
  write.csv(oSegments,file=paste0(getwd(),"/output/","Segments_",nameSpecies,".csv"))
  
  ##CLEAN ENVIRONMENT
  rm("data","nRows","nMin","nMax","nSimulation","nSimulationIterations","aSimulationOutput","aSimulationNumbers","aSegmentBucket",
     "aRandomNumbers","y",
     "nRandom","aSegmentFound","nSegment","oSegmentData","oSegments",
     "dfSimulation","aResults")
  
}



####################################################

pelc.significant.func <- function(mat=Bmat, rvec=realvec, siglevel=c(0.9, 0.95, 0.99)){
  ## pelc.significant.func
  ## Looks for segments (1 to 92) with a significant number of breaks among the 21 species of interest.
  ##
  ## EXAMPLE:
  ## pelc.significant.func()
  
  ## Check dimensions and entries.  I'm going to hard-code #segments=92 and #species=21 here for checks.
  if(nrow(mat)!=92) stop("I'm expecting a matrix of bootstrap replicates with 92 rows, one for each segment.")
  if(any(mat!=trunc(mat)))stop("All the entries of the bootstrap matrix should be integers.  Something's wrong?")
  if(!all(mat>=0 & mat <=21)) stop("All entries of bootstrap matrix should be between 0 and 21 species!")
  if(any(rvec!=trunc(rvec)))stop("All the entries of rvec should be integers.  Something's wrong?")
  if(!all(rvec>=0 & rvec <=21)) stop("All entries of rvec should be between 0 and 21 species!")
  
  ## First of all, find the 90%, 95%, and 99% quantiles of the bootstrap replicates for each location.
  ## We want to get out a matrix with 92 rows and 3 columns, where the 3 columns give the 90, 95, and 99%
  ## values of the bootstrap replicates.  Call this quant.mat.
  quant.mat <- t(apply(mat, 1, function(x)quantile(x, prob=siglevel)))
  
  ## We now need to compare the real data with the numbers in quant.mat.
  ## For each of the 92 sites, we want to know if it is significant at the 10%, 5%, and 1% levels.
  ## In other words, whether the real data value for that segment exceeds the 90%, 95%, and 99% level
  ## that we've found in quant.mat.
  ## We can make a new matrix called sig.mat which has 92 rows and 3 columns.
  ## The first column of sig.mat gives all the segments that are significant at the 10% level
  ## (weakest evidence of breaks): that's the lightest grey shading in Pelc et al Fig 2.
  ## The second column gives the segments significant at the 5% level (mid-grey in Pelc et al) and the third
  ## column gives the segments that are significant at the 1% level (strongest evidence of a break).
  sig.mat <- sapply(1:length(siglevel), function(i) as.numeric(realvec > quant.mat[,i]))
  ## Change sig.mat into a data frame with meaningful titles:
  sig.mat <- data.frame(Segment=1:92, sig.mat)
  names(sig.mat) <- c("Segment", paste0("Sig", (1-siglevel)*100, "%"))
  
  list(quant=quant.mat, sig=sig.mat)
}

##################################################

sig.plot <- function(sigres, segdat=seg.dat){
  ## sig.plot
  ## EXAMPLE:
  ## res.out <- pelc.significant.func()
  ## sig.plot(res.out$sig, seg.dat)
  
  ## Number of segments:
  nseg <- nrow(segdat)
  
  ## Append the first row onto the bottom of segdat:
  segdat <- rbind(segdat, segdat[1,])
  
  ## Set up a vector of colours for the nseg segments:
  ## grey is not significant;
  ## yellow is 10%; orange is 5%; red is 1%.
  colvec <- rep("grey", nseg)
  ## Replace segments significant at the 10% level with yellow:
  colvec[sigres[["Sig10%"]]==1] <- "yellow"
  ## Replace segments significant at the 5% level with orange:
  colvec[sigres[["Sig5%"]]==1] <- "orange"
  ## Replace segments significant at the 1% level with red:
  colvec[sigres[["Sig1%"]]==1] <- "red"
  
  ## Plot the map:
  plot(-1, xlab="Longitude", ylab="Latitude", type="n", xlim=range(segdat$Longitude)+c(-0.1, 0.1),
       ylim=range(segdat$Latitude)+c(-0.1, 0.1))
  ## Legend
  legend("bottomright", pt.bg=c("yellow", "orange", "red"), pch=22, 
         col = par("col"), legend=c("10%", "5%", "1%"), cex=2, pt.cex=3)
  ## Go through segments one at a time and plot:
  for(i in 1:nseg){
    lines(segdat$Longitude[c(i, i+1)], segdat$Latitude[c(i, i+1)], col=1, lwd=9)
    lines(segdat$Longitude[c(i, i+1)], segdat$Latitude[c(i, i+1)], col=colvec[i], lwd=8)
    
  }
  
  return(data.frame(segdat[1:nseg,], colour=colvec))
}


###############################################

sig.plot2 <- function(sigres, segdat=seg.dat){
  ## sig.plot
  ## EXAMPLE:
  ## res.out <- pelc.significant.func()
  ## sig.plot(res.out$sig, seg.dat)
  
  ## Number of segments:
  nseg <- nrow(segdat)
  
  ## Append the first row onto the bottom of segdat:
  segdat <- rbind(segdat, segdat[1,])
  
  ## Set up a vector of colours for the nseg segments:
  ## grey is not significant;
  ## yellow is 10%; orange is 5%; red is 1%.
  colvec <- rep("grey", nseg)
  ## Replace segments significant at the 10% level with yellow:
  colvec[sigres[["Sig20%"]]==1] <- "yellow"
  ## Replace segments significant at the 5% level with orange:
  colvec[sigres[["Sig10%"]]==1] <- "orange"
  ## Replace segments significant at the 1% level with red:
  colvec[sigres[["Sig5%"]]==1] <- "red"
  
  ## Plot the map:
  plot(-1, xlab="Longitude", ylab="Latitude", type="n", xlim=range(segdat$Longitude)+c(-0.1, 0.1),
       ylim=range(segdat$Latitude)+c(-0.1, 0.1))
  ## Legend
  legend("bottomright", pt.bg=c("yellow", "orange", "red"), pch=22, 
         col = par("col"), legend=c("20%", "10%", "5%"), cex=2, pt.cex=3)
  ## Go through segments one at a time and plot:
  for(i in 1:nseg){
    lines(segdat$Longitude[c(i, i+1)], segdat$Latitude[c(i, i+1)], col=1, lwd=9)
    lines(segdat$Longitude[c(i, i+1)], segdat$Latitude[c(i, i+1)], col=colvec[i], lwd=8)
    
  }
  
  return(data.frame(segdat[1:nseg,], colour=colvec))
}
