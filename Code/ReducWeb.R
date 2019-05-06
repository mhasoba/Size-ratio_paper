# Analyze size-ratios in coexisting CR pair data. This script averages the
# CR size ratios across all the resources of each consumer. Replaces Resource name with number of resources

rm(list=ls()) # clear objects
graphics.off() #close all open figures and graphics objects

IDir <- "../Data/"
ODir <- "../Results/"
DataName <- "AllWebCRsHiRes.csv"

Data <- as.matrix(read.table(paste(IDir,DataName,sep = ""), header = TRUE, fill = TRUE,sep = ","))
OutCols = list(NULL,c("Community","General.habitat","Description","Source","Consumer","Resource","D","Cmass","GeoAvRmass","Log10Cmass","Log10Rmass",
                      "log10AvkRC","SDLog10RMass","SELog10RMass","MinLog10RMass" ,"MaxLog10RMass","RanLog10RMass", "MedLog10RMass", "log10MedkRC"))
ReducData <- matrix(data = "",nrow(Data),sapply(OutCols,length)[2],dimnames = OutCols)#Preallocate empty array for output
WebNames <- unique(Data[,"Community"])
ct <- 0
AxFnt <- 2 #axis font type
AxFntSz <- 1;
graphics.off() #close all open figures and graphics objects
for (i in 1:length(WebNames)){
  TmpInds <- which((Data[,"Community"] == WebNames[i]))
  TmpData <- Data[TmpInds,]
  TmpCnames <- as.matrix(paste(TmpData[,"Consumer"],TmpData[,"Cmass"],TmpData[,"D"]))# eparate each consumer by dimensionality 
  TmpConsLst <- unique(TmpCnames)
  ReducData[(ct+1):(ct+length(TmpConsLst)),"Community"] <- TmpData[1,"Community"] #first populate all the fields that apply across consumers 
  ReducData[(ct+1):(ct+length(TmpConsLst)),"General.habitat"] <- TmpData[1,"General.habitat"]
  ReducData[(ct+1):(ct+length(TmpConsLst)),"Description"] <- TmpData[1,"Description"] 
  ReducData[(ct+1):(ct+length(TmpConsLst)),"Source"] <- TmpData[1,"Source"]
  for (j in 1:nrow(TmpConsLst)){ #now for each unique consumer x D combo
    TmpInds <- which(TmpCnames == TmpConsLst[j])
    TmpTmpData <- TmpData[TmpInds,,drop = FALSE] #Extract data for just that consumer
    ReducData[ct+j,"Consumer"] <- TmpTmpData[1,"Consumer"]
    ReducData[ct+j,"Resource"] <- nrow(TmpTmpData)
    ReducData[ct+j,"D"] <- TmpTmpData[1,"D"]
    ReducData[ct+j,"Cmass"] <- TmpTmpData[1,"Cmass"]
    ReducData[ct+j,"GeoAvRmass"] <- exp(mean(log(as.numeric(TmpTmpData[,"Rmass"])))) #geometric mean
    ReducData[ct+j,"Log10Cmass"] <- TmpTmpData[1,"Log10Cmass"]
    ReducData[ct+j,"Log10Rmass"] <- log10(as.numeric(ReducData[ct+j,"GeoAvRmass"]))
    ReducData[ct+j,"log10AvkRC"] <- log10(as.numeric(ReducData[ct+j,"GeoAvRmass"])/as.numeric(ReducData[ct+j,"Cmass"]))
    ReducData[ct+j,"SDLog10RMass"] <- sd(log10(as.numeric(TmpTmpData[,"Rmass"])))
    ReducData[ct+j,"SELog10RMass"] <- as.numeric(ReducData[ct+j,"SDLog10RMass"])/sqrt(as.numeric(ReducData[ct+j,"Resource"]))
    ReducData[ct+j,"MinLog10RMass"] <- log10(min(as.numeric(TmpTmpData[,"Rmass"])))
    ReducData[ct+j,"MaxLog10RMass"] <- log10(max(as.numeric(TmpTmpData[,"Rmass"])))
    ReducData[ct+j,"RanLog10RMass"] <- as.numeric(ReducData[ct+j,"MaxLog10RMass"]) - as.numeric(ReducData[ct+j,"MinLog10RMass"])
    ReducData[ct+j,"MedLog10RMass"] <- median(log10(as.numeric(TmpTmpData[,"Rmass"]))) #median
    ReducData[ct+j,"log10MedkRC"] <- as.numeric(ReducData[ct+j,"MedLog10RMass"]) - as.numeric(ReducData[ct+j,"Log10Cmass"])
  }
  #plotting
  y <- as.numeric(ReducData[(ct+1):(ct+length(TmpConsLst)),"SELog10RMass"])
  y1 <- as.numeric(ReducData[(ct+1):(ct+length(TmpConsLst)),"RanLog10RMass"])
  x <- as.numeric(ReducData[(ct+1):(ct+length(TmpConsLst)),"Log10Cmass"])
  fit <- summary(lm(y ~ x)) #OLS estimate of slope    
  fit1 <- summary(lm(y1 ~ x)) 
  pdf(paste(ODir,DataName,"_ConSiz_Vs_ResSiz_SE_",WebNames[i],".pdf",sep = ""),4, 4) #initialize plot for printing
  plot(x,y); abline(lm(y ~ x))
  mtext(paste("slope =",as.character(format(fit$coefficients[2],digits = 3)), "(",as.character(format(fit$coefficients[2],digits = 3)),")",
              "; R2 = ",as.character(format(fit$r.squared,digits = 2))), side = 3, line = 0,font = AxFnt,cex = AxFntSz*0.75)
  title(WebNames[i])
  dev.off()
  ct <- ct+j
}
ReducData <- ReducData[1:ct,]
write.csv(ReducData,paste(ODir,DataName,"_Reduc.csv",sep = ""),row.names = FALSE)
