# Modify based on https://github.com/rpeckner-broad/Specter/blob/master/SpecterQuant.R

FillCoeffs = function(CoeffList, MaxFillNum){
  if (length(CoeffList) <= 2)
    return (CoeffList)
  for (i in 2:length(CoeffList)){
    if (CoeffList[i]==0 & CoeffList[i-1]>0){
      j = i+1
      while(j<=length(CoeffList) & CoeffList[j]==0 & j-i+1<=MaxFillNum)
        j = j+1
      
      if (j<=length(CoeffList) & CoeffList[j]>0 & j-i<=MaxFillNum){
        growth = (CoeffList[j] - CoeffList[i-1]) / (j-i+1)
        for (k in i:(j-1))
          CoeffList[k] = CoeffList[k-1] + growth
      }
    }
  }
  return (CoeffList)
}

SparkCoeffs = function(rawCoeffsFilepath,ExperimentHeader=NULL, lib=NULL,MS1=FALSE)
{
  DIA_RefSpectraCoeffs = read.csv(rawCoeffsFilepath,header=FALSE,stringsAsFactors = FALSE)
  names(DIA_RefSpectraCoeffs) = c("Coeff","Scan","Sequence", "Charge","PrecursorMZ","RetentionTime",
                                  "Level","Correlation") 
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[c(2,1,3,4,5,6,8)]
  
  DIA_RefSpectraCoeffs$Coeff = as.numeric(DIA_RefSpectraCoeffs$Coeff)
  DIA_RefSpectraCoeffs$Scan = strtoi(DIA_RefSpectraCoeffs$Scan) + 1
  DIA_RefSpectraCoeffs$Sequence = as.character(DIA_RefSpectraCoeffs$Sequence)
  DIA_RefSpectraCoeffs$Charge = strtoi(DIA_RefSpectraCoeffs$Charge)
  DIA_RefSpectraCoeffs$PrecursorMZ = as.numeric(DIA_RefSpectraCoeffs$PrecursorMZ)
  DIA_RefSpectraCoeffs$RetentionTime = as.numeric(DIA_RefSpectraCoeffs$RetentionTime)
  DIA_RefSpectraCoeffs$Correlation = as.numeric(DIA_RefSpectraCoeffs$Correlation)
  
  DIA_RefSpectraCoeffs = na.omit(DIA_RefSpectraCoeffs)
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[with(DIA_RefSpectraCoeffs,order(Scan, Sequence, Charge)),]
  
  if (max(DIA_RefSpectraCoeffs$RetentionTime) < 600)
    DIA_RefSpectraCoeffs$RetentionTime = 60*DIA_RefSpectraCoeffs$RetentionTime
  
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[which(DIA_RefSpectraCoeffs$Coeff != 0),]
  DIA_RefSpectraCoeffs = data.table(DIA_RefSpectraCoeffs)
  setkey(DIA_RefSpectraCoeffs,Sequence,Charge)
  
  return(DIA_RefSpectraCoeffs)
}

QuantifyAllFromCoeffs = function(FIGSData,header,FillCoeffsFlag) {
  IDs = data.frame(unique(FIGSData[,c('Sequence','Charge'),with=FALSE]))
  IsThereEnoughData = unlist(Map(function(i) EnoughData(FIGSData,IDs,i,RTWindow=FALSE),
                                 seq(1,nrow(IDs))))
  
  f = function(k) {if (IsThereEnoughData[k]) {
    QuantifyPeptides(FIGSData,IDs,k,header,FillCoeffsFlag,RTWindow=FALSE)} 
    else {list(0,1,0,0,0,0,1,0,0,0,0,0,0)}}
  
  Quant = Map(f,seq(1,nrow(IDs)))
  Quants = IDs
  
  Quants$Quantity = as.numeric(unlist(lapply(Quant,function(l) l[[1]][1])))
  Quants$LBPval = as.numeric(unlist(lapply(Quant,function(l) l[[2]][1])))
  Quants$SNR = as.numeric(unlist(lapply(Quant,function(l) l[[3]][1])))
  Quants$MaxCoeffTime = as.numeric(unlist(lapply(Quant,function(l) l[[4]][1])))
  Quants$MaxCoeff = as.numeric(unlist(lapply(Quant,function(l) l[[5]][1])))
  Quants$PeakWidth = as.numeric(unlist(lapply(Quant,function(l) l[[6]][1])))
  Quants$LBPvalOnUniformGrid = as.numeric(unlist(lapply(Quant,function(l) l[[7]][1])))
  Quants$PeakVariance = as.numeric(unlist(lapply(Quant,function(l) l[[8]][1])))
  Quants$PeakSkewness = as.numeric(unlist(lapply(Quant,function(l) l[[9]][1])))
  Quants$PeakKurtosis = as.numeric(unlist(lapply(Quant,function(l) l[[10]][1])))
  Quants$PeakStart = as.numeric(unlist(lapply(Quant,function(l) l[[11]][1])))
  Quants$PeakEnd = as.numeric(unlist(lapply(Quant,function(l) l[[12]][1])))
  Quants$Correlation = as.numeric(unlist(lapply(Quant,function(l) l[[13]][1])))
  
  return(Quants)
}

FindPeaksInData = function(Data,Identifiers,i,header,FillCoeffsFlag,maxFillNum=1,Multiplexed=FALSE,IntensityCutoff=0,
                           QuantileCutoff=0,RTWindow=TRUE,smooth="rollmean",FilterWindow= 3,KZiters=3) {
  PeptideData = data.frame(Data[list(Identifiers$Sequence[i],Identifiers$Charge[i])])
  
  if (RTWindow)
    PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
  
  PeptideData = PeptideData[c("Scan","RetentionTime","Coeff","PrecursorMZ")]
  PeptideData = PeptideData[!duplicated(PeptideData$RetentionTime),]
  PeptideData$Coeff[which(PeptideData$Coeff < 1)] = 0
  
  PeptidePeaks = NULL
  if (length(which(PeptideData$Coeff > 0)) > 3)
  {
    MinScanIndex = min(PeptideData$Scan,na.rm=TRUE)
    MaxScanIndex = max(PeptideData$Scan,na.rm=TRUE)
    
    HeaderRange = data.frame(header[J(MinScanIndex:MaxScanIndex)])
    HeaderRange = subset(HeaderRange,msLevel == 2 & floor(precursorMZ) %in% floor(unique(PeptideData$PrecursorMZ)))
    
    DataOnUniformGrid = HeaderRange[c("seqNum","retentionTime")]
    DataOnUniformGrid$Coeff = numeric(nrow(DataOnUniformGrid))
    DataOnUniformGrid$Coeff[which(DataOnUniformGrid$seqNum %in% PeptideData$Scan)] = PeptideData$Coeff  
    
    DataOnUniformGrid = DataOnUniformGrid$Coeff
    
    if (FillCoeffsFlag)
    {
      DataOnUniformGrid = FillCoeffs(DataOnUniformGrid, maxFillNum)
    }
    PeptidePeaks = findpeaks(DataOnUniformGrid, nups = 2, ndowns = 2, npeaks= 10,sortstr = TRUE)
    
    if (smooth == "rollmean")
    {
      KZData = kz(DataOnUniformGrid,m=FilterWindow,k=KZiters)
      KZData[which(DataOnUniformGrid < 1)] = 0
      
      if (FillCoeffsFlag)
      {
        KZData = FillCoeffs(KZData, maxFillNum)
      }
      PeptidePeaks = findpeaks(KZData, nups = 2, ndowns = 2, npeaks= 10,sortstr = TRUE)
    }
    
  }
  return(PeptidePeaks)
}

QuantifyPeptides = function(Data,Identifiers,i,header,FillCoeffsFlag,maxFillNum=1,IntensityCutoff=0,QuantileCutoff=0,RTWindow=TRUE,smooth="rollmean",FilterWindow= 3,KZiters=3) {
  PeptideData = data.frame(Data[list(Identifiers$Sequence[i],Identifiers$Charge[i])]) 
  if (RTWindow)
    PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
  
  PeptideData = PeptideData[c("Scan","RetentionTime","Coeff","PrecursorMZ","Correlation")]
  PeptideData = PeptideData[!duplicated(PeptideData$RetentionTime),]
  PeptideData$Coeff[which(PeptideData$Coeff < 1)] = 0  
  
  area = 0
  
  RawPeaks = FindPeaksInData(Data,Identifiers,i,header,FillCoeffsFlag,maxFillNum=maxFillNum,IntensityCutoff=IntensityCutoff,
                             QuantileCutoff=QuantileCutoff,RTWindow=RTWindow,smooth="none",FilterWindow=FilterWindow)
  
  SmoothPeaks = FindPeaksInData(Data,Identifiers,i,header,FillCoeffsFlag,maxFillNum=maxFillNum,IntensityCutoff=IntensityCutoff,
                                QuantileCutoff=QuantileCutoff,RTWindow=RTWindow,smooth=smooth,FilterWindow=FilterWindow,KZiters=KZiters)
  
  BoxTestPval = 1
  BoxTestPvalOnGrid = 1
  MaxCoeff= 0
  TimeAtTopOfPeak = 0
  snr = 0
  PeakCentralMoments = numeric(3)
  result = list(Area=area,BoxTestPval=BoxTestPval,SNR=snr,RT=TimeAtTopOfPeak,MaxCoeff=0,PeakWidth=0,BoxTestPvalOnGrid=BoxTestPvalOnGrid,
                Variance=0,Skewness=0,Kurtosis=0,PeakStart=0,PeakEnd=0,
                Correlation=0)
  
  if ((length(which(PeptideData$Coeff > 1)) > 3) & !is.null(SmoothPeaks)) {
    MinScanIndex = min(PeptideData$Scan,na.rm=TRUE)
    MaxScanIndex = max(PeptideData$Scan,na.rm=TRUE)
    
    HeaderRange = data.frame(header[J(MinScanIndex:MaxScanIndex)]) 
    HeaderRange = subset(HeaderRange,msLevel == 2 & floor(precursorMZ) %in% floor(unique(PeptideData$PrecursorMZ)))
    
    DataOnUniformGrid = HeaderRange[c("seqNum","retentionTime")]
    DataOnUniformGrid$Coeff = numeric(nrow(DataOnUniformGrid))
    DataOnUniformGrid$Coeff[which(DataOnUniformGrid$seqNum %in% PeptideData$Scan)] = PeptideData$Coeff
    
    if (FillCoeffsFlag)
    {
      DataOnUniformGrid$Coeff = FillCoeffs(DataOnUniformGrid$Coeff, maxFillNum)
    }
    
    start = SmoothPeaks[1,3]
    end = SmoothPeaks[1,4]
    
    PeakDataOnUniformGrid = DataOnUniformGrid[start:end,]
    
    area = trapz(PeakDataOnUniformGrid$retentionTime, PeakDataOnUniformGrid$Coeff) 
    
    TimeAtTopOfPeak = DataOnUniformGrid$retentionTime[SmoothPeaks[1,2]]
    MaxCoeff=DataOnUniformGrid$Coeff[SmoothPeaks[1,2]]
    
    PeakCentralMoments = all.moments(scale(PeakDataOnUniformGrid$Coeff,center=FALSE),order.max = 4,central = TRUE)[3:5]
    
    BoxTest = Box.test(PeptideData$Coeff,type="Ljung-Box")
    BoxTestPval = round(BoxTest$p.value,digits=5)
    
    BoxTestOnGrid = Box.test(DataOnUniformGrid$Coeff,type="Ljung-Box")
    BoxTestPvalOnGrid = round(BoxTestOnGrid$p.value,digits=5)
    snr=area/sd(DataOnUniformGrid$Coeff[-seq(start,end)])
    
    peakWidth = DataOnUniformGrid$retentionTime[end] - DataOnUniformGrid$retentionTime[start]
    
    Correlation = mean(PeptideData$Correlation)
    
    result = list(Area=area,BoxTestPval=BoxTestPval,SNR=snr,RT=TimeAtTopOfPeak,MaxCoeff=MaxCoeff,PeakWidth=peakWidth,
                  BoxTestPvalOnGrid=BoxTestPvalOnGrid,Variance=PeakCentralMoments[1],
                  Skewness=PeakCentralMoments[2],Kurtosis=PeakCentralMoments[3],
                  PeakStart=DataOnUniformGrid$retentionTime[start],PeakEnd=DataOnUniformGrid$retentionTime[end],
                  Correlation=Correlation)
  }
  
  return(result)
} 

EnoughData = function(Data,Identifiers,i,RTWindow=TRUE) {
  PeptideData = Data[list(Identifiers$Sequence[i],Identifiers$Charge[i])]
  
  if (RTWindow)
    PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
  
  if (nrow(PeptideData) > 0)
  {return(TRUE)} else {
    return(FALSE)
  }
}

library(kza)
library(pracma)
library(zoo)
library(moments)
library(data.table)
library(MASS)

args=commandArgs(TRUE)
dirPath = args[1]     # The path to "Coeffs.csv"
FillCoeffsFlag = FALSE

headerPath = strsplit(dirPath,split="/")[[1]]
headerPath = headerPath[-length(headerPath)]
headerPath = paste(headerPath,collapse="/")
headerPath = file.path(headerPath, "header.csv")

resultsPath = file.path(dirPath, "Coeffs.csv")
decoyResultsPath = file.path(dirPath, "DecoyCoeffs.csv")

h = read.csv(headerPath,header=FALSE,stringsAsFactors=FALSE)
names(h) = c("precursorMZ","retentionTime","seqNum")
h$retentionTime = 60*h$retentionTime
h$seqNum = h$seqNum + 1
h$msLevel = 2
h = data.table(h)
setkey(h,seqNum)

results = SparkCoeffs(resultsPath)
resultsWithDecoys = SparkCoeffs(decoyResultsPath)

quants = QuantifyAllFromCoeffs(results,h,FillCoeffsFlag)
quantsWithDecoys = QuantifyAllFromCoeffs(resultsWithDecoys,h,FillCoeffsFlag)

quantsWithDecoys = rbind(quants,quantsWithDecoys)

quantsWithDecoys$Type = ifelse(grepl("DECOY",quantsWithDecoys$Sequence),"Decoy","Target")
quantsWithDecoys = quantsWithDecoys[which(quantsWithDecoys$Quantity > 0),]

if (length(which(quantsWithDecoys$Type == "Decoy")) > 0) {                                                 
  FIGSLDA = lda(Type~ PeakVariance + PeakSkewness + PeakKurtosis + Correlation,
                data = quantsWithDecoys)
  
  plda = predict(object = FIGSLDA, newdata = quantsWithDecoys)
  
  D = data.frame(Sequence = quantsWithDecoys$Sequence, 
                 Charge = quantsWithDecoys$Charge, 
                 Type=quantsWithDecoys$Type,
                 Score=as.numeric(unlist(plda$x)),
                 
                 Quantity = quantsWithDecoys$Quantity,
                 PeakStart = quantsWithDecoys$PeakStart,
                 PeakEnd = quantsWithDecoys$PeakEnd
  )
  
  FDR = function(score)  {
    return(length(which(D$Type == "Decoy" & D$Score > score))/length(which(D$Score > score)))
  }
  
  fdr = as.numeric(sapply(D$Score,FDR))
  
  cutoff = min(D$Score[which(fdr < 0.01)],na.rm=TRUE)
  
  D = D[which(D$Score >= cutoff),]
  quants = quants[which(paste(quants$Sequence,quants$Charge) %in% paste(D$Sequence,D$Charge)),]
  quants$Score = D$Score[match(paste(quants$Sequence,quants$Charge),paste(D$Sequence,D$Charge))]
} else {
  quants$Score = 1
}

outputPath = gsub("Coeffs","Quants",resultsPath)
write.csv(quants[c("Sequence","Charge","Quantity","Score","PeakStart","PeakEnd")],file=outputPath,quote=FALSE,row.names=FALSE)
