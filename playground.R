library("DNAcopy")
args<-commandArgs(TRUE)
inputName<-args[1]
pixelWidth<-args[2]
pnts <- read.table(inputName,header=F)

outputPlotName <- paste(inputName,".png", sep='')
png(outputPlotName, width=as.numeric(pixelWidth))


plot(pnts[,1], pnts[,2])

#write.table(segs2[,2:6], file=outputName, row.names=T, col.names=T, quote=F, sep="\t")


