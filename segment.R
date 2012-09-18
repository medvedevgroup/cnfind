library("DNAcopy")
args            <-commandArgs(TRUE)
normse          <- as.numeric(args[1])
pvalCutoff      <- as.numeric(args[2])
inputName       <- args[3]
outputName      <- args[4]
chr             <- args[5]
ratioColumn     <- args[6]
minCutoff       <- as.numeric(args[7])
maxCutoff       <- as.numeric(args[8])

win             <- read.table(inputName,header=T)
CNA.object      <- CNA( genomdat = win[,ratioColumn], chrom = win$Chr, maploc = win$Start, data.type = 'logratio')
CNA.smoothed    <- smooth.CNA(CNA.object)
segsRaw         <- segment(CNA.smoothed, verbose=0, min.width=2)
segs            <- segsRaw$output
#segs2          <- subset(segs$output, segs$output$seg.mean != "12.345")
names(segs)[names(segs) == "chrom"] <- "Chr"
names(segs)[names(segs) == "loc.start"] <- "Start"
names(segs)[names(segs) == "loc.end"] <- "End"

#fix loc.end -- this is telling the starting location of the last segment, but we really want the ending location of the last segment
RealEnd <- sapply(segs$End, function(x) {win$End[win$Start == x]})
segs$End <- RealEnd

if (normse != 0) {
	normmean         <- 1.0
	segs$seused     <- normse / sqrt(segs$num.mark) 
	#segs$seused     <- normse / sqrt((segs$End - segs$Start) / 1000000)
    segs$Wald       <- ((2^segs$seg.mean) - normmean) / segs$seused
    segs$pval       <- 2 * pnorm(-abs(segs$Wald))
    segs$Call       <- segs$Wald / abs(segs$Wald)
	segs$Call[segs$pval > pvalCutoff] = 0
	segs$seused     <- format(segs$seused, digits=3)
	segs$pval       <- format(segs$pval, digits=3)
	segs$Wald       <- format(segs$Wald, digits=3)
} else {
	segs$seused     <- "NA"
	segs$Wald       <- "NA"
	segs$pval       <- "NA"
	segs$Call       <- 0
	segs$Call[segs$seg.mean > maxCutoff] <- 1 
	segs$Call[segs$seg.mean < minCutoff] <- -1
}

segs$Chr        <- chr
#options(scipen=3) 
write.table(segs, outputName, sep='\t', quote=F, row.names=F)


outputPlotName <- paste(outputName, ".pdf", sep='')
pdf(outputPlotName)
plot(segsRaw, plot.type="w", ylim=c(-2,2))
dev.off()






