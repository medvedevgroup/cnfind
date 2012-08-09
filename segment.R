#!/home/pmedvedev/R-2.15.0/bin/Rscript
library("DNAcopy")
args            <-commandArgs(TRUE)
normse          <- as.numeric(args[1])
pvalCutoff      <- as.numeric(args[2])
inputName       <- args[3]
outputName      <- args[4]
chr             <- args[5]
ratioColumn     <- args[6]


win             <- read.table(inputName,header=T)
CNA.object      <- CNA( genomdat = win[,ratioColumn], chrom = win$Chr, maploc = win$Start, data.type = 'logratio')
CNA.smoothed    <- smooth.CNA(CNA.object)
segs            <- segment(CNA.smoothed, verbose=0, min.width=2)
#segs            <- segment(CNA.object, verbose=0, min.width=2)
segs2           <- subset(segs$output, segs$output$seg.mean != "12.345")
names(segs2)[names(segs2) == "chrom"] <- "Chr"
names(segs2)[names(segs2) == "loc.start"] <- "Start"
names(segs2)[names(segs2) == "loc.end"] <- "End"



if (normse != 0) {
	normmean         <- 1.0
	segs2$seused     <- normse / sqrt((segs2$End - segs2$Start) / 1000000)
    segs2$Wald       <- ((2^segs2$seg.mean) - normmean) / segs2$seused
    segs2$pval       <- 2 * pnorm(-abs(segs2$Wald))
    segs2$Call       <- segs2$Wald / abs(segs2$Wald)
	segs2$Call[segs2$pval > pvalCutoff] = 0
	segs2$Chr        <- chr
}

write.table(segs2, outputName, sep='\t', quote=F, row.names = F)


outputPlotName <- paste(outputName, ".pdf", sep='')
pdf(outputPlotName)
plot(segs, plot.type="w")
dev.off()






