#Copyright 2012, Paul Medvedev
#
#This file is part of Cnfind
#
#Cnfind is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Cnfind is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with cnfind (see file gpl.txt). If not, see <http://www.gnu.org/licenses/>.
#

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
merge_segs      <- as.numeric(args[9])

win             <- read.table(inputName,header=T)


CNA.gdat        <- win[,ratioColumn]
CNA.object      <- CNA( genomdat = CNA.gdat, chrom = win$Chr, maploc = win$Start, data.type = 'logratio')
CNA.smoothed    <- smooth.CNA(CNA.object)

if (merge_segs != 0) {
	CNA.sd          <- sd(CNA.gdat - append(CNA.gdat[2:length(CNA.gdat)], 0))
	segsRaw         <- segment(CNA.smoothed, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=0.3*sqrt(2.0)/CNA.sd)
} else {
	segsRaw         <- segment(CNA.smoothed, verbose=0, min.width=2)

}

segs            <- segsRaw$output
names(segs)[names(segs) == "chrom"] <- "Chr"
names(segs)[names(segs) == "loc.start"] <- "Start"
names(segs)[names(segs) == "loc.end"] <- "End"

#fix loc.end -- this is telling the starting location of the last segment, but we really want the ending location of the last segment
LastMarker      <- segs$End
RealEnd         <- sapply(segs$End, function(x) {win$End[win$Start == x]})
segs$End        <- RealEnd

if (normse != 0) {
	normmean        <- 1.0
	segs$seused     <- normse / sqrt(segs$num.mark) 
	#segs$seused     <- normse / sqrt((segs$End - segs$Start) / 1000000)
	segs$Wald       <- ((2^segs$seg.mean) - normmean) / segs$seused
	segs$pval       <- 2 * pnorm(-abs(segs$Wald))
	segs$Call       <- segs$Wald / abs(segs$Wald)
	segs$Call[segs$pval > pvalCutoff] <- 0
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

#add rmse collumn
temp_rmse <- c()
for (i in 1:nrow(segs)) {
	st <- segs[i,]$Start
	en   <- segs[i,]$End
	ratmean <- 2 ^ segs[i,]$seg.mean
	nummark <- segs[i,]$num.mark
	this_rmse <- sqrt(sum(((win$Ratio[win$End <= en & win$Start >= st]) - ratmean) ^ 2) / nummark)
	temp_rmse <- c(temp_rmse, this_rmse)
}
segs$rmse <- format(temp_rmse, digits=3)
segs$LastMarker <- LastMarker



write.table(segs, outputName, sep='\t', quote=F, row.names=F)


outputPlotName <- paste(outputName, ".pdf", sep='')
pdf(outputPlotName)
plot(segsRaw, plot.type="w", ylim=c(-2,2))
dev.off()






