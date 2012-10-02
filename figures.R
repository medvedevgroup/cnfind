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

args            <-commandArgs(TRUE)
work_dir        <- args[1]
outputBase      <- args[2]
onpc=FALSE

#chrNames=read.table("~/data/hg19comp/chrom_lengths")$V1
#chrLens = read.table("~/data/hg19comp/chrom_lengths")$V2
#chrNames = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
#chrLens  = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
if (onpc == TRUE){
	priSegmentsFile <- "~/../Desktop/backup/Primary/all.segments"
	metSegmentsFile <- "~/../Desktop/backup/Met/all.segments"
	segSegmentsFile <- "~/../Desktop/backup/Serum/all.segments"
	options(width=130)
} else {
	priSegmentsFile <- paste(work_dir, "/Primary/all.segments", sep="")
	metSegmentsFile <- paste(work_dir, "/Met/all.segments", sep="")
	segSegmentsFile <- paste(work_dir, "/Serum/all.segments", sep="")
	options(width=160)
}

chrNames = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
chrLens  = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
resolution <- 500000
reducedChrLens <- floor(chrLens / resolution)
partialSumChrLens = c(0) 
for (i in 2:(length(chrLens) + 1)) { partialSumChrLens[i] <- partialSumChrLens[i-1] + reducedChrLens[i-1] }
chrName2Num <- function(chr) { if (chr == "chrX") {  23 } else { as.numeric(substr(chr, 4,5)) } }
pos2axis <- function(chr, pos) { partialSumChrLens[chrName2Num(chr)] + floor(pos/resolution) }
pos2axisV <- Vectorize(function(chr, pos) { partialSumChrLens[chrName2Num(chr)] + floor(pos/resolution) })

addPos2Segs <- function (segs) {
	segs$dstart <- pos2axisV(segs$Chr, segs$Start) + 1
	segs$dend   <- pos2axisV(segs$Chr, segs$End)
	segs
}


segFunGain <- function (pos, segs) {
	res =  max(0,segs$CorrRatio[segs$dend >= pos & segs$dstart <= pos]) -1
	if (length(res) == 0) { 
		res = 0	
	} else if (res[1] < 0) { #insertions only
		res = 0
	} else {
		res[1]
	} 
}
segFunLoss <- function (pos, segs) {
	res = segs$CorrRatio[segs$dend >= pos & segs$dstart <= pos]
	if (length(res) == 0) { 
		res = 0	
	} else if (res[1] > 1) { #losses only
		res = 0 
	} else {
		 1 - res[1]
	} 
}

getTripleZero = function (vals) (c(0,vals[1:length(vals) - 1]) == 0) &  (c(vals[2:length(vals)], 0) == 0)

segFunQval <- function (pos, segs, type) {
	hits = segs$dend >= pos & segs$dstart <= pos
	if (length(hits[hits]) == 0) { #if there were no hits
		1
	} else {
		rat = segs$CorrRatio[hits] 
		if (type == "gain" & rat[1] <= 1) {
			1
		} else if (type == "loss" & rat[1] >= 1) {
			1
		} else {
			qval = segs$qval[hits]
			qval[1]
		}
	}
}

#computeQvals <- function(mypval, allpvals) { myqvals <- sum(allpvals[allpvals <= mypval]) }

getSegmentsQval = function(segmentsFile) {
	yvals <- seq(1, pos2axis("chr23",0))
	segs = read.table(segmentsFile, header=T)
	segs = subset(segs, Chr != "chrX")
	segs <- addPos2Segs (segs)
	segs$qval = sapply(segs$pval, function(curpval) { disc=segs$pval[segs$pval <= curpval]; sum(disc) / length(disc) } )
	Gain <- unlist(sapply(yvals, segFunQval, segs, type="gain"))
	Loss <- unlist(sapply(yvals, segFunQval, segs, type="loss"))
	Loss[length(yvals)] = 1
	Gain[length(yvals)] = 1
	list(Gain=Gain, Loss=Loss, Segs=segs)
}


getSegments = function(segmentsFile) {
	yvals <- seq(1, pos2axis("chr23",0))
	segs = read.table(segmentsFile, header=T)
	segs = subset(segs, Chr != "chrX")
	segs <- addPos2Segs (segs)
	Gain <- unlist(sapply(yvals, segFunGain, segs))
	Loss <- unlist(sapply(yvals, segFunLoss, segs))
	Gain[getTripleZero(Gain)] = NA
	Loss[getTripleZero(Loss)] = NA
	Loss[length(yvals)] = 0 * Loss[length(yvals) - 1] #either 0 or NA
	Gain[length(yvals)] = 0 * Gain[length(yvals) - 1] #either 0 or NA
	list(Gain=Gain, Loss=Loss, Segs=segs)
}

setUpPlot = function (yrange) {
	shadeColor <- rgb(173, 216, 230, alpha=80, maxColorValue=255) #light blue, transparent
	plot(c(-1,1), yrange,  axes=FALSE, bty = 'o', type='n', xlab="CN Abs Diff", ylab="")
	abline(v=0); abline(v=0.5); abline(v=-0.5); 
	axis(4, at=partialSumChrLens[1:length(chrNames)], labels=substr(chrNames, 4, 5), las=2)
	axis(1, at=seq(-1, 1, 0.5), labels = c(2, 1, 0, 1, 2))
	axis(2); #useful for debugging
	boxLims = par("usr")
	rect(boxLims[1], partialSumChrLens[seq(1, length(partialSumChrLens), 2) ], boxLims[2], partialSumChrLens[seq(2, length(partialSumChrLens), 2)], col=shadeColor, lty="blank")
	box()
}

setUpQPlot = function (yrange) {
	shadeColor <- rgb(173, 216, 230, alpha=80, maxColorValue=255) #light blue, transparent
	plot(c(-1,1), yrange,  axes=FALSE, bty = 'o', type='n',xlab="", ylab="")
	abline(v=0); abline(v=0.5); 
	axis(4, at=partialSumChrLens[1:length(chrNames)], labels=substr(chrNames, 4, 5), las=2)
	axis(1, at=c(-1, -0.5, 0, 0.5, 1), labels = c(minq, minq/2, 0, 1, 2))
	axis(2); #useful for debugging
	title(xlab="log(q-val)", adj=0.25)
	title(xlab="CN Abs Diff", adj=0.75)
	boxLims = par("usr")
	rect(boxLims[1], partialSumChrLens[seq(1, length(partialSumChrLens), 2) ], boxLims[2], partialSumChrLens[seq(2, length(partialSumChrLens), 2)], col=shadeColor, lty="blank")
	box()
}


yvals <- seq(1, pos2axis("chr23",0))
yrange=rev(range(yvals))
scalef=10

xPri = getSegments(priSegmentsFile)
xMet = getSegments(metSegmentsFile)
xSer = getSegmentsQval(segSegmentsFile)
mergedPri = xPri
mergedPri$Gain[xMet$Gain >= xPri$Gain] = 0
mergedPri$Loss[xMet$Loss >= xPri$Loss] = 0

#for debugging qvals
#so = order(xSer$Segs$pval); plot(1:length(xSer$Segs$pval), xSer$Segs$pval[so], col="blue"); points(1:length(xSer$Segs$pval), xSer$Segs$qval[so], col="red");

if (!onpc) {
	outName <- paste(outputBase, "fig_blood_big.pdf", sep="")
	pdf(outName, paper="US", width=0, height=0)
}
minq = -40
xvalsSerGain = pmax(minq, log10(xSer$Gain)) / (-1 * minq)
xvalsSerLoss = pmax(minq, log10(xSer$Loss)) / (-1 * minq)

setUpQPlot(yrange)
lines(x=xvalsSerGain, y=yvals, col="blue", lwd=2)
lines(x=xvalsSerLoss, y=yvals, col="red", lwd=2)
lines(x=xMet$Gain, y=yvals, col="blue", lwd=2)
lines(x=xMet$Loss, y=yvals, col="red", lwd=2)
lines(x=mergedPri$Gain, y=yvals, col="blue", lwd=2)
lines(x=mergedPri$Loss, y=yvals, col="red", lwd=2)
mtext("Met+Pri", adj=0.75, side=3, line=0, cex=1.2); mtext("Serum", adj=0.25, side=3, line=0, cex=1.2)
if  (!onpc) dev.off()

if (!onpc) {
	outName <- paste(outputBase, "fig_blood_tripanel.pdf", sep="")
	pdf(outName, paper="USr", width=0, height=0)
}

layout( matrix(c(1,2,3), 1)) #side-by-side panels
setUpQPlot(yrange)
lines(x=xvalsSerGain, y=yvals, col="blue", lwd=2)
lines(x=xvalsSerLoss, y=yvals, col="red", lwd=2)
lines(x=xMet$Gain, y=yvals, col="blue", lwd=2)
lines(x=xMet$Loss, y=yvals, col="red", lwd=2)
lines(x=mergedPri$Gain, y=yvals, col="blue", lwd=2)
lines(x=mergedPri$Loss, y=yvals, col="red", lwd=2)
mtext("Met+Pri", adj=0.75, side=3, line=0, cex=1.2); mtext("Serum", adj=0.25, side=3, line=0, cex=1.2)

setUpQPlot(yrange)
lines(x=xvalsSerGain, y=yvals, col="blue", lwd=2)
lines(x=xvalsSerLoss, y=yvals, col="red", lwd=2)
lines(x=xPri$Gain, y=yvals, col="blue", lwd=2)
lines(x=xPri$Loss, y=yvals, col="red", lwd=2)
mtext("Primary", adj=0.75, side=3, line=0, cex=1.2); mtext("Serum", adj=0.25, side=3, line=0, cex=1.2)

setUpQPlot(yrange)
lines(x=xvalsSerGain, y=yvals, col="blue", lwd=2)
lines(x=xvalsSerLoss, y=yvals, col="red", lwd=2)
lines(x=xMet$Gain, y=yvals, col="blue", lwd=2)
lines(x=xMet$Loss, y=yvals, col="red", lwd=2)
mtext("Metastatic", adj=0.75, side=3, line=0, cex=1.2); mtext("Serum", adj=0.25, side=3, line=0, cex=1.2)

dev.off()

if (!onpc) {
	outName <- paste(outputBase, "fig_tumor.pdf", sep="")
	pdf(outName, paper="USr", width=0, height=0)
}
layout( matrix(c(1,2,3), 1)) #side-by-side panels

setUpPlot(yrange)
lines(x=-1 * xMet$Gain, y=yvals, col="blue", lwd=2)
lines(x=xPri$Gain, y=yvals, col="blue", lwd=2)
lines(x=-1 * xMet$Loss, y=yvals, col="red", lwd=2)
lines(x=xPri$Loss, y=yvals, col="red", lwd=2)
mtext("Primary", adj=0.75, side=3, line=0, cex=1.2); mtext("Metastatic", adj=0.25, side=3, line=0, cex=1.2)
title("Gain/Losses"); 

setUpPlot(yrange)
lines(x=-1 *xMet$Gain, y=yvals, col="blue", lwd=2)
lines(x=xPri$Gain, y=yvals, col="blue", lwd=2)
mtext("Primary", adj=0.75, side=3, line=0, cex=1.2); mtext("Metastatic", adj=0.25, side=3, line=0, cex=1.2)
title("Amplifications")

setUpPlot(yrange)
lines(x=-1 * xMet$Loss, y=yvals, col="red", lwd=2)
lines(x=xPri$Loss, y=yvals, col="red", lwd=2)
mtext("Primary", adj=0.75, side=3, line=0, cex=1.2); mtext("Metastatic", adj=0.25, side=3, line=0, cex=1.2)
title("Losses")


dev.off()
