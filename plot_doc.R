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

args           <- commandArgs(TRUE)
chr            <- args[1]
pixel_width    <- args[2]
outBase        <- args[3]
snpFile        <- args[4] #set to nosnps to not plot snps
segmentsFile   <- args[5] #set to nocalls to not plot calls
firstbuffy     <- args[6] #set to use lrOsEs for first sample 
callsFile10    <- args[7] #set to nocalls to not plot calls
callsFile50    <- args[8] #set to nocalls to not plot calls
ColName        <- args[9]
basicArgs      <- 9

numFiles       <- (length(args) - basicArgs ) / 3

outputPlotName <- paste(outBase, ".png", sep='')
png(outputPlotName, width=as.integer(pixel_width), height=400)

#plot first doc 
inputName      <- args[basicArgs + 1]
inputSample    <- args[basicArgs + 2 ]
inputCol       <- args[basicArgs + 3 ]
names          <- c(inputSample)
cols           <- c(inputCol)
tbl            <- read.table(inputName,header=T)
firstColName = ColName
if (firstbuffy == 1) {
	firstColName = "lrOsEs"
} 

plot(tbl$Start, tbl[,firstColName], col=inputCol, xlab = chr, ylab = "LogRatio", ylim = c(-2,2))
curFile <- 2

#plot segment track
if (segmentsFile != "nosegments") {
	segs           <- read.table(segmentsFile, header=T);
	segs$seg.mean  <- segs$seg.mean
	gains          <- segs$Call > 0
	losses         <- segs$Call < 0
	neutral        <- segs$Call == 0
	#segs$seg.mean  <- segs$seg.mean - 1
	#points(tbl$Start, tbl[,ColName] - 1, col="gray") #plots the gray points of the first segment, dup
	segments(segs$Start[gains],   segs$seg.mean[gains],   segs$End[gains],   segs$seg.mean[gains],   col="blue",  lwd=3)
	segments(segs$Start[losses],  segs$seg.mean[losses],  segs$End[losses],  segs$seg.mean[losses],  col="red",   lwd=3)
	segments(segs$Start[neutral], segs$seg.mean[neutral], segs$End[neutral], segs$seg.mean[neutral], col="black", lwd=3)
}




#plot horizontal lines
#abline(a=-1, b=0)
abline(a=0, b=0)
#abline(a=0.5850, b=0)

#plot snps
if (snpFile != "nosnps") {
	snps <- read.table(snpFile, header=T)
	attach(snps) #we get contig, position, bfreq
	ourChr <- contig == chr #index set that identifies points of interest to our chr
	points(position[ourChr], bfreq[ourChr]  - 2, col="black")
}

if (callsFile10 != "nocalls") {
	calls     <- read.table(callsFile10, header=T)
	lossCalls <- calls$Call < 0
	gainCalls <- calls$Call > 0

	yvals     <- rep(-1.9, length(calls$Start[lossCalls]))
	segments(calls$Start[lossCalls], yvals, calls$End[lossCalls], yvals, col="red", lwd=2)

	yvals     <- rep(-1.9, length(calls$Start[gainCalls]))
	segments(calls$Start[gainCalls], yvals, calls$End[gainCalls], yvals, col="blue", lwd=2)
}

if (callsFile50 != "nocalls") {
	calls     <- read.table(callsFile50, header=T)
	lossCalls <- calls$Call < 0 & calls$Chr == chr
	gainCalls <- calls$Call > 0 & calls$Chr == chr

	yvals     <- rep(-2.0, length(calls$Start[lossCalls]))
	segments(calls$Start[lossCalls], yvals, calls$End[lossCalls], yvals, col="red", lwd=2)

	yvals     <- rep(-2.0, length(calls$Start[gainCalls]))
	segments(calls$Start[gainCalls], yvals, calls$End[gainCalls], yvals, col="blue", lwd=2)
}


while (curFile <= numFiles) {
	inputName    <- args[basicArgs - 2 + curFile * 3 ]
	inputSample  <- args[basicArgs - 1 + curFile * 3 ]
	inputCol     <- args[basicArgs - 0 + curFile * 3 ]
	names        <- c(names, inputSample)
	cols         <- c(cols, inputCol)
	tbl          <- read.table(inputName,header=T)
	points(tbl$Start, tbl[,ColName], col=inputCol)

	curFile <- curFile + 1
}


legend(x="topright", legend=names, col=cols, text.col=cols, pch=1 )
