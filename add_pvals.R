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

args            <- commandArgs(TRUE)
f               <- file("stdin")
#normFile        <- args[1]
normmean        <- 1.0
normse          <- as.numeric(args[1])
#callFile        <- args[2] 
#pvalCutoff      <- as.numeric(args[3])

win             <- read.table(f, header = T)
#sumExpected     <- sum(win$Expected)
#sumObserved     <- sum(win$Observed)
#win$ExpectedRaw <- win$Expected
#win$Expected    <- win$Expected * (sumObserved / sumExpected)
#win$RawRatio    <- win$Ratio
#win$RawLogRatio <- win$LogRatio
#win$Ratio       <- win$Observed / win$Expected
#win$LogRatio    <- log2(win$Ratio)


if (normse != 0) {
	Wald            <- (win$Ratio - normmean) / normse
	win$pval        <- 2 * pnorm(-abs(Wald))
	win$Call        <- Wald / abs(Wald)
	win$sdUsed      <- rep(normse, length(win$Call))
	#calls           <- subset(win, win$pval < pvalCutoff)[c("Chr", "Start", "End", "Call")]
	#write.table(calls, callFile, sep='\t', quote=F, row.names = F)
} else {
	win$pval        <- rep("NA", length(win$Ratio))
	win$Call        <- rep("NA", length(win$Ratio))
	win$sdUsed      <- rep("NA", length(win$Ratio))
}

write.table(win, stdout(), sep='\t', quote=F, row.names = F)

#if (normFile != "nonorm") {
#	normwin         <- read.table(normFile, header=T)
#	normmean        <- mean(normwin$Ratio)
#	normse          <- sd(normwin$Ratio)
#	Wald            <- (win$Ratio - normmean) / normse
#	win$pval        <- 2 * pnorm(-abs(Wald))
#	win$Call        <- Wald / abs(Wald)
#	win$sdUsed      <- rep(normse, length(win$Call))
#	write(paste("FYI:", normFile, "has mean", normmean, "and standard error", normse), stderr())
#	calls <- subset(win, win$pval < pvalCutoff)[c("Chr", "Start", "End", "Call")]
#	write.table(calls, callFile, sep='\t', quote=F, row.names = F)
#}







