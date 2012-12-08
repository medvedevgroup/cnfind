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
f               <- file("stdin")
normalWinName   <- args[1]
sam2normScaling <- as.numeric(args[2])
outputName      <- args[3]

win             <- read.table(f, header = T)

win$ObsN  = NA
win$rOnEn = NA
win$lrOsOn = NA
win$rOsOn = NA
win$lrRsRn = NA
win$rRsRn = NA


#calculate covRatio between normal and sample
if (normalWinName != "nonorm") {
	chr = levels(win$Chr[1])[1] #help me dear g-d
	normalWin       <- read.table(normalWinName,header=T)
	win$ObsN = normalWin$Observed[normalWin$Chr == chr]
	win$rOnEn = normalWin$rOsEs[normalWin$Chr == chr]
	win$rOsOn = format(win$Observed / (sam2normScaling * win$ObsN), digits=4)
	win$lrOsOn = format(log2(as.numeric(win$rOsOn)), digits=4)
	win$rRsRn = format(win$rOsEs / win$rOnEn, digits=4)
	win$lrRsRn = format(log2(as.numeric(win$rRsRn)), digits=4)
}
 
write.table(win, outputName, sep='\t', quote=F, row.names=F)








