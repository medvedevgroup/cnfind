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
outputName      <- args[1]

winin           <- read.table(f, header = F)
Chr      = winin$V1
Start    = winin$V2
End      = winin$V3
rOsEs    = format( winin$V5 , digits=4)
Observed = winin$V7
Expected = winin$V9
Masked   = winin$V11
ObservedFull = winin$V13
GC = winin$V14
lrOsEs    = format(log2(as.numeric(rOsEs)), digits=4)



winout = data.frame(Chr = Chr, Start = Start, End = End, lrOsEs = lrOsEs, rOsEs = rOsEs, Observed = Observed, Expected = Expected, Masked = Masked, ObservedFull = ObservedFull, GC = GC)
write.table(winout, outputName, sep='\t', quote=F, row.names=F)








