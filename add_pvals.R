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
normmean        <- 1.0
normse          <- as.numeric(args[1])

win             <- read.table(f, header = T)

if (normse != 0) {
	Wald            <- (win$Ratio - normmean) / normse
	win$pval        <- 2 * pnorm(-abs(Wald))
	win$Call        <- Wald / abs(Wald)
	win$sdUsed      <- rep(normse, length(win$Call))
} else {
	win$pval        <- rep("NA", length(win$Ratio))
	win$Call        <- rep("NA", length(win$Ratio))
	win$sdUsed      <- rep("NA", length(win$Ratio))
}

write.table(win, stdout(), sep='\t', quote=F, row.names = F)


