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
inputName       <- args[1]
win             <- read.table(inputName, header = T)
observed        <- win$Ratio

myerr           <- function(purity) { 
	truth <- (observed - 1 + purity ) / purity
	perr1 <- (1 - truth) ^ 2
	perr2 <- (0.5 - truth) ^ 2 
	perr3 <- (1.5 - truth) ^ 2
	perr  <- pmin(perr1, perr2, perr3)
	terr  <- sqrt(sum(perr))
	return(terr)
}

xvals <- seq(1,100) * 0.01
plot(xvals, sapply(xvals, myerr), ylim = c(0,..5))







