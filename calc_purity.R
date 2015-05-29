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
outBase         <- args[2]
win             <- read.table(inputName, header = T)
observed        <- win$rOsEs

#filter out NA's
observed <- observed[!is.na(observed)]

myerr           <- function(purity) { 
    truth <- (observed - 1 + purity ) / purity
    perr1 <- (1 - truth) ^ 2
    perr2 <- (0.5 - truth) ^ 2 
    perr3 <- (0 - truth) ^ 2 
    perr4 <- (1.5 - truth) ^ 2
    perr5 <- (2.0 - truth) ^ 2
    perr  <- pmin(perr1, perr2, perr3, perr4, perr5)
    terr  <- sqrt(sum(perr))
    return(terr)
}

purvals <- seq(1,100) * 0.01
errvals <- sapply(purvals, myerr)

bestp   <- 0
besterr <- 10000000

for (p in 1:100) {
    er <- myerr(p * 0.01)
    #cat("At ", p, " we have error ", er, "\n")
    if (er < besterr) {
        besterr <- er
        bestp   <- p * 0.01
    }
}

cat("purity_estimate ", bestp, "\n" )

outputPlotName <- paste(outBase, ".purity_error.pdf", sep='')
pdf(outputPlotName)
plot(purvals, errvals, ylim = c(besterr,besterr* 2))
dev.off()

outputPlotName <- paste(outBase, ".purity_evidence.pdf", sep='')
pdf(outputPlotName)
plot(seq(1,length(observed)), (observed - 1 + bestp) / bestp, ylim = c(0,4))
dev.off()

