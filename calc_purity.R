args            <- commandArgs(TRUE)
#f               <- file("stdin")
inputName<-args[1]
win             <- read.table(inputName, header = T)
observed <- win$RawRatio

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







