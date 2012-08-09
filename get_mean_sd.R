#!/home/pmedvedev/R-2.15.0/bin/Rscript
args            <- commandArgs(TRUE)
f               <- file("stdin")
column          <- args[1]
win             <- read.table(f, header = T)
mymean          <- mean(win[,column])
mysd            <- sd  (win[,column])
write(paste(mymean, mysd), stdout())

