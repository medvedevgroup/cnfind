
#plots all the gc bins
flatchroms <-  c("chr5", "chr6", "chr8", "chr9", "chr10","chr12","chr13","chr14","chr19","chr22")
autosomes  <-  c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13","chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")



png("bins.png", width=1000)
gc <- read.table("gcsummary", header = T)
xvals <- seq(1,40) * 2.5
plot(xvals,gc$chr1, ylim = c(0,.13))
sapply(autosomes, function(chr) points(x=xvals, y=gc[,chr]))
dev.off()

