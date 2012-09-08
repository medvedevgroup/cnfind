bu <- read.table("/home/pmedvedev/blood/data/KIBurun/all.win", header=T)
se <- read.table("/home/pmedvedev/blood/data/KISerun/all.win", header=T)
flatchroms <-  c("chr5", "chr6", "chr8", "chr9", "chr10","chr12","chr13","chr14","chr19","chr22")
allchroms <- levels(se$Chr)[-24]
autosomes <- levels(se$Chr)[c(-23,-24)]


chr <- function(tbl, contig) tbl$Chr == contig
mymean <- function(table, contig) { sum(table$Observed[chr(table, contig)]) / sum(table$Expected[chr(table,contig)]) }
means <- function(tbl, chroms) {sapply(chroms, mymean, table=tbl)}

png("means.png", width=1000)
#par(mfrow = c(1,2))
plot(means(se,autosomes[-1]), means(bu,autosomes[-1]), asp=1)
points(means(se, flatchroms), means(bu,flatchroms), col="red")
abline(0.99,0); abline(1.01,0); abline(1,0); 
abline(v=1.0); abline(v=1.02); abline(v=1.01)
dev.off()

