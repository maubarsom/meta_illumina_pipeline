#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
library(ggplot2)

if(length(args) == 2){
	kmers <- read.table(args[1],header=F,sep=' ',col.names=c("kmer_freq","count"))
	pdf(args[2])
	ggplot(kmers, aes(x=kmer_freq, y=count)) + geom_line() + xlim(c(0,100))
	dev.off()
}
