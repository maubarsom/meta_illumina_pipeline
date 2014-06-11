#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) == 2){
	kmers <- read.table(args[1],header=F,sep=' ',col.names=c("kmer_freq","count"))
	pdf(args[2])
	#plot(kmers$kmer_freq, kmers$count,xlim=c(0,100),ylim=c(0,1000000))
	plot(kmers$kmer_freq, kmers$count,xlim=c(0,100))
	dev.off()
}